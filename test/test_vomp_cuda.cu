#include <boost/lexical_cast.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <gme_GaussianMixtureEstimator.h>

#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <rofl/common/macros.h>
#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

// #include <gme_GaussianMixtureEstimator.h>
#include <bin_utils.h>
#include <dsd_utils.h>
#include <gme_gaussian_metric.h>

// CUDA
#include <dsd_utils.cuh>
#include <device_launch_parameters.h>
#include <thrust/host_vector.h>

int main(int argc, char **argv)
{
    std::string filenameCfg, filePath;
    std::string filenameIn, filenameOutIse, filenameOutRtc;
    std::set<fs::path> sortedByName;

    double sigmaIn;
    std::vector<double> sigmasIn;
    std::vector<double> weightsIn;
    bool enableGmmHier;
    bool plotEllipsesArrows;

    double sigmaMin, res, iseThresh;

    double angleWin;
    int nSamples;
    bool enablePclVisualization;

    int numThreads;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", filenameIn,
        std::string("/Geode/Urban_tunnel_01/LiDAR/bin/"));
    params.getParam<std::string>("out_ise", filenameOutIse,
                                 std::string("ise.csv"));
    params.getParam<std::string>("out_rtc", filenameOutRtc,
                                 std::string("rtc.csv"));
    params.getParam<double>("sigmaIn", sigmaIn, 0.10);
    params.getParam<double>("sigmaMin", sigmaMin, 0.05);
    params.getParam<double>("res", res, 1);
    params.getParam<double>("iseThresh", iseThresh, 0.2);
    params.getParam<double>("angleWin", angleWin, M_PI / 180.0 * 5.0);
    params.getParam<int>("nSamples", nSamples, 180);
    params.getParam<bool>("enableGmmHier", enableGmmHier, true);
    params.getParam<bool>("plotEllipsesArrows", plotEllipsesArrows, false);
    params.getParam<int>("numThreads", numThreads, 256);
    params.getParam<bool>("enablePclVisualization", enablePclVisualization, false);

    params.adaptTildeInPaths();
    params.getParam<std::string>("in", filenameIn, "sample.csv");

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n"
              << std::endl;

    BinReader binReader;
    binReader.setVehiclePtsMinMax(dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0));
    binReader.readCloudBin(filenameIn, BinReader::LidarType::VELODYNE);

    dsd::VectorVector3 musIn = binReader.getCloud();

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
        new pcl::PointCloud<pcl::PointXYZ>);
    dsd::binToPcl(musIn, cloud); // just for visualization purposes

    std::cout << "Read " << musIn.size() << " points" << std::endl;

    size_t n = musIn.size();
    for (size_t i = 0; i < n; ++i)
    {
        // ROFL_VAR2(i, musIn[i].transpose())
    }

    sigmasIn.resize(n);
    std::fill(sigmasIn.begin(), sigmasIn.end(), sigmaIn * sigmaIn);
    weightsIn.resize(n);
    std::fill(weightsIn.begin(), weightsIn.end(), 1.0 / n);

    ////////////////////////////////////////////////////////////////////////

    double execTimeI = 0.0;

    dsd::VectorVector3 musOut;
    dsd::VectorMatrix3 sigmasOut;
    std::vector<double> weightsOut;

    /* Calling GMM Hierarchical Estimation */
    if (enableGmmHier)
    {
        rofl::ScopedTimer gmeTimer("gme hier timer");

        gme::GaussianMixtureEstimatorHierarchical3d gme;
        gme.setSigmaMin(sigmaMin);
        // gme.setCovarWidth(sigmaMin); // not used
        gme.setIseThreshold(iseThresh);
        gme.setCellSizeMax(res);
        gme.compute(musIn);

        // gme.initIsotropic(musIn);
        gme.exportGaussians(musOut, sigmasOut, weightsOut);
        ROFL_VAR3(musOut.size(), sigmasOut.size(), weightsOut.size());
        execTimeI += gmeTimer.elapsedTimeMs();
        ROFL_VAR1(gmeTimer.elapsedTimeMs())
    }
    rofl::Profiler::getProfiler().printStats(std::cout);

    // std::cout << "Loading point cloud from \"" << filenameIn << "\""
    //           << std::endl;
    // if (pcl::io::loadPCDFile(filenameIn, *cloud) < 0)
    // {
    //     std::cerr << "Cannot load point cloud from \"" << filenameIn
    //               << "\"" << std::endl;
    //     return 1;
    // }

    size_t cloudSz = cloud->size();
    std::cout << "PCL cloud has " << cloudSz << " points" << std::endl;

    std::for_each(weightsOut.begin(), weightsOut.end(),
                  [cloudSz](double &n)
                  { n *= cloudSz; });

    //
    dsd::VectorVector3 musVomp;
    std::vector<double> kappasVomp;
    std::vector<double> weightsVomp;
    int szPadded = dsd::computeSzPadded(musOut.size(), nSamples, numThreads);
    ROFL_VAR1("von Mises stats")
    {
        rofl::ScopedTimer vonMisesStats("von Mises stats");

        dsd::vonMisesStats3dCuda(musVomp, kappasVomp, weightsVomp, szPadded, musOut, sigmasOut, weightsOut); // TODO: AUTOMATE szPadded

        execTimeI += vonMisesStats.elapsedTimeMs();
        ROFL_VAR1(vonMisesStats.elapsedTimeMs());
    }
    rofl::Profiler::getProfiler().printStats(std::cout);

    /********************START OF CUDA-RELATED PART************************/
    // CUDA MALLOC! -> vomp
    int mukwSz = musVomp.size(); //= szPadded
    ROFL_VAR1(mukwSz);           // !! MANDATORY PADDING!

    ROFL_ASSERT(mukwSz == szPadded && mukwSz == kappasVomp.size() && mukwSz == weightsVomp.size())

    int totalVompSz = musVomp.size() * 2 * nSamples * nSamples; // musOut.size() * 2
    double *vompDevice;
    cudaMalloc((void **)&vompDevice, totalVompSz * sizeof(double));
    ROFL_VAR1(totalVompSz); // !! MANDATORY PADDING!
    // cudaMemcpy(kernelInput, dataChunk.data(), (dataChunk.size()) * sizeof (cuars::Vec2d), cudaMemcpyHostToDevice);
    // cudaMemset()

    // return 0;

    // CUDA MALLOC! -> mus
    double3 *musDevice;
    cudaMalloc((void **)&musDevice, mukwSz * sizeof(double3));
    cudaMemcpy(musDevice, musVomp.data(), (musVomp.size()) * sizeof(double3), cudaMemcpyHostToDevice); // cudaMemcpy on double3 ??
    // CUDA MALLOC! -> k
    double *kDevice;
    cudaMalloc((void **)&kDevice, mukwSz * sizeof(double));
    cudaMemcpy(kDevice, kappasVomp.data(), (kappasVomp.size()) * sizeof(double), cudaMemcpyHostToDevice);

    // kernel input params:
    //  double *vomp, int nSamples, //!! mu size needs to be 3 times the size of k, w
    //  const double3 *mu, const double *k, const double *w, const int mukwSz
    int numBlocks = szPadded * 2 * nSamples * nSamples / numThreads; // TODO: automate numBlocks computation
    ROFL_VAR1(numBlocks);
    cudaEvent_t startKernelAll, stopKernelAll; // timing using CUDA events
    cudaEventCreate(&startKernelAll);
    cudaEventCreate(&stopKernelAll);
    cudaEventRecord(startKernelAll);
    vomp3d_kernel<<<numBlocks, numThreads>>>(vompDevice, nSamples, musDevice, kDevice, szPadded);
    cudaEventRecord(stopKernelAll);

    // std::vector<double> vompHost(totalVompSz, 0.0);
    // cudaMemcpy(vompHost.data(), vompDevice, (totalVompSz) * sizeof(double), cudaMemcpyDeviceToHost);

    cudaError_t cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
        printf("kernel launch failed with error \"%s\".\n", cudaGetErrorString(cudaerr));

    cudaEventSynchronize(stopKernelAll);
    float millisecondsKernelAll = 0.0f;
    cudaEventElapsedTime(&millisecondsKernelAll, startKernelAll, stopKernelAll);
    std::cout << "millisecondsKernelAll " << millisecondsKernelAll << " ms" << std::endl;

    cudaEventDestroy(startKernelAll);
    cudaEventDestroy(stopKernelAll);

    cudaFree(kDevice);
    cudaFree(musDevice);

    // maxVal = vompHost.empty() ? -1 : *std::max_element(vompHost.begin(), vompHost.end());
    // minVal = vompHost.empty() ? -1 : *std::min_element(vompHost.begin(), vompHost.end());
    // int maxId = vompHost.empty()
    //                 ? -1
    //                 : std::distance(vompHost.begin(),
    //                                 std::max_element(vompHost.begin(), vompHost.end()));

    // ROFL_VAR3(maxId, minVal, maxVal);

    // SUMMATION
    // CUDA MALLOC! -> w
    double *wDevice;
    cudaMalloc((void **)&wDevice, mukwSz * sizeof(double));
    // cuars::VecVec2d dataChunk(points.begin() + indicesStartEnd.first, points.begin() + (indicesStartEnd.first + currChunkSz));
    cudaMemcpy(wDevice, weightsVomp.data(), (weightsVomp.size()) * sizeof(double), cudaMemcpyHostToDevice);
    // CUDA MALLOC! -> vompSums
    double *vompSumsDevice;
    cudaMalloc((void **)&vompSumsDevice, nSamples * 2 * nSamples * sizeof(double));
    cudaMemset(vompSumsDevice, 0.0, nSamples * 2 * nSamples);

    cudaEvent_t startKernelDevice, stopKernelDevice; // timing using CUDA events
    cudaEventCreate(&startKernelDevice);
    cudaEventCreate(&stopKernelDevice);
    cudaEventRecord(startKernelDevice);
    cudaEventRecord(stopKernelDevice);
    vomp3d_summation_kernel<<<2 * nSamples * nSamples, 1>>>(vompSumsDevice, nSamples, //!! mu size needs to be 3 times the size of k, w
                                                            vompDevice, wDevice, mukwSz);

    std::vector<double> vomp(nSamples * 2 * nSamples, 0.0);
    cudaMemcpy(vomp.data(), vompSumsDevice, (nSamples * 2 * nSamples) * sizeof(double), cudaMemcpyDeviceToHost);

    cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
        printf("kernel launch failed with error \"%s\".\n", cudaGetErrorString(cudaerr));

    cudaEventSynchronize(stopKernelDevice);
    float millisecondsKernelDevice = 0.0f;
    cudaEventElapsedTime(&millisecondsKernelDevice, startKernelDevice, stopKernelDevice);
    std::cout << "millisecondsKernelDevice " << millisecondsKernelDevice << " ms" << std::endl;

    cudaEventDestroy(startKernelDevice);
    cudaEventDestroy(stopKernelDevice);

    // int i = 0;
    // for (auto &v : vomp)
    // {
    //     ROFL_VAR2(i, v)
    //     ++i;
    // }

    cudaFree(wDevice);
    cudaFree(vompDevice);
    /********************END OF CUDA-RELATED PART************************/

    // for (size_t i = 0; i < vomp.size(); ++i)
    // {
    // ROFL_VAR3(i, 360.0 * i / nSamples, vomp[i]);
    // }

    /**
     * Plot VMM distribution values
     */
    // dsd::plotVomp(vsmm, minVal, maxVal, viewer);

    double maxVal = vomp.empty() ? -1 : *std::max_element(vomp.begin(), vomp.end());
    double minVal = vomp.empty() ? -1 : *std::min_element(vomp.begin(), vomp.end());
    int maxId = vomp.empty()
                    ? -1
                    : std::distance(vomp.begin(),
                                    std::max_element(vomp.begin(), vomp.end()));
    ROFL_VAR3(maxId, minVal, maxVal)

    /**
     * Find peaks of VMM
     */
    ROFL_VAR1("von Mises max")
    std::vector<int> vompMaxima;
    std::vector<double> maximaValues;
    std::vector<std::pair<double, double>> thetaPhiMaxima;
    {
        rofl::ScopedTimer vonMisesMax("von Mises max");
        dsd::vompMax(vompMaxima, maximaValues, thetaPhiMaxima, vomp, nSamples, angleWin);
        execTimeI += vonMisesMax.elapsedTimeMs();
        ROFL_VAR1(vonMisesMax.elapsedTimeMs())
    }

    int maxIdx = 0;
    for (auto &tpm : thetaPhiMaxima)
    {
        // ROFL_VAR3(maximaValues[maxIdx], tpm.first, tpm.second)
        maxIdx++;
    }

    // plot von Mises mixture vector(matrix)
    double dtheta = M_PI / nSamples;

    int sizePhis = 2 * nSamples;
    int sizeThetas = nSamples;
    Eigen::MatrixXd vompEigenMat(sizePhis, sizeThetas);
    for (int jTheta = 0; jTheta < sizeThetas; ++jTheta)
    {
        // double theta = dtheta * jTheta;

        for (int jPhi = 0; jPhi < sizePhis; ++jPhi)
        {
            // int j = jTheta * sizePhis + jPhi;
            // double phi = dtheta * jPhi;
            vompEigenMat(jPhi, jTheta) = vomp[jTheta * sizePhis + jPhi];
        }
    }
    std::ofstream fileOutVomp("vomp.csv");
    fileOutVomp << vompEigenMat << std::endl;
    fileOutVomp.close();

    if (enablePclVisualization)
    {
        pcl::visualization::PCLVisualizer::Ptr viewer(
            new pcl::visualization::PCLVisualizer("3D Viewer"));
        viewer->setBackgroundColor(0.9, 0.9, 0.9);
        viewer->addCoordinateSystem(2.0);
        viewer->initCameraParameters();

        viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");

        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "cloud");
        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_COLOR, 0.05, 0.05, 0.05, "cloud");

        // pass to viewer

        // viewer->removeAllPointClouds();
        // viewer->removeAllShapes();

        std::cout << "Plotting result of gme" << std::endl;

        // Plot Ellipses for GMM covariance matrices and relative max eigenvector
        // direction on PCL visualizer
        if (plotEllipsesArrows)
            dsd::plotEllipsesArrows3d(viewer, musOut, sigmasOut, weightsOut);

        while (!viewer->wasStopped())
        {
            viewer->spinOnce(100);
            // ROFL_VAR1("viewer->spinOnce(100)")
        }
    }

    return 0;
}