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

    params.adaptTildeInPaths();



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

    pcl::visualization::PCLVisualizer::Ptr viewer(
        new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0.9, 0.9, 0.9);
    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();

    viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");

    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "cloud");
    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_COLOR, 0.05, 0.05, 0.05, "cloud");

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

    // pass to viewer

    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");

    std::cout << "Plotting result of gme" << std::endl;

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

    // Plot Ellipses for GMM covariance matrices and relative max eigenvector
    // direction on PCL visualizer
    if (plotEllipsesArrows)
        dsd::plotEllipsesArrows3d(viewer, musOut, sigmasOut, weightsOut);

    //
    std::vector<double> vmm;
    double maxVal, minVal;
    // musOut.erase(musOut.begin() + 1, musOut.end());
    // sigmasOut.erase(sigmasOut.begin() + 1, sigmasOut.end());
    // weightsOut.erase(weightsOut.begin() + 1, weightsOut.end());
    int szPadded = dsd::computeSzPadded(musOut.size(), nSamples, numThreads);
    ROFL_VAR1("von Mises stats")
    {
        rofl::ScopedTimer vonMisesStats("von Mises stats");
        dsd::vonMisesStats3d(vmm, nSamples, musOut, sigmasOut, weightsOut, szPadded);
        execTimeI += vonMisesStats.elapsedTimeMs();
        ROFL_VAR1(vonMisesStats.elapsedTimeMs());
    }
    rofl::Profiler::getProfiler().printStats(std::cout);

    //plot von Mises mixture vector(matrix)
    double dtheta = M_PI / nSamples;

    int sizePhis = 2 * nSamples;
    int sizeThetas = nSamples;
    Eigen::MatrixXd vmmEigenMat(sizePhis, sizeThetas);
    for (int jTheta = 0; jTheta < sizeThetas; ++jTheta)
    {
        // double theta = dtheta * jTheta;

        for (int jPhi = 0; jPhi < sizePhis; ++jPhi)
        {
            // int j = jTheta * sizePhis + jPhi;
            // double phi = dtheta * jPhi;
            vmmEigenMat(jPhi, jTheta) = vmm[jTheta * sizePhis + jPhi];
        }
    }
    std::ofstream fileOutVmm("vmm.csv");
    fileOutVmm << vmmEigenMat << std::endl;
    fileOutVmm.close();

    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");
    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_COLOR, 0.05, 0.05, 0.05, "cloud");

    // for (size_t i = 0; i < vmm.size(); ++i)
    // {
    // ROFL_VAR3(i, 360.0 * i / nSamples, vmm[i]);
    // }
    maxVal = vmm.empty() ? -1 : *std::max_element(vmm.begin(), vmm.end());
    minVal = vmm.empty() ? -1 : *std::min_element(vmm.begin(), vmm.end());
    int maxId = vmm.empty()
                    ? -1
                    : std::distance(vmm.begin(),
                                    std::max_element(vmm.begin(), vmm.end()));

    ROFL_VAR3(maxId, minVal, maxVal);

    int i = 0;
    for (auto &v : vmm)
    {
        ROFL_VAR2(i, v)
        ++i;
    }

    /**
     * Plot VMM distribution values
     */
    // dsd::plotVmm(vmm, minVal, maxVal, viewer);

    /**
     * Find peaks of VMM
     */
    ROFL_VAR1("von Mises max")
    std::vector<int> vmmMaxima;
    std::vector<std::pair<double, double>> thetaPhiMaxima;
    std::vector<double> maximaValues;
    {
        rofl::ScopedTimer vonMisesMax("von Mises max");
        dsd::fvmMax(vmmMaxima, maximaValues, thetaPhiMaxima, vmm, nSamples, angleWin);
        execTimeI += vonMisesMax.elapsedTimeMs();
        ROFL_VAR1(vonMisesMax.elapsedTimeMs());
    }

    int maxIdx = 0;
    for (auto &tpm : thetaPhiMaxima)
    {
        ROFL_VAR3(maximaValues[maxIdx], tpm.first, tpm.second)
        maxIdx++;
    }

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
    }

    return 0;
}