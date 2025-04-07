#include <iostream>

#include <rofl/common/macros.h>
#include <rofl/common/param_map.h>

#include <gme_GaussianMixtureEstimator.h>
#include <dsd_utils.h>
#include <csv_utils.h>

#include <rofl/common/profiler.h>

#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

int main(int argc, char **argv)
{
    std::string filenameCfg, filePath;

    bool enableGmmHier, removeVehiclePts;

    int nSamples;

    dsd::VectorVector2 musIn;
    double sigmaMin, res, iseThresh;
    dsd::VectorVector2 musOut;
    dsd::VectorMatrix2 sigmasOut;
    std::vector<double> weightsOut;
    double angleWin;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", filePath,
        std::string("sample.csv"));

    // hier GME
    params.getParam<double>("sigmaMin", sigmaMin, 0.05);
    params.getParam<double>("res", res, 1);
    params.getParam<double>("iseThresh", iseThresh, 0.2);
    params.getParam<double>("angleWin", angleWin, M_PI / 180.0 * 5.0);

    params.getParam<int>("nSamples", nSamples, 360);

    params.getParam<bool>("enableGmmHier", enableGmmHier, true);

    params.getParam<bool>("removeVehiclePts", removeVehiclePts, false);

    // adapt tildes
    params.adaptTildeInPaths();
    params.getParam<std::string>("in", filePath, std::string(""));

    std::cout << "-------\n"
              << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    CsvScan csvScanMain;
    readCsvScan(filePath.c_str(), csvScanMain);

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
        new pcl::PointCloud<pcl::PointXYZ>);
    musIn.clear();

    for (int i = 0; i < csvScanMain.pts.size(); ++i)
    {
        // copy scan into PCL cloud (2D)
        if (std::isfinite(csvScanMain.pts[i](0)) &&
            std::isfinite(csvScanMain.pts[i](1)))
        {
            pcl::PointXYZ pI(pcl::PointXYZ(csvScanMain.pts[i](0),
                                           csvScanMain.pts[i](1), 0.0f));
            Eigen::Vector2d pIVehicleCoord; // pI in Vehicle coordinates;
            Eigen::Vector2d pIEig;      // pI in Eigen
            pIEig << pI.x, pI.y;
            pIVehicleCoord =
                pIEig; // initializing it here with no real purpose
                       // laserToVehicle(pIEig, pIVehicleCoord, laser2VehicleX, laser2VehicleY,
                       // laser2VehicleT); if (pI.getVector3fMap().norm() <
                       // fimRange && !dsd::isVehicleShapePt(pIVehicleCoord))
            // ROFL_VAR1(pIEig.transpose());

            if (dsd::isVehicleShapePt(pIEig) && removeVehiclePts)
                continue;

            ROFL_VAR3(i, csvScanMain.pts[i](0), csvScanMain.pts[i](1));
            cloud->push_back(pI);
            // std::cout << "i " << i << " pI " <<
            // pI.getVector3fMap().transpose() << " norm " <<
            // pI.getVector3fMap().norm() << std::endl;

            musIn.push_back(pIEig);
            // std::cout << "i " << i << " pIEig " << pIEig.transpose() <<
            // std::endl;
        }
    }
    ROFL_VAR1(musIn.size());

    /* Calling GMM Hierarchical Estimation */
    if (enableGmmHier)
    {
        rofl::ScopedTimer gmeTimer("gme hier timer");

        gme::GaussianMixtureEstimatorHierarchical gme;
        gme.setSigmaMin(sigmaMin);
        // gme.setCovarWidth(sigmaMin); // not used
        gme.setIseThreshold(iseThresh);
        gme.setCellSizeMax(res);
        gme.compute(musIn);

        gme.exportGaussians(musOut, sigmasOut, weightsOut);
        ROFL_VAR3(musOut.size(), sigmasOut.size(), weightsOut.size());
    };
    rofl::Profiler::getProfiler().printStats(std::cout);

    // pass to viewer

    pcl::visualization::PCLVisualizer::Ptr viewer(
        new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0.9, 0.9, 0.9);
    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();

    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");
    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");
    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_COLOR, 0.05, 0.05, 0.05,
        "cloud");

    bool plotGauss = true;
    if (plotGauss)
    {
        std::cout << "Plotting result of gme" << std::endl;

        // std::cout << "Loading point cloud from \"" << srcFilename << "\"" <<
        // std::endl; if (pcl::io::loadPCDFile(srcFilename, *cloud) < 0) {
        //     std::cerr << "Cannot load point cloud from \"" << srcFilename <<
        //     "\"" << std::endl; return 1;
        // }

        size_t cloudSz = cloud->size();
        std::cout << "Loaded cloud with " << cloudSz << " points" << std::endl;

        std::for_each(weightsOut.begin(), weightsOut.end(),
                      [cloudSz](double &n)
                      { n *= cloudSz; });

        // dsd::plotEllipsesArrows(viewer, musOut, sigmasOut, weightsOut);

        std::vector<double> vmm;
        double maxVal, minVal;
        // musOut.erase(musOut.begin() + 1, musOut.end());
        // sigmasOut.erase(sigmasOut.begin() + 1, sigmasOut.end());
        // weightsOut.erase(weightsOut.begin() + 1, weightsOut.end());
        {
            rofl::ScopedTimer("von Mises stats");
            dsd::vonMisesStats(vmm, nSamples, musOut, sigmasOut, weightsOut);
        }
        rofl::Profiler::getProfiler().printStats(std::cout);

        for (size_t i = 0; i < vmm.size(); ++i)
        {
            // ROFL_VAR3(i, 360.0 * i / nSamples, vmm[i]);
        }
        maxVal = *std::max_element(vmm.begin(), vmm.end());
        minVal = *std::min_element(vmm.begin(), vmm.end());
        int maxId = std::distance(vmm.begin(),
                                  std::max_element(vmm.begin(), vmm.end()));

        ROFL_VAR3(maxId, minVal, maxVal);

        dsd::plotEllipsesArrows(viewer, musOut, sigmasOut, weightsOut);

        dsd::plotVmm(vmm, minVal, maxVal, viewer);

        // plot von Mises mixture vector(matrix)

        Eigen::MatrixXd vmmEigenMat(nSamples, 1);

        for (int jAngle = 0; jAngle < nSamples; ++jAngle)
        {
            // int j = jTheta * sizePhis + jAngle;
            // double phi = dtheta * jAngle;
            vmmEigenMat(jAngle, 0) = vmm[jAngle];
        }

        std::ofstream fileOutVmm("vmm2d.csv");
        fileOutVmm << vmmEigenMat << std::endl;
        fileOutVmm.close();

        //////////

        std::vector<int> vmmMaxima;
        dsd::vonMisesMax(vmmMaxima, vmm, angleWin);
        for (int k = 0; k < vmmMaxima.size(); ++k)
        {
            ROFL_VAR3(k, vmmMaxima[k], vmm[vmmMaxima[k]]);
        }
    }

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
        std::this_thread::sleep_for(100ms);
    }

    return 0;
}
