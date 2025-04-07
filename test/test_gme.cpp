#include <iostream>
#include <fstream>
#include <filesystem>
#include <set>
#include <string>
#include <boost/lexical_cast.hpp>
#include <thread>

#include <rofl/common/param_map.h>

#include <csv_utils.h>

#include <gme_GaussianMixtureEstimator.h>

#include <dsd_utils.h>

#include <gme_gaussian_metric.h>

#include <rofl/common/profiler.h>

namespace fs = std::filesystem;

using namespace std::chrono_literals;

int main(int argc, char **argv)
{

    std::string filenameCfg, filePath, folderPath;

    /* Declare variables */
    dsd::VectorVector2 muIn;
    double sigmaMin, res, iseThresh;
    dsd::VectorVector2 muOut;
    dsd::VectorMatrix2 sigmaOut;
    std::vector<double> weightOut;

    double distGap, distSplit;

    bool removeVehiclePoints;

    bool enableGmmScan, enableGmmHier;

    int jump;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>("in", filePath,
                                 std::string("sample.csv"));
    params.getParam<std::string>("folder", folderPath,
                                 "");

    // hier GME
    params.getParam<double>("sigmaMin", sigmaMin, 0.05);
    params.getParam<double>("res", res, 1);
    params.getParam<double>("iseThresh", iseThresh, 0.2);

    // scan GME
    params.getParam<double>("distGap", distGap, 0.4);
    params.getParam<double>("distSplit", distSplit, 0.05);

    // methods enabling booleans
    params.getParam<bool>("enableGmmScan", enableGmmScan, true);
    params.getParam<bool>("enableGmmHier", enableGmmHier, false);

    params.getParam<bool>("removeVehiclePoints", removeVehiclePoints, false);

    params.getParam<int>("jump", jump, 80);

    // adapt tildes
    params.adaptTildeInPaths();
    params.getParam<std::string>("in", filePath, std::string(""));
    params.getParam<std::string>("folder", folderPath,
                                 "");

    std::cout << "-------\n"
              << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    // //--- filenames are unique so we can use a set
    // std::set<fs::path> sortedByName;

    pcl::visualization::PCLVisualizer::Ptr viewer(
        new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0, 0, 0);
    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();

    //--- filenames are unique so we can use a set
    std::set<fs::path> sortedByName;

    int counter = 0;

    for (auto &entry : fs::directory_iterator(folderPath))
        sortedByName.insert(entry.path());
    for (const auto &entry : sortedByName)
    {
        counter++;
        if (counter % jump != 0)
            continue;

        ROFL_VAR1(entry);

        CsvScan csvScanMain;
        readCsvScan(entry.c_str(), csvScanMain);

        // now copy scan into PCL point cloud

        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        muIn.clear();

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
                pIVehicleCoord = pIEig; // initializing it here with no real purpose
                                    // laserToVehicle(pIEig, pIVehicleCoord, laser2VehicleX, laser2VehicleY, laser2VehicleT);
                                    // if (pI.getVector3fMap().norm() < fimRange && !dsd::isVehicleShapePt(pIVehicleCoord))
                // ROFL_VAR1(pIEig.transpose());

                if (dsd::isVehicleShapePt(pIEig))
                    continue;

                cloud->push_back(pI);
                // std::cout << "i " << i << " pI " <<
                // pI.getVector3fMap().transpose() << " norm " <<
                // pI.getVector3fMap().norm() << std::endl;

                muIn.push_back(pIEig);
                // std::cout << "i " << i << " pIEig " << pIEig.transpose() << std::endl;
            }
        }

        /* Calling GMM Scan Estimation */
        if (enableGmmScan)
        {
            rofl::ScopedTimer gmeTimer("gme scan timer");

            gme::GaussianMixtureEstimatorScan gme;
            gme.setSigmaMin(sigmaMin);
            // gme.setCovarWidth(sigmaMin); // not used
            gme.setDistanceGap(distGap);
            gme.setDistanceSplit(distSplit);
            gme.compute(muIn);

            gme.exportGaussians(muOut, sigmaOut, weightOut);
            ROFL_VAR3(muOut.size(), sigmaOut.size(), weightOut.size());
        };
        rofl::Profiler::getProfiler().printStats(std::cout);

        /* Calling GMM Hierarchical Estimation */
        if (enableGmmHier)
        {
            rofl::ScopedTimer gmeTimer("gme hier timer");

            gme::GaussianMixtureEstimatorHierarchical gme;
            gme.setSigmaMin(sigmaMin);
            // gme.setCovarWidth(sigmaMin); // not used
            gme.setIseThreshold(iseThresh);
            gme.setCellSizeMax(res);
            gme.compute(muIn);

            gme.exportGaussians(muOut, sigmaOut, weightOut);
            ROFL_VAR3(muOut.size(), sigmaOut.size(), weightOut.size());
        };
        rofl::Profiler::getProfiler().printStats(std::cout);

        // pass to viewer
        viewer->removeAllPointClouds();
        viewer->removeAllShapes();

        viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");

        bool plotGauss = true;
        if (plotGauss)
        {
            std::cout << "Plotting result of gme" << std::endl;

            // std::cout << "Loading point cloud from \"" << srcFilename << "\"" << std::endl;
            // if (pcl::io::loadPCDFile(srcFilename, *cloud) < 0) {
            //     std::cerr << "Cannot load point cloud from \"" << srcFilename << "\"" << std::endl;
            //     return 1;
            // }

            size_t cloudSz = cloud->size();
            std::cout << "Loaded cloud with " << cloudSz << " points" << std::endl;

            std::for_each(weightOut.begin(), weightOut.end(), [cloudSz](double &n)
                          { n *= cloudSz; });

            dsd::plotEllipsesArrows(viewer, muOut, sigmaOut, weightOut);

            viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud");
            viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 1.0, "cloud");
        }

        std::chrono::milliseconds sleepclockMs = std::chrono::milliseconds(2500);
        std::this_thread::sleep_for(sleepclockMs);
        viewer->spinOnce(2500);
    }

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
        std::this_thread::sleep_for(100ms);
    }

    return 0;
}
