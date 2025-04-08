#include <fstream>
#include <iostream>
#include <string>
#include <thread>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <rofl/common/macros.h>
#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

#include <bin_utils.h>
#include <dsd_utils.h>
#include <gme_GaussianMixtureEstimator.h>

using namespace std::chrono_literals;

int main(int argc, char **argv)
{
    std::string filenameCfg, fileIn;
    std::string lidarTypeStr;
    BinReader::LidarType lidarType;
    dsd::VectorVector3 musIn;
    std::vector<double> sigmasIn;
    std::vector<double> weightsIn;
    // GME parameters
    gme::GaussianMixtureEstimatorHierarchical3d gme3d;
    double sigmaIn, niseThr;
    dsd::VectorVector3 musOut;
    dsd::VectorMatrix3 covarsOut;
    std::vector<double> weightsOut;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", fileIn,
        std::string("/Geode/Urban_tunnel_01/LiDAR/bin/"));
    params.getParam<double>("sigmaIn", sigmaIn, 0.10);
    params.getParam<double>("niseThr", niseThr, 0.20);
    params.getParam<std::string>("lidar", lidarTypeStr, "Velodyne");

    std::cout << "-------\n"
              << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    if (lidarTypeStr == "Velodyne")
    {
        ROFL_MSG("LIDAR Type: Velodyne");
        lidarType = BinReader::LidarType::VELODYNE;
    }
    else if (lidarTypeStr == "Ouster")
    {
        ROFL_MSG("LIDAR Type: Ouster");
        lidarType = BinReader::LidarType::OUSTER;
    }
    else if (lidarTypeStr == "Livox")
    {
        ROFL_MSG("LIDAR Type: Livox");
        lidarType = BinReader::LidarType::LIVOX;
    }
    else
    {
        ROFL_MSG("LIDAR Type: INVALID!\n  selected default type Velodyne");
        lidarType = BinReader::LidarType::VELODYNE;
    }

    // Reads input file
    BinReader binReader;
    binReader.setVehiclePtsMinMax(dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0));
    binReader.readCloudBin(fileIn, lidarType);
    std::cout << "read file \"" << fileIn << "\"" << std::endl;

    musIn = binReader.getCloud();
    std::cout << "read " << musIn.size() << " input points" << std::endl;

    gme3d.setIseThreshold(niseThr);
    gme3d.setSigmaMin(sigmaIn);
    gme3d.setCellSizeMax(16.0 * sigmaIn);
    gme3d.compute(musIn);

    gme3d.exportGaussians(musOut, covarsOut, weightsOut);

    std::cout << "computed " << gme3d.gaussians().size() << "  Gaussian kernels"
              << std::endl;
    for (int k = 0; k < musOut.size() && k < covarsOut.size(); ++k)
    {
        ROFL_MSG("k " << k << ": mu " << musOut[k].transpose() << "\n"
                      << covarsOut[k]);
    }

    std::cout << "plotting point cloud and Gaussians" << std::endl;

    pcl::visualization::PCLVisualizer::Ptr viewer(
        new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0, 0, 0);
    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();

    // viewer->removeAllPointClouds();
    // viewer->removeAllShapes();

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
        new pcl::PointCloud<pcl::PointXYZ>);
    for (int i = 0; i < musIn.size(); ++i)
    {
        // copy scan into PCL cloud (2D)
        if (std::isfinite(musIn[i](0)) && std::isfinite(musIn[i](1)) &&
            std::isfinite(musIn[i](2)))
        {
            pcl::PointXYZ pI;
            pI.x = musIn[i](0);
            pI.y = musIn[i](1);
            pI.z = musIn[i](2);
            cloud->push_back(pI);
        }
    }

    viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");
    dsd::plotEllipsesArrows3d(viewer, musOut, covarsOut, weightsOut);

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
        std::this_thread::sleep_for(100ms);
    }

    return 0;
}