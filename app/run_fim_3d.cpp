#include <boost/lexical_cast.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <rofl/common/macros.h>
#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

// #include <gme_GaussianMixtureEstimator.h>
#include <bin_utils.h>
#include <dsd_utils.h>
#include <fim3d.h>

void convertEigenToPcl(const dsd::VectorVector3 &points,
                       pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
                       pcl::PointCloud<pcl::Normal>::Ptr cloudNormals);

int main(int argc, char **argv)
{
    std::string filenameCfg, filePath;
    std::string folderInPath;
    std::set<fs::path> sortedByName;
    std::string lidarTypeStr;
    BinReader::LidarType lidarType;
    std::string resultsBasePath;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudIn(
        new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::Normal>::Ptr normalsIn(
        new pcl::PointCloud<pcl::Normal>);

    double sigmaIn;
    std::vector<double> sigmasIn;
    std::vector<double> weightsIn;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", folderInPath,
        std::string("/Geode/Urban_tunnel_01/LiDAR/bin/"));
    params.getParam<double>("sigmaIn", sigmaIn, 0.10);
    params.getParam<std::string>("resultsBasePath", resultsBasePath,
                                 "../results_fim_3d/");
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

    if (!fs::exists(resultsBasePath))
        fs::create_directory(resultsBasePath);

    BinReader binReader;
    binReader.setVehiclePtsMinMax(dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0));

    for (auto &entry : fs::directory_iterator(folderInPath))
    {
        sortedByName.insert(entry.path());
    }

    std::string folderAppendNameStamped = dsd::generateStampedString("", "");
    std::string folderId =
        fs::path(folderInPath).parent_path().parent_path().filename().string();

    std::string filenameOutFim3d = resultsBasePath + folderAppendNameStamped +
                                   "_" + folderId + "_fim3d.out";

    std::string filenameOutTime = resultsBasePath + folderAppendNameStamped +
                                  "_" + folderId + "_fim3d_time.out";

    std::ofstream fileFim(filenameOutFim3d);
    ROFL_ASSERT(fileFim);

    // Read files from the input directory
    int counter = 0;
    dsd::Matrix6 fim3D;
    for (const auto &entry : sortedByName)
    {
        binReader.readCloudBin(entry.c_str(), lidarType);

        {
            rofl::ScopedTimer timer("FIM_3D_WITH NORMALS");
            convertEigenToPcl(binReader.getCloud(), cloudIn, normalsIn);
            computeFIM3D(fim3D, cloudIn, normalsIn, sigmaIn);
        }

        Eigen::SelfAdjointEigenSolver<dsd::Matrix6> eigsolFim(fim3D);
        double eigMinFim = fabs(eigsolFim.eigenvalues()(0));
        double eigMaxFim = fabs(eigsolFim.eigenvalues()(5));
        if (eigMinFim > eigMaxFim)
        {
            std::swap(eigMinFim, eigMaxFim);
        }
        // std::cout << "FIM\n"
        //           << fim3D << "\n"
        //           << "  eigenvalues " << eigsolFim.eigenvalues().transpose()
        //           << ", min " << eigMinFim << " max " << eigMaxFim << ", cond
        //           "
        //           << (eigMaxFim / eigMinFim) << ", eigMin/eigMax "
        //           << (eigMinFim / eigMaxFim) << std::endl;
        std::cout << counter << ", " << entry.stem() << ", "
                  << eigsolFim.eigenvalues()(0) << ", "
                  << eigsolFim.eigenvalues()(1) << ", "
                  << eigsolFim.eigenvalues()(2) << ", "
                  << eigsolFim.eigenvalues()(3) << ", "
                  << eigsolFim.eigenvalues()(4) << ", "
                  << eigsolFim.eigenvalues()(5) << ", "
                  << (eigMaxFim / eigMinFim) << ", " << (eigMinFim / eigMaxFim)
                  << "\n";
        // fileFim << "FIM\n"
        //         << fim3D << "\n"
        //         << "  eigenvalues " << eigsolFim.eigenvalues().transpose()
        //         << ", min " << eigMinFim << " max " << eigMaxFim << ", cond "
        //         << (eigMaxFim / eigMinFim) << ", eigMin/eigMax "
        //         << (eigMinFim / eigMaxFim) << std::endl;

        if (fileFim.is_open())
        {
            fileFim << counter << ", " << entry.stem() << ", "
                    << eigsolFim.eigenvalues()(0) << ", "
                    << eigsolFim.eigenvalues()(1) << ", "
                    << eigsolFim.eigenvalues()(2) << ", "
                    << eigsolFim.eigenvalues()(3) << ", "
                    << eigsolFim.eigenvalues()(4) << ", "
                    << eigsolFim.eigenvalues()(5) << ", "
                    << (eigMaxFim / eigMinFim) << ", "
                    << (eigMinFim / eigMaxFim) << "\n";
            fileFim.flush();
            // fileFim.close();
        }
        else
        {
            std::cerr << "fileFim NOT open" << std::endl;
            ROFL_ASSERT(0)
        }

        rofl::Profiler::getProfiler().printStats(std::cout);
        std::ofstream fileTime(filenameOutTime);
        if (fileTime.is_open())
        {
            rofl::Profiler::getProfiler().printStats(fileTime);
            fileTime.close();
        }

        counter++;
    }

    fileFim.close();

    return 0;
}

void convertEigenToPcl(const dsd::VectorVector3 &points,
                       pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
                       pcl::PointCloud<pcl::Normal>::Ptr cloudNormals)
{
    cloud->clear();
    cloudNormals->clear();
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    pcl::PointXYZ pI;

    cloud->points.reserve(points.size());
    for (int i = 0; i < points.size(); ++i)
    {
        // copy scan into PCL cloud (2D)
        if (std::isfinite(points[i](0)) && std::isfinite(points[i](1)) &&
            std::isfinite(points[i](2)))
        {
            pI.x = points[i](0);
            pI.y = points[i](1);
            pI.z = points[i](2);
            cloud->push_back(pI);
        }
    }

    // Computes the normals
    ne.setInputCloud(cloud);
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(
        new pcl::search::KdTree<pcl::PointXYZ>());
    ne.setSearchMethod(tree);
    ne.setRadiusSearch(0.03);
    ne.compute(*cloudNormals);
}