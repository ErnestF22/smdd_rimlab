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
#include <gme_gaussian_metric.h>

int main(int argc, char** argv) {
    std::string filenameCfg, filePath;
    std::string folderInPath, filenameOutIse, filenameOutRtc, filenameOutTime;
    std::set<fs::path> sortedByName;
    std::string lidarTypeStr;
    BinReader::LidarType lidarType;

    double sigmaIn;
    std::vector<double> sigmasIn;
    std::vector<double> weightsIn;
    std::string resultsBasePath;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", folderInPath,
        std::string("/Geode/Urban_tunnel_01/LiDAR/bin/"));
    params.getParam<std::string>("out_ise", filenameOutIse,
                                 std::string("ise.csv"));
    params.getParam<std::string>("out_rtc", filenameOutRtc,
                                 std::string("rtc.csv"));
    params.getParam<std::string>("out_time", filenameOutTime,
                                 std::string("times.txt"));
    params.getParam<double>("sigmaIn", sigmaIn, 0.10);
    params.getParam<std::string>("resultsBasePath", resultsBasePath,
                                 "../results_ise_rtc_3d/");
    params.getParam<std::string>("lidar", lidarTypeStr, "Velodyne");

    std::cout << "-------\n" << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    if (lidarTypeStr == "Velodyne") {
        ROFL_MSG("LIDAR Type: Velodyne");
        lidarType = BinReader::LidarType::VELODYNE;
    } else if (lidarTypeStr == "Ouster") {
        ROFL_MSG("LIDAR Type: Ouster");
        lidarType = BinReader::LidarType::OUSTER;
    } else if (lidarTypeStr == "Livox") {
        ROFL_MSG("LIDAR Type: Livox");
        lidarType = BinReader::LidarType::LIVOX;
    } else {
        ROFL_MSG("LIDAR Type: INVALID!\n  selected default type Velodyne");
        lidarType = BinReader::LidarType::VELODYNE;
    }

    if (!fs::exists(resultsBasePath))
        fs::create_directory(resultsBasePath);

    BinReader binReader;
    binReader.setVehiclePtsMinMax(dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0));

    // Reads file from the input directory
    for (auto& entry : fs::directory_iterator(folderInPath))
        sortedByName.insert(entry.path());

    std::string folderAppendNameStamped = dsd::generateStampedString("", "");
    std::string folderId =
        fs::path(folderInPath).parent_path().parent_path().filename().string();

    filenameOutIse = resultsBasePath + folderAppendNameStamped + "_" +
                     folderId + "_ise3d.out";
    filenameOutRtc = resultsBasePath + folderAppendNameStamped + "_" +
                     folderId + "_rtc3d.out";
    filenameOutTime = resultsBasePath + folderAppendNameStamped + "_" +
                      folderId + "_time3d.out";
    std::cout << "filenameOutIse: " << std::endl;
    std::cout << "filenameOutIse: " << std::endl;

    std::ofstream fileIse(filenameOutIse);
    ROFL_ASSERT(fileIse);
    std::ofstream fileRtc(filenameOutRtc);
    ROFL_ASSERT(fileRtc);

    int counter = 0;
    for (const auto& entry : sortedByName) {
        binReader.readCloudBin(entry.c_str(), lidarType);

        dsd::VectorVector3 musIn = binReader.getCloud();

        // pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
        //     new pcl::PointCloud<pcl::PointXYZ>);
        // dsd::binToPcl(musIn, cloud);

        // pcl::visualization::PCLVisualizer::Ptr viewer(
        //     new pcl::visualization::PCLVisualizer("3D Viewer"));
        // viewer->setBackgroundColor(0.9, 0.9, 0.9);
        // viewer->addCoordinateSystem(1.0);
        // viewer->initCameraParameters();

        // viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");

        // viewer->setPointCloudRenderingProperties(
        //     pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "cloud");
        // viewer->setPointCloudRenderingProperties(
        //     pcl::visualization::PCL_VISUALIZER_COLOR, 0.05, 0.05, 0.05,
        //     "cloud");

        std::cout << "Read " << musIn.size() << " points" << std::endl;

        size_t n = musIn.size();
        // for (size_t i = 0; i < n; ++i)
        // {
        //     std::cout << "  " << i << ": " << musIn[i].transpose() <<
        //     std::endl;
        // }
        // std::cout << "Read " << n << " points" << std::endl;

        sigmasIn.resize(n);
        std::fill(sigmasIn.begin(), sigmasIn.end(), sigmaIn * sigmaIn);
        weightsIn.resize(n);
        std::fill(weightsIn.begin(), weightsIn.end(), 1.0 / n);

        dsd::Vector6 gradIse, gradRtc;
        dsd::Matrix6 hessianIse, hessianRtc;
        ROFL_VAR2(gradIse.transpose(), gradRtc.transpose());
        gme::computeGmmIseRtcGradHessian3D(musIn, sigmasIn, weightsIn, musIn,
                                           sigmasIn, weightsIn, gradIse,
                                           hessianIse, gradRtc, hessianRtc);

        Eigen::SelfAdjointEigenSolver<dsd::Matrix6> eigsolIse(hessianIse);
        double eigMinIse = fabs(eigsolIse.eigenvalues()(0));
        double eigMaxIse = fabs(eigsolIse.eigenvalues()(5));
        if (eigMaxIse < eigMinIse) {
            std::swap(eigMinIse, eigMaxIse);
        }
        std::cout << "ISE gradient: " << gradIse.transpose() << "\n"
                  << "ISE Hessian\n"
                  << hessianIse << "\n"
                  << "  eigenvalues " << eigsolIse.eigenvalues().transpose()
                  << ", min " << eigMinIse << " max " << eigMaxIse << ", cond "
                  << (eigMaxIse / eigMinIse) << ", eigMin/eigMax "
                  << (eigMinIse / eigMaxIse) << std::endl;
        fileIse << counter << ", " << entry.stem() << ", "
                << eigsolIse.eigenvalues()(0) << ", "
                << eigsolIse.eigenvalues()(1) << ", "
                << eigsolIse.eigenvalues()(2) << ", "
                << eigsolIse.eigenvalues()(3) << ", "
                << eigsolIse.eigenvalues()(4) << ", "
                << eigsolIse.eigenvalues()(5) << ", " << (eigMaxIse / eigMinIse)
                << ", " << (eigMinIse / eigMaxIse) << "\n";

        Eigen::SelfAdjointEigenSolver<dsd::Matrix6> eigsolRtc(hessianRtc);
        double eigMinRtc = fabs(eigsolRtc.eigenvalues()(0));
        double eigMaxRtc = fabs(eigsolRtc.eigenvalues()(5));
        if (eigMaxRtc < eigMinRtc) {
            std::swap(eigMinRtc, eigMaxRtc);
        }
        std::cout << "RTC gradient: " << gradRtc.transpose() << "\n"
                  << "RTC Hessian\n"
                  << hessianRtc << "\n"
                  << "  eigenvalues " << eigsolRtc.eigenvalues().transpose()
                  << ", min " << eigMinRtc << " max " << eigMaxRtc << ", cond "
                  << (eigMaxRtc / eigMinRtc) << ", eigMin/eigMax "
                  << (eigMinRtc / eigMaxRtc) << std::endl;
        fileRtc << counter << ", " << entry.stem() << ", "
                << eigsolRtc.eigenvalues()(0) << ", "
                << eigsolRtc.eigenvalues()(1) << ", "
                << eigsolRtc.eigenvalues()(2) << ", "
                << eigsolRtc.eigenvalues()(3) << ", "
                << eigsolRtc.eigenvalues()(4) << ", "
                << eigsolRtc.eigenvalues()(5) << ", " << (eigMaxRtc / eigMinRtc)
                << ", " << (eigMinRtc / eigMaxRtc) << "\n";

        rofl::Profiler::getProfiler().printStats(std::cout);
        std::ofstream fileTime(filenameOutTime);
        if (fileTime.is_open()) {
            rofl::Profiler::getProfiler().printStats(fileTime);
            fileTime.close();
        }

        counter++;
    }

    fileIse.close();
    fileRtc.close();
    // while (!viewer->wasStopped()) {
    //     viewer->spinOnce(100);
    // }

    return 0;
}