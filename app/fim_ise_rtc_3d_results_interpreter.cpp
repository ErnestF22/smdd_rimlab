#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>

#include <bin_utils.h>
#include <dsd_utils.h>

#include <pcl/visualization/pcl_visualizer.h>

#include <rofl/common/param_map.h>

namespace fs = std::filesystem;

struct Line {
    // counter
    int counter;
    // entry.stem()
    std::string entryName;
    // eigsolIse.eigenvalues()(0)
    // eigsolIse.eigenvalues()(1)
    // eigsolIse.eigenvalues()(2)
    // eigsolIse.eigenvalues()(3)
    // eigsolIse.eigenvalues()(4)
    // eigsolIse.eigenvalues()(5)
    std::vector<double> eigenvalues;

    // (eigMaxIse / eigMinIse)
    double conditionNumber;
    // (eigMinIse / eigMaxIse)
    double invConditionNumber;

    Line()
        : counter(-1),
          entryName(""),
          eigenvalues(6, 0.0),
          conditionNumber(-1),
          invConditionNumber(-1) {}

    Line(int counter,
         std::string entryName,
         const std::vector<double>& eigenvalues,
         double conditionNumber,
         double invConditionNumber)
        : counter(counter),
          entryName(entryName),
          eigenvalues(eigenvalues),
          conditionNumber(conditionNumber),
          invConditionNumber(invConditionNumber) {
        ROFL_ASSERT(eigenvalues.size() == 6);
    }

    ~Line() { eigenvalues.clear(); }

    friend std::ostream& operator<<(std::ostream& os, const Line& m) {
        os << "counter: " << m.counter << ", entryName: " << m.entryName
           << std::endl
           << "eigenvalues: " << m.eigenvalues[0] << " " << m.eigenvalues[1]
           << " " << m.eigenvalues[2] << " " << m.eigenvalues[3] << " "
           << m.eigenvalues[4] << " " << m.eigenvalues[5] << std::endl
           << "cond: " << m.conditionNumber
           << " invcond: " << m.invConditionNumber << std::endl;
        return os;
    }
};

bool isDegenerate(const Line& l, double thr, bool lesserTrue);

int main(int argc, char** argv) {
    rofl::ParamMap params;

    std::string filenameIn;
    std::string filenameCfg;
    std::string gnuplotBasePath;
    double thr;
    bool lesserTrue;
    int numCloudToVisualize;
    std::string folderForVisualization;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    params.getParam<std::string>(
        "in", filenameIn,
        "file.out");
    params.getParam<std::string>("gnuplotBasePath", gnuplotBasePath,
                                 "../gnuplot/");
    // params.getParam<std::string>("method", method, "fim3d");
    params.getParam<double>("thr", thr, -1.0);
    params.getParam<bool>("lesserTrue", lesserTrue, true);
    params.getParam<int>("numCloudToVisualize", numCloudToVisualize, 0);
    params.getParam<std::string>(
        "folderForVisualization", folderForVisualization,
        "/Geode/sensor_data/Urban_tunnel/"
        "Urban_Tunnel01/LiDAR/bin/");

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n" << std::endl;

    // sort by name input files
    std::set<fs::path> sortedByName;
    // setup output
    if (!fs::exists(gnuplotBasePath))
        fs::create_directory(gnuplotBasePath);
    std::string folderAppendNameStamped = dsd::generateStampedString("", "");
    std::string folderMethodId =
        (fs::path(filenameIn).filename().string())
            .substr(folderAppendNameStamped.size() + 1);
    ROFL_VAR1(fs::path(folderMethodId).stem().string())
    std::string filenameOut = gnuplotBasePath + folderAppendNameStamped + "_" +
                              fs::path(folderMethodId).stem().string() +
                              +".gnuplotdata";

    std::string method =
        fs::path(folderMethodId)
            .stem()
            .string()
            .substr(fs::path(folderMethodId).stem().string().size() -
                    5);  // fim3d, rtc3d or ise3d
    ROFL_VAR1(method)

    std::ifstream fileIn(filenameIn);

    std::ofstream fileOut(filenameOut);

    if (!fileIn.is_open()) {
        std::cerr << "Error: Could not open fileIn " << filenameIn << std::endl;
        return 1;
    }

    if (!fileOut.is_open()) {
        std::cerr << "Error: Could not open fileOut " << filenameOut
                  << std::endl;
        return 1;
    }

    ROFL_VAR1(filenameOut)

    std::string line;
    int lineCtr = 0;
    std::vector<Line> lines;
    while (std::getline(fileIn, line)) {
        Line l;
        std::string s;
        std::istringstream iss(line);
        int inlineCtr = 0;
        while (std::getline(iss, s, ',')) {
            // counter
            if (inlineCtr == 0) {
                l.counter = std::stoi(s);
                ROFL_ASSERT(lineCtr == l.counter)
            } else if (inlineCtr == 1) {
                // entry.stem()
                std::string entryName = s;
                // ROFL_VAR1(entryName)
                entryName.erase(
                    std::remove_copy(entryName.begin(), entryName.end(),
                                     entryName.begin(), '\"'),
                    entryName.end());
                // ROFL_VAR1(entryName)

                l.entryName = entryName;
            } else if (inlineCtr == 2) {
                // eigsolIse.eigenvalues()(0)
                l.eigenvalues[0] = std::stod(s);
            } else if (inlineCtr == 3) {
                // eigsolIse.eigenvalues()(1)

                l.eigenvalues[1] = std::stod(s);
            } else if (inlineCtr == 4) {
                // eigsolIse.eigenvalues()(2)

                l.eigenvalues[2] = std::stod(s);
            } else if (inlineCtr == 5) {
                // eigsolIse.eigenvalues()(3)

                l.eigenvalues[3] = std::stod(s);
            } else if (inlineCtr == 6) {
                // eigsolIse.eigenvalues()(4)

                l.eigenvalues[4] = std::stod(s);
            } else if (inlineCtr == 7) {
                // eigsolIse.eigenvalues()(5)

                l.eigenvalues[5] = std::stod(s);
            } else if (inlineCtr == 8) {
                // (eigMaxIse / eigMinIse)

                l.conditionNumber = std::stod(s);
            } else if (inlineCtr == 9) {
                // (eigMinIse / eigMaxIse)

                l.invConditionNumber = std::stod(s);
            } else {
                ROFL_VAR1("Bad inlineCtr")
            }
            inlineCtr++;
        }
        // ROFL_VAR2(lineCtr, l)
        lines.push_back(l);
        lineCtr++;
    }

    for (int i = 0; i < lines.size(); ++i) {
        fileOut << isDegenerate(lines[i], thr, lesserTrue) << " "
                << lines[i].conditionNumber << std::endl;
    }

    // Visualizing desired cloud
    auto cloudName = lines[numCloudToVisualize].entryName;

    BinReader binReader;
    binReader.setVehiclePtsMinMax(dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0));
    std::string cloudFilename = folderForVisualization + cloudName + ".bin";
    std::string::iterator end_pos =
        std::remove(cloudFilename.begin(), cloudFilename.end(), ' ');
    cloudFilename.erase(end_pos, cloudFilename.end());
    ROFL_VAR1(cloudFilename)
    binReader.readCloudBin(cloudFilename, BinReader::LidarType::VELODYNE);

    dsd::VectorVector3 musIn = binReader.getCloud();

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
        new pcl::PointCloud<pcl::PointXYZ>);
    dsd::binToPcl(musIn, cloud);  // just for visualization purposes

    // viewer->removeAllPointClouds();
    // viewer->removeAllShapes();

    // pcl::visualization::PCLVisualizer::Ptr viewer(
    //     new pcl::visualization::PCLVisualizer("3D Viewer"));
    // viewer->setBackgroundColor(0.9, 0.9, 0.9);
    // viewer->addCoordinateSystem(1.0);
    // viewer->initCameraParameters();

    // viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");

    // viewer->setPointCloudRenderingProperties(
    //     pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");
    // viewer->setPointCloudRenderingProperties(
    //     pcl::visualization::PCL_VISUALIZER_COLOR, 0.05, 0.05, 0.05, "cloud");

    fileIn.close();
    fileOut.close();

    // while (!viewer->wasStopped()) {
    //     viewer->spinOnce(100);
    // }

    return 0;
}

bool isDegenerate(const Line& l, double thr, bool lesserTrue) {
    if (l.conditionNumber < thr) {
        if (lesserTrue)
            return true;
        else
            return false;
    } else {
        if (lesserTrue)
            return false;
        else
            return true;
    }
}