#include <boost/lexical_cast.hpp>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <thread>

#include <eigen3/Eigen/Dense>

#include <rofl/common/param_map.h>

#include "dsd_utils.h"

namespace fs = std::filesystem;

using namespace std::chrono_literals;

struct CsvScan
{
    std::string header;
    std::vector<Eigen::Vector2f> pts;
};

bool isValid(double r)
{
    return !std::isinf(r) && !std::isnan(r);
}

void deserializeRow(const std::string &row, Eigen::Vector2f &pt)
{
    std::stringstream ss(row);
    std::vector<std::string> ptCoordStrings;

    // ROFL_VAR1("\n");
    int coordId = 0;
    for (std::string strI; ss >> strI;)
    {
        ptCoordStrings.push_back(strI);

        if (ss.peek() == ',')
            ss.ignore();

        strI.erase(std::remove(strI.begin(), strI.end(), ','), strI.end());
        // ROFL_VAR1(strI);

        double ptCoord = std::stod(strI);

        // ROFL_VAR1(ptCoord);
        pt(coordId) = ptCoord;

        coordId++;
    }
    // ROFL_VAR1(pt.transpose());
}

void readCsvScan(std::string fname, CsvScan &csvsc)
{
    std::fstream fout;

    fout.open(fname, std::ios::in);

    std::string line;
    getline(fout, csvsc.header, '\n');
    ROFL_VAR1(csvsc.header);
    while (getline(fout, line, '\n'))
    {
        // ROFL_VAR1(line);

        // add all the column data
        // of a row to a vector

        Eigen::Vector2f pt;
        deserializeRow(line, pt);
        // ROFL_VAR1(line);

        csvsc.pts.push_back(pt);
    }

    fout.close();
}

struct point
{
    float x;
    float y;
    float phi;
};

struct Scan
{
    uint64_t id;

    // velocity
    float vx;
    float vy;
    float vOmega;

    // pose ground truth
    float vx_gt;
    float vy_gt;
    float vOmega_gt;

    // ranges and points
    std::vector<float> ranges;
    std::vector<point> points;
};

void saveCsvLog(const Scan &scan,
                const std::string &scanIdStr,
                const std::string &logName,
                bool firstTime)
{
    std::fstream fout;

    std::string dirName = "./csv/" + logName + "/";

    if (!(fs::exists(dirName)))
    {
        std::cout << "Doesn't Exist" << std::endl;

        if (fs::create_directory(dirName))
            std::cout << "....Successfully Created !" << std::endl;
    }
    else if (fs::exists(dirName) && firstTime)
    {
        fs::remove_all(dirName);
    }

    fout.open(dirName + logName + "_" + scanIdStr + ".csv",
              std::ios::out | std::ios::app);

    std::string header = "x, y";
    fout << header << std::endl;
    // Insert the data to file
    for (size_t i = 0; i < scan.points.size(); ++i)
    {
        fout << scan.points[i].x << ", " << scan.points[i].y << "\n";
    }
    fout << "\n";
    fout.close();
}
