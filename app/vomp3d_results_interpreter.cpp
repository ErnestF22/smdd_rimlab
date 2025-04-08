#include <iostream>
#include <fstream>
#include <filesystem>
#include <set>

#include <dsd_utils.h>

#include <rofl/common/param_map.h>

struct Maximum
{
    double theta;
    double phi;
    double val;
    // double exectime;

    Maximum() : theta(0.0), phi(0.0), val(-1.0) {}

    Maximum(double theta, double phi, double val) : theta(theta), phi(phi), val(val) {}

    ~Maximum() {}

    friend std::ostream &operator<<(std::ostream &os, const Maximum &m)
    {
        os << "theta: " << m.theta << ", phi: " << m.phi << ", val: " << m.val;
        return os;
    }
};

std::string firstIntInString(std::string const &str);

bool isDegenerate(const std::vector<Maximum> &maxima);

namespace fs = std::filesystem;

int main(int argc, char **argv)
{
    rofl::ParamMap params;

    std::string folderInPath;
    std::string filenameCfg;
    std::string gnuplotBasePath;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // params.getParam<std::string>("in", folderInPath, "/home/rimlab/rclone/dsd_paper_results/results_vomp3d/20250214_1511_55_Urban_Tunnel01/1693022008.716670513.out");
    params.getParam<std::string>("in", folderInPath, "../results_vomp3d/20250214_1511_55_Urban_Tunnel01/");
    params.getParam<std::string>("gnuplotBasePath", gnuplotBasePath, "../gnuplot/");

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n"
              << std::endl;

    // sort by name input files
    std::set<fs::path> sortedByName;
    for (auto &entry : fs::directory_iterator(folderInPath))
        sortedByName.insert(entry.path());

    // setup output

    std::string folderAppendNameStamped = dsd::generateStampedString("", "");
    std::string folderId = (fs::path(folderInPath).parent_path().filename().string()).substr(folderAppendNameStamped.size() + 1);
    ROFL_VAR1(fs::path(folderId).stem().string())

    std::string tmp = fs::path(folderInPath).parent_path().parent_path().filename().string();
    std::string method = tmp.substr(tmp.size() - 5); // fim3d, rtc3d or ise3d
    ROFL_VAR1(method)

    if (!fs::exists(gnuplotBasePath + "/" + method + "_gnuplotdata/"))
        fs::create_directories(gnuplotBasePath + "/" + method + "_gnuplotdata/");

    std::string filenameOut = gnuplotBasePath + "/" + method + "_gnuplotdata/" + folderAppendNameStamped + "_" +
                              fs::path(folderId).stem().string() + "_" + method + ".gnuplotdata";

    ROFL_VAR1(filenameOut)
    // return 0;

    std::ofstream fileOut(filenameOut);

    std::ofstream meanExecTimeFileOut(gnuplotBasePath + "/" + method + "_gnuplotdata/" + folderAppendNameStamped + "_" +
                                      fs::path(folderId).stem().string() + "_" + method + "_meanExecTimes.gnuplotdata");

    if (!fileOut.is_open())
    {
        std::cerr << "Error: Could not open fileOut " << filenameOut << std::endl;
        return 1;
    }

    double totalExecTime = 0.0;
    int execTimeCtr = 0;

    for (const auto &entry : sortedByName)
    {
        std::cout << entry << std::endl;

        if (entry.extension().string() != ".out")
            continue;

        std::ifstream fileIn(entry);

        if (!fileIn.is_open())
        {
            std::cerr << "Error: Could not open fileIn " << entry << std::endl;
            return 1;
        }

        std::vector<Maximum> maxima;

        std::string line;
        int lineCtr = 0;
        int numMaxima;
        std::getline(fileIn, line);

        // ROFL_VAR2(lineCtr, line);
        numMaxima = stoi(firstIntInString(line));
        ROFL_VAR1(numMaxima);

        lineCtr++;

        double exectime;
        for (int m = 0; m < numMaxima; ++m)
        {
            double theta, phi, val;
            // maximum # 0
            std::getline(fileIn, line);
            ROFL_ASSERT(stoi(firstIntInString(line)) == m)
            // theta 0.0698132
            std::getline(fileIn, line);
            {
                std::string s;
                std::istringstream iss(line);
                while (std::getline(iss, s, ' '))
                {
                    // printf("`%s'\n", s.c_str());
                    if (s != "theta")
                    {
                        theta = stod(s);
                    }
                }
            }
            // phi 4.36332
            std::getline(fileIn, line);
            {
                std::string s;
                std::istringstream iss(line);
                while (std::getline(iss, s, ' '))
                {
                    // printf("`%s'\n", s.c_str());
                    if (s != "phi")
                    {
                        phi = stod(s);
                    }
                }
            }
            // val 0.0159072
            std::getline(fileIn, line);
            {
                std::string s;
                std::istringstream iss(line);
                while (std::getline(iss, s, ' '))
                {
                    // printf("`%s'\n", s.c_str());
                    if (s != "val")
                    {
                        val = stod(s);
                    }
                }
            }
            // 142.906 [ms]
            std::getline(fileIn, line);
            {
                std::string s;
                std::istringstream iss(line);
                while (std::getline(iss, s, ' '))
                {
                    // printf("`%s'\n", s.c_str());
                    if (s != "[ms]")
                    {
                        exectime = stod(s); // OBS. Written by mistake in output once per peak...
                    }
                }
            }
            // ROFL_VAR5(m, theta, phi, val, exectime);
            Maximum maximum(theta, phi, val);
            maxima.push_back(maximum);
            ROFL_VAR1(maximum)
            ROFL_VAR2(entry, m)
            lineCtr++;
        }

        totalExecTime += exectime;
        execTimeCtr++;

        bool isDegSeqI = isDegenerate(maxima);
        fileOut << isDegSeqI << std::endl;

        fileIn.close();
    }

    if (execTimeCtr > 0)
    {
        double meanExecTime = totalExecTime / execTimeCtr;
        meanExecTimeFileOut << meanExecTime << std::endl;
    }
    else
    {
        meanExecTimeFileOut << "no cases to be evaluated" << std::endl;
    }

    meanExecTimeFileOut.close();

    fileOut.close();

    return 0;
}

std::string firstIntInString(std::string const &str)
{
    char const *digits = "0123456789";
    std::size_t const n = str.find_first_of(digits);
    if (n != std::string::npos)
    {
        std::size_t const m = str.find_first_not_of(digits, n);
        return str.substr(n, m != std::string::npos ? m - n : m);
    }
    return std::string();
}

bool isDegenerate(const std::vector<Maximum> &maxima)
{

    if (maxima.size() < 4 || maxima.size() > 16)
    {
        return true;
    }

    if (maxima.size() >= 4 && maxima.size() < 9)
    {
        std::vector<double> ds;
        for (int i = 0; i < maxima.size(); ++i)
        {
            ds.push_back(maxima[i].val);
        }
        std::sort(ds.begin(), ds.end(), std::greater<double>());
        if (ds[0] > 10.0 * ds[4])
        {
            return true;
        }
    }
    return false;
}
