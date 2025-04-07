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
#include <csv_utils.h>
#include <dsd_utils.h>
#include <gme_gaussian_metric.h>

int main(int argc, char** argv) {
    std::string filenameCfg, filePath;
    std::string folderInPath, filenameOutIse, filenameOutRtc, filenameOutTime;
    std::set<fs::path> sortedByName;

    dsd::VectorVector2 musIn;
    double sigmaIn;
    std::vector<double> sigmasIn;
    std::vector<double> weightsIn;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>("in", folderInPath,
                                 std::string("sample.csv"));
    params.getParam<std::string>("out_ise", filenameOutIse,
                                 std::string("ise.csv"));
    params.getParam<std::string>("out_rtc", filenameOutRtc,
                                 std::string("rtc.csv"));
    params.getParam<std::string>("out_time", filenameOutTime,
                                 std::string("time.txt"));
    params.getParam<double>("sigmaIn", sigmaIn, 0.05);

    std::cout << "-------\n" << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::ofstream fileIse(filenameOutIse);
    ROFL_ASSERT(fileIse);
    std::ofstream fileRtc(filenameOutRtc);
    ROFL_ASSERT(fileRtc);

    for (auto& entry : fs::directory_iterator(folderInPath))
        sortedByName.insert(entry.path());

    int counter = 0;
    for (const auto& entry : sortedByName) {
        // Reads the CSV file containing the point coordinates
        std::cout << "\n-----\n" << counter << " " << entry.stem() << std::endl;
        CsvScan csvScanMain;
        readCsvScan(entry.c_str(), csvScanMain);

        // Reads the scan, filters the invalid points or those belonging to the
        // Vehicle, fills the vector of points/mean values muIn
        musIn.clear();
        for (int i = 0; i < csvScanMain.pts.size(); ++i) {
            // copy scan into PCL cloud (2D)
            if (std::isfinite(csvScanMain.pts[i](0)) &&
                std::isfinite(csvScanMain.pts[i](1))) {
                Eigen::Vector2d pIVehicleCoord;  // pI in Vehicle coordinates;
                Eigen::Vector2d pIEig;       // pI in Eigen
                pIEig << csvScanMain.pts[i](0), csvScanMain.pts[i](1);
                pIVehicleCoord = pIEig;  // initializing it here with no real
                                     // purpose laserToVehicle(pIEig, pIVehicleCoord,
                                     // laser2VehicleX, laser2VehicleY, laser2VehicleT); if
                                     // (pI.getVector3fMap().norm() < fimRange
                                     // && !dsd::isVehicleShapePt(pIVehicleCoord))
                // ROFL_VAR1(pIEig.transpose());

                if (dsd::isVehicleShapePt(pIEig))
                    continue;

                musIn.push_back(pIEig);
            }
        }

        size_t n = musIn.size();
        // ROFL_VAR1(n);
        sigmasIn.resize(n);
        std::fill(sigmasIn.begin(), sigmasIn.end(), sigmaIn * sigmaIn);
        weightsIn.resize(n);
        std::fill(weightsIn.begin(), weightsIn.end(), 1.0 / n);

        dsd::Vector3 gradIse, gradRtc;
        dsd::Matrix3 hessianIse, hessianRtc;
        gme::computeGmmIseRtcGradHessian(musIn, sigmasIn, weightsIn, musIn,
                                         sigmasIn, weightsIn, gradIse,
                                         hessianIse, gradRtc, hessianRtc);

        Eigen::SelfAdjointEigenSolver<dsd::Matrix3> eigsolIse(hessianIse);
        double eigMinIse = fabs(eigsolIse.eigenvalues()(0));
        double eigMaxIse = fabs(eigsolIse.eigenvalues()(2));
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
                << eigsolIse.eigenvalues()(2) << ", " << (eigMaxIse / eigMinIse)
                << ", " << (eigMinIse / eigMaxIse) << "\n";

        Eigen::SelfAdjointEigenSolver<dsd::Matrix3> eigsolRtc(hessianRtc);
        double eigMinRtc = fabs(eigsolRtc.eigenvalues()(0));
        double eigMaxRtc = fabs(eigsolRtc.eigenvalues()(2));
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
                << eigsolRtc.eigenvalues()(2) << ", " << (eigMaxRtc / eigMinRtc)
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

    return 0;
}