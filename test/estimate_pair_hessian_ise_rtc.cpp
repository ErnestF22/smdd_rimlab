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

int main(int argc, char **argv)
{
    std::string filenameCfg, filePath;
    std::string folderInPath, filenameOutIse, filenameOutRtc;
    std::set<fs::path> sortedByName;

    dsd::VectorVector2 musDst, musSrc;
    double sigmaIn;
    std::vector<double> sigmasDst, sigmasSrc;
    std::vector<double> weightsDst, weightsSrc;

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
    params.getParam<double>("sigmaIn", sigmaIn, 0.05);

    std::cout << "-------\n"
              << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::ofstream fileIse(filenameOutIse);
    ROFL_ASSERT(fileIse);
    std::ofstream fileRtc(filenameOutRtc);
    ROFL_ASSERT(fileRtc);

    for (auto &entry : fs::directory_iterator(folderInPath))
        sortedByName.insert(entry.path());

    int counter = 0;
    musSrc.clear();
    sigmasSrc.clear();
    weightsSrc.clear();
    for (const auto &entry : sortedByName)
    {
        // Reads the CSV file containing the point coordinates
        std::cout << "\n-----\n"
                  << counter << " " << entry.stem() << std::endl;
        CsvScan csvScanMain;
        readCsvScan(entry.c_str(), csvScanMain);

        // Reads the scan, filters the invalid points or those belonging to the
        // Vehicle, fills the vector of points/mean values muIn
        musDst.clear();
        for (int i = 0; i < csvScanMain.pts.size(); ++i)
        {
            // copy scan into PCL cloud (2D)
            if (std::isfinite(csvScanMain.pts[i](0)) &&
                std::isfinite(csvScanMain.pts[i](1)))
            {
                Eigen::Vector2d pIVehicleCoord; // pI in Vehicle coordinates;
                Eigen::Vector2d pIEig;          // pI in Eigen
                pIEig << csvScanMain.pts[i](0), csvScanMain.pts[i](1);
                pIVehicleCoord = pIEig; // initializing it here with no real
                                        // purpose laserToVehicle(pIEig, pIVehicleCoord,
                                        // laser2VehicleX, laser2VehicleY, laser2VehicleT); if
                                        // (pI.getVector3fMap().norm() < fimRange
                                        // && !dsd::isVehicleShapePt(pIVehicleCoord))
                // ROFL_VAR1(pIEig.transpose());

                if (dsd::isVehicleShapePt(pIEig))
                    continue;

                musDst.push_back(pIEig);
            }
        }

        size_t n = musDst.size();
        // ROFL_VAR1(n);
        sigmasDst.resize(n);
        std::fill(sigmasDst.begin(), sigmasDst.end(), sigmaIn * sigmaIn);
        weightsDst.resize(n);
        std::fill(weightsDst.begin(), weightsDst.end(), 1.0 / n);
        ROFL_VAR3(musDst.size(), sigmasDst.size(), weightsDst.size());

        if (!musSrc.empty())
        {
            dsd::Vector3 gradIse, gradRtc;
            dsd::Matrix3 hessianIse, hessianRtc;
            dsd::Transform2 transformSrcDst;

            transformSrcDst = dsd::Transform2::Identity();
            gme::estimateTransformIse(transformSrcDst, musSrc, sigmasSrc,
                                      weightsSrc, musDst, sigmasDst,
                                      weightsDst);

            gme::computeGmmIseRtcGradHessian(
                musSrc, sigmasSrc, weightsSrc, musDst, sigmasDst, weightsDst,
                gradIse, hessianIse, gradRtc, hessianRtc);

            std::cout << "comparing src scan with " << musSrc.size()
                      << " points and dst scan with " << musDst.size()
                      << " points" << std::endl;

            Eigen::SelfAdjointEigenSolver<dsd::Matrix3> eigsolIse(hessianIse);
            double eigMinIse = fabs(eigsolIse.eigenvalues()(0));
            double eigMaxIse = fabs(eigsolIse.eigenvalues()(2));
            if (eigMaxIse < eigMinIse)
            {
                std::swap(eigMinIse, eigMaxIse);
            }
            std::cout << "ISE gradient: " << gradIse.transpose() << "\n"
                      << "ISE Hessian\n"
                      << hessianIse << "\n"
                      << "  eigenvalues " << eigsolIse.eigenvalues().transpose()
                      << ", min " << eigMinIse << " max " << eigMaxIse
                      << ", cond " << (eigMaxIse / eigMinIse) << std::endl;
            fileIse << counter << ", " << entry.stem() << ", "
                    << eigsolIse.eigenvalues()(0) << ", "
                    << eigsolIse.eigenvalues()(1) << ", "
                    << eigsolIse.eigenvalues()(2) << ", "
                    << (eigMaxIse / eigMinIse) << "\n";

            Eigen::SelfAdjointEigenSolver<dsd::Matrix3> eigsolRtc(hessianRtc);
            double eigMinRtc = fabs(eigsolRtc.eigenvalues()(0));
            double eigMaxRtc = fabs(eigsolRtc.eigenvalues()(2));
            if (eigMaxRtc < eigMinRtc)
            {
                std::swap(eigMinRtc, eigMaxRtc);
            }
            std::cout << "RTC gradient: " << gradRtc.transpose() << "\n"
                      << "RTC Hessian\n"
                      << hessianRtc << "\n"
                      << "  eigenvalues " << eigsolRtc.eigenvalues().transpose()
                      << ", min " << eigMinRtc << " max " << eigMaxRtc
                      << ", cond " << (eigMaxRtc / eigMinRtc) << std::endl;
            fileRtc << counter << ", " << entry.stem() << ", "
                    << eigsolRtc.eigenvalues()(0) << ", "
                    << eigsolRtc.eigenvalues()(1) << ", "
                    << eigsolRtc.eigenvalues()(2) << ", "
                    << (eigMaxRtc / eigMinRtc) << "\n";
        }

        musSrc = musDst;
        sigmasSrc = sigmasDst;
        weightsSrc = weightsDst;

        counter++;
    }

    fileIse.close();
    fileRtc.close();

    return 0;
}