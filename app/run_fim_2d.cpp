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

#include <fim2d.h>
// #include <fim3d.h>
#include <orientation_rimlab.h>

int main(int argc, char **argv)
{
    std::string filenameCfg, filePath;
    std::string folderInPath, filenameOutFim, filenameOutTime;
    std::set<fs::path> sortedByName;

    dsd::VectorVector2 musIn;
    double sigmaIn;
    std::vector<double> sigmasIn;
    std::vector<double> weightsIn;
    float tx, ty, thetaRobot;
    float sigma;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>("in", folderInPath,
                                 std::string(""));
    params.getParam<std::string>("out_fim", filenameOutFim,
                                 std::string("fim.csv"));
    params.getParam<float>("sigma", sigma, 1.0f);
    params.getParam<float>("tx", tx, 0.0f);
    params.getParam<float>("ty", ty, 0.0f);
    params.getParam<float>("thetaRobot", thetaRobot, 0.0f);
    params.getParam<std::string>("out_time", filenameOutTime,
                                 std::string("time.txt"));

    std::cout << "-------\n"
              << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::ofstream fileFim(filenameOutFim);
    ROFL_ASSERT(fileFim);

    for (auto &entry : fs::directory_iterator(folderInPath))
        sortedByName.insert(entry.path());

    int counter = 0;
    for (const auto &entry : sortedByName)
    {
        // Reads the CSV file containing the point coordinates
        std::cout << "\n-----\n"
                  << counter << " " << entry.stem() << std::endl;
        CsvScan csvScanMain;
        readCsvScan(entry.c_str(), csvScanMain);

        std::vector<float> ranges;
        std::vector<float> phis;
        float r, phi;

        // Reads the scan, filters the invalid points or those belonging to the
        // Vehicle, fills the vector of points/mean values muIn
        musIn.clear();
        for (int i = 0; i < csvScanMain.pts.size(); ++i)
        {
            // copy scan into PCL cloud (2D)
            if (std::isfinite(csvScanMain.pts[i](0)) &&
                std::isfinite(csvScanMain.pts[i](1)))
            {
                Eigen::Vector2d pIVehicleCoord; // pI in Vehicle coordinates;
                Eigen::Vector2d pIEig;      // pI in Eigen
                pIEig << csvScanMain.pts[i](0), csvScanMain.pts[i](1);
                pIVehicleCoord = pIEig; // initializing it here with no real
                                    // purpose laserToVehicle(pIEig, pIVehicleCoord,
                                    // laser2VehicleX, laser2VehicleY, laser2VehicleT); if
                                    // (pI.getVector3fMap().norm() < fimRange
                                    // && !dsd::isVehicleShapePt(pIVehicleCoord))
                // ROFL_VAR1(pIEig.transpose());

                if (dsd::isVehicleShapePt(pIEig))
                    continue;

                musIn.push_back(pIEig);
                r = pIEig.norm();
                phi = atan2(pIEig(1), pIEig(0));
                ranges.push_back(r);
                phis.push_back(phi);
            }
        }

        size_t n = ranges.size();
        assert(phis.size() == n);
        std::vector<float> alphas(n, std::numeric_limits<double>::quiet_NaN()),
            cov0Alphas(n, 0.0f);
        int whsize = 5;
        // double theta0 = thetaRobot;
        for (int i = 0; i < n; i++)
        {
            std::vector<float> rangesPart;
            std::vector<float> phisPart; //(phis.begin(), phis.begin() + 100);
            int jmin = std::max<int>(0, i - whsize);
            int jmax = std::min<int>(n - 1, i + whsize);
            for (int j = jmin; j <= jmax; ++j)
            {
                if (j != i && isValid(ranges[i]))
                {
                    // ROFL_VAR2(i, j);
                    rangesPart.push_back(ranges[j]);
                    phisPart.push_back(phis[j]);
                }
            }
            double theta0 = phis[i];
            double rho0 = ranges[i];
            float alpha = 42, cov0Alpha = 32;
            if (!isValid(rho0))
            {
                continue;
            }
            int wsize = rangesPart.size();
            filterOrientationEigen(theta0, rho0, wsize, phisPart, rangesPart,
                                   alpha, cov0Alpha);
            // std::cout << "alpha " << alpha << std::endl;
            // std::cout << "cov0Alpha " << cov0Alpha << std::endl;
            // std::cout << "angleInc " << angleIncMain << std::endl;
            alphas[i] = alpha;
            cov0Alphas[i] = cov0Alpha;
        }

        Eigen::Matrix3f fim2D;
        computeFIM2D(fim2D, sigma, thetaRobot, alphas, phis, ranges);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigsolFim(fim2D);
        double eigMinFim = fabs(eigsolFim.eigenvalues()(0));
        double eigMaxFim = fabs(eigsolFim.eigenvalues()(2));
        if (eigMaxFim < eigMinFim)
        {
            std::swap(eigMinFim, eigMaxFim);
        }
        std::cout << "  eigenvalues " << eigsolFim.eigenvalues().transpose()
                  << ", min " << eigMinFim << " max " << eigMaxFim
                  << ", eigMax/eigMin " << (eigMaxFim / eigMinFim)
                  << ", eigMin/eigMax " << (eigMinFim / eigMaxFim) << std::endl;
        fileFim << counter << ", " << entry.stem() << ", "
                << eigsolFim.eigenvalues()(0) << ", "
                << eigsolFim.eigenvalues()(1) << ", "
                << eigsolFim.eigenvalues()(2) << ", " << (eigMaxFim / eigMinFim)
                << ", " << (eigMinFim / eigMaxFim) << "\n";

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