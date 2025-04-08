#include <iostream>

#include <rofl/common/macros.h>
#include <rofl/common/param_map.h>

#include <gme_GaussianMixtureEstimator.h>
#include <dsd_utils.h>

#include <csv_utils.h>

#include <rofl/common/profiler.h>

#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

int main(int argc, char **argv)
{
    std::string filenameCfg;
    std::string folderInPath;

    bool enableGmmHier;

    int nSamples;

    dsd::VectorVector2 musIn;
    double sigmaMin, res, iseThresh;
    dsd::VectorVector2 musOut;
    dsd::VectorMatrix2 sigmasOut;
    std::vector<double> weightsOut;
    double angleWin;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>("in", folderInPath,
                                 std::string("sample.csv"));

    // hier GME
    params.getParam<double>("sigmaMin", sigmaMin, 0.05);
    params.getParam<double>("res", res, 1);
    params.getParam<double>("iseThresh", iseThresh, 0.2);
    params.getParam<double>("angleWin", angleWin, M_PI / 180.0 * 5.0);

    params.getParam<int>("nSamples", nSamples, 360);

    params.getParam<bool>("enableGmmHier", enableGmmHier, true);

    // adapt tildes
    params.adaptTildeInPaths();
    params.getParam<std::string>("in", folderInPath, std::string("sample.csv"));

    std::cout << "-------\n"
              << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    // end of param reading section

    //--- filenames are unique so we can use a set
    std::set<fs::path> sortedByName;

    pcl::visualization::PCLVisualizer::Ptr viewer(
        new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0.9, 0.9, 0.9);
    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();

    int numZeroFails = 0;
    int numEqualFails = 0;
    int numBiggerFails = 0;

    double exectime = 0.0;

    bool firstScan = true;
    for (auto &entry : fs::directory_iterator(folderInPath))
        sortedByName.insert(entry.path());

    int numScans = sortedByName.size();

    for (const auto &entry : sortedByName)
    {
        numScans++;
        std::cout << entry.stem() << std::endl;

        CsvScan csvScanMain;
        readCsvScan(entry.c_str(), csvScanMain);

        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
            new pcl::PointCloud<pcl::PointXYZ>);
        musIn.clear();

        for (int i = 0; i < csvScanMain.pts.size(); ++i)
        {
            // copy scan into PCL cloud (2D)
            if (std::isfinite(csvScanMain.pts[i](0)) &&
                std::isfinite(csvScanMain.pts[i](1)))
            {
                pcl::PointXYZ pI(pcl::PointXYZ(csvScanMain.pts[i](0),
                                               csvScanMain.pts[i](1), 0.0f));
                Eigen::Vector2d pIVehicleCoord; // pI in Vehicle coordinates;
                Eigen::Vector2d pIEig;          // pI in Eigen
                pIEig << pI.x, pI.y;
                pIVehicleCoord = pIEig; // initializing it here with no real
                                        // purpose laserToVehicle(pIEig, pIVehicleCoord,
                                        // laser2VehicleX, laser2VehicleY, laser2VehicleT); if
                                        // (pI.getVector3fMap().norm() < fimRange
                                        // && !dsd::isVehicleShapePt(pIVehicleCoord))
                // ROFL_VAR1(pIEig.transpose());

                if (dsd::isVehicleShapePt(pIEig))
                    continue;

                cloud->push_back(pI);
                // std::cout << "i " << i << " pI " <<
                // pI.getVector3fMap().transpose() << " norm " <<
                // pI.getVector3fMap().norm() << std::endl;

                musIn.push_back(pIEig);
                // std::cout << "i " << i << " pIEig " << pIEig.transpose() <<
                // std::endl;
            }
        }
        ROFL_VAR1(musIn.size());

        double execTimeI = 0.0;

        /* Calling GMM Hierarchical Estimation */
        if (enableGmmHier)
        {
            rofl::ScopedTimer gmeTimer("gme hier timer");

            gme::GaussianMixtureEstimatorHierarchical gme;
            gme.setSigmaMin(sigmaMin);
            // gme.setCovarWidth(sigmaMin); // not used
            gme.setIseThreshold(iseThresh);
            gme.setCellSizeMax(res);
            gme.compute(musIn);

            gme.exportGaussians(musOut, sigmasOut, weightsOut);
            ROFL_VAR3(musOut.size(), sigmasOut.size(), weightsOut.size());
            execTimeI += gmeTimer.elapsedTimeMs();
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

            // std::cout << "Loading point cloud from \"" << srcFilename << "\""
            // << std::endl; if (pcl::io::loadPCDFile(srcFilename, *cloud) < 0)
            // {
            //     std::cerr << "Cannot load point cloud from \"" << srcFilename
            //     <<
            //     "\"" << std::endl; return 1;
            // }

            size_t cloudSz = cloud->size();
            std::cout << "Loaded cloud with " << cloudSz << " points"
                      << std::endl;

            std::for_each(weightsOut.begin(), weightsOut.end(),
                          [cloudSz](double &n)
                          { n *= cloudSz; });

            // Plot Ellipses for GMM covariance matrices and relative max eigenvector direction on PCL visualizer
            dsd::plotEllipsesArrows(viewer, musOut, sigmasOut, weightsOut);

            //
            std::vector<double> vomp;
            double maxVal, minVal;
            // musOut.erase(musOut.begin() + 1, musOut.end());
            // sigmasOut.erase(sigmasOut.begin() + 1, sigmasOut.end());
            // weightsOut.erase(weightsOut.begin() + 1, weightsOut.end());
            {
                rofl::ScopedTimer vonMisesStats("von Mises stats");
                dsd::vonMisesStats(vomp, nSamples, musOut, sigmasOut,
                                   weightsOut);
                execTimeI += vonMisesStats.elapsedTimeMs();
            }
            rofl::Profiler::getProfiler().printStats(std::cout);

            viewer->setPointCloudRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");
            viewer->setPointCloudRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_COLOR, 0.05, 0.05, 0.05,
                "cloud");

            // for (size_t i = 0; i < vomp.size(); ++i)
            // {
            // ROFL_VAR3(i, 360.0 * i / nSamples, vomp[i]);
            // }
            maxVal = *std::max_element(vomp.begin(), vomp.end());
            minVal = *std::min_element(vomp.begin(), vomp.end());
            int maxId = std::distance(vomp.begin(),
                                      std::max_element(vomp.begin(), vomp.end()));

            ROFL_VAR3(maxId, minVal, maxVal);

            /**
             * Plot VMM distribution values
             */
            dsd::plotVomp(vomp, minVal, maxVal, viewer);

            /**
             * Find peaks of VMM
             */
            std::vector<int> vompMaxima;
            {
                rofl::ScopedTimer vonMisesMax("von Mises max");
                dsd::vonMisesMax(vompMaxima, vomp, angleWin);
                execTimeI += vonMisesMax.elapsedTimeMs();
            }

            // compute error metrics (numZeroFails, numEqualFails, numWayBigger)
            if (vompMaxima.size() == 0)
            {
                numZeroFails++;
            }

            bool isEqualFl = false;
            for (int k = 0; k < vompMaxima.size(); ++k)
            {
                ROFL_VAR3(k, vompMaxima[k], vomp[vompMaxima[k]]);

                ROFL_ASSERT(vomp[vompMaxima[k]] > 0)

                // numEqualFails check
                if (k > 0 && vompMaxima.size() == 2 || k > 1) // this should be ok since they need to be equal in couples
                {
                    for (int k2 = 0; k2 < vompMaxima.size(); ++k2)
                    {
                        if (k2 == k)
                            continue;
                        if (dsd::isEqualFloats(vomp[vompMaxima[k]], vomp[vompMaxima[k2]]))
                        {
                            isEqualFl = true;
                            break;
                        }
                    }
                }
            }
            if (!isEqualFl)
            {
                ROFL_VAR1("Increasing numEqualFails")
                numEqualFails++;
            }

            // isWayBigger check
            bool isWayBigger = true;
            if (vompMaxima.size() > 2)
            {
                for (int k = 0; k < vompMaxima.size(); ++k)
                    for (int k2 = 0; k2 < vompMaxima.size(); ++k2)
                    {
                        if (k2 == k)
                            continue;

                        if (!dsd::isEqualFloats(vomp[vompMaxima[k]], vomp[vompMaxima[k2]]))
                        {
                            if (vomp[vompMaxima[k]] < vomp[vompMaxima[k2]])
                            {
                                if (2 * vomp[vompMaxima[k]] >= vomp[vompMaxima[k2]])
                                    isWayBigger = false;
                            }
                            else // if (vomp[vompMaxima[k]] > vomp[vompMaxima[k2]])
                                if (vomp[vompMaxima[k]] <= 2 * vomp[vompMaxima[k2]])
                                    isWayBigger = false;
                        }
                    }
            }
            if (!isWayBigger && !isEqualFl)
            {
                ROFL_VAR1("Increasing wayBigger")
                numBiggerFails++;
            }

            // save outputs to file:
            // vomp peaks
            std::string logName = std::filesystem::path(folderInPath)
                                      .parent_path()
                                      .filename()
                                      .string();
            std::string csvFilenameOut = "./" + logName +
                                         "_peaks"
                                         ".csv";
            ROFL_VAR1(csvFilenameOut);
            std::string entryStr = entry.string();
            if (firstScan)
            {
                std::ofstream fout;
                // std::string scanIdStr = boost::lexical_cast<std::string,
                // uint64_t>(scan.id);
                if (fs::exists(fs::path(csvFilenameOut)))
                    fout.open(csvFilenameOut, std::ios::out | std::ios::trunc);
                else
                    fout.open(csvFilenameOut, std::ios::out);

                std::string header = "& scan#, peaks#";

                fout << header << std::endl;
                unsigned first = entryStr.find_last_of('_') + 1;
                unsigned last = entryStr.find_last_of('.');
                std::string scanIdStr = entryStr.substr(first, last - first);
                if (vompMaxima.size() > 2)
                {
                    std::vector<double> vompSorted;
                    for (int i = 0; i < vompMaxima.size(); ++i)
                        vompSorted.push_back(vomp[vompMaxima[i]]);
                    std::sort(vompSorted.begin(), vompSorted.end(), std::greater<double>());
                    if (vompSorted[0] < 3 * vompSorted[2])
                    {
                        ROFL_VAR2(vompSorted[0], vompSorted[2])
                        fout << scanIdStr << ", " << 2.5 << std::endl;
                    }
                    else
                        fout << scanIdStr << ", " << vompMaxima.size() << std::endl;
                }
                else
                    fout << scanIdStr << ", " << vompMaxima.size() << std::endl;

                fout.close();
                firstScan = false;
            }
            else
            {
                std::ofstream fout;
                // std::string scanIdStr = boost::lexical_cast<std::string,
                // uint64_t>(scan.id);
                if (fs::path(csvFilenameOut).empty())
                {
                    ROFL_ASSERT_VAR1(0, fs::path(csvFilenameOut).empty());
                }
                else
                {
                    fout.open(csvFilenameOut, std::ios::out | std::ios::app);
                }

                unsigned first = entryStr.find_last_of('_') + 1;
                unsigned last = entryStr.find_last_of('.');
                std::string scanIdStr = entryStr.substr(first, last - first);
                if (vompMaxima.size() > 2)
                {
                    std::vector<double> vompSorted;
                    for (int i = 0; i < vompMaxima.size(); ++i)
                        vompSorted.push_back(vomp[vompMaxima[i]]);
                    std::sort(vompSorted.begin(), vompSorted.end(), std::greater<double>());
                    if (vompSorted[0] < 3 * vompSorted[2])
                    {
                        ROFL_VAR2(vompSorted[0], vompSorted[2])
                        fout << scanIdStr << ", " << 2.5 << std::endl;
                    }
                    else
                        fout << scanIdStr << ", " << vompMaxima.size() << std::endl;
                }
                else
                    fout << scanIdStr << ", " << vompMaxima.size() << std::endl; //!! TOFIX: bug that at times does not allow only scan id to be printed to file (but rather a longer string)
                fout.close();
            }

            // spin viewer/sleep thread
            viewer->spinOnce(500);
            //     std::this_thread::sleep_for(100ms);
        }
        exectime += execTimeI;
    }

    ROFL_VAR3(numZeroFails, numEqualFails, numBiggerFails);

    ROFL_VAR2(numScans, exectime / numScans);

    std::cout << "Mean exectime is in ms" << std::endl;

    // while (!viewer->wasStopped())
    // {
    //     viewer->spinOnce(100);
    //     std::this_thread::sleep_for(100ms);
    // }

    return 0;
}