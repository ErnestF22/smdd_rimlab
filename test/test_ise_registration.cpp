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

using PointT = pcl::PointXYZ;
using PointCloudT = pcl::PointCloud<PointT>;

void readScanCsv(const std::string& filename, dsd::VectorVector2& points);

void convertToGmm(const dsd::VectorVector2& points,
                  gme::IsotropicGaussianMixtureModel2::Ptr gmm,
                  double angleSigma);

void convertToPcl(const dsd::VectorVector2 points,
                  PointCloudT::Ptr cloud,
                  const dsd::Transform2& transf);

void plotGmm(pcl::visualization::PCLVisualizer::Ptr viewer,
             gme::IsotropicGaussianMixtureModel2::Ptr gmm,
             const dsd::Transform2& transf,
             std::string cloudLabel);

int main(int argc, char** argv) {
    std::string filenameCfg, filenameSrc, filenameDst;
    dsd::VectorVector2 pointsSrc, pointsDst;
    gme::IsotropicGaussianMixtureModel2::Ptr gmmSrc(
        new gme::IsotropicGaussianMixtureModel2);
    gme::IsotropicGaussianMixtureModel2::Ptr gmmDst(
        new gme::IsotropicGaussianMixtureModel2);
    gme::GmmRegistrationIse gmmIse;
    double angleSigma;
    dsd::Transform2 transf;
    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "src", filenameSrc, std::string("csv/degenere1/degenere1_000100.csv"));
    params.getParam<std::string>(
        "dst", filenameDst, std::string("csv/degenere1/degenere1_000101.csv"));
    params.getParam<double>("angleSigma", angleSigma, 0.017453);

    std::cout << "-------\n" << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    readScanCsv(filenameSrc, pointsSrc);
    std::cout << "read scan from file \"" << filenameSrc << "\": found "
              << pointsSrc.size() << " points" << std::endl;
    convertToGmm(pointsSrc, gmmSrc, angleSigma);
    std::cout << "  gmmSrc: means " << gmmSrc->means().size() << " sigmas "
              << gmmSrc->sigmas().size() << " weights "
              << gmmSrc->weights().size() << std::endl;

    readScanCsv(filenameDst, pointsDst);
    std::cout << "read scan from file \"" << filenameDst << "\": found "
              << pointsDst.size() << " points" << std::endl;
    convertToGmm(pointsDst, gmmDst, angleSigma);
    std::cout << "  gmmDst: means " << gmmDst->means().size() << " sigmas "
              << gmmDst->sigmas().size() << " weights "
              << gmmDst->weights().size() << std::endl;

    gmmIse.setSource(gmmSrc);
    gmmIse.setDestination(gmmDst);
    transf = dsd::Transform2::Identity();
    gmmIse.setInitGuess(transf);

    // Initializes viewer
    pcl::visualization::PCLVisualizer::Ptr viewer(
        new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0.9, 0.9, 0.9);
    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();
    viewer->setCameraPosition(1.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    transf = dsd::Transform2::Identity();
    gmmIse.init();
    for (int iter = 0; iter < 50; ++iter) {
        std::cout << "\n---\niteration " << iter << std::endl;
        gmmIse.iterateOnce();

        viewer->removeAllPointClouds();
        viewer->removeAllShapes();

        transf = gmmIse.getTransform();
        plotGmm(viewer, gmmSrc, transf, "cloud_src");
        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0,
            "cloud_src");

        plotGmm(viewer, gmmDst, dsd::Transform2::Identity(), "cloud_dst");
        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 1.0,
            "cloud_dst");
        viewer->spinOnce(100);
        std::this_thread::sleep_for(10ms);
    }

    while (!viewer->wasStopped()) {
        viewer->spinOnce(100);
        std::this_thread::sleep_for(100ms);
    }
    return 0;
}

void readScanCsv(const std::string& filename, dsd::VectorVector2& points) {
    CsvScan csvScanMain;
    readCsvScan(filename, csvScanMain);

    // Reads the scan, filters the invalid points or those belonging to the
    // Vehicle, fills the vector of points/mean values muIn
    points.clear();
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

            points.push_back(pIEig);
        }
    }
}

void convertToGmm(const dsd::VectorVector2& points,
                  gme::IsotropicGaussianMixtureModel2::Ptr gmm,
                  double angleSigma) {
    size_t n = points.size();
    gme::IsotropicGaussianMixtureModel2::VectorScalar sigmas(n);
    gme::IsotropicGaussianMixtureModel2::VectorScalar weights(n, 1.0 / n);
    for (int i = 0; i < n; ++i) {
        sigmas[i] = angleSigma * points[i].norm();
        ROFL_VAR3(angleSigma, points[i].norm(), sigmas[i]);
    }
    gmm->insert(points, sigmas, weights);
}

void convertToPcl(const dsd::VectorVector2 points,
                  PointCloudT::Ptr cloud,
                  const dsd::Transform2& transf) {
    PointT p;
    if (cloud.get() == nullptr) {
        cloud.reset(new PointCloudT);
    }
    cloud->resize(points.size());
    for (int i = 0; i < points.size(); ++i) {
        dsd::Vector2 ptrans = transf * points[i];
        cloud->points[i].x = ptrans(0);
        cloud->points[i].y = ptrans(1);
        cloud->points[i].z = 0.0;
    }
}

void plotGmm(pcl::visualization::PCLVisualizer::Ptr viewer,
             gme::IsotropicGaussianMixtureModel2::Ptr gmm,
             const dsd::Transform2& transf,
             std::string cloudLabel) {
    PointCloudT::Ptr cloud(new PointCloudT);
    convertToPcl(gmm->means(), cloud, transf);

    dsd::plotEllipsesArrows(viewer, gmm->means(), gmm->getCovariances(),
                       gmm->weights());

    viewer->addPointCloud<pcl::PointXYZ>(cloud, cloudLabel);
    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, cloudLabel);
}
