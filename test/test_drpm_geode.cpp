#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <set>
#include <string>
#include <filesystem>
#include <vector>

#include "drpm_degeneracy.h"

#include <rofl/common/macros.h>
#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <bin_utils.h>
#include <dsd_utils.h>

const double stdev_points = 0.100000000000;
const double stdev_normals = 0.050000000000;

namespace fs = std::filesystem;

int main(int argc, char **argv)
{
    std::string filenameCfg, filePath;
    std::string filenameIn, filenameOutIse, filenameOutRtc;
    std::set<fs::path> sortedByName;

    double sigmaIn;
    std::vector<double> sigmasIn;
    std::vector<double> weightsIn;
    bool enableGmmHier;
    bool plotEllipsesArrows;

    double sigmaMin, res, iseThresh;

    double angleWin;
    int nSamples;

    int numThreads;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", filenameIn,
        std::string("/home/rimlab/Datasets/Geode/1693022008.716670513.bin"));
    params.getParam<std::string>("out_ise", filenameOutIse,
                                 std::string("ise.csv"));
    params.getParam<std::string>("out_rtc", filenameOutRtc,
                                 std::string("rtc.csv"));
    params.getParam<double>("sigmaIn", sigmaIn, 0.10);
    params.getParam<double>("sigmaMin", sigmaMin, 0.05);
    params.getParam<double>("res", res, 1);
    params.getParam<double>("iseThresh", iseThresh, 0.2);
    params.getParam<double>("angleWin", angleWin, M_PI / 180.0 * 5.0);
    params.getParam<int>("nSamples", nSamples, 180);
    params.getParam<bool>("enableGmmHier", enableGmmHier, true);
    params.getParam<bool>("plotEllipsesArrows", plotEllipsesArrows, false);
    params.getParam<int>("numThreads", numThreads, 256);

    params.adaptTildeInPaths();

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n"
              << std::endl;

    // ROFL param reading end

    // Reading GEODE .bin file

    BinReader binReader;
    binReader.setVehiclePtsMinMax(dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0),
                                  dsd::Vector2(-3.0, 3.0));
    binReader.readCloudBin(filenameIn, BinReader::LidarType::VELODYNE);

    dsd::VectorVector3 musIn = binReader.getCloud();

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
        new pcl::PointCloud<pcl::PointXYZ>);
    dsd::binToPcl(musIn, cloud); // just for visualization purposes

    pcl::visualization::PCLVisualizer::Ptr viewer(
        new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0.9, 0.9, 0.9);
    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();

    viewer->addPointCloud<pcl::PointXYZ>(cloud, "cloud");

    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "cloud");
    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_COLOR, 0.05, 0.05, 0.05, "cloud");

    std::cout << "Read " << musIn.size() << " points" << std::endl;

    size_t n = musIn.size();
    for (size_t i = 0; i < n; ++i)
    {
        // ROFL_VAR2(i, musIn[i].transpose())
    }

    sigmasIn.resize(n);
    std::fill(sigmasIn.begin(), sigmasIn.end(), sigmaIn * sigmaIn);
    weightsIn.resize(n);
    std::fill(weightsIn.begin(), weightsIn.end(), 1.0 / n);

    // Points, normals and covariances must be expressed in the same frame of reference
    // For the conditioning of the Hessian, it is preferable to use the LiDAR frame (and not the world frame)
    const std::vector<Eigen::Vector3d> points = musIn;

    // Create the normal estimation class, and pass the input dataset to it
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    ne.setInputCloud(cloud);
    // Create an empty kdtree representation, and pass it to the normal estimation object.
    // Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
    ne.setSearchMethod(tree);
    // Output datasets
    pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>);
    // Use all neighbors in a sphere of radius 3cm
    ne.setRadiusSearch(0.5);
    // Compute the features
    ne.compute(*cloud_normals);
    // cloud_normals->size () should have the same size as the input cloud->size ()*

    std::vector<Eigen::Vector3d> normals;
    int ptsNanCtr = 0;
    for (int i = 0; i < cloud_normals->size(); ++i)
    {
        if (std::isnan(cloud->points[i].x) || std::isnan(cloud->points[i].y) || std::isnan(cloud->points[i].z))
            ptsNanCtr++;
        // else
        // {
        //     ROFL_VAR4(i, cloud->points[i].x, cloud->points[i].y, cloud->points[i].z);
        // }
        //
        Eigen::Vector3d normal(cloud_normals->points[i].normal_x,
                               cloud_normals->points[i].normal_y,
                               cloud_normals->points[i].normal_z);
        // if (!std::isnan(cloud_normals->points[i].curvature) && fabs(cloud_normals->points[i].curvature) > 1e-6)
        //     ROFL_VAR2(i, cloud_normals->points[i].curvature);
        normal.normalize();
        normals.push_back(normal);
    }

    ROFL_VAR1(ptsNanCtr)

    const std::vector<double> weights_squared(cloud->size(), 1.0); // does it need to be normalized??

    const auto normal_covariances = drpm_degeneracy::GetIsotropicCovariances(normals.size(), stdev_normals);

    // for (size_t i = 0; i < normal_covariances.size(); ++i)
    // {
    //     ROFL_VAR2(i, normal_covariances[i].transpose())
    // }

    ROFL_VAR3(points.size(), normals.size(), weights_squared.size());
    const auto H = drpm_degeneracy::ComputeHessian(points, normals, weights_squared);
    ROFL_VAR3(H.rows(), H.cols(), H);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> eigensolver(H);

    const auto eigenvectors = eigensolver.eigenvectors();
    const auto eigenvalues = eigensolver.eigenvalues();

    ROFL_VAR3(eigenvectors.rows(), eigenvectors.cols(), eigenvectors);

    Eigen::Matrix<double, 6, 6> noise_mean;
    Eigen::Matrix<double, 6, 1> noise_variance;
    const double snr_factor = 10.0;

    std::tie(noise_mean, noise_variance) = drpm_degeneracy::ComputeNoiseEstimate<double, double>(points, normals, weights_squared, normal_covariances, eigenvectors, stdev_points);
    Eigen::Matrix<double, 6, 1> non_degeneracy_probabilities = drpm_degeneracy::ComputeSignalToNoiseProbabilities<double>(H, noise_mean, noise_variance, eigenvectors, snr_factor);

    std::cout << "The non-degeneracy probabilities are: " << std::endl;
    std::cout << non_degeneracy_probabilities.transpose() << std::endl;

    std::cout << "For the eigenvectors of the Hessian: " << std::endl;
    std::cout << eigenvectors << std::endl;

    // // The following exemplifies how to solve the system of equations using the probabilities
    // // Dummy right hand side rhs = Jtb
    // const Eigen::Matrix<double, 6, 1> rhs = Eigen::Matrix<double, 6, 1>::Zero(6, 1);
    // const auto estimate = degeneracy::SolveWithSnrProbabilities(eigenvectors, eigenvalues, rhs, non_degeneracy_probabilities);

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
    }

    return 0;
}