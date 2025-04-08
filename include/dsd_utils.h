/**
 * ARS - Angular Radon Spectrum
 * Copyright (C) 2017-2020 Dario Lodi Rizzini.
 *
 * ARS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ARS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ARS.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef DSD_UTILS_H_
#define DSD_UTILS_H_

#include <rofl/common/macros.h>
#include <cmath>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <eigen3/Eigen/Dense>

#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <rofl/common/macros.h>

#include <boost/lexical_cast.hpp>

#include <find_peaks_matrix.h>

namespace dsd
{
    using Vector2 = Eigen::Vector2d;

    using VectorVector2 = std::vector<Vector2>;

    using Matrix2 = Eigen::Matrix2d;

    using VectorMatrix2 = std::vector<Matrix2>;

    using Transform2 = Eigen::Affine2d;

    using Vector3 = Eigen::Vector3d;

    using VectorVector3 = std::vector<Vector3>;

    using Matrix3 = Eigen::Matrix3d;

    using VectorMatrix3 = std::vector<Matrix3>;

    using Transform3 = Eigen::Affine3d;

    using Vector4 = Eigen::Matrix<double, 4, 1>;

    using Vector5 = Eigen::Matrix<double, 5, 1>;

    using Vector6 = Eigen::Matrix<double, 6, 1>;

    using Matrix6 = Eigen::Matrix<double, 6, 6>;

    using Matrix4x6 = Eigen::Matrix<double, 4, 6>;

    std::string generateStampedString(const std::string prefix, const std::string postfix);

    /**
     * Compute szPadded
     */
    int computeSzPadded(int numGaussians, int nSamples, int numThreads);

    /**
     * Computes the diagonalization of the given positive definite matrix m.
     *
     *   m = rot(theta) * diag([lmax lmin]) rot(-theta)
     *
     * @param m the input positive definite matrix
     * @param lmin the minimum eigenvalue
     * @param lmax the maximum eigenvalue
     * @param theta the angle of the eigenvector corresponding to lmax w.r.t. axis x
     */
    void diagonalize(const Matrix2 &m, double &lmin, double &lmax, double &theta);

    /**
     * Computes the diagonalization of the given positive definite matrix m.
     * The relation among the matrices:
     *
     *   m = v * l * v.transpose()
     *
     * @param m the input positive definite matrix
     * @param l the matrix of eigenvalues (the maximum eigenvalue first)
     * @param v the matrix with eigenvectors on columns
     */
    void diagonalize(const Matrix2 &m, Matrix2 &l, Matrix2 &v);

    /**
     * Computes the diagonalization of the given positive definite matrix m.
     * The relation among the matrices:
     *
     *   m = v * l * v.transpose()
     *
     * @param m the input positive definite matrix
     * @param l the matrix of eigenvalues (the maximum eigenvalue first)
     * @param v the matrix with eigenvectors on columns
     *
     * Return false if
     * (m - v * l * v.transpose()).cwiseAbs().maxCoeff() < 1e-4
     * and return true otherwise.
     */
    bool diagonalize3d(const Matrix3 &m, Matrix3 &l, Matrix3 &v);

    /**
     * Saturates the eigenvalues of the input covariance matrix.
     * @param covar
     * @param sigmaMinSquare
     */
    void saturateEigenvalues(Matrix2 &covar, double sigmaMinSquare);

    /**
     * @brief Create a Rotation Matrix object, saving it in @param rotM, with params
     * @param ax
     * @param ay
     * @param az
     * as the rotations along corresponding axes (expressed in radians)
     */
    void createRotationMatrix(Eigen::Affine3d &rotM,
                              double ax,
                              double ay,
                              double az);

    /**
     * Von Mises distribution with scalar inputs
     * f(x | \mu, \kappa) = \frac{\exp(\kappa \cos(x - \mu))}{2 \pi I_0(\kappa)}
     */
    double vonMises(double k, double x, double mu);

    double vonMisesFisher(const Vector3 &x, const Vector3 &mu, double k);

    /**
     * @brief Compute von Mises distribution on a GMM
     */
    void vonMisesMixture(std::vector<double> &vomp,
                         int n,
                         const std::vector<double> &mu,
                         const std::vector<double> &k,
                         const std::vector<double> &w);

    void vonMisesMixture3d(std::vector<double> &vomp,
                           int nSamples,
                           const std::vector<Vector3> &mu,
                           const std::vector<double> &k,
                           const std::vector<double> &w,
                           int szPadded);

    /**
     * @brief Evaluate von Mises distribution with params @param mus, @param sigmas,
     * @param weights saving computed values in @param vomp which contains @param n
     * angles sample windows
     */
    void vonMisesStats(std::vector<double> &vomp,
                       int n,
                       const VectorVector2 &mus,
                       const VectorMatrix2 &sigmas,
                       const std::vector<double> &weights);

    /**
     * vonMisesStats3d BUT without last part of vonMisesMixture3d (which is performed in CUDA)
     */
    void vonMisesStats3dCuda(VectorVector3 &musVonMises,
                             std::vector<double> &kappas,
                             std::vector<double> &weightsVomp,
                             int szPadded,
                             const VectorVector3 &mus,
                             const VectorMatrix3 &sigmas,
                             const std::vector<double> &weights);

    void vonMisesStats3d(std::vector<double> &vomp,
                         int nSamples,
                         const VectorVector3 &mus,
                         const VectorMatrix3 &sigmas,
                         const std::vector<double> &weights,
                         int szPadded);

    /**
     * @brief Compute maxima of @param vomp
     * returning them in @param maximaIdx;
     * @param angleWin is the peakFinder's window
     */
    void vonMisesMax(std::vector<int> &maximaIdx,
                     const std::vector<double> &vomp,
                     double angleWin,
                     int findPeaksWindow = 0);

    void vompMax(std::vector<int> &maximaIdx,
                 std::vector<double> &maximaValues,
                 std::vector<std::pair<double, double>> &thetaPhiMaxima,
                 const std::vector<double> &vomp,
                 int nSamples,
                 double angleWin,
                 int findPeaksWindow = 0);

    /**
     * @brief Returns whether two floating point numbers are equal up to a threshold
     * thr
     */
    bool isEqualFloats(double a, double b, double thr = 1e-5);

    bool isFiniteVector3(const Vector3 &v);

    /**
     * @brief Struct used in CSV I/O
     */
    struct SaveCsvOptions
    {
        bool save;
        std::string filename;

        SaveCsvOptions() : save(false), filename("") {};

        SaveCsvOptions(bool s, std::string f)
        {
            save = s;
            filename = f;
        };

        ~SaveCsvOptions() {};
    };

    /**
     * @brief Find smallest eigenvalue and associated eigenvector (returning them in
     * reference params) returning also idx (depends on Eigen's default ordering
     * during EigenSolver)
     */
    void findSmallestEigs(float &eigvalMin,
                          Eigen::Vector2f &eigvec,
                          Eigen::Index &idx,
                          const Eigen::MatrixXf &covmat);

    /** @brief In Vehicle coordinates, the shape of the Vehicle points contained in interest
     * scans is: xMin= -2.0, xMax= 2.5, yMin=-1.0, yMax= 1.0. This method is useful
     * when we want to remove Vehicle points from the scan: given calibration values
     * (laser2VehicleX,laser2VehicleY,laser2VehicleT), we need to express scan points in Vehicle
     * coordinate ed eliminate the points that return true when passed into this
     * function FLOAT version*/
    bool isVehicleShapePt(const Eigen::Vector2f &ptVehicleCoord,
                          float xMin = -2.0,
                          float xMax = 2.5,
                          float yMin = -1.0,
                          float yMax = 1.0);

    /** @brief In Vehicle coordinates, the shape of the Vehicle points contained in interest
     * scans is: xMin= -2.0, xMax= 2.5, yMin=-1.0, yMax= 1.0. This method is useful
     * when we want to remove Vehicle points from the scan: given calibration values
     * (laser2VehicleX,laser2VehicleY,laser2VehicleT), we need to express scan points in Vehicle
     * coordinate ed eliminate the points that return true when passed into this
     * function DOUBLE version*/
    bool isVehicleShapePt(const Eigen::Vector2d &ptVehicleCoord,
                          double xMin = -2.0,
                          double xMax = 2.5,
                          double yMin = -1.0,
                          double yMax = 1.0);

    /**
     * @brief Plot an ellipse in PCL visualizer, with size depending on
     * input @param eigvals, @param eigvecs and @param weights
     */
    void plotEllipseArrowI(pcl::visualization::PCLVisualizer::Ptr viz,
                           int idx,
                           const Vector2 &mean,
                           const Matrix2 &covar,
                           const Matrix2 &eigvals,
                           const Matrix2 &eigvecs,
                           const double weight,
                           SaveCsvOptions saveOutput = SaveCsvOptions());

    /**
     * @brief Plot an ellipse in PCL visualizer, with size depending on
     * input @param eigvals, @param eigvecs and @param weights
     * (3D version of plotEllipseArrowI)
     */
    void plotSingleEllipseArrows3d(pcl::visualization::PCLVisualizer::Ptr viz,
                                   int idx,
                                   const Vector3 &mean,
                                   const Matrix3 &covar,
                                   const Matrix3 &eigvals,
                                   const Matrix3 &eigvecs,
                                   const double weight,
                                   SaveCsvOptions saveOutput);

    /**
     * @brief Compute axes of a Gaussian covariance matrix
     */
    void computeGaussianAxes(dsd::VectorMatrix2 &eigvals,
                             dsd::VectorMatrix2 &eigvecs,
                             const dsd::VectorMatrix2 &covar,
                             const std::vector<double> &weights);

    /**
     * @brief Compute axes of a Gaussian covariance matrix (3d case)
     */
    void computeGaussianAxes3d(VectorMatrix3 &eigvals,
                               VectorMatrix3 &eigvecs,
                               const VectorMatrix3 &covar,
                               const std::vector<double> &weights);

    /**
     * @brief Same as plotEllipseArrow but this works on an array of inputs
     */
    void plotEllipsesArrows(pcl::visualization::PCLVisualizer::Ptr viz,
                            const dsd::VectorVector2 &means,
                            const dsd::VectorMatrix2 &covars,
                            const std::vector<double> &weights,
                            dsd::SaveCsvOptions saveOutput = dsd::SaveCsvOptions());

    /**
     * @brief Same as plotEllipseArrows but this works on 3d inputs
     */
    void plotEllipsesArrows3d(pcl::visualization::PCLVisualizer::Ptr viz,
                              const dsd::VectorVector3 &means,
                              const dsd::VectorMatrix3 &covars,
                              const std::vector<double> &weights,
                              dsd::SaveCsvOptions saveOutput = dsd::SaveCsvOptions());

    /**
     * @brief Plot vomp distribution values (eight-shape when degenerate, shuriken
     * otherwise) as a weighted (color-wise and distance-from-center-wise) set of
     * circles
     */
    void plotVomp(const std::vector<double> &vomp,
                  double vompMin,
                  double vompMax,
                  pcl::visualization::PCLVisualizer::Ptr viewer,
                  dsd::SaveCsvOptions saveOutput = dsd::SaveCsvOptions());

    /**
     * Convert Velodyne point cloud from bin (GEODE dataset) format to PCL
     * mostly for visualization purposes
     */
    void binToPcl(const dsd::VectorVector3 &bin, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

} // namespace dsd

#endif /* DSD_UTILS_H_ */
