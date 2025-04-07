#ifndef GAUSSIAN_METRIC_H_
#define GAUSSIAN_METRIC_H_

#include <cmath>

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>

#include <dsd_utils.h>
#include <rofl/common/macros.h>

namespace gme {

/**
 * @brief Computes the (left) gradient and Hessian of the square distance
 * function in the identity between the transformed point mu1 and the point mu2:
 *
 *    dist(T) = \| T * mu1 - mu2 \|^2
 *
 * @param mu1
 * @param mu2
 * @param distSq
 * @param grad
 * @param hessian
 */
void computeDistanceGradHessian(const dsd::Vector2& mu1,
                                const dsd::Vector2& mu2,
                                double& distSq,
                                dsd::Vector3& grad,
                                dsd::Matrix3& hessian);

void computeDistanceGradHessian3D(const dsd::Vector3& mu1,
                                  const dsd::Vector3& mu2,
                                  double& distSq,
                                  dsd::Vector6& grad,
                                  dsd::Matrix6& hessian);

/**
 * @brief Compute the Integral Square Error between two Gaussian
 * distributions.
 *
 * ISE(g1,g2) = \int_{R^d} (g(x - mu1, Sigma1) - g(x - mu2, Sigma2))^2 * dx
 *
 * @param mean1
 * @param covar1
 * @param mean2
 * @param covar2
 * @return double
 */
double computeGaussianIse(const dsd::Vector2& mean1,
                          const dsd::Matrix2& covar1,
                          const dsd::Vector2& mean2,
                          const dsd::Matrix2& covar2);

/**
 * @brief Compute the Integral Square Error using brute-force procedure, i.e.
 * by summing the ISE of each Gaussian:
 *
 * ISE(gmm1,gmm2) = \sum_{i,j} w1[i] * w2[j] * ISE[g1[i], g2[j]]
 *
 * where gmm1 = \sum_{i} w1[i] * g1[i], gmm2 = \sum_{i} w2[i] * g2[i]
 *
 */
double computeGmmIseBrute(const dsd::VectorVector2& means1,
                          const dsd::VectorMatrix2& covars1,
                          const std::vector<double>& weights1,
                          const dsd::VectorVector2& means2,
                          const dsd::VectorMatrix2& covars2,
                          const std::vector<double>& weights2);

/**
 * Computes the gradient and the Hessian of the ISE in the IDENTITY
 * between the transformed isotropic Gaussian (mean1, sigma1) and
 * the Gaussian (mean2, sigma2).
 */
void computeGaussianIseGradHessian(const dsd::Vector2& mean1,
                                   double sigma1,
                                   const dsd::Vector2& mean2,
                                   double sigma2,
                                   double& ise,
                                   dsd::Vector3& grad,
                                   dsd::Matrix3& hessian);

/**
 * Computes the gradient and the Hessian of the Randon Transform Correlation in
 * the IDENTITY between the transformed isotropic Gaussian (mean1, sigma1) and
 * the Gaussian (mean2, sigma2).
 */
void computeRadonCorrGradHessian(const dsd::Vector2& mean1,
                                 double sigma1,
                                 const dsd::Vector2& mean2,
                                 double sigma2,
                                 double& ise,
                                 dsd::Vector3& grad,
                                 dsd::Matrix3& hessian);

/**
 * Computes the gradient and Hessian of ISE and RTC in the identity
 * between two isotropic GMMs.
 */
void computeGmmIseRtcGradHessian(const dsd::VectorVector2& means1,
                                 const std::vector<double>& sigmas1,
                                 const std::vector<double>& weights1,
                                 const dsd::VectorVector2& means2,
                                 const std::vector<double>& sigmas2,
                                 const std::vector<double>& weights2,
                                 dsd::Vector3& gradIse,
                                 dsd::Matrix3& hessianIse,
                                 dsd::Vector3& gradRtc,
                                 dsd::Matrix3& hessianRtc);

/**
 * Computes the gradient and Hessian of ISE and RTC in the identity
 * between two isotropic GMMs in 3D space.
 */
void computeGmmIseRtcGradHessian3D(const dsd::VectorVector3& means1,
                                   const std::vector<double>& sigmas1,
                                   const std::vector<double>& weights1,
                                   const dsd::VectorVector3& means2,
                                   const std::vector<double>& sigmas2,
                                   const std::vector<double>& weights2,
                                   dsd::Vector6& gradIse,
                                   dsd::Matrix6& hessianIse,
                                   dsd::Vector6& gradRtc,
                                   dsd::Matrix6& hessianRtc);

/**
 * Computes the transformation using a local minimization of ISE.
 */
void estimateTransformIse(dsd::Transform2& transform,
                          const dsd::VectorVector2& meansSrc,
                          const std::vector<double>& sigmasSrc,
                          const std::vector<double>& weightsSrc,
                          const dsd::VectorVector2& meansDst,
                          const std::vector<double>& sigmasDst,
                          const std::vector<double>& weightsDst);

// ------------------------------------------------------
// ISOTROPIC GMM
// ------------------------------------------------------

/**
 * @brief class IsotropicGaussianMixtureModel2
 *
 */
class IsotropicGaussianMixtureModel2 {
   public:
    static const size_t Dim = 2;
    using Scalar = double;
    using Type = IsotropicGaussianMixtureModel2;
    using Ptr = std::shared_ptr<Type>;
    using VectorD = Eigen::Matrix<Scalar, Dim, 1>;
    using MatrixD = Eigen::Matrix<Scalar, Dim, Dim>;
    using VectorScalar = std::vector<Scalar>;
    using VectorVectorD = std::vector<VectorD>;
    using VectorMatrixD = std::vector<MatrixD>;

    using PointT = pcl::PointXY;
    using PointCloudT = pcl::PointCloud<PointT>;
    using SearchT = pcl::KdTreeFLANN<PointT>;

    IsotropicGaussianMixtureModel2();

    virtual ~IsotropicGaussianMixtureModel2();

    const VectorVectorD& means() const;

    const VectorScalar& sigmas() const;

    VectorMatrixD getCovariances() const;

    const VectorScalar& weights() const;

    size_t size() const;

    void insert(const VectorVectorD& means,
                const VectorScalar& sigmas,
                const VectorScalar& weights);

    void insert(const VectorVectorD& means, const VectorScalar& sigmas);

    void insert(const VectorVectorD& means, const Scalar& sigma);

    void initSearch();

    void findRange(const VectorD& query, Scalar r, std::vector<int>& indices);

   private:
    VectorVectorD means_;
    VectorScalar sigmas_;
    VectorScalar weights_;
    PointCloudT::Ptr cloud_;
    SearchT::Ptr search_;
};

// ------------------------------------------------------
// GMM ISE REGISTRATION
// ------------------------------------------------------

class GmmRegistrationIse {
   public:
    static const size_t Dim = 2;
    static const size_t L = 3;
    static constexpr double ALPHA = 1e-3;

    using Scalar = double;
    using GmmType = IsotropicGaussianMixtureModel2;
    using GmmPtr = std::shared_ptr<GmmType>;
    using VectorL = Eigen::Matrix<Scalar, L, 1>;
    using MatrixL = Eigen::Matrix<Scalar, L, L>;
    using TransformT = Eigen::Transform<Scalar, Dim, Eigen::Affine>;

    GmmRegistrationIse();

    virtual ~GmmRegistrationIse();

    void setSource(GmmPtr gmm);

    void setDestination(GmmPtr gmm);

    void setInitGuess(const TransformT& t);

    const TransformT& getTransform() const;

    void init();

    void iterateOnce();

   private:
    GmmPtr gmmSrc_;
    GmmPtr gmmDst_;
    double iseSelfSrc_;
    double iseSelfDst_;
    TransformT transform_;
    double iseLast_;
    VectorL gradLast_;
    MatrixL hessianLast_;
    int iterationMax_;

    double computeSelfIse(GmmPtr gmm);
};

}  // namespace gme

#endif