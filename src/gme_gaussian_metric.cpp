#include <gme_gaussian_metric.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>

#include <rofl/common/functions.h>
#include <rofl/common/profiler.h>

namespace gme
{

    void computeDistanceGradHessian(const dsd::Vector2 &mu1,
                                    const dsd::Vector2 &mu2,
                                    double &distSq,
                                    dsd::Vector3 &grad,
                                    dsd::Matrix3 &hessian)
    {
        dsd::Vector3 mu12;
        dsd::Matrix3 M1, V1, M12, V12;
        dsd::Vector3 grad2;
        dsd::Matrix3 hessian2;

        mu12(0) = mu1(0) - mu2(0);
        mu12(1) = mu1(1) - mu2(1);
        mu12(2) = 0.0;
        M1 << 1.0, 0.0, -mu1(1), 0.0, 1.0, mu1(0), 0.0, 0.0, 0.0;
        V1 << 0.0, 0.0, mu1(1), 0.0, 0.0, -mu1(0), mu1(0), mu1(1), 0.0;
        M12 << 1.0, 0.0, -mu12(1), 0.0, 1.0, mu12(0), 0.0, 0.0, 0.0;
        V12 << 0.0, 0.0, mu12(1), 0.0, 0.0, -mu12(0), mu12(0), mu12(1), 0.0;

        // Square distance in the identity transform T = I
        distSq = mu12.dot(mu12);

        // By definition grad is a 1x3 row vector:
        //   grad = 2 mu12^T M_1
        // But here we represent it as a column vector, i.e. its transpose.
        grad = 2.0 * M1.transpose() * mu12;
        grad2 << 2.0 * mu12(0), 2.0 * mu12(1),
            2.0 * (mu12(1) * mu1(0) - mu12(0) * mu1(1));

        // Hessian
        hessian = M1.transpose() * M1 +
                  0.5 * (V12.transpose() * M1 + M1.transpose() * V12);
        hessian2 << 1.0, 0.0, -0.5 * mu12(1), 0.0, 1.0, 0.5 * mu12(0),
            -0.5 * mu12(1), 0.5 * mu12(0), 0.5 * mu1.dot(mu1 + mu2);
        // ROFL_VAR4(grad.transpose(), grad2.transpose(), hessian, hessian2);
    }

    void computeDistanceGradHessian3D(const dsd::Vector3 &mu1,
                                      const dsd::Vector3 &mu2,
                                      double &distSq,
                                      dsd::Vector6 &grad,
                                      dsd::Matrix6 &hessian)
    {
        dsd::Vector4 mu12;
        dsd::Matrix4x6 M1, V12;

        mu12(0) = mu1(0) - mu2(0);
        mu12(1) = mu1(1) - mu2(1);
        mu12(2) = mu1(2) - mu2(2);
        mu12(3) = 0.0;
        M1 = dsd::Matrix4x6::Zero();
        M1(0, 1) = mu1(2);
        M1(1, 0) = -mu1(2);
        M1(0, 2) = -mu1(1);
        M1(2, 0) = mu1(1);
        M1(1, 2) = mu1(0);
        M1(2, 1) = -mu1(0);
        M1(0, 3) = 1.0;
        M1(1, 4) = 1.0;
        M1(2, 5) = 1.0;
        V12 = dsd::Matrix4x6::Zero();
        V12(0, 1) = -mu12(2);
        V12(1, 0) = mu12(2);
        V12(0, 2) = mu12(1);
        V12(2, 0) = -mu12(1);
        V12(1, 2) = -mu12(0);
        V12(2, 1) = mu12(0);
        V12(3, 3) = mu12(0);
        V12(3, 4) = mu12(1);
        V12(3, 5) = mu12(2);

        // Square distance in the identity transform T = I
        distSq = mu12.dot(mu12);

        // By definition grad is a 1x3 row vector:
        //   grad = 2 mu12^T M_1
        // But here we represent it as a column vector, i.e. its transpose.
        grad = 2.0 * M1.transpose() * mu12;

        // Hessian
        hessian = M1.transpose() * M1 +
                  0.5 * (V12.transpose() * M1 + M1.transpose() * V12);
    }

    double computeGaussianIse(const dsd::Vector2 &mean1,
                              const dsd::Matrix2 &covar1,
                              const dsd::Vector2 &mean2,
                              const dsd::Matrix2 &covar2)
    {
        // Product of gaussians:
        //    g(x - mu1, Sigma1) * g(x - mu2, Sigma2)
        //  = g(mu2 - mu1, Sigma1 + Sigma2) * g(x - mu12, Sigma12)
        //
        //  where Sigma12 = (Sigma1^{-1} + Sigma2^{-1})^{-1},
        //        mu12 = Sigma12 * (Sigma1^{-1} * mu1 + Sigma2^{-1} * mu2)
        //
        //  Integral of product of Gaussians (with dimension d=2 in our case):
        //  PI[g1,g2]
        //  = \int_{R^d} g(x - mu1, Sigma1) * g(x - mu2, Sigma2) * dx
        //  = g(mu2 - mu1, Sigma1 + Sigma2) * \int_{R^d} g(x - mu12, Sigma12) * dx
        //  = g(mu2 - mu1, Sigma1 + Sigma2)
        //
        //  ISE(g1,g2)
        //  = \int_{R^d} (g1(x) - g2(x))^2 * dx
        //  = \int_{R^d} (g1(x)^2 + g2(x)^2 - 2 * g1(x) * g2(x))^2 * dx
        //  = PI[g1,g1] + PI[g2,g2] - 2 * PI[g1,g2]
        //  = (2 * pi * det(2 * Sigma1))^{-d/2} + (2 * pi * det(2 * Sigma2))^{-d/2}
        //    -2 * g(mu2 - mu1, Sigma1 + Sigma2)
        double square1 = 1.0 / (8.0 * M_PI * covar1.determinant());
        double square2 = 1.0 / (8.0 * M_PI * covar2.determinant());

        dsd::Matrix2 covar12 = covar1 + covar2;
        double norm12 = 1.0 / (2.0 * M_PI * covar12.determinant());
        double prod12 = norm12 * exp(-0.5 * (mean1 - mean2).transpose() *
                                     covar12.inverse() * (mean1 - mean2));
        return (square1 + square2 - 2.0 * prod12);
    }

    double computeGmmIseBrute(const dsd::VectorVector2 &means1,
                              const dsd::VectorMatrix2 &covars1,
                              const std::vector<double> &weights1,
                              const dsd::VectorVector2 &means2,
                              const dsd::VectorMatrix2 &covars2,
                              const std::vector<double> &weights2)
    {
        size_t n1 = means1.size();
        ROFL_ASSERT_VAR3(covars1.size() == n1 && weights1.size() == n1, n1,
                         covars1.size(), weights1.size());
        size_t n2 = means2.size();
        ROFL_ASSERT_VAR3(covars2.size() == n2 && weights2.size() == n2, n2,
                         covars2.size(), weights2.size());

        double ise = 0.0;
        for (size_t i = 0; i < n1; ++i)
        {
            for (size_t j = 0; j < n2; ++j)
            {
                ise += weights1[i] * weights2[j] *
                       computeGaussianIse(means1[i], covars1[i], means2[j],
                                          covars2[j]);
            }
        }
        return ise;
    }

    void computeGaussianIseGradHessian(const dsd::Vector2 &mean1,
                                       double sigma1,
                                       const dsd::Vector2 &mean2,
                                       double sigma2,
                                       double &ise,
                                       dsd::Vector3 &grad,
                                       dsd::Matrix3 &hessian)
    {
        double g1, g2, distSq, normK, sigmaSq12;
        dsd::Vector3 gradDist;
        dsd::Matrix3 hessianDist;

        rofl::ScopedTimer timer("ISE_HESSIAN_2D");

        computeDistanceGradHessian(mean1, mean2, distSq, gradDist, hessianDist);

        Eigen::SelfAdjointEigenSolver<dsd::Matrix3> eigsol(hessianDist);
        // ROFL_VAR2(hessianDist, eigsol.eigenvalues().transpose());

        // Gaussian function (with x = distSq):
        //   K = 1 / (2 * M_PI * (sigma1^2 + sigma2^2))
        //   g(x)  = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //   g1(x) = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //             * (-0.5 / (sigma1^2 + sigma2^2))
        //   g2(x) = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //             * (0.25 / (sigma1^2 + sigma2^2)^2)
        //
        // ISE[Gaussian1, Gaussian2]
        //    = K(sigma1) * g(0) + K(sigma2) * g(0) - 2 * g(distSq12)
        // Hence, gradient and Hessian must be multiplied by -2
        sigmaSq12 = sigma1 * sigma1 + sigma2 * sigma2;
        normK = 1.0 / (2.0 * M_PI * sigmaSq12);
        g1 = (-0.5 * normK / sigmaSq12) * exp(-0.5 * distSq / sigmaSq12);
        // ROFL_VAR3(0.25 * normK, (sigmaSq12 * sigmaSq12),
        //           exp(-0.5 * distSq / sigmaSq12))
        g2 = (0.25 * normK / (sigmaSq12 * sigmaSq12)) *
             exp(-0.5 * distSq / sigmaSq12);

        ise = normK * exp(-0.5 * distSq / sigmaSq12);
        grad = g1 * gradDist;
        hessian = (g2 * gradDist * gradDist.transpose() + g1 * hessianDist);
    }

    void computeGaussianIseGradHessian3D(const dsd::Vector3 &mean1,
                                         double sigma1,
                                         const dsd::Vector3 &mean2,
                                         double sigma2,
                                         double &ise,
                                         dsd::Vector6 &grad,
                                         dsd::Matrix6 &hessian)
    {
        double g1, g2, distSq, normK, sigmaSq12;
        dsd::Vector6 gradDist;
        dsd::Matrix6 hessianDist;

        rofl::ScopedTimer timer("ISE_HESSIAN_3D");

        computeDistanceGradHessian3D(mean1, mean2, distSq, gradDist, hessianDist);

        Eigen::SelfAdjointEigenSolver<dsd::Matrix6> eigsol(hessianDist);
        // ROFL_VAR2(hessianDist, eigsol.eigenvalues().transpose());

        // Gaussian function (with x = distSq):
        //   K = 1 / (2 * M_PI * (sigma1^2 + sigma2^2))
        //   g(x)  = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //   g1(x) = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //             * (-0.5 / (sigma1^2 + sigma2^2))
        //   g2(x) = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //             * (0.25 / (sigma1^2 + sigma2^2)^2)
        //
        // ISE[Gaussian1, Gaussian2]
        //    = K(sigma1) * g(0) + K(sigma2) * g(0) - 2 * g(distSq12)
        // Hence, gradient and Hessian must be multiplied by -2
        sigmaSq12 = sigma1 * sigma1 + sigma2 * sigma2;
        normK = 1.0 / (2.0 * M_PI * sigmaSq12);
        g1 = (-0.5 * normK / sigmaSq12) * exp(-0.5 * distSq / sigmaSq12);
        // ROFL_VAR3(0.25 * normK, (sigmaSq12 * sigmaSq12),
        //           exp(-0.5 * distSq / sigmaSq12))
        g2 = (0.25 * normK / (sigmaSq12 * sigmaSq12)) *
             exp(-0.5 * distSq / sigmaSq12);

        ise = normK * exp(-0.5 * distSq / sigmaSq12);
        grad = g1 * gradDist;
        hessian = (g2 * gradDist * gradDist.transpose() + g1 * hessianDist);

        // ROFL_VAR7(g1, g2, gradDist.transpose(), "\n  ", hessianDist, "\n  ",
        //           hessian);
    }

    void computeRadonCorrGradHessian(const dsd::Vector2 &mean1,
                                     double sigma1,
                                     const dsd::Vector2 &mean2,
                                     double sigma2,
                                     dsd::Vector3 &grad,
                                     dsd::Matrix3 &hessian)
    {
        double g1, g2, distSq, normK, sigmaSq12;
        dsd::Vector3 gradDist;
        dsd::Matrix3 hessianDist;
        std::vector<double> pnebis;

        rofl::ScopedTimer timer("RTC_HESSIAN_2D");

        sigmaSq12 = sigma1 * sigma1 + sigma2 * sigma2;
        computeDistanceGradHessian(mean1, mean2, distSq, gradDist, hessianDist);
        rofl::evaluatePnebiVector(2, distSq / (4.0 * sigmaSq12), pnebis);

        // Radon Transform Correlation (with x = distSq):
        //   K = 1 / (2 * M_PI * (sigma1^2 + sigma2^2))
        //   g(x)  = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //   g1(x) = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //             * (-0.5 / (sigma1^2 + sigma2^2))
        //   g2(x) = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //             * (0.25 / (sigma1^2 + sigma2^2)^2)
        g1 = pnebis[1] - pnebis[0];
        g2 = 1.5 * pnebis[0] - 2.0 * pnebis[1] + 0.5 * pnebis[2];

        grad = g1 * gradDist;
        hessian = g2 * gradDist * gradDist.transpose() + g1 * hessianDist;
    }

    void computeRadonCorrGradHessian3D(const dsd::Vector3 &mean1,
                                       double sigma1,
                                       const dsd::Vector3 &mean2,
                                       double sigma2,
                                       dsd::Vector6 &grad,
                                       dsd::Matrix6 &hessian)
    {
        double g1, g2, distSq, normK, sigmaSq12;
        dsd::Vector6 gradDist;
        dsd::Matrix6 hessianDist;
        std::vector<double> pnebis;

        rofl::ScopedTimer timer("RTC_HESSIAN_3D");

        sigmaSq12 = sigma1 * sigma1 + sigma2 * sigma2;
        computeDistanceGradHessian3D(mean1, mean2, distSq, gradDist, hessianDist);
        rofl::evaluatePnebiVector(2, distSq / (4.0 * sigmaSq12), pnebis);

        // Radon Transform Correlation (with x = distSq):
        //   K = 1 / (2 * M_PI * (sigma1^2 + sigma2^2))
        //   g(x)  = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //   g1(x) = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //             * (-0.5 / (sigma1^2 + sigma2^2))
        //   g2(x) = K * exp(-0.5 * x / (sigma1^2 + sigma2^2))
        //             * (0.25 / (sigma1^2 + sigma2^2)^2)
        g1 = pnebis[1] - pnebis[0];
        g2 = 1.5 * pnebis[0] - 2.0 * pnebis[1] + 0.5 * pnebis[2];

        grad = g1 * gradDist;
        hessian = g2 * gradDist * gradDist.transpose() + g1 * hessianDist;
    }

    void computeGmmIseRtcGradHessian(const dsd::VectorVector2 &means1,
                                     const std::vector<double> &sigmas1,
                                     const std::vector<double> &weights1,
                                     const dsd::VectorVector2 &means2,
                                     const std::vector<double> &sigmas2,
                                     const std::vector<double> &weights2,
                                     dsd::Vector3 &gradIse,
                                     dsd::Matrix3 &hessianIse,
                                     dsd::Vector3 &gradRtc,
                                     dsd::Matrix3 &hessianRtc)
    {
        dsd::Vector3 gradIse12, gradRtc12;
        dsd::Matrix3 hessianIse12, hessianRtc12;
        double ise12, iseTot;

        iseTot = 0.0;
        gradIse = dsd::Vector3::Zero();
        hessianIse = dsd::Matrix3::Zero();
        gradRtc = dsd::Vector3::Zero();
        hessianRtc = dsd::Matrix3::Zero();

        rofl::ScopedTimer timerIseRtcPair("ISE_RTC_PAIR_2D");

        // Creates a KDTree for indexing points
        using PointT = pcl::PointXY;

        pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
        for (int i = 0; i < means1.size(); ++i)
        {
            PointT p;
            p.x = means1[i](0);
            p.y = means1[i](1);
            cloud->push_back(p);
        }
        pcl::KdTreeFLANN<PointT>::Ptr tree(new pcl::KdTreeFLANN<PointT>);
        tree->setInputCloud(cloud);

        // Brute force procedure
        // for (int i1 = 0; i1 < means1.size(); ++i1) {
        //     for (int i2 = 0; i2 < means2.size(); ++i2) {
        //         if ((means1[i1] - means2[i2]).norm() >
        //             4.0 * sigmas1[i1] * sigmas2[i2])
        //             continue;
        //         computeGaussianIseGradHessian(means1[i1], sigmas1[i1],
        //         means2[i2],
        //                                       sigmas2[i2], grad12, hessian12);
        //         grad = grad + weights1[i1] * weights2[i2] * grad12;
        //         hessian = hessian + weights1[i1] * weights2[i2] * hessian12;
        //     }
        // }

        // Nearest neighbor computation.
        // Motivation: if | means1[i1] - means1[i2] |^2 << eps, then the
        // contributions of the pair to the gradient and the Hessian are negligible.
        std::vector<int> pointIdxRadiusSearch;
        std::vector<float> pointRadiusSquaredDistance;
        for (int i1 = 0; i1 < cloud->size(); ++i1)
        {
            double radius = 4.0 * sigmas1[i1];
            {
                rofl::ScopedTimer timerIseRtcSearch("ISE_RTC_SEARCH_2D");
                tree->radiusSearch(cloud->points[i1], radius, pointIdxRadiusSearch,
                                   pointRadiusSquaredDistance);
            }
            for (auto &i2 : pointIdxRadiusSearch)
            {
                computeGaussianIseGradHessian(means1[i1], sigmas1[i1], means2[i2],
                                              sigmas2[i2], ise12, gradIse12,
                                              hessianIse12);
                iseTot += ise12;
                gradIse -= 2.0 * weights1[i1] * weights2[i2] * gradIse12;
                hessianIse -= 2.0 * weights1[i1] * weights2[i2] * hessianIse12;

                computeRadonCorrGradHessian(means1[i1], sigmas1[i1], means2[i2],
                                            sigmas2[i2], gradRtc12, hessianRtc12);
                gradRtc = gradRtc + weights1[i1] * weights2[i2] * gradRtc12;
                hessianRtc =
                    hessianRtc + weights1[i1] * weights2[i2] * hessianRtc12;
            }
        }
        hessianIse *= 1.0 / iseTot;
    }

    void computeGmmIseRtcGradHessian3D(const dsd::VectorVector3 &means1,
                                       const std::vector<double> &sigmas1,
                                       const std::vector<double> &weights1,
                                       const dsd::VectorVector3 &means2,
                                       const std::vector<double> &sigmas2,
                                       const std::vector<double> &weights2,
                                       dsd::Vector6 &gradIse,
                                       dsd::Matrix6 &hessianIse,
                                       dsd::Vector6 &gradRtc,
                                       dsd::Matrix6 &hessianRtc)
    {
        dsd::Vector6 gradIse12, gradRtc12;
        dsd::Matrix6 hessianIse12, hessianRtc12;
        double ise12, iseTot;

        iseTot = 0.0;
        gradIse = dsd::Vector6::Zero();
        hessianIse = dsd::Matrix6::Zero();
        gradRtc = dsd::Vector6::Zero();
        hessianRtc = dsd::Matrix6::Zero();

        rofl::ScopedTimer timerIseRtcAll("ISE_RTC_ALL_3D");

        // Creates a KDTree for indexing points
        using PointT = pcl::PointXYZ;
        pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
        for (int i = 0; i < means2.size(); ++i)
        {
            PointT p;
            p.x = means2[i](0);
            p.y = means2[i](1);
            p.z = means2[i](2);
            cloud->push_back(p);
        }
        pcl::KdTreeFLANN<PointT>::Ptr tree(new pcl::KdTreeFLANN<PointT>);
        tree->setInputCloud(cloud);

        // Nearest neighbor computation.
        // Motivation: if | means1[i1] - means1[i2] |^2 << eps, then the
        // contributions of the pair to the gradient and the Hessian are negligible.
        std::vector<int> pointIdxRadiusSearch;
        std::vector<float> pointRadiusSquaredDistance;
        for (int i1 = 0; i1 < cloud->size(); ++i1)
        {
            double radius = 4.0 * sigmas1[i1];
            PointT q;
            q.x = means1[i1](0);
            q.y = means1[i1](1);
            q.z = means1[i1](2);
            {
                rofl::ScopedTimer timerIseRtcSearch("ISE_RTC_SEARCH_3D");
                tree->radiusSearch(q, radius, pointIdxRadiusSearch,
                                   pointRadiusSquaredDistance);
            }
            for (auto &i2 : pointIdxRadiusSearch)
            {
                computeGaussianIseGradHessian3D(means1[i1], sigmas1[i1], means2[i2],
                                                sigmas2[i2], ise12, gradIse12,
                                                hessianIse12);
                iseTot += ise12;
                gradIse -= 2.0 * weights1[i1] * weights2[i2] * gradIse12;
                hessianIse -= 2.0 * weights1[i1] * weights2[i2] * hessianIse12;
                // ROFL_VAR3(i1, i2, hessianIse.transpose());

                computeRadonCorrGradHessian3D(means1[i1], sigmas1[i1], means2[i2],
                                              sigmas2[i2], gradRtc12, hessianRtc12);
                gradRtc = gradRtc + weights1[i1] * weights2[i2] * gradRtc12;
                hessianRtc =
                    hessianRtc + weights1[i1] * weights2[i2] * hessianRtc12;
            }
        }
        hessianIse *= 1.0 / iseTot;
    }

    void estimateTransformIse(dsd::Transform2 &transform,
                              const dsd::VectorVector2 &meansSrc,
                              const std::vector<double> &sigmasSrc,
                              const std::vector<double> &weightsSrc,
                              const dsd::VectorVector2 &meansDst,
                              const std::vector<double> &sigmasDst,
                              const std::vector<double> &weightsDst)
    {
        dsd::Vector2 meanTransf;
        dsd::Vector3 gradIse, gradIsePair, x;
        dsd::Matrix3 hessianIse, hessianIsePair;
        double valIse, isePair;
        dsd::Transform2 transformInc;

        // Creates a KDTree for indexing points
        using PointT = pcl::PointXY;
        pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
        for (int i = 0; i < meansDst.size(); ++i)
        {
            PointT p;
            p.x = meansDst[i](0);
            p.y = meansDst[i](1);
            cloud->push_back(p);
        }
        pcl::KdTreeFLANN<PointT>::Ptr tree(new pcl::KdTreeFLANN<PointT>);
        tree->setInputCloud(cloud);

        //
        for (int k = 0; k < 60; ++k)
        {
            std::vector<int> pointIdxRadiusSearch;
            std::vector<float> pointRadiusSquaredDistance;
            valIse = 0.0;
            gradIse = dsd::Vector3::Zero();
            hessianIse = dsd::Matrix3::Zero();
            // Computes ISE as well as its gradient and Hessian at current pose
            for (int is = 0; is < meansSrc.size(); ++is)
            {
                meanTransf = transform * meansSrc[is];
                PointT pT;
                pT.x = meanTransf(0);
                pT.y = meanTransf(1);
                double radius = 8.0 * sigmasSrc[is];
                tree->radiusSearch(cloud->points[is], radius, pointIdxRadiusSearch,
                                   pointRadiusSquaredDistance);
                for (auto &id : pointIdxRadiusSearch)
                {
                    computeGaussianIseGradHessian(
                        meanTransf, sigmasSrc[is], meansDst[id], sigmasDst[id],
                        isePair, gradIsePair, hessianIsePair);
                    valIse += isePair;
                    gradIse =
                        gradIse + weightsSrc[is] * weightsDst[id] * gradIsePair;
                    hessianIse = hessianIse +
                                 weightsSrc[is] * weightsDst[id] * hessianIsePair;
                }
            }
            // Computes the new transform
            //   ISE[x] ~= const + gradIse^T * x + x^T * hessianIse * x
            // Local minimum estimated by solving linear system:
            //   2 * hessianIse * x = -gradIse
            // where x = [tx, ty, theta]
            x = hessianIse.colPivHouseholderQr().solve(-0.5 * gradIse);
            transformInc = Eigen::Rotation2D(x(2));
            transformInc.pretranslate(x.head<2>());
            transform = transformInc * transform;
            std::cout << "iteration " << k << ": ISE " << valIse << ": transform\n"
                      << transform.matrix() << std::endl;
        }
    }

    // ------------------------------------------------------
    // ISOTROPIC GMM
    // ------------------------------------------------------

    IsotropicGaussianMixtureModel2::IsotropicGaussianMixtureModel2()
        : means_(),
          sigmas_(),
          weights_(),
          cloud_(new PointCloudT),
          search_(new SearchT) {}

    IsotropicGaussianMixtureModel2::~IsotropicGaussianMixtureModel2() {}

    const IsotropicGaussianMixtureModel2::VectorVectorD &
    IsotropicGaussianMixtureModel2::means() const
    {
        return means_;
    }

    const IsotropicGaussianMixtureModel2::VectorScalar &
    IsotropicGaussianMixtureModel2::sigmas() const
    {
        return sigmas_;
    }

    const IsotropicGaussianMixtureModel2::VectorScalar &
    IsotropicGaussianMixtureModel2::weights() const
    {
        return weights_;
    }

    dsd::VectorMatrix2 IsotropicGaussianMixtureModel2::getCovariances() const
    {
        size_t n = sigmas_.size();
        dsd::VectorMatrix2 covars(n);
        for (int i = 0; i < n; ++i)
        {
            double sq = sigmas_[i] * sigmas_[i];
            covars[i] << sq, 0.0, 0.0, sq;
        }
        return covars;
    }

    size_t IsotropicGaussianMixtureModel2::size() const
    {
        return means_.size();
    }

    void IsotropicGaussianMixtureModel2::insert(const VectorVectorD &means,
                                                const VectorScalar &sigmas,
                                                const VectorScalar &weights)
    {
        ROFL_ASSERT_VAR3(
            means.size() == sigmas.size() && means.size() == weights.size(),
            means.size(), sigmas.size(), weights.size());
        means_ = means;
        sigmas_ = sigmas;
        weights_ = weights;
    }

    void IsotropicGaussianMixtureModel2::insert(const VectorVectorD &means,
                                                const VectorScalar &sigmas)
    {
        ROFL_ASSERT_VAR2(means.size() == sigmas.size(), means.size(),
                         sigmas.size());
        means_ = means;
        sigmas_ = sigmas;
        size_t n = means.size();
        if (n > 0)
        {
            weights_.resize(n);
            weights_.assign(n, Scalar(1.0 / n));
        }
    }

    void IsotropicGaussianMixtureModel2::insert(const VectorVectorD &means,
                                                const Scalar &sigma)
    {
        means_ = means;
        size_t n = means.size();

        if (n > 0)
        {
            sigmas_.resize(n);
            sigmas_.assign(n, sigma);
            weights_.resize(n);
            weights_.assign(n, Scalar(1.0 / n));
        }
    }

    void IsotropicGaussianMixtureModel2::initSearch()
    {
        PointT p;
        cloud_->clear();
        for (int i = 0; i < means_.size(); ++i)
        {
            p.x = means_[i](0);
            p.y = means_[i](1);
            cloud_->push_back(p);
        }
        search_->setInputCloud(cloud_);
    }

    void IsotropicGaussianMixtureModel2::findRange(const VectorD &query,
                                                   Scalar r,
                                                   std::vector<int> &indices)
    {
        PointT q;
        q.x = query(0);
        q.y = query(1);
        std::vector<float> distances;
        search_->radiusSearch(q, (float)r, indices, distances);
    }

    // ------------------------------------------------------
    // GMM ISE REGISTRATION
    // ------------------------------------------------------

    GmmRegistrationIse::GmmRegistrationIse()
        : gmmSrc_(new GmmType),
          gmmDst_(new GmmType),
          iseSelfSrc_(0.0),
          iseSelfDst_(0.0),
          transform_(TransformT::Identity()),
          iterationMax_(100) {}

    GmmRegistrationIse::~GmmRegistrationIse() {}

    void GmmRegistrationIse::setSource(GmmPtr gmm)
    {
        gmmSrc_ = gmm;
    }

    void GmmRegistrationIse::setDestination(GmmPtr gmm)
    {
        gmmDst_ = gmm;
    }

    void GmmRegistrationIse::setInitGuess(const TransformT &t)
    {
        transform_ = t;
    }

    const GmmRegistrationIse::TransformT &GmmRegistrationIse::getTransform() const
    {
        return transform_;
    }

    void GmmRegistrationIse::init()
    {
        gmmSrc_->initSearch();
        gmmDst_->initSearch();
        // Computes ISE constant terms
        iseSelfSrc_ = computeSelfIse(gmmSrc_);
        iseSelfDst_ = computeSelfIse(gmmDst_);
        ROFL_VAR2(iseSelfSrc_, iseSelfDst_);
    }

    void GmmRegistrationIse::iterateOnce()
    {
        size_t numSrc, numDst;
        std::vector<int> indicesDst;
        Scalar radius, sigma, iseSD, iseTot;
        VectorL gradSD, gradTot, x;
        MatrixL hessianSD, hessianTot;
        TransformT delta;

        numSrc = gmmSrc_->means().size();
        numDst = gmmDst_->means().size();
        iseTot = 0.0;
        gradTot = VectorL::Zero();
        hessianTot = MatrixL::Zero();
        for (size_t is = 0; is < numSrc; ++is)
        {
            indicesDst.clear();
            // 1/sqrt(2 * PI *sigma^2) * exp(d^2 / (2 sigma^2)) = alpha
            // exp(-d^2 / (2 sigma^2)) =
            // d^2 = 2 * sigma^2 * log(sqrt(2 * PI) * sigma * alpha)
            // d = sqrt(2) * sigma * sqrt(log(sqrt(2 * PI) * sigma * alpha))
            sigma = gmmSrc_->sigmas()[is];
            radius =
                sqrt(2.0) * sigma * sqrt(log(sqrt(2.0 * M_PI) * sigma * ALPHA));
            radius = 0.05;
            gmmDst_->findRange(gmmSrc_->means()[is], radius, indicesDst);
            // ROFL_VAR4(is, sigma, radius, indicesDst.size());

            // Computes cross-kernel ISE
            for (auto &id : indicesDst)
            {
                computeGaussianIseGradHessian(
                    gmmSrc_->means()[is], gmmSrc_->sigmas()[is],
                    transform_.inverse() * gmmDst_->means()[id],
                    gmmDst_->sigmas()[id], iseSD, gradSD, hessianSD);
                iseTot = iseTot +
                         gmmSrc_->weights()[is] * gmmDst_->weights()[id] * iseSD;
                gradTot = gradTot +
                          gmmSrc_->weights()[is] * gmmDst_->weights()[id] * gradSD;
                hessianTot = hessianTot + gmmSrc_->weights()[is] *
                                              gmmDst_->weights()[id] * hessianSD;
            }
        }
        iseTot = iseSelfSrc_ + iseSelfDst_ - 2.0 * iseTot;

        // Computes the new transform
        //   f(x) ~= const + gradIse^T * x + x^T * hessianIse * x
        // Local minimum estimated by solving linear system:
        //   2 * hessianIse * x = -gradIse
        // where x = [tx, ty, theta]
        x = hessianTot.colPivHouseholderQr().solve(-0.5 * gradTot);
        delta = Eigen::Rotation2D(x(2));
        delta.pretranslate(x.head<2>());
        transform_ = transform_ * delta;

        Eigen::SelfAdjointEigenSolver<dsd::Matrix3> eigsol(hessianTot);

        ROFL_MSG("ise " << iseTot << "\ngrad " << gradTot.transpose()
                        << "\nHessian\n"
                        << hessianTot << "\nHessian eigenvalues: "
                        << eigsol.eigenvalues().transpose() << "\ntransform\n"
                        << transform_.matrix());
    }

    double GmmRegistrationIse::computeSelfIse(GmmPtr gmm)
    {
        double ise, distSq, sigmaSq, normK, radius;
        std::vector<int> indices;

        gmm->initSearch();
        ise = 0.0;
        for (size_t is = 0; is < gmm->size(); ++is)
        {
            indices.clear();
            radius = 0.05;
            gmm->findRange(gmm->means()[is], radius, indices);

            // Computes cross-kernel ISE
            for (auto &id : indices)
            {
                sigmaSq = gmm->sigmas()[is] * gmm->sigmas()[is] +
                          gmm->sigmas()[id] * gmm->sigmas()[id];
                normK = 1.0 / (2.0 * M_PI * sigmaSq);
                distSq = (gmm->means()[is] - gmm->means()[id]).norm();
                distSq = distSq * distSq;
                ise += gmm->weights()[is] * gmm->weights()[id] * normK *
                       exp(-0.5 * distSq / sigmaSq);
            }
        }
        return ise;
    }

} // namespace gme