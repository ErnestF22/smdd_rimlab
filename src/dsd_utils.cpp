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
#include <dsd_utils.h>
#include <find_peak.h>
#include <rofl/common/macros.h>

namespace dsd
{
    std::string generateStampedString(const std::string prefix, const std::string postfix)
    {
        boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
        std::ostringstream formatter;
        std::string formatstring = prefix + "%Y%m%d_%H%M_%S" + postfix;
        formatter.imbue(std::locale(std::cout.getloc(), new boost::posix_time::time_facet(formatstring.c_str())));
        formatter << now;
        return formatter.str();
    }

    void diagonalize(const Matrix2 &m, double &lmin, double &lmax, double &theta)
    {
        double a, b, c, s;

        //   [ ct  st] * [m00 m01] * [ ct -st] = [ ct  st] * [m00*ct+m01*st,
        //   -m00*st+m01*ct]
        //   [-st  ct]   [m10 m11]   [ st  ct]   [-st  ct]   [m10*ct+m11*st,
        //   -m10*st+m11*ct]
        // = [ m00*ct*ct + m01*ct*st + m10*ct*st + m11*st*st, -m00*ct*st + m01*ct*ct
        // - m10*st*st + m11*ct*st ]
        //   [ -m00*ct*st + m01*ct*ct -m10*st*st + m11*ct*st,  m00*st*st - m01*ct*st
        //   - m10*st*ct + m11*ct*ct ]
        // non_diag = -m00*ct*st + m01*ct*ct - m10*st*st + m11*ct*st
        //          = ct * st * (m11 - m00) + ct^2 * m01 - st^2 * m10
        //          = ct * st * (m11 - m00) + (1 + cos(2*t)) / 2 * m01 - (1 -
        //          cos(2*t)) / 2 * m10 = ct * st * (m11 - m00) + (m01 - m10) / 2 +
        //          cos(2*t) * (m01 + m10) / 2 = sin(2*t) * a + cos(2*t) * b + (m01
        //          - m10) / 2

        // Diagonalizes sigma12
        a = 0.5 * (m(1, 1) - m(0, 0));
        b = 0.5 * (m(0, 1) + m(1, 0));
        // ROFL_VARIABLE2(a, b);

        theta = 0.5 * atan2(-b, a);

        c = cos(theta);
        s = sin(theta);
        lmax = m(0, 0) * c * c + m(1, 1) * s * s + (m(0, 1) + m(1, 0)) * c * s;
        lmin = m(0, 0) * s * s + m(1, 1) * c * c - (m(0, 1) + m(1, 0)) * c * s;
        // ROFL_VARIABLE3(theta, lmax, lmin);

        if (lmax < lmin)
        {
            theta += 0.5 * M_PI;
            std::swap(lmax, lmin);
            // ROFL_VAR1("after swap: lmin " << lmin << " < lmax " << lmax << ",
            // theta + PI/2: " << theta);
        }
    }

    void diagonalize(const Matrix2 &m, Matrix2 &l, Matrix2 &v)
    {
        double lmin, lmax, theta;

        diagonalize(m, lmin, lmax, theta);
        l = Matrix2::Zero();
        l.diagonal() << lmax, lmin;
        v = Eigen::Rotation2Dd(-theta);
    }

    void saturateEigenvalues(Matrix2 &covar, double sigmaMinSquare)
    {
        Matrix2 v;
        double lmin, lmax, theta;

        diagonalize(covar, lmin, lmax, theta);
        if (lmin < sigmaMinSquare)
        {
            lmin = sigmaMinSquare;
        }
        if (lmax < sigmaMinSquare)
        {
            lmax = sigmaMinSquare;
        }
        covar << lmax, 0.0, 0.0, lmin;
        v = Eigen::Rotation2Dd(theta);
        covar = v * covar * v.transpose();
    }

    double vonMises(double x, double mu, double k)
    {
        double bess = std::cyl_bessel_i(0, k);
        double denom = 1 / (2 * M_PI * bess);
        return denom * exp(k * cos(x - mu));
    }

    void vonMisesMixture(std::vector<double> &vomp,
                         int nSamples,
                         const std::vector<double> &mu,
                         const std::vector<double> &k,
                         const std::vector<double> &w)
    {
        double dtheta = 2.0 * M_PI / nSamples;

        vomp.resize(nSamples);
        vomp.assign(nSamples, 0.0f);
        for (int i = 0; i < nSamples; ++i)
        {
            double theta = dtheta * i;
            // vomp[i] = 0.0;
            for (int j = 0; j < mu.size(); ++j)
            {
                vomp[i] += w[j] * vonMises(theta, mu[j], k[j]);
                // ROFL_VAR1(vomp[i]);
            }
        }
    }

    void vonMisesStats(std::vector<double> &vomp,
                       int nSamples,
                       const VectorVector2 &mus,
                       const VectorMatrix2 &sigmas,
                       const std::vector<double> &weights)
    {
        std::vector<double> phis;
        std::vector<double> kappas;
        std::vector<double> weightsVomp;
        double phiTmp;

        size_t numGaussians = mus.size();
        phis.resize(2 * numGaussians);
        kappas.resize(2 * numGaussians);
        weightsVomp.resize(2 * numGaussians);

        for (int i = 0; i < numGaussians; ++i)
        {
            double lmin, lmax;
            dsd::diagonalize(sigmas[i], lmin, lmax, phiTmp);
            if (fabs(lmin) > 1e-5 && fabs(lmax) > 1e-5 && std::isfinite(lmin) &&
                std::isfinite(lmax))
            {
                if (lmin / lmax < 0)
                    continue;
            }
            else
                continue;
            phis[2 * i] = phiTmp;
            kappas[2 * i] = lmax / lmin;
            phis[2 * i + 1] = phiTmp + M_PI;
            kappas[2 * i + 1] = kappas[2 * i];
            weightsVomp[2 * i] = weights[i];
            weightsVomp[2 * i + 1] = weights[i];

            // ROFL_VAR3(i, phis[i], kappas[i]);
        }

        vonMisesMixture(vomp, nSamples, phis, kappas, weightsVomp);
    }

    void vonMisesMax(std::vector<int> &maximaIdx,
                     const std::vector<double> &vomp,
                     double angleWin,
                     int findPeaksWindow)
    {
        auto func = [&](int idx) -> double
        { return vomp[idx]; };
        double angleRes = 2.0 * M_PI / vomp.size();

        if (findPeaksWindow == 0)
            findPeaksWindow = (int)ceil(angleWin / angleRes);

        ROFL_VAR1(findPeaksWindow)
        double minval = 0.0;
        find_peak_circular(func, 0, vomp.size(), findPeaksWindow, minval,
                           std::back_inserter(maximaIdx));

        std::vector<int> maximaIdx2;
        find_peak(func, 0, vomp.size(), findPeaksWindow, minval,
                  std::back_inserter(maximaIdx2));

        if (maximaIdx.size() % 2 != 0)
        {
            for (int i = 0; i < vomp.size(); ++i)
            {
                ROFL_VAR3(i, angleRes * i, vomp[i]);
            }
            ROFL_VAR1("maxima of find_peak_circular");
            for (auto &m : maximaIdx)
                ROFL_VAR2(m, vomp[m]);
            ROFL_VAR1("maxima of find_peak");
            for (auto &m : maximaIdx2)
                ROFL_VAR2(m, vomp[m]);
            ROFL_ASSERT_VAR1(0, maximaIdx.size())
        }
    }

    bool isEqualFloats(double a, double b, double thr)
    {
        double val = fabs(a - b);
        if (val > thr)
            return false;
        return true;
    }

    void findSmallestEigs(float &eigvalMin,
                          Eigen::Vector2f &eigvec,
                          Eigen::Index &idx,
                          const Eigen::MatrixXf &covmat)
    {
        // ROFL_VAR1(covmat);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(covmat);
        // std::cout << "The eigenvalues of A are:" << std::endl
        //           << es.eigenvalues() << std::endl;
        // std::cout << "The matrix of eigenvectors, V, is:" << std::endl
        //           << es.eigenvectors() << std::endl
        //           << std::endl;

        Eigen::MatrixXf D = es.eigenvalues().asDiagonal();
        Eigen::MatrixXf V = es.eigenvectors();
        // std::cout << "Finally, V * D * V^(-1) = " << std::endl
        //           << V * D * V.inverse() << std::endl;

        eigvalMin = es.eigenvalues().minCoeff();
        // std::cout << "eigvalMin " << eigvalMin << std::endl;

        std::vector<Eigen::Index> idxs;
        // Searching for minimum eigenvalue index
        for (Eigen::Index i = 0; i < es.eigenvalues().size(); ++i)
            if (es.eigenvalues()(i) == eigvalMin)
                idxs.push_back(i);

        // ROFL_VAR1(idxs.size());
        idx = idxs[0]; // Maybe TODO here when same eigenvalue appears multiple
                       // times?
        eigvec = V.col(idx);
    }

    bool isVehicleShapePt(const Eigen::Vector2f &ptVehicleCoord,
                          float xMin,
                          float xMax,
                          float yMin,
                          float yMax)
    {
        if (ptVehicleCoord(0) > xMin && ptVehicleCoord(0) < xMax && ptVehicleCoord(1) > yMin &&
            ptVehicleCoord(1) < yMax)
            return true;

        return false;
    }

    bool isVehicleShapePt(const Eigen::Vector2d &ptVehicleCoord,
                          double xMin,
                          double xMax,
                          double yMin,
                          double yMax)
    {
        if (ptVehicleCoord(0) > xMin && ptVehicleCoord(0) < xMax && ptVehicleCoord(1) > yMin &&
            ptVehicleCoord(1) < yMax)
            return true;

        return false;
    }

    void plotEllipseArrowI(pcl::visualization::PCLVisualizer::Ptr viz,
                           int idx,
                           const Vector2 &mean,
                           const Matrix2 &covar,
                           const Matrix2 &eigvals,
                           const Matrix2 &eigvecs,
                           const double weight,
                           SaveCsvOptions saveOutput)
    {
        //    /** \brief Add an ellipsoid from the given parameters
        //          * \param[in] transform a transformation to apply to the
        //          ellipsoid from 0,0,0
        //          * \param[in] radius_x the ellipsoid's radius along its local
        //          x-axis
        //          * \param[in] radius_y the ellipsoid's radius along its local
        //          y-axis
        //          * \param[in] radius_z the ellipsoid's radius along its local
        //          z-axis
        //          * \param[in] id the ellipsoid id/name (default: "ellipsoid")
        //          * \param[in] viewport (optional) the id of the new viewport
        //          (default: 0)
        //          */
        //        bool
        //        addEllipsoid (const Eigen::Isometry3d &transform,
        //                      double radius_x, double radius_y, double radius_z,
        //                      const std::string &id = "ellipsoid",
        //                      int viewport = 0);

        double lmin, lmax, angle;
        //    set object 1 ellipse center 1.5, 1  size 6, 12  angle 60 front fs
        //    empty bo 3 plot '-' with points
        // Confidence 0.95 -> chi2 5.991 -> axis

        dsd::diagonalize(covar, lmin, lmax, angle);

        double translX = mean(0);
        double translY = mean(1);
        //    double rot = angle;

        Eigen::Isometry3d transfIso;
        Eigen::AffineCompact3d rotMat =
            Eigen::Affine3d(Eigen::AngleAxisd(angle, Eigen::Vector3d(0, 0, 1)));
        Eigen::AffineCompact3d translMat =
            Eigen::Affine3d(Eigen::Translation3d(mean(0), mean(1), 0.0));
        Eigen::AffineCompact3d transfMat = translMat * rotMat;
        transfIso = transfMat.matrix();

        //    double radiusX = 5.991 * sqrt(lmax); // -> change confidence interval
        double radiusX = 2 * sqrt(lmax);
        double radiusY = 2 * sqrt(lmin);
        double radiusZ = 0.05; // sigmaMin

        const std::string idEll =
            "ellipsoid" + boost::lexical_cast<std::string, int>(idx);
        if (!viz->addEllipsoid(transfIso, radiusX, radiusY, radiusZ, idEll))
        {
            std::cerr << "Failed adding following ellipsoid: " << idEll
                      << std::endl;
        }
        else
        {
            viz->setShapeRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, idEll);
            viz->setShapeRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_OPACITY, 0.2, idEll);

            // std::cout << "Added ellipsoid " << idEll
            //           << " corresponding to:" << std::endl;
            // std::cout << "set object " << idx << " ellipse center " << mean(0) <<
            // ", " << mean(1)
            //           << " size " << radiusX << ", " << radiusY << " angle " <<
            //           (180.0 / M_PI * angle)
            //           << " front fs empty bo 3\n";

            if (saveOutput.save)
            {
                // std::string scanIdStr = boost::lexical_cast<std::string,
                // uint64_t>(scan.id); std::string logName =
                // std::filesystem::path(filePath).stem().string(); if
                // (!fs::path("./ellipsoids/" + logName + "/").empty())
                //     isFirstTime = false;
                // saveCsvLog(scan, scanIdStr, logName, isFirstTime);
            }
        }

        size_t idxMax = 0;
        if (eigvals(0, 0) >= eigvals(1, 1))
            idxMax = 0;
        else
            idxMax = 1;

        const std::string idArrRight =
            "arrow_right_" + boost::lexical_cast<std::string, int>(idx);
        pcl::PointXYZ p1(pcl::PointXYZ(mean(0), mean(1), 0.0f));
        pcl::PointXYZ p2(pcl::PointXYZ(mean(0) + 0.01 * weight * eigvecs(idxMax, 0),
                                       mean(1) + 0.01 * weight * eigvecs(idxMax, 1),
                                       0.0f));
        if (!viz->addLine(p1, p2, 0.0, 1.0, 0.0, idArrRight))
        {
            std::cerr << "Failed adding following arrow: " << idArrRight
                      << std::endl;
        }
        else
        {
            // viz->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
            // 1.0, 0.0, 0.0, idArrRight);
            viz->setShapeRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, idArrRight);

            // std::cout << "Added arrow " << idArrRight
            //           << " corresponding to: " << std::endl;
            // ROFL_VAR1(weight);
            // std::cout << "set object " << idx << " arrow center " << mean(0) <<
            // ", " << mean(1)
            //           << " size " << radiusX << ", " << radiusY << " angle " <<
            //           (180.0 / M_PI * angle)
            //           << " front fs empty bo 3\n";
        }

        const std::string idArrLeft =
            "arrow_left_" + boost::lexical_cast<std::string, int>(idx);
        pcl::PointXYZ p1l(pcl::PointXYZ(mean(0), mean(1), 0.0f));
        pcl::PointXYZ p2l(pcl::PointXYZ(mean(0) - 0.01 * weight * eigvecs(0, 0),
                                        mean(1) - 0.01 * weight * eigvecs(0, 1),
                                        0.0f));
        if (!viz->addLine(p1l, p2l, 0.0, 1.0, 0.0, idArrLeft))
        {
            std::cerr << "Failed adding following arrow: " << idArrLeft
                      << std::endl;
        }
        else
        {
            // viz->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
            // 1.0, 0.0, 0.0, idArrRight);
            viz->setShapeRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, idArrLeft);

            // std::cout << "Added arrow " << idArrLeft
            //           << " corresponding to: " << std::endl;
            // ROFL_VAR1(weight);
            // std::cout << "set object " << idx << " arrow center " << mean(0) <<
            // ", " << mean(1)
            //           << " size " << radiusX << ", " << radiusY << " angle " <<
            //           (180.0 / M_PI * angle)
            //           << " front fs empty bo 3\n";
        }
    }

    void computeGaussianAxes(VectorMatrix2 &eigvals,
                             VectorMatrix2 &eigvecs,
                             const VectorMatrix2 &covar,
                             const std::vector<double> &weights)
    {
        eigvals.clear();
        eigvecs.clear();
        size_t sz = covar.size();
        for (size_t i = 0; i < sz; ++i)
        {
            const Matrix2 sigma = covar[i];
            Matrix2 eigval, eigvec;
            diagonalize(sigma, eigval, eigvec);
            eigvals.push_back(eigval);
            eigvecs.push_back(eigvec);
        }
    }

    void plotEllipsesArrows(pcl::visualization::PCLVisualizer::Ptr viz,
                            const dsd::VectorVector2 &means,
                            const dsd::VectorMatrix2 &covars,
                            const std::vector<double> &weights,
                            dsd::SaveCsvOptions saveOutput)
    {
        VectorMatrix2 eigvals, eigvecs;
        computeGaussianAxes(eigvals, eigvecs, covars, weights);
        // viz->removeAllShapes();
        for (int i = 0; i < means.size() && i < covars.size(); ++i)
        {
            plotEllipseArrowI(viz, i, means[i], covars[i], eigvals[i], eigvecs[i],
                              weights[i], saveOutput);
        }
    }

    void plotVomp(const std::vector<double> &vomp,
                  double vompMin,
                  double vompMax,
                  pcl::visualization::PCLVisualizer::Ptr viewer,
                  dsd::SaveCsvOptions saveOutput)
    {
        double radiusBig = 2.0;
        for (int i = 0; i < vomp.size(); ++i)
        {
            pcl::ModelCoefficients circleCoeff;
            double radius = 0.1, angle;
            angle = 2.0 * M_PI / vomp.size() * i;
            double ratio = (vomp[i] - vompMin) / (vompMax - vompMin);

            circleCoeff.values.resize(3); // We need 3 values
            circleCoeff.values[0] = (vomp[i] / vompMax) * radiusBig * cos(angle);
            circleCoeff.values[1] = (vomp[i] / vompMax) * radiusBig * sin(angle);
            circleCoeff.values[2] = radius;
            std::string circleId =
                "circle" + boost::lexical_cast<std::string, int>(i);
            viewer->addCircle(circleCoeff, circleId);

            if (saveOutput.save)
            {
                // std::string scanIdStr = boost::lexical_cast<std::string,
                // uint64_t>(scan.id); std::string logName =
                // std::filesystem::path(filePath).stem().string(); if
                // (!fs::path("./ellipsoids/" + logName + "/").empty())
                //     isFirstTime = false;
                // saveCsvLog(scan, scanIdStr, logName, isFirstTime);
            }

            viewer->setShapeRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 1.0 * ratio,
                circleId);
            viewer->setShapeRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_OPACITY, 0.9 * ratio, circleId);
            viewer->setShapeRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_SHADING_FLAT, 1, circleId);
        }
    }

    /**
     * Convert Velodyne point cloud from bin (GEODE dataset) format to PCL
     * mostly for visualization purposes
     */
    void binToPcl(const dsd::VectorVector3 &bin,
                  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
    {
        cloud->clear();
        for (int i = 0; i < bin.size(); ++i)
        {
            pcl::PointXYZ ptPcl(pcl::PointXYZ(bin[i](0), bin[i](1), bin[i](2)));
            // if (!removeVehicleShapePts)
            //     cloud->push_back(ptPcl);
            // else
            // {
            //     if (!isVehicleShapePt(ptPcl))
            //         cloud->push_back(ptPcl);
            // }
            cloud->push_back(ptPcl);
        }
    }

} // namespace dsd
