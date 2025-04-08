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
    int computeSzPadded(int numGaussians, int nSamples, int numThreads)
    {
        int szPadded = 2 * numGaussians - 1; //!! -1
        int mod = -1;
        int totNumSamples = 2 * nSamples * nSamples;
        while (mod != 0)
        {
            szPadded++;
            mod = (szPadded * totNumSamples) % numThreads;
        }
        ROFL_VAR1(szPadded);
        return szPadded;
    }

    void createRotationMatrix(Eigen::Affine3d &rotM,
                              double ax,
                              double ay,
                              double az)
    {
        Eigen::Affine3d rx =
            Eigen::Affine3d(Eigen::AngleAxisd(ax, Eigen::Vector3d(1, 0, 0)));
        Eigen::Affine3d ry =
            Eigen::Affine3d(Eigen::AngleAxisd(ay, Eigen::Vector3d(0, 1, 0)));
        Eigen::Affine3d rz =
            Eigen::Affine3d(Eigen::AngleAxisd(az, Eigen::Vector3d(0, 0, 1)));
        rotM = rz * ry * rx;
    }

    bool diagonalize3d(const Matrix3 &m, Matrix3 &l, Matrix3 &v)
    {
        // double lmin, lmax, theta;
        // diagonalize(m, lmin, lmax, theta);
        // l = Matrix2::Zero();
        // l.diagonal() << lmax, lmin;
        // v = Eigen::Rotation2Dd(-theta);

        Eigen::EigenSolver<Matrix3> eigensolver;
        eigensolver.compute(m);
        l.setZero();
        l.diagonal() = eigensolver.eigenvalues().real();
        v = eigensolver.eigenvectors().real();

        if ((m - v * l * v.transpose()).cwiseAbs().maxCoeff() < 1e-4)
            return false;
        return true;
    }

    double vonMisesFisher(const Vector3 &x, const Vector3 &mu, double k)
    {
        // in this context k is the weight
        double c3k = k / (4 * M_PI * sinh(k));
        double e = exp(k * mu.dot(x));
        return c3k * e;
    }

    void vonMisesMixture3d(std::vector<double> &vomp,
                           int nSamples,
                           const std::vector<Vector3> &mu,
                           const std::vector<double> &k,
                           const std::vector<double> &w,
                           int szPadded)
    {
        double dtheta = M_PI / nSamples;

        int numCols = 2 * nSamples * nSamples;

        int vompSz = numCols;
        vomp.resize(vompSz);
        vomp.assign(vompSz, 0.0f);
        Eigen::MatrixXi validsCtr(Eigen::MatrixXi::Zero(2 * nSamples, nSamples));
        int index = 0;
        int sizePhis = 2 * nSamples;
        int sizeThetas = nSamples;

        for (int jTheta = 0; jTheta < sizeThetas; ++jTheta)
        {
            double theta = dtheta * jTheta;

            for (int jPhi = 0; jPhi < sizePhis; ++jPhi)
            {
                // ROFL_VAR2(jTheta, jPhi)
                int j = jTheta * sizePhis + jPhi;
                double phi = dtheta * jPhi;

                // vomp[i] = 0.0;
                for (int i = 0; i < szPadded; ++i)
                {
                    if (i >= mu.size())
                    {
                        ++index;
                        continue;
                    }
                    // ROFL_VAR1(mukwSz)
                    Vector3 xThetaPhi(sin(theta) * cos(phi), sin(theta) * sin(phi),
                                      cos(theta));
                    if (isFiniteVector3(xThetaPhi))
                    {
                        auto vmf = vonMisesFisher(xThetaPhi, mu[i], k[i]);
                        if (std::isfinite(vmf))
                        {
                            validsCtr(jPhi, jTheta)++;
                            vomp[j] += w[i] * vmf; //!!!!
                            // DEBUGGING CUDA
                            // int i = index % mukwSz;
                            int jCheck = (index - i) / szPadded;

                            if (i >= szPadded)
                                printf("i >= mukwSz!! %d \t %d\n", i, szPadded);

                            if (jCheck >= numCols)
                                printf("j >= numCols!! %d \t %d\n", j, numCols);

                            ROFL_ASSERT_VAR2(j == jCheck, j, jCheck)

                            // printf("index %d \t i %d \t j %d \t jTheta %d \t jPhi %d \t theta %f \t phi %f \t vmf %f\n",
                            //        index, i, j, jTheta, jPhi, theta, phi, vmf);
                        }
                        else
                        {
                            // ROFL_VAR2(j, "vmfJ NOT finite")
                        }
                    }
                    else
                    {
                        // ROFL_VAR2(j, "xThetaPhi NOT finite")
                    }
                    // ROFL_VAR1(vomp[j])
                    index++;
                }
            }
        }
        ROFL_VAR1(validsCtr)
        ROFL_VAR2(validsCtr.sum(), nSamples * 2 * nSamples)
    }

    void vonMisesStats3dCuda(VectorVector3 &musVomp,
                             std::vector<double> &kappasVomp,
                             std::vector<double> &weightsVomp,
                             int szPadded,
                             const VectorVector3 &mus,
                             const VectorMatrix3 &sigmas,
                             const std::vector<double> &weights)
    {
        double phiTmp;

        size_t numGaussians = mus.size();

        musVomp.resize(szPadded, Vector3::Zero());
        weightsVomp.resize(szPadded, 0.0);
        kappasVomp.resize(szPadded, 0.0);

        int validDiagCtr = 0;
        for (int i = 0; i < numGaussians; ++i)
        {
            // ROFL_VAR1(i)
            Matrix3 eigvalsMat, eigvecs;
            bool isDiVehiclealid = dsd::diagonalize3d(sigmas[i], eigvalsMat, eigvecs);
            if (!isDiVehiclealid)
            {
                musVomp[2 * i] = Vector3::Zero();
                kappasVomp[2 * i] = 0;
                musVomp[2 * i + 1] = Vector3::Zero();
                kappasVomp[2 * i + 1] = kappasVomp[2 * i];
                weightsVomp[2 * i] = 0;
                weightsVomp[2 * i + 1] = 0;
                continue;
            }
            auto eigvals = eigvalsMat.diagonal();
            Eigen::Index minIdx, maxIdx;
            double lmin = eigvals.minCoeff(&minIdx);
            double lmax = eigvals.maxCoeff(&maxIdx);
            if (fabs(lmin) > 1e-5 && fabs(lmax) > 1e-5 && std::isfinite(lmin) &&
                std::isfinite(lmax))
            {
                ROFL_ASSERT_VAR2((lmin / lmax) >= 0, lmin, lmax);
            }
            else
                continue;

            validDiagCtr++;
            // ROFL_VAR1(validDiagCtr);
            musVomp[2 * i] = eigvecs.col(minIdx);
            if (lmax > 10.0 * lmin)
            {
                kappasVomp[2 * i] = 10.0;
            }
            else
            {
                kappasVomp[2 * i] = lmax / lmin;
            }

            musVomp[2 * i + 1] = -eigvecs.col(minIdx);
            kappasVomp[2 * i + 1] = kappasVomp[2 * i];
            weightsVomp[2 * i] = weights[i];
            weightsVomp[2 * i + 1] = weights[i];

            // ROFL_VAR3(i, phis[i], kappas[i]);
        }
    }

    void vonMisesStats3d(std::vector<double> &vomp,
                         int nSamples,
                         const VectorVector3 &mus,
                         const VectorMatrix3 &sigmas,
                         const std::vector<double> &weights,
                         int szPadded)
    {
        VectorVector3 musVonMises;
        std::vector<double> kappas;
        std::vector<double> weightsVomp;
        double phiTmp;

        size_t numGaussians = mus.size();
        musVonMises.resize(2 * numGaussians);
        kappas.resize(2 * numGaussians);
        weightsVomp.resize(2 * numGaussians);

        int validDiagCtr = 0;
        for (int i = 0; i < numGaussians; ++i)
        {
            ROFL_VAR1(i)
            Matrix3 eigvalsMat, eigvecs;
            bool isDiVehiclealid = dsd::diagonalize3d(sigmas[i], eigvalsMat, eigvecs);
            if (!isDiVehiclealid)
            {
                musVonMises[2 * i] = Vector3::Zero();
                kappas[2 * i] = 0;
                musVonMises[2 * i + 1] = Vector3::Zero();
                kappas[2 * i + 1] = kappas[2 * i];
                weightsVomp[2 * i] = 0;
                weightsVomp[2 * i + 1] = 0;
                continue;
            }
            auto eigvals = eigvalsMat.diagonal();
            Eigen::Index minIdx, maxIdx;
            double lmin = eigvals.minCoeff(&minIdx);
            double lmax = eigvals.maxCoeff(&maxIdx);
            if (fabs(lmin) > 1e-5 && fabs(lmax) > 1e-5 && std::isfinite(lmin) &&
                std::isfinite(lmax))
            {
                ROFL_ASSERT_VAR2((lmin / lmax) >= 0, lmin, lmax);
            }
            else
                continue;

            validDiagCtr++;
            ROFL_VAR1(validDiagCtr);
            musVonMises[2 * i] = eigvecs.col(minIdx);
            if (lmax > 10.0 * lmin)
            {
                kappas[2 * i] = 10.0;
            }
            else
            {
                kappas[2 * i] = lmax / lmin;
            }

            musVonMises[2 * i + 1] = -eigvecs.col(minIdx);
            kappas[2 * i + 1] = kappas[2 * i];
            weightsVomp[2 * i] = weights[i];
            weightsVomp[2 * i + 1] = weights[i];

            // ROFL_VAR3(i, phis[i], kappas[i]);
        }
        // ROFL_ASSERT(0)
        ROFL_VAR1(numGaussians)
        vonMisesMixture3d(vomp, nSamples, musVonMises, kappas, weightsVomp, szPadded);
    }

    void vompMax(std::vector<int> &maximaIdx,
                 std::vector<double> &maximaValues,
                 std::vector<std::pair<double, double>> &thetaPhiMaxima,
                 const std::vector<double> &vomp,
                 int nSamples,
                 double angleWin,
                 int findPeaksWindow)
    {
        // TODO: improve I/O (very inefficient atm)

        Eigen::MatrixXd vompM(2 * nSamples, nSamples);
        for (int itheta = 0; itheta < vompM.cols(); ++itheta)
            for (int iphi = 0; iphi < vompM.rows(); ++iphi)
            {
                int idx = itheta * nSamples + iphi;
                vompM(iphi, itheta) = vomp[idx];
            }

        int neighborhoodSz = 2; // l start
        int totGridSz = 2 * nSamples * nSamples;

        // return;

        Eigen::MatrixXi maximaMat(Eigen::MatrixXi::Ones(2 * nSamples, nSamples));
        int numMaxima = maximaMat.sum();

        bool multipleFpmMaxima = false;
        while (numMaxima > 0.1 * totGridSz && !multipleFpmMaxima)
        {
            ROFL_VAR2("Running fpm with l neighborhoodSz", neighborhoodSz)
            FPM fpm(vompM, neighborhoodSz); // TODO: growing window size
            fpm.run();
            multipleFpmMaxima = fpm.multipleMaxima_;
            maximaMat = fpm.getIndicesMax();

            neighborhoodSz++;
            numMaxima = maximaMat.sum();
        }

        // ROFL_VAR1(vompM)
        // ROFL_VAR1(maximaMat)

        maximaIdx.clear();
        maximaIdx.resize(vomp.size(), 0);

        double dtheta = M_PI / nSamples;

        for (int itheta = 0; itheta < vompM.cols(); ++itheta)
        {
            double theta = dtheta * itheta;

            for (int iphi = 0; iphi < vompM.rows(); ++iphi)
            {
                int idx = itheta * nSamples + iphi;
                auto mPhiTheta = maximaMat(iphi, itheta);
                maximaIdx[idx] = mPhiTheta;
                //
                if (mPhiTheta == 1)
                {
                    double phi = dtheta * iphi;
                    thetaPhiMaxima.push_back(std::pair<double, double>(theta, phi));
                    maximaValues.push_back(vomp[idx]);
                }
            }
        }
    }

    bool isFiniteVector3(const Vector3 &v)
    {
        if (!std::isfinite(v(0, 0)) ||
            !std::isfinite(v(1, 0) || !std::isfinite(v(2, 0))))
            return false;
        return true;
    }

    void plotSingleEllipseArrows3d(pcl::visualization::PCLVisualizer::Ptr viz,
                                   int idx,
                                   const Vector3 &mean,
                                   const Matrix3 &covar,
                                   const Matrix3 &eigvals,
                                   const Matrix3 &eigvecs,
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

        Eigen::Isometry3d transfIso;
        Eigen::AffineCompact3d translMat =
            Eigen::Affine3d(Eigen::Translation3d(mean(0), mean(1), mean(2)));
        Eigen::AffineCompact3d transfMat = translMat;
        transfIso = transfMat.matrix();

        double radiusX =
            0.5 * sqrt(eigvals.diagonal()(0)); // TODO: 2 * sqrt() actually needed?
        double radiusY = 0.5 * sqrt(eigvals.diagonal()(1));
        double radiusZ = 0.5 * sqrt(eigvals.diagonal()(2));

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
                pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 1.0, idEll);
            viz->setShapeRenderingProperties(
                pcl::visualization::PCL_VISUALIZER_OPACITY, 0.2, idEll);

            // std::cout << "Added ellipsoid " << idEll
            //           << " corresponding to:" << std::endl;
            // std::cout << "set object " << idx
            //           << " ellipse center " << mean(0) << ", " << mean(1) << ", "
            //           << mean(2)
            //           << " size " << radiusX << ", " << radiusY << ", " <<
            //           radiusZ << std::endl;

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

        std::vector<double> eigvalsVector;
        eigvalsVector.push_back(eigvals.diagonal()(0));
        eigvalsVector.push_back(eigvals.diagonal()(1));
        eigvalsVector.push_back(eigvals.diagonal()(2));
        std::sort(std::begin(eigvalsVector), std::end(eigvalsVector),
                  std::less<double>{}); // sorting in descending order

        std::vector<double> redPalette{1.0f, 0.0f, 0.0f};
        std::vector<double> greenPalette{0.0f, 1.0f, 0.0f};
        std::vector<double> bluePalette{0.0f, 0.0f, 1.0f};

        for (int i = 0; i < 3; ++i)
        {
            const std::string idArrRightI =
                "arrow_right_" + boost::lexical_cast<std::string, int>(idx) + "_" +
                boost::lexical_cast<std::string, int>(i);
            pcl::PointXYZ p1r(pcl::PointXYZ(mean(0), mean(1), mean(2)));
            pcl::PointXYZ p2r(
                pcl::PointXYZ(mean(0) + 0.001 * weight * eigvecs(i, 0),
                              mean(1) + 0.001 * weight * eigvecs(i, 1),
                              mean(2) + 0.001 * weight * eigvecs(i, 2)));
            if (!viz->addLine(p1r, p2r, redPalette[i], greenPalette[i],
                              bluePalette[i], idArrRightI))
            {
                std::cerr << "Failed adding following arrow: " << idArrRightI
                          << std::endl;
            }
            else
            {
                // viz->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
                // 1.0, 0.0, 0.0, idArrRight);
                viz->setShapeRenderingProperties(
                    pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, idArrRightI);

                // std::cout << "Added arrow " << idArrRight
                //           << " corresponding to: " << std::endl;
                // ROFL_VAR1(weight);
                // std::cout << "set object " << idx << " arrow center " << mean(0)
                // <<
                // ", " << mean(1)
                //           << " size " << radiusX << ", " << radiusY << " angle "
                //           << (180.0 / M_PI * angle)
                //           << " front fs empty bo 3\n";
            }

            const std::string idArrLeftI =
                "arrow_left_" + boost::lexical_cast<std::string, int>(idx) + "_" +
                boost::lexical_cast<std::string, int>(i);
            pcl::PointXYZ p1l(pcl::PointXYZ(mean(0), mean(1), mean(2)));
            pcl::PointXYZ p2l(
                pcl::PointXYZ(mean(0) - 0.001 * weight * eigvecs(i, 0),
                              mean(1) - 0.001 * weight * eigvecs(i, 1),
                              mean(2) - 0.001 * weight * eigvecs(i, 2)));
            if (!viz->addLine(p1l, p2l, redPalette[i], greenPalette[i],
                              bluePalette[i], idArrLeftI))
            {
                std::cerr << "Failed adding following arrow: " << idArrLeftI
                          << std::endl;
            }
            else
            {
                // viz->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
                // 1.0, 0.0, 0.0, idArrRight);
                viz->setShapeRenderingProperties(
                    pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, idArrLeftI);

                // std::cout << "Added arrow " << idArrLeft
                //           << " corresponding to: " << std::endl;
                // ROFL_VAR1(weight);
                // std::cout << "set object " << idx << " arrow center " << mean(0)
                // <<
                // ", " << mean(1)
                //           << " size " << radiusX << ", " << radiusY << " angle "
                //           << (180.0 / M_PI * angle)
                //           << " front fs empty bo 3\n";
            }
        }
    }

    void computeGaussianAxes3d(VectorMatrix3 &eigvals,
                               VectorMatrix3 &eigvecs,
                               const VectorMatrix3 &covar,
                               const std::vector<double> &weights)
    {
        ROFL_VAR1("Running computeGaussianAxes3d()")
        eigvals.clear();
        eigvecs.clear();
        size_t sz = covar.size();
        for (size_t i = 0; i < sz; ++i)
        {
            ROFL_VAR1(i)
            const Matrix3 sigma = covar[i];
            Matrix3 eigval, eigvec;
            if (diagonalize3d(sigma, eigval, eigvec))
            {
                eigvals.push_back(eigval);
                eigvecs.push_back(eigvec);
            }
            // TODO: check how this affects future sizes
        }
        ROFL_VAR1("End of computeGaussianAxes3d()")
    }

    void plotEllipsesArrows3d(pcl::visualization::PCLVisualizer::Ptr viz,
                              const dsd::VectorVector3 &means,
                              const dsd::VectorMatrix3 &covars,
                              const std::vector<double> &weights,
                              dsd::SaveCsvOptions saveOutput)
    {
        VectorMatrix3 eigvals, eigvecs;
        computeGaussianAxes3d(eigvals, eigvecs, covars, weights);
        // viz->removeAllShapes();
        for (int i = 0; i < means.size() && i < covars.size(); ++i)
        {
            ROFL_VAR3("plotSingleEllipseArrows3d", i, weights[i])
            if (weights[i] > 8.0)
                plotSingleEllipseArrows3d(viz, i, means[i], covars[i], eigvals[i],
                                          eigvecs[i], weights[i], saveOutput);
        }
    }

} // namespace dsd
