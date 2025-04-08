/**
 * ARS - Angular Radon Spectrum
 * Copyright (C) 2017 Dario Lodi Rizzini.
 *           (C) 2021 Dario Lodi Rizzini, Ernesto Fontana.
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
#include <gme_GaussianMixtureEstimator.h>
#include <signal.h>
#include <deque>
#include <iterator>
#include <vector>

#include <dsd_utils.h>

namespace gme
{

    //-----------------------------------------------------
    // GaussianMixtureEstimator
    //-----------------------------------------------------

    GaussianMixtureEstimator::GaussianMixtureEstimator()
        : gaussians_()
    { //: means_(), covROFL_(), weights_() {
    }

    GaussianMixtureEstimator::~GaussianMixtureEstimator() {}

    void GaussianMixtureEstimator::clearGaussians()
    {
        if (gaussians_.size() > 0)
            gaussians_.clear();
    }

    size_t GaussianMixtureEstimator::size() const
    {
        // return weights_.size();
        return gaussians_.size();
    }

    const dsd::Vector2 &GaussianMixtureEstimator::mean(int i) const
    {
        //            ROFL_ASSERT(0 <= i && i < means_.size());
        //            return means_.at(i);
        ROFL_ASSERT(0 <= i && i < gaussians_.size());
        return gaussians_[i].mean;
    }

    const dsd::Matrix2 &GaussianMixtureEstimator::covariance(int i) const
    {
        //            ROFL_ASSERT(0 <= i && i < covROFL_.size());
        //            return covROFL_.at(i);
        ROFL_ASSERT(0 <= i && i < gaussians_.size());
        return gaussians_[i].covar;
    }

    double GaussianMixtureEstimator::weight(int i) const
    {
        //            ROFL_ASSERT(0 <= i && i < weights_.size());
        //            return weights_.at(i);
        ROFL_ASSERT(0 <= i && i < gaussians_.size());
        return gaussians_[i].weight;
    }

    const GaussianMixtureEstimator::VectorGaussian &
    GaussianMixtureEstimator::gaussians() const
    {
        return gaussians_;
    }

    void GaussianMixtureEstimator::exportGaussians(
        dsd::VectorVector2 &means,
        dsd::VectorMatrix2 &covariances,
        std::vector<double> &weights) const
    {
        int n = gaussians_.size();
        means.clear();
        covariances.clear();
        weights.clear();
        for (int i = 0; i < n; ++i)
        {
            // if (!checkIsGaussianNull(gaussians_[i].covar, gaussians_[i].weight))
            // {
            means.push_back(gaussians_[i].mean);
            covariances.push_back(gaussians_[i].covar);
            weights.push_back(gaussians_[i].weight);
            //}
        }
    }

    void GaussianMixtureEstimator::exportGaussians3d(
        dsd::VectorVector3 &means,
        dsd::VectorMatrix3 &covariances,
        std::vector<double> &weights) const
    {
        int n = gaussians3_.size();
        means.clear();
        covariances.clear();
        weights.clear();
        for (int i = 0; i < n; ++i)
        {
            if (!checkIsGaussianNull(gaussians3_[i].covar, gaussians3_[i].weight) &&
                i % 100 == 0)
            {
                means.push_back(gaussians3_[i].mean);
                covariances.push_back(gaussians3_[i].covar);
                weights.push_back(gaussians3_[i].weight);
            }
        }
    }

    void GaussianMixtureEstimator::executeEM(const dsd::VectorVector2 &samples,
                                             int stepNum)
    {
        Eigen::MatrixXd c(samples.size(), gaussians_.size());
        // Eigen::VectorXd sumC(gaussians_.size());
        double sumC;

        for (int step = 0; step < stepNum; ++step)
        {
            // Expectation
            for (int i = 0; i < samples.size(); ++i)
            {
                sumC = 0.0;
                for (int j = 0; j < gaussians_.size(); ++j)
                {
                    c(i, j) = gaussians_[j].weight * gaussians_[j].eval(samples[i]);
                    sumC += c(i, j);
                }
                for (int j = 0; j < gaussians_.size(); ++j)
                {
                    c(i, j) = c(i, j) / sumC;
                }
            }

            // Maximization
            for (int j = 0; j < gaussians_.size(); ++j)
            {
                gaussians_[j].mean = dsd::Vector2::Zero();
                gaussians_[j].covar = dsd::Matrix2::Zero();
                gaussians_[j].weight = 0.0;
                sumC = 0.0;

                // Computes the mean value
                for (int i = 0; i < samples.size(); ++i)
                {
                    gaussians_[j].mean += samples[i] * c(i, j);
                    sumC += c(i, j);
                }
                gaussians_[j].mean = gaussians_[j].mean / sumC;
                gaussians_[j].weight = sumC / samples.size();

                // Computes the covariance
                for (int i = 0; i < samples.size(); ++i)
                {
                    gaussians_[j].covar +=
                        (samples[i] - gaussians_[j].mean) *
                        (samples[i] - gaussians_[j].mean).transpose() * c(i, j);
                }
                gaussians_[j].covar = gaussians_[j].covar / sumC;
            }
        }
    }

    //-----------------------------------------------------
    // GaussianMixtureEstimatorScan
    //-----------------------------------------------------

    GaussianMixtureEstimatorScan::GaussianMixtureEstimatorScan()
        : GaussianMixtureEstimator(),
          intervals_(),
          distanceGap_(0.8),
          distanceSplit_(0.5),
          sigmaMin_(0.1) {}

    GaussianMixtureEstimatorScan::~GaussianMixtureEstimatorScan() {}

    void GaussianMixtureEstimatorScan::compute(const dsd::VectorVector2 &samples)
    {
        dsd::Vector2 mean;
        dsd::Matrix2 covar;
        std::deque<IndexInterval> intervals;
        IndexInterval interv, interv1, interv2;
        double dist, distMax, w;
        int farthest;

        int sum = 0;

        // Splits the scan points into intervals when a gap between consecutive
        // points is found
        interv.first = 0;
        for (int i = 1; i < samples.size(); ++i)
        {
            dist = (samples[i] - samples[i - 1]).norm();
            if (dist > distanceGap_)
            {
                interv.last = i - 1;
                interv.num = interv.last - interv.first + 1;
                intervals.push_back(interv);
                // ROFL_VAR1("interv [" << interv.first << ", " << interv.last << "]
                // num " << interv.num);
                interv.first = i;
            }
        }
        interv.last = samples.size() - 1;
        interv.num = interv.last - interv.first + 1;
        intervals.push_back(interv);

        std::cout << "\n----\n"
                  << std::endl;

        // Searches for aligned points in interval
        while (!intervals.empty())
        {
            // Extracts the first interval
            interv = intervals.front();
            intervals.pop_front();
            // Checks if the interval is split according to a policy based on
            // distance of farthest point from segment
            findFarthest(samples, interv.first, interv.last, farthest, distMax);
            if (distMax > distanceSplit_)
            {
                // Interval is split at farthest point. Formally:
                // - interv1: [interv.first, farthest-1] (farthest NOT included**)
                // - interv2: [farthest, interv.last]
                // ** the fathest is not included in the first interval, but it's
                // used
                //    in the computation of the gaussian!!! So we have an overlap:
                // - interv1: [interv.first, farthest] (farthest NOT included**)
                interv1.first = interv.first;
                interv1.last = farthest;
                interv1.num = farthest - interv.first;
                interv2.first = farthest;
                interv2.last = interv.last;
                interv2.num = interv.num - interv1.num;
                intervals.push_front(interv1);
                intervals.push_front(interv2);
            }
            else
            {
                // ROFL_VAR1("interv [" << interv.first << ", " << interv.last << "]
                // num " << interv.num);
                Gaussian g;
                // estimateGaussianFromPoints(samples, interv.first, interv.last,
                // g.mean, g.covar);
                estimateGaussianFromSegment(samples, interv.first, interv.last,
                                            g.mean, g.covar);
                w = interv.num * 1.0 / samples.size();
                g.weight = w;
                gaussians_.push_front(g);
                // means_.push_back(mean);
                // covROFL_.push_back(covar);
                // weights_.push_back(w);
                intervals_.push_front(interv);
                sum += interv.num;
            }
        }

        ROFL_VAR1(sum);
    }

    void GaussianMixtureEstimatorScan::findFarthest(
        const dsd::VectorVector2 &points,
        int first,
        int last,
        int &farthest,
        double &distMax) const
    {
        dsd::Vector2 dp;
        double t, ct, st, r, dist;

        ROFL_ASSERT(0 <= first && last < points.size() && first <= last);
        dp = points[last] - points[first];
        t = atan2(dp(0), -dp(1));
        ct = cos(t);
        st = sin(t);
        r = points[first](0) * ct + points[first](1) * st;

        distMax = 0.0;
        for (int i = first; i <= last; ++i)
        {
            dist = fabs(points[i](0) * ct + points[i](1) * st - r);
            if (dist > distMax)
            {
                distMax = dist;
                farthest = i;
            }
        }
    }

    void GaussianMixtureEstimatorScan::estimateGaussianFromPoints(
        const dsd::VectorVector2 &points,
        int first,
        int last,
        dsd::Vector2 &mean,
        dsd::Matrix2 &covar) const
    {
        dsd::Matrix2 l, v;
        dsd::Vector2 tmp;
        double sigmaMinSquare = sigmaMin_ * sigmaMin_;

        ROFL_ASSERT(first >= 0 && last < points.size());

        // Computes the mean value vector
        mean = dsd::Vector2::Zero();
        for (int i = first; i <= last; ++i)
        {
            mean += points[i];
        }
        mean = mean / (last - first + 1);

        // Computes the covariance
        covar = dsd::Matrix2::Zero();
        for (int i = first; i <= last; ++i)
        {
            tmp = (points[i] - mean);
            covar = tmp * tmp.transpose();
        }

        if (first == last)
        {
            // Only one point: use the point uncertainty
            covar << sigmaMinSquare, 0.0, 0.0, sigmaMinSquare;
        }
        else
        {
            covar = covar / (last - first);
            dsd::diagonalize(covar, l, v);
            if (l(0, 0) < sigmaMinSquare)
                l(0, 0) = sigmaMinSquare;
            if (l(1, 1) < sigmaMinSquare)
                l(1, 1) = sigmaMinSquare;
            covar = v * l * v.transpose();
        }
    }

    void GaussianMixtureEstimatorScan::estimateGaussianFromSegment(
        const dsd::VectorVector2 &points,
        int first,
        int last,
        dsd::Vector2 &mean,
        dsd::Matrix2 &covar) const
    {
        dsd::Matrix2 v;
        dsd::Vector2 tmp;
        double lmin, lmax, theta, dfirst, dlast;
        double sigmaMinSquare = sigmaMin_ * sigmaMin_;

        ROFL_ASSERT(first >= 0 && last < points.size());

        // Computes the mean value vector
        mean = dsd::Vector2::Zero();
        for (int i = first; i <= last; ++i)
        {
            mean += points[i];
        }
        mean = mean / (last - first + 1);

        // Computes the covariance
        covar = dsd::Matrix2::Zero();
        for (int i = first; i <= last; ++i)
        {
            tmp = (points[i] - mean);
            covar = tmp * tmp.transpose();
        }

        if (first == last)
        {
            // Only one point: use the point uncertainty
            covar << sigmaMinSquare, 0.0, 0.0, sigmaMinSquare;
        }
        else
        {
            covar = covar / (last - first);
            dsd::diagonalize(covar, lmin, lmax, theta);
            dfirst = cos(theta) * points[first](0) + sin(theta) * points[first](1);
            dlast = cos(theta) * points[last](0) + sin(theta) * points[last](1);
            lmax = 0.2 * (dlast - dfirst) * (dlast - dfirst);
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
            //            ROFL_VAR1("[" << first << "," << last << "]: lmin " <<
            //            lmin << ", lmax " << lmax << ", theta[deg] " << (180.0 /
            //            M_PI * theta)
            //                    << "\ncovar\n" << covar);
            //            dsd::diagonalize(covar, l, v);
            //            if (l(0, 0) < sigmaMinSquare)
            //                l(0, 0) = sigmaMinSquare;
            //            if (l(1, 1) < sigmaMinSquare)
            //                l(1, 1) = sigmaMinSquare;
            //            covar = v * l * v.transpose();
        }
    }

    //-----------------------------------------------------
    // GaussianMixtureEstimatorHierarchical
    //-----------------------------------------------------

    GaussianMixtureEstimatorHierarchical::GaussianMixtureEstimatorHierarchical()
        : data_(),
          sigmaMin_(0.1),
          covarWidth_(0.2),
          chi2Thres_(5.99146),
          iseThres_(0.3),
          inlierPerc_(0.60),
          levelMax_(32) {}

    GaussianMixtureEstimatorHierarchical::~GaussianMixtureEstimatorHierarchical() {}

    double GaussianMixtureEstimatorHierarchical::getSigmaMin() const
    {
        return sigmaMin_;
    }

    double GaussianMixtureEstimatorHierarchical::getIseThreshold() const
    {
        return iseThres_;
    }

    void GaussianMixtureEstimatorHierarchical::setSigmaMin(double sm)
    {
        sigmaMin_ = sm;
        data_.setRes(sigmaMin_);
    }

    void GaussianMixtureEstimatorHierarchical::setCovarWidth(double cw)
    {
        covarWidth_ = cw;
    }

    void GaussianMixtureEstimatorHierarchical::setInlierPerc(double ip)
    {
        inlierPerc_ = ip;
        if (inlierPerc_ < 0.0)
            inlierPerc_ = 0.0;
        else if (inlierPerc_ >= 1.0)
            inlierPerc_ = 1.0;
    }

    void GaussianMixtureEstimatorHierarchical::setChiConfidence(double conf)
    {
        static const int CHI2_DOF = 2;
        ROFL_ASSERT(0.0 <= conf && conf <= 1.0);
        boost::math::chi_squared mydist(CHI2_DOF);
        chi2Thres_ = boost::math::quantile(mydist, conf);
    }

    void GaussianMixtureEstimatorHierarchical::setIseThreshold(double iseTh)
    {
        iseThres_ = iseTh;
    }

    void GaussianMixtureEstimatorHierarchical::setCellSizeMax(double s)
    {
        levelMax_ = log2Mod((int)ceil(s / sigmaMin_));
    }

    void GaussianMixtureEstimatorHierarchical::initIsotropic(
        const dsd::VectorVector3 &samples)
    {
        for (int i = 0; i < samples.size(); ++i)
        {
            Gaussian3d g(samples[i], Eigen::Matrix3d::Identity(), 1.0f);

            gaussians3_.push_back(g);
        }
    }

    void GaussianMixtureEstimatorHierarchical::compute(
        const dsd::VectorVector2 &samples)
    {
        std::deque<ConstInterval> intervals;
        ConstInterval intervalCur;
        ConstIterator mid;
        Gaussian g;
        double num;
        int level;

        num = (double)samples.size();
        if (num == 0)
        {
            return;
        }

        gaussians_.clear();
        data_.insert(samples);
        intervals.push_back(std::make_pair(std::begin(data_), std::end(data_)));
        while (!intervals.empty())
        {
            intervalCur = intervals.front();
            intervals.pop_front();
            level = data_.computeLevel(intervalCur.first, intervalCur.second);
            if (intervalCur.second != intervalCur.first)
            {
                //                ROFL_VAR4(intervalCur.first->index.transpose(),
                //                (intervalCur.second - 1)->index.transpose(),
                //                level, levelMax_);
                //                ROFL_VAR1(std::distance(intervalCur.first,
                //                intervalCur.second)); for (auto it =
                //                intervalCur.first; it != intervalCur.second; ++it)
                //                {
                //                  std::cout << "  [" << it->value.transpose() <<
                //                  "] -> index [" << it->index.transpose() << "]
                //                  pos " << it->pos << "\n";
                //                }
            }
            if (level <= levelMax_ &&
                (level == 0 ||
                 estimateGaussianISE(intervalCur.first, intervalCur.second, g.mean,
                                     g.covar, g.weight)))
            {
                //                ROFL_VAR1("insert gaussian in [" <<
                //                g.mean.transpose() << "]");
                gaussians_.push_front(g);
            }
            else
            {
                mid = data_.findSplit(intervalCur.first, intervalCur.second);
                //                ROFL_VAR1("splitting into\n" << "  ([" <<
                //                intervalCur.first->index.transpose() << "], [" <<
                //                mid->index.transpose() << "])\n"
                //                   << "  ([" << mid->index.transpose() << "], ["
                //                   << (intervalCur.second-1)->index.transpose() <<
                //                   "])");
                if (mid != intervalCur.first && mid != intervalCur.second)
                {
                    ROFL_ASSERT(std::distance(intervalCur.first, mid) >= 0);
                    intervals.push_back(std::make_pair(intervalCur.first, mid));
                    ROFL_ASSERT(std::distance(mid, intervalCur.second) >= 0);
                    intervals.push_back(std::make_pair(mid, intervalCur.second));
                }
            }
        }
    }

    bool GaussianMixtureEstimator::checkIsGaussianNull(const dsd::Matrix2 &covar,
                                                       const double &w) const
    {
        if (!dsd::isEqualFloats(w, 0.0f))
            return false;
        if (!dsd::isEqualFloats(covar.norm(), 0.0f))
            return false;

        return true;
    }

    bool GaussianMixtureEstimator::checkIsGaussianNull(const dsd::Matrix3 &covar,
                                                       const double &w) const
    {
        if (!dsd::isEqualFloats(w, 0.0f))
            return false;
        if (!dsd::isEqualFloats(covar.norm(), 0.0f))
            return false;

        return true;
    }

    bool GaussianMixtureEstimatorHierarchical::estimateGaussianFromPoints(
        const ConstIterator &beg,
        const ConstIterator &end,
        dsd::Vector2 &mean,
        dsd::Matrix2 &covar,
        double &w) const
    {
        dsd::Matrix2 l, v, infoMat;
        dsd::Vector2 tmp;
        double sigmaMinSquare = sigmaMin_ * sigmaMin_;
        double distSqr;
        int num, inlier;

        // Computes the mean value vector
        mean = dsd::Vector2::Zero();
        num = 0;
        for (auto it = beg; it != end; ++it)
        {
            mean += it->value;
            num++;
        }
        mean = mean / num;

        //	ROFL_VARIABLE2(num, mean.transpose());

        // Computes the covariance
        covar = dsd::Matrix2::Zero();
        for (auto it = beg; it != end; ++it)
        {
            tmp = (it->value - mean);
            covar += tmp * tmp.transpose();
        }

        if (num <= 1)
        {
            // Only one point: use the point uncertainty
            covar << sigmaMinSquare, 0.0, 0.0, sigmaMinSquare;
        }
        else
        {
            covar = covar / (num - 1);
            dsd::diagonalize(covar, l, v);
            if (l(0, 0) < sigmaMinSquare)
                l(0, 0) = sigmaMinSquare;
            if (l(1, 1) < sigmaMinSquare)
                l(1, 1) = sigmaMinSquare;
            covar = v * l * v.transpose();
        }

        inlier = 0;
        infoMat = covar.inverse();
        // ROFL_VAR1("covar\n" << covar << "\neigenvalues:\n" << l.transpose() <<
        // "\ninfoMat\n" << infoMat);
        for (auto it = beg; it != end; ++it)
        {
            tmp = (it->value - mean);
            distSqr = tmp.transpose() * infoMat * tmp;
            if (distSqr < chi2Thres_)
            {
                inlier++;
            }
        }
        //	ROFL_VARIABLE2(inlier, chi2Thres_);
        w = num / data_.size();
        return (inlier >= inlierPerc_ * num);
    }

    bool GaussianMixtureEstimatorHierarchical::estimateGaussianFromSegment(
        const ConstIterator &beg,
        const ConstIterator &end,
        dsd::Vector2 &mean,
        dsd::Matrix2 &covar,
        double &w) const
    {
        dsd::Matrix2 l, v, infoMat;
        dsd::Vector2 tmp;
        double sigmaMinSquare = sigmaMin_ * sigmaMin_;
        int num, inlier;
        double lmin, lmax, theta, ct, st, distSqr, d, dfirst, dlast;

        // Computes the mean value vector
        mean = dsd::Vector2::Zero();
        num = 0;
        for (auto it = beg; it != end; ++it)
        {
            mean += it->value;
            num++;
        }
        mean = mean / num;

        // Computes the covariance
        covar = dsd::Matrix2::Zero();
        for (auto it = beg; it != end; ++it)
        {
            tmp = (it->value - mean);
            covar += tmp * tmp.transpose();
        }

        if (num <= 1)
        {
            // Only one point: use the point uncertainty
            covar << sigmaMinSquare, 0.0, 0.0, sigmaMinSquare;
        }
        else
        {
            covar = covar / (num - 1);
            dsd::diagonalize(covar, lmin, lmax, theta);
            ct = cos(theta);
            st = sin(theta);
            dfirst = 1e+6;
            dlast = -1e+6;
            for (auto it = beg; it != end; ++it)
            {
                d = ct * it->value(0) + st * it->value(1);
                if (d < dfirst)
                    dfirst = d;
                if (d > dlast)
                    dlast = d;
            }
            lmax = covarWidth_ * (dlast - dfirst) * (dlast - dfirst);
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

        inlier = 0;
        infoMat = covar.inverse();
        // ROFL_VAR1("covar\n" << covar << "\neigenvalues:\n" << l.transpose() <<
        // "\ninfoMat\n" << infoMat);
        for (auto it = beg; it != end; ++it)
        {
            tmp = (it->value - mean);
            distSqr = tmp.transpose() * infoMat * tmp;
            if (distSqr < chi2Thres_)
            {
                inlier++;
            }
        }
        //	ROFL_VARIABLE2(inlier, chi2Thres_);
        w = num / data_.size();

        return (inlier >= inlierPerc_ * num);
    }

    bool GaussianMixtureEstimatorHierarchical::estimateGaussianISE(
        const ConstIterator &beg,
        const ConstIterator &end,
        dsd::Vector2 &mean,
        dsd::Matrix2 &covar,
        double &wMerged) const
    {
        ROFL_ASSERT(beg != end);

        double sigmaMinSqrd = sigmaMin_ * sigmaMin_;
        double wOrig = 1.0 / data_.size(); // w_i

        int num = 0;
        mean = dsd::Vector2::Zero();
        for (auto it = beg; it != end; ++it)
        {
            mean += it->value;
            num++;
        }
        ROFL_ASSERT(num > 0);
        mean /= num;

        // std::cout << "num " << num << std::endl;
        // std::cout << "mean " << mean << std::endl;

        wMerged = wOrig * num; // w_L

        covar = dsd::Matrix2::Zero();
        for (auto it = beg; it != end; ++it)
        {
            dsd::Vector2 tmp = it->value - mean;
            covar += tmp * tmp.transpose();
        }
        covar = covar / num;
        covar(0, 0) += sigmaMinSqrd;
        covar(1, 1) += sigmaMinSqrd;

        // std::cout << "covar\n" << covar << std::endl;

        double jNN = 0.0, jLL = 0.0, jNL = 0.0;

        jLL = wMerged * wMerged / (2 * M_PI * sqrt((2 * covar).determinant()));

        double normOrig = 1.0 / (4 * M_PI * sigmaMinSqrd); // normNN
        double gaussConstOrig = -0.25 / sigmaMinSqrd;      // const2

        dsd::Matrix2 covarNL = covar + sigmaMinSqrd * dsd::Matrix2::Identity();
        dsd::Matrix2 infoNL = covarNL.inverse();
        double normNL = 1.0 / (2 * M_PI * sqrt(covarNL.determinant()));

        for (auto it = beg; it != end; ++it)
        {
            jNN += 1.0;
            for (auto it2 = it + 1; it2 != end; ++it2)
            {
                jNN +=
                    2 * exp(gaussConstOrig *
                            (it->value - it2->value).dot(it->value - it2->value));
            }

            jNL += exp(-0.5 * (it->value - mean).transpose() * infoNL *
                       (it->value - mean));
        }
        jNN *= wOrig * wOrig * normOrig;

        jNL *= wOrig * wMerged * normNL;

        double nise = 1 - 2 * jNL / (jNN + jLL);

        if (num == 1)
        {
            // std::cout << "nise " << nise << std::endl;
            //             std::cout << "wMerged " << wMerged << std::endl;
            //                         std::cout << "wOrig " << wOrig << std::endl;
            // std::cout << "mean\n" << mean.transpose() << std::endl;
            // std::cout << "beg " << (beg->value).transpose() << std::endl;
            // std::cout << "normOrig " << normOrig << std::endl;
            // std::cout << "normNL " << normNL << std::endl;
            // std::cout << "jNN " << jNN << std::endl;
            // std::cout << "jNL " << jNL << std::endl;
            // std::cout << "jLL " << jLL << std::endl;
            return true;
            //            std::cin.get();
        }
        // ROFL_VAR4(num, mean.transpose(), nise, iseThres_);

        //        std::cout << "--------------\n\n\n" << std::endl;

        return nise < iseThres_;
    }

    //-----------------------------------------------------
    // GaussianMixtureEstimatorHierarchical
    //-----------------------------------------------------

    GaussianMixtureEstimatorHierarchical3d::GaussianMixtureEstimatorHierarchical3d()
        : data_(),
          sigmaMin_(0.1),
          covarWidth_(0.2),
          chi2Thres_(5.99146),
          iseThres_(0.3),
          inlierPerc_(0.60),
          levelMax_(32) {}

    GaussianMixtureEstimatorHierarchical3d::
        ~GaussianMixtureEstimatorHierarchical3d() {}

    double GaussianMixtureEstimatorHierarchical3d::getSigmaMin() const
    {
        return sigmaMin_;
    }

    double GaussianMixtureEstimatorHierarchical3d::getIseThreshold() const
    {
        return iseThres_;
    }

    void GaussianMixtureEstimatorHierarchical3d::setSigmaMin(double sm)
    {
        sigmaMin_ = sm;
        data_.setRes(sigmaMin_);
    }

    void GaussianMixtureEstimatorHierarchical3d::setCovarWidth(double cw)
    {
        covarWidth_ = cw;
    }

    void GaussianMixtureEstimatorHierarchical3d::setInlierPerc(double ip)
    {
        inlierPerc_ = ip;
        if (inlierPerc_ < 0.0)
            inlierPerc_ = 0.0;
        else if (inlierPerc_ >= 1.0)
            inlierPerc_ = 1.0;
    }

    void GaussianMixtureEstimatorHierarchical3d::setChiConfidence(double conf)
    {
        static const int CHI2_DOF = 2;
        ROFL_ASSERT(0.0 <= conf && conf <= 1.0);
        boost::math::chi_squared mydist(CHI2_DOF);
        chi2Thres_ = boost::math::quantile(mydist, conf);
    }

    void GaussianMixtureEstimatorHierarchical3d::setIseThreshold(double iseTh)
    {
        iseThres_ = iseTh;
    }

    void GaussianMixtureEstimatorHierarchical3d::setCellSizeMax(double s)
    {
        levelMax_ = log2Mod((int)ceil(s / sigmaMin_));
    }

    const GaussianMixtureEstimatorHierarchical3d::VectorGaussian &
    GaussianMixtureEstimatorHierarchical3d::gaussians() const
    {
        return gaussians_;
    }

    void GaussianMixtureEstimatorHierarchical3d::exportGaussians(
        dsd::VectorVector3 &means,
        dsd::VectorMatrix3 &covariances,
        std::vector<double> &weights) const
    {
        int n = gaussians_.size();
        means.clear();
        covariances.clear();
        weights.clear();
        for (int i = 0; i < n; ++i)
        {
            if (!checkIsGaussianNull(gaussians_[i].covar, gaussians_[i].weight))
            {
                means.push_back(gaussians_[i].mean);
                covariances.push_back(gaussians_[i].covar);
                weights.push_back(gaussians_[i].weight);
            }
        }
    }

    void GaussianMixtureEstimatorHierarchical3d::compute(
        const dsd::VectorVector3 &samples)
    {
        std::deque<ConstInterval> intervals;
        ConstInterval intervalCur;
        ConstIterator mid;
        Gaussian g;
        double num;
        int level;

        num = (double)samples.size();
        if (num == 0)
        {
            return;
        }

        gaussians_.clear();
        data_.insert(samples);
        ROFL_VAR1(data_.size());
        intervals.push_back(std::make_pair(std::begin(data_), std::end(data_)));
        while (!intervals.empty())
        {
            intervalCur = intervals.front();
            intervals.pop_front();
            level = data_.computeLevel(intervalCur.first, intervalCur.second);
            // ROFL_VAR4(intervalCur.first - data_.begin(),
            //           intervalCur.second - data_.begin(), level, levelMax_);
            if (level <= levelMax_ &&
                (level == 0 ||
                 estimateGaussianISE(intervalCur.first, intervalCur.second, g.mean,
                                     g.covar, g.weight)))
            {
                // ROFL_MSG("insert gaussian in [" << g.mean.transpose() << "]");
                gaussians_.push_front(g);
            }
            else
            {
                mid = data_.findSplit(intervalCur.first, intervalCur.second);
                // ROFL_MSG("splitting into\n"
                //          << "  ([" << intervalCur.first->index.transpose() << "],
                //          ["
                //          << mid->index.transpose() << "])\n"
                //          << "  ([" << mid->index.transpose() << "], ["
                //          << (intervalCur.second - 1)->index.transpose() << "])");
                if (mid != intervalCur.first && mid != intervalCur.second)
                {
                    ROFL_ASSERT(std::distance(intervalCur.first, mid) >= 0);
                    intervals.push_back(std::make_pair(intervalCur.first, mid));
                    ROFL_ASSERT(std::distance(mid, intervalCur.second) >= 0);
                    intervals.push_back(std::make_pair(mid, intervalCur.second));
                }
            }
        }
    }

    bool GaussianMixtureEstimatorHierarchical3d::estimateGaussianISE(
        const ConstIterator &beg,
        const ConstIterator &end,
        dsd::Vector3 &mean,
        dsd::Matrix3 &covar,
        double &wMerged) const
    {
        ROFL_ASSERT(beg != end);

        double sigmaMinSqrd = sigmaMin_ * sigmaMin_;
        double wOrig = 1.0 / data_.size(); // w_i

        int num = 0;
        mean = dsd::Vector3::Zero();
        for (auto it = beg; it != end; ++it)
        {
            mean += it->value;
            num++;
        }
        ROFL_ASSERT(num > 0);
        mean /= num;

        // std::cout << "num " << num << std::endl;
        // std::cout << "mean " << mean << std::endl;

        wMerged = wOrig * num; // w_L

        covar = dsd::Matrix3::Zero();
        for (auto it = beg; it != end; ++it)
        {
            dsd::Vector3 tmp = it->value - mean;
            covar += tmp * tmp.transpose();
        }
        covar = covar / num;
        for (int i = 0; i < 3; ++i)
        {
            covar(i, i) += sigmaMinSqrd;
        }

        // std::cout << "covar\n" << covar << std::endl;

        double jNN = 0.0, jLL = 0.0, jNL = 0.0;

        jLL = wMerged * wMerged / (2 * M_PI * sqrt((2 * covar).determinant()));

        double normOrig = 1.0 / (4 * M_PI * sigmaMinSqrd); // normNN
        double gaussConstOrig = -0.25 / sigmaMinSqrd;      // const2

        dsd::Matrix3 covarNL = covar + sigmaMinSqrd * dsd::Matrix3::Identity();
        dsd::Matrix3 infoNL = covarNL.inverse();
        double normNL = 1.0 / (2 * M_PI * sqrt(covarNL.determinant()));

        for (auto it = beg; it != end; ++it)
        {
            jNN += 1.0;
            for (auto it2 = it + 1; it2 != end; ++it2)
            {
                jNN +=
                    2 * exp(gaussConstOrig *
                            (it->value - it2->value).dot(it->value - it2->value));
                // ROFL_VAR4(jNN, it - beg, it2 - beg, end - beg);
            }

            jNL += exp(-0.5 * (it->value - mean).transpose() * infoNL *
                       (it->value - mean));
        }
        jNN *= wOrig * wOrig * normOrig;

        jNL *= wOrig * wMerged * normNL;

        double nise = 1 - 2 * jNL / (jNN + jLL);

        if (num == 1)
        {
            // std::cout << "nise " << nise << std::endl;
            //             std::cout << "wMerged " << wMerged << std::endl;
            //                         std::cout << "wOrig " << wOrig << std::endl;
            // std::cout << "mean\n" << mean.transpose() << std::endl;
            // std::cout << "beg " << (beg->value).transpose() << std::endl;
            // std::cout << "normOrig " << normOrig << std::endl;
            // std::cout << "normNL " << normNL << std::endl;
            // std::cout << "jNN " << jNN << std::endl;
            // std::cout << "jNL " << jNL << std::endl;
            // std::cout << "jLL " << jLL << std::endl;
            return true;
            //            std::cin.get();
        }
        // ROFL_VAR4(num, mean.transpose(), nise, iseThres_);

        //        std::cout << "--------------\n\n\n" << std::endl;

        return nise < iseThres_;
    }

    bool GaussianMixtureEstimatorHierarchical3d::checkIsGaussianNull(
        const dsd::Matrix3 &covar,
        const double &w) const
    {
        if (!dsd::isEqualFloats(w, 0.0f))
            return false;
        if (!dsd::isEqualFloats(covar.norm(), 0.0f))
            return false;

        return true;
    }

} // namespace gme
