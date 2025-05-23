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
#ifndef GAUSSIANMIXTUREESTIMATOR_H
#define GAUSSIANMIXTUREESTIMATOR_H

#include <Eigen/Dense>
#include <boost/math/distributions/chi_squared.hpp>
#include <deque>
#include <iostream>
#include <vector>

#include <dsd_utils.h>
#include <gme_MortonOctree.h>
#include <rofl/common/macros.h>

namespace gme
{

    //-----------------------------------------------------
    // GaussianMixtureEstimator
    //-----------------------------------------------------

    /**
     * Class GaussianMixtureEstimator provides a general interface for the
     * estimators of Gaussian Mixture Models (GMMs) that compute the Gaussian
     * parameters (mean vectors, covariance matrices, weights, etc.) from observed
     * samples.
     */
    class GaussianMixtureEstimator
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        struct Gaussian
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            dsd::Vector2 mean;
            dsd::Matrix2 covar;
            double weight;

            Gaussian() {}

            Gaussian(const dsd::Vector2 &m, const dsd::Matrix2 &c, double w)
            {
                mean = m;
                covar = c;
                weight = w;
            }

            ~Gaussian() {}

            double eval(const dsd::Vector2 &v) const
            {
                double k = 1.0 / sqrt(2.0 * M_PI * covar.determinant());
                double arg = (v - mean).transpose() * covar.inverse() * (v - mean);
                return k * exp(-0.5 * arg);
            }
        };

        struct Gaussian3d
        {
            dsd::Vector3 mean;
            dsd::Matrix3 covar;
            double weight;

            Gaussian3d() : mean(), covar(dsd::Matrix3::Identity()), weight(1.0) {}

            Gaussian3d(const dsd::Vector3 &m, const dsd::Matrix3 &c, double w)
            {
                mean = m;
                covar = c;
                weight = w;
            }

            ~Gaussian3d() {}

            double eval(const dsd::Vector3 &v) const
            {
                double k = 1.0 / sqrt(2.0 * M_PI * covar.determinant());
                double arg = (v - mean).transpose() * covar.inverse() * (v - mean);
                return k * exp(-0.5 * arg);
            }
        };

#if __cplusplus < 201703L
        using VectorGaussian =
            std::deque<Gaussian, Eigen::aligned_allocator<Gaussian>>;
#else
        using VectorGaussian = std::deque<Gaussian>;
#endif

        /**
         * Default constructor.
         */
        GaussianMixtureEstimator();

        /**
         * Destructor.
         */
        virtual ~GaussianMixtureEstimator();

        /**
         * Clear gaussians_ vector
         */
        void clearGaussians();

        /**
         * Computes the Gaussian parameters from the given samples.
         * @param samples
         */
        virtual void compute(const dsd::VectorVector2 &samples) = 0;

        /**
         * Returns the number of components/hypotheses of the mixture.
         * @return
         */
        size_t size() const;

        /**
         * Returns the estimated mean value of i-th Gaussian distribution in the
         * mixture.
         * @param i the index of the distribution/hypothesis
         * @return the mean vector
         */
        const dsd::Vector2 &mean(int i) const;

        /**
         * Returns the estimated covariance of i-th Gaussian distribution in the
         * mixture.
         * @param i the index of the distribution/hypothesis
         * @return the covariance matrix
         */
        const dsd::Matrix2 &covariance(int i) const;

        /**
         * Returns the estimated weight of i-th Gaussian distribution in the
         * mixture, i.e. the probability that the i-th component/hypothesis is
         * drawn.
         * @param i the index of the distribution/hypothesis
         * @return the weight
         */
        double weight(int i) const;

        /**
         * Returns a const reference to the vector of gaussians.
         * @return
         */
        const VectorGaussian &gaussians() const;

        /**
         * Returns a const reference to the vector of gaussians3.
         */
        const std::vector<Gaussian3d> &gaussians3() const;

        /**
         * Check whether covariance has norm null or weight is close to 0
         */
        bool checkIsGaussianNull(const dsd::Matrix2 &covar, const double &w) const;

        /**
         * Check whether covariance has norm null or weight is close to 0 (3D
         * version)
         */
        bool checkIsGaussianNull(const dsd::Matrix3 &covar, const double &w) const;

        /**
         * Exports the Gaussian mixture parameters, i.e. means, covariances and
         * weights, into separate vectors.
         * @param means std::vector of mean vectors
         * @param covariances std::vector of covariance matrices
         * @param weights std::vector of weights
         */
        void exportGaussians(dsd::VectorVector2 &means,
                             dsd::VectorMatrix2 &covariances,
                             std::vector<double> &weights) const;

        /**
         * 3D version of exportGaussians
         */
        void exportGaussians3d(dsd::VectorVector3 &means,
                               dsd::VectorMatrix3 &covariances,
                               std::vector<double> &weights) const;

        /**
         * Executes Expectation Maximization (EM) updating the Gaussian
         * mixture stored in the class.
         * @param samples vector of samples
         * @param stepNum number of iteration of EM
         */
        void executeEM(const dsd::VectorVector2 &samples, int stepNum = 1);

    protected:
        //        dsd::VectorVector2 means_;
        //        dsd::VectorMatrix2 covROFL_;
        //        std::vector<double> weights_;
        VectorGaussian gaussians_;

        std::vector<Gaussian3d> gaussians3_;
    };

    //-----------------------------------------------------
    // GaussianMixtureEstimatorScan
    //-----------------------------------------------------

    class GaussianMixtureEstimatorScan : public GaussianMixtureEstimator
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        // using IndexInterval = std::pair<int, int>;
        struct IndexInterval
        {
            int first;
            int last;
            int num;
        };

        /**
         * Default constructor.
         */
        GaussianMixtureEstimatorScan();

        /**
         * Destructor.
         */
        virtual ~GaussianMixtureEstimatorScan();

        /**
         * Sets the threshold above which a gap between consecutive points is
         * detected.
         * @param dg the distance gap threshold
         */
        void setDistanceGap(double dg) { distanceGap_ = dg; }

        /**
         * Sets the splitting distance for segment detection
         * @param ds the distance threshold to split
         */
        void setDistanceSplit(double ds) { distanceSplit_ = ds; }

        /**
         * Sets the minimum value of standard deviation of Gaussians.
         * @param sm the minimum standard deviation
         */
        void setSigmaMin(double sm) { sigmaMin_ = sm; }

        /**
         * Computes the Gaussian parameters from the given samples.
         * @param samples sorted in counter-clockwise order
         */
        virtual void compute(const dsd::VectorVector2 &samples);

        /**
         * Returns the i-th interval.
         * @param i
         */
        const IndexInterval &interval(int i) const
        {
            ROFL_ASSERT(0 <= i && i < intervals_.size());
            return intervals_.at(i);
        }

    private:
        std::deque<IndexInterval> intervals_; // used for debug
        double distanceGap_;
        double distanceSplit_;
        double sigmaMin_;

        /**
         * Finds the farthest point in the set from the line through points first
         * and last, i.e. points[first] and points[last]. The farthest index belongs
         * to interval first and last.
         * @param points the complete set of points
         * @param first the index of the first point in the segment interval
         * @param last the index of the last point in the segment interval
         * @param farthest the index of the farthest point from the line
         * @param distMax the distance of the farthest point from the line
         */
        void findFarthest(const dsd::VectorVector2 &points,
                          int first,
                          int last,
                          int &farthest,
                          double &distMax) const;

        /**
         * Computes the Gaussian mean and covariance matrix of points in interval
         * between first and last.
         * @param points the set of points
         * @param first the index of the first point
         * @param last the intex of the last point
         * @param mean the mean vector
         * @param covar the covariance matrix
         */
        void estimateGaussianFromPoints(const dsd::VectorVector2 &points,
                                        int first,
                                        int last,
                                        dsd::Vector2 &mean,
                                        dsd::Matrix2 &covar) const;

        /**
         * Computes the Gaussian distribution, i.e. its parameters, assuming the
         * input points are not samples of the Gaussian, but rather approximately
         * uniform sampled points on the segment.
         * @param points the set of points
         * @param first the index of the first point
         * @param last the intex of the last point
         * @param mean the mean vector
         * @param covar the covariance matrix
         */
        void estimateGaussianFromSegment(const dsd::VectorVector2 &points,
                                         int first,
                                         int last,
                                         dsd::Vector2 &mean,
                                         dsd::Matrix2 &covar) const;
    };

    //-----------------------------------------------------
    // GaussianMixtureEstimatorHierarchical
    //-----------------------------------------------------

    /**
     * Computes a Gaussian Mixture approximately exploiting some of the ideas in
     *
     * Benjamin Eckart, Kihwan Kim Jan Kau,
     * "HGMR: Hierarchical Gaussian Mixtures for Adaptive 3D Registration",
     * ECCV 2018.
     *
     */
    class GaussianMixtureEstimatorHierarchical : public GaussianMixtureEstimator
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // Private
        using PointContainer = MortonOctree<2, double, int32_t>;
        using Iterator = PointContainer::Iterator;
        using ConstIterator = PointContainer::ConstIterator;
        using ConstInterval = std::pair<ConstIterator, ConstIterator>;

        /**
         * Default constructor.
         */
        GaussianMixtureEstimatorHierarchical();

        /**
         * Destructor.
         */
        virtual ~GaussianMixtureEstimatorHierarchical();

        double getSigmaMin() const;

        double getIseThreshold() const;

        /**
         * Sets the minimum value of standard deviation of Gaussians.
         * @param sm the minimum standard deviation
         */
        void setSigmaMin(double sm);

        void setCovarWidth(double cw);

        void setInlierPerc(double ip);

        void setChiConfidence(double conf);

        void setIseThreshold(double iseTh);

        void setCellSizeMax(double s);

        /**
         * Init as 3D Isotropic GME
         */
        void initIsotropic(const dsd::VectorVector3 &samples);

        /**
         * Computes the Gaussian parameters from the given samples.
         * @param samples sorted in counter-clockwise order
         */
        virtual void compute(const dsd::VectorVector2 &samples);

    private:
        PointContainer data_;
        double sigmaMin_;
        double covarWidth_;
        double chi2Thres_;
        double iseThres_;
        double inlierPerc_;
        int levelMax_;

        bool estimateGaussianFromPoints(const ConstIterator &beg,
                                        const ConstIterator &end,
                                        dsd::Vector2 &mean,
                                        dsd::Matrix2 &covar,
                                        double &w) const;

        bool estimateGaussianFromSegment(const ConstIterator &beg,
                                         const ConstIterator &end,
                                         dsd::Vector2 &mean,
                                         dsd::Matrix2 &covar,
                                         double &w) const;

        bool estimateGaussianISE(const ConstIterator &beg,
                                 const ConstIterator &end,
                                 dsd::Vector2 &mean,
                                 dsd::Matrix2 &covar,
                                 double &wMerged) const;
    };

    /**
     * Class for computing Hierachical Gaussian optimization.
     * It is not elegant to write a duplicate of the 2D algotihm designed
     * for 3D case.
     */
    class GaussianMixtureEstimatorHierarchical3d
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // Private
        using PointContainer = MortonOctree<3, double, int32_t>;
        using Iterator = PointContainer::Iterator;
        using ConstIterator = PointContainer::ConstIterator;
        using ConstInterval = std::pair<ConstIterator, ConstIterator>;

        // TODO: remove Gaussian3d definition from
        using Gaussian = GaussianMixtureEstimator::Gaussian3d;
        using VectorGaussian = std::deque<Gaussian>;

        /**
         * Default constructor.
         */
        GaussianMixtureEstimatorHierarchical3d();

        /**
         * Destructor.
         */
        virtual ~GaussianMixtureEstimatorHierarchical3d();

        double getSigmaMin() const;

        double getIseThreshold() const;

        /**
         * Sets the minimum value of standard deviation of Gaussians.
         * @param sm the minimum standard deviation
         */
        void setSigmaMin(double sm);

        void setCovarWidth(double cw);

        void setInlierPerc(double ip);

        void setChiConfidence(double conf);

        void setIseThreshold(double iseTh);

        void setCellSizeMax(double s);

        const VectorGaussian &gaussians() const;

        void exportGaussians(dsd::VectorVector3 &means,
                             dsd::VectorMatrix3 &covariances,
                             std::vector<double> &weights) const;

        /**
         * Computes the Gaussian parameters from the given samples.
         * @param samples sorted in counter-clockwise order
         */
        virtual void compute(const dsd::VectorVector3 &samples);

    private:
        PointContainer data_;
        VectorGaussian gaussians_;
        double sigmaMin_;
        double covarWidth_;
        double chi2Thres_;
        double iseThres_;
        double inlierPerc_;
        int levelMax_;

        bool estimateGaussianISE(const ConstIterator &beg,
                                 const ConstIterator &end,
                                 dsd::Vector3 &mean,
                                 dsd::Matrix3 &covar,
                                 double &wMerged) const;

        bool checkIsGaussianNull(const dsd::Matrix3 &covar, const double &w) const;
    };

} // namespace gme

#endif /* GAUSSIANMIXTUREESTIMATOR_H */
