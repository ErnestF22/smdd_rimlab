#ifndef FPM_H_
#define FPM_H_

#include <eigen3/Eigen/Dense>

#include <rofl/common/macros.h>

namespace dsd
{
    class FPM
    {
    public:
        int r_;
        int c_;

        Eigen::MatrixXd M_;

        Eigen::MatrixXi Mbools_;

        int l_;
        int sideSz_;
        int neighborhoodSz_;

        bool multipleMaxima_;

        FPM() {}

        FPM(const Eigen::MatrixXd &M, int l)
        {
            M_.resizeLike(M);

            l_ = l;

            M_ = M;
            r_ = M.rows();
            c_ = M.cols();

            Mbools_.resizeLike(M);
            Mbools_.setOnes();

            sideSz_ = 2 * l + 1;
            neighborhoodSz_ = sideSz_ * sideSz_ - 1;

            multipleMaxima_ = false;
        }

        ~FPM() {}

        void is8neighborhoodMax(int i, int j, Eigen::MatrixXi &maskIJ)
        {
            double centerVal = M_(i, j);

            auto eightneighb = M_.block(i - l_, j - l_, sideSz_, sideSz_);

            ROFL_ASSERT_VAR1(centerVal - eightneighb(l_, l_) < 1e-8, centerVal - eightneighb(l_, l_))

            Eigen::Index maxRow, maxCol;
            float max = eightneighb.maxCoeff(&maxRow, &maxCol);

            maskIJ(maxRow, maxCol) = 1;
        }

        void is8neighborhoodMaxEdge(int i, int j, Eigen::MatrixXi &maskIJ)
        {
            double centerVal = M_(i, j);

            Eigen::MatrixXd eightneighbEdge(sideSz_, sideSz_);

            // fill eightneighbEdge
            for (int i2 = -l_; i2 <= l_; ++i2)
                for (int j2 = -l_; j2 <= l_; ++j2)
                {
                    int ie = i + i2;
                    int je = j + j2;
                    if (ie >= M_.rows())
                        ie -= M_.rows();
                    else if (ie < 0)
                        ie += M_.rows();
                    if (je >= M_.cols())
                        je -= M_.cols();
                    else if (je < 0)
                        je += M_.cols();

                    eightneighbEdge(i2 + l_, j2 + l_) = M_(ie, je);
                }

            ROFL_ASSERT(centerVal - eightneighbEdge(l_, l_) < 1e-8)

            // ROFL_VAR4(i, j, centerVal, eightneighbEdge);

            Eigen::Index maxRow, maxCol;
            float max = eightneighbEdge.maxCoeff(&maxRow, &maxCol);

            // ROFL_VAR3(max, maxRow, maxCol)

            maskIJ(maxRow, maxCol) = 1;
        }

        void applyMaskIJ(const Eigen::MatrixXi &maskIJ, int i, int j)
        {
            for (int i2 = -l_; i2 <= l_; ++i2)
                for (int j2 = -l_; j2 <= l_; ++j2)
                {
                    Mbools_(i + i2, j + j2) *= maskIJ(i2 + l_, j2 + l_);
                }
        }

        void applyMaskIJedge(const Eigen::MatrixXi &maskIJedge, int i, int j)
        {
            for (int i2 = -l_; i2 <= l_; ++i2)
                for (int j2 = -l_; j2 <= l_; ++j2)
                {
                    int ie = i + i2;
                    int je = j + j2;
                    if (ie >= Mbools_.rows())
                        ie -= Mbools_.rows();
                    else if (ie < 0)
                        ie += Mbools_.rows();
                    if (je >= Mbools_.cols())
                        je -= Mbools_.cols();
                    else if (je < 0)
                        je += Mbools_.cols();

                    // eightneighbEdge(i2 + 1, j2 + 1) = M_(ie, je);

                    Mbools_(ie, je) *= maskIJedge(i2 + l_, j2 + l_);
                }
        }

        void run()
        {
            ROFL_VAR2(r_, c_);
            ROFL_VAR3(l_, sideSz_, neighborhoodSz_);

            // First, run on elements that are NOT edges

            for (int i = l_; i < r_ - l_; ++i)
                for (int j = l_; j < c_ - l_; ++j)
                {
                    if (Mbools_(i, j) == 0)
                        continue;

                    Eigen::MatrixXi maskIJ(Eigen::MatrixXi::Zero(sideSz_, sideSz_));
                    is8neighborhoodMax(i, j, maskIJ);
                    applyMaskIJ(maskIJ, i, j);
                }

            // std::cout << "Mbools_ NO EDGES" << std::endl
            //           << Mbools_ << std::endl;

            // AAA
            // D*B
            // CCB

            /* A */
            // ROFL_VAR1("edge A")
            for (int i = 0; i < l_; ++i)
                for (int j = 0; j < M_.cols(); ++j)
                {
                    // ROFL_VAR3("edgeA", i, j)

                    if (Mbools_(i, j) == 0)
                        continue;

                    Eigen::MatrixXi maskIJ(Eigen::MatrixXi::Zero(sideSz_, sideSz_));
                    is8neighborhoodMaxEdge(i, j, maskIJ);
                    applyMaskIJedge(maskIJ, i, j);
                }

            // std::cout << "Mbools_ after edges A" << std::endl
            //           << Mbools_ << std::endl;

            /* B */
            // ROFL_VAR1("edge B")
            for (int i = 0; i < l_; ++i)
                for (int j = l_; j < M_.rows(); ++j)
                {
                    // ROFL_VAR3("edgeB", j, M_.cols() - i - 1)
                    if (Mbools_(j, M_.cols() - i - 1) == 0)
                        continue;

                    Eigen::MatrixXi maskIJ(Eigen::MatrixXi::Zero(sideSz_, sideSz_));
                    is8neighborhoodMaxEdge(j, M_.cols() - i - 1, maskIJ);
                    applyMaskIJedge(maskIJ, j, M_.cols() - i - 1);
                }

            // std::cout << "Mbools_  after edges B" << std::endl
            //           << Mbools_ << std::endl;

            /* C */
            // ROFL_VAR1("edge C")
            for (int i = 0; i < l_; ++i)
                for (int j = M_.cols() - 2; j >= 0; --j)
                {
                    // ROFL_VAR3("edgeC", i, j)

                    if (Mbools_(M_.rows() - i - 1, j) == 0)
                        continue;

                    Eigen::MatrixXi maskIJ(Eigen::MatrixXi::Zero(sideSz_, sideSz_));
                    is8neighborhoodMaxEdge(M_.rows() - i - 1, j, maskIJ);
                    applyMaskIJedge(maskIJ, M_.rows() - i - 1, j);
                }

            // std::cout << "Mbools_  after edges C" << std::endl
            //           << Mbools_ << std::endl;

            /* D */
            // ROFL_VAR1("edge D")
            for (int i = 0; i < l_; ++i)
                for (int j = M_.rows() - 2; j > 0; --j)
                {
                    // ROFL_VAR3("edgeD", i, j)

                    if (Mbools_(j, 0) == 0)
                        continue;

                    Eigen::MatrixXi maskIJ(Eigen::MatrixXi::Zero(sideSz_, sideSz_));
                    is8neighborhoodMaxEdge(j, i, maskIJ);
                    applyMaskIJedge(maskIJ, j, i);
                }

            // if no maxima are output (i.e., there are multiple maxima, sufficiently close to one another)
            if (Mbools_.sum() == 0)
            {
                multipleMaxima_ = true;
                Eigen::Index maxRow, maxCol;
                double max = M_.maxCoeff(&maxRow, &maxCol);
                auto MboolsTmp = M_.array() > (max - 1e-3);
                // ROFL_VAR1(MboolsTmp)
                Mbools_ = MboolsTmp.cast<int>();
                // ROFL_VAR1(Mbools_)
            }
        }

        Eigen::MatrixXi getIndicesMax() const
        {
            return Mbools_;
        }
    };

} // end of namespace DSD

#endif /*FPM_H_*/