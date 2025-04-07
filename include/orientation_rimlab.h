#include <eigen3/Eigen/Dense>

#include <vector>

#include <rofl/common/macros.h>

/** A very cool algorithm for finding the orientation */

void filterOrientationEigen(double theta0,
                            double rho0,
                            int n,
                            const std::vector<float> &thetas,
                            const std::vector<float> &rhos,
                            float &alpha,
                            float &cov0_alpha)
{
    // ROFL_VAR3(theta0, rho0, n);

    /* Y = L x + R epsilon */

    size_t i;

    std::vector<bool> considerRow(n, true);
    int numRowsBad = 0;
    for (i = 0; i < n; i++)
    {
        // ROFL_VAR3(i, thetas[i] - theta0, theta0);

        if (fabs(thetas[i] - theta0) < 1e-9)
        {
            considerRow[i] = false;
            numRowsBad++;
            continue;
            // ROFL_ASSERT_VAR3(0, i, thetas[i], theta0);
        }
    }

    Eigen::MatrixXf Y(n - numRowsBad, 1);
    Y.setZero();
    Eigen::MatrixXf L(n - numRowsBad, 1);
    L.setOnes();
    Eigen::MatrixXf R(n - numRowsBad, n - numRowsBad + 1);
    R.setZero();

    int rowsConsidered = 0;
    for (i = 0; i < n; i++)
    {
        if (considerRow[i])
        {
            Y(rowsConsidered, 0) = (rhos[i] - rho0) / (thetas[i] - theta0);
            R(rowsConsidered, 0) = -1 / (thetas[i] - theta0);
            R(rowsConsidered, rowsConsidered + 1) = +1 / (thetas[i] - theta0);
            rowsConsidered++;
        }
    }

    // ROFL_VAR1("Passed first for cycle in filter_orientation()");

    Eigen::MatrixXf eRinv = (R * R.transpose()).inverse();
    Eigen::MatrixXf vcov_f1 = (L.transpose() * eRinv * L).inverse();
    Eigen::MatrixXf vf1 = vcov_f1 * L.transpose() * eRinv * Y;

    // ROFL_VAR3(eRinv, vcov_f1, vf1);

    double cov_f1 = vcov_f1(0, 0);
    double f1 = vf1(0, 0);

    alpha = theta0 - atan(f1 / rho0);

    // ROFL_VAR2(rho0, atan(f1 / rho0));

    if (cos(alpha) * cos(theta0) + sin(alpha) * sin(theta0) > 0)
        alpha = alpha + M_PI;

    double dalpha_df1 = rho0 / (rho0 * rho0 + f1 * f1);
    double dalpha_drho = -f1 / (rho0 * rho0 + f1 * f1);

    // *cov0_alpha = square(dalpha_df1) * cov_f1 + square(dalpha_drho);
    cov0_alpha = dalpha_df1 * dalpha_df1 * cov_f1 + dalpha_drho * dalpha_drho;

    if (std::isnan(alpha) || std::isnan(cov0_alpha))
    {
        // ROFL_MSG("alpha OR cov0_alpha NAN");
        // ROFL_VAR2(alpha, cov0_alpha);

        // std::cout << "   f1 = " << f1 << " cov = " << cov_f1 << std::endl;
        // std::cout << "   f1/rho0 = " << f1 / rho0 << std::endl;
        // std::cout << "   atan = " << atan(f1 / rho0) << std::endl;
        // std::cout << "   theta0= " << theta0 << std::endl;
    }

    // ROFL_VAR1("Exiting filter_orientation()");
}
