#include <gsl/gsl_nan.h>
#include <math.h>

#include <csm/csm.h>

#include <rofl/common/macros.h>

void filter_orientation(double theta0,
                        double rho0,
                        size_t n,
                        const double *thetas,
                        const double *rhos,
                        double &alpha,
                        double &cov0_alpha);

void find_neighbours(LDP ld,
                     int i,
                     int max_num,
                     int *indexes,
                     size_t *num_found);

/** Requires the "cluster" field to be set */
void ld_compute_orientation(LDP ld, int size_neighbourhood, double sigma);

/** A very cool algorithm for finding the orientation */

void filter_orientation(double theta0,
                        double rho0,
                        size_t n,
                        const double *thetas,
                        const double *rhos,
                        double &alpha,
                        double &cov0_alpha)
{
    ROFL_VAR3(theta0, rho0, n);

    // egsl_push();
    /* Y = L x + R epsilon */
    Eigen::MatrixXf Y(n, 1);
    Y.setZero();
    Eigen::MatrixXf L(n, 1);
    L.setOnes();
    Eigen::MatrixXf R(n, n + 1);
    R.setZero();

    size_t i;
    for (i = 0; i < n; i++)
    {
        ROFL_VAR3(i, thetas[i] - theta0, theta0);

        if (fabs(thetas[i] - theta0) < 1e-6)
        {
            ROFL_ASSERT_VAR3(0, i, thetas[i], theta0);
        }

        Y(i, 0) = (rhos[i] - rho0) / (thetas[i] - theta0);
        R(i, 0) = -1 / (thetas[i] - theta0);
        R(i, i + 1) = +1 / (thetas[i] - theta0);
    }

    ROFL_VAR1("Passed first for cycle in filter_orientation()");

    Eigen::MatrixXf eRinv = (R * R.transpose()).inverse();
    ;
    Eigen::MatrixXf vcov_f1 = (L.transpose() * eRinv * L).inverse();
    Eigen::MatrixXf vf1 = vcov_f1 * L.transpose() * eRinv * Y;

    // ROFL_VAR3(eRinv.gslm->data, vcov_f1.gslm->data, vf1.gslm->data);

    double cov_f1 = vcov_f1(0, 0);
    double f1 = vf1(0, 0);

    alpha = theta0 - atan(f1 / rho0);

    ROFL_VAR2(rho0, atan(f1 / rho0));

    if (cos(alpha) * cos(theta0) + sin(alpha) * sin(theta0) > 0)
        alpha = alpha + M_PI;

    double dalpha_df1 = rho0 / (rho0 * rho0 + f1 * f1);
    double dalpha_drho = -f1 / (rho0 * rho0 + f1 * f1);

    cov0_alpha = (dalpha_df1 * dalpha_df1) * cov_f1 + dalpha_drho * dalpha_drho;

    if (std::isnan(alpha))
    {
        std::cout << "Y " << Y << std::endl;
        std::cout << "L " << L << std::endl;
        std::cout << "R " << R << std::endl;
        std::cout << "eRinv " << eRinv << std::endl;
        std::cout << "vcov_f1 " << vcov_f1 << std::endl;

        printf("   f1 = %f cov =%f \n", f1, cov_f1);
        printf("   f1/rho = %f \n", f1 / rho0);
        printf("   atan = %f \n", atan(f1 / rho0));
        printf("   theta0= %f \n", theta0);
    }

    // egsl_pop();

    ROFL_VAR1("Exiting filter_orientation()");
}

/* indexes: an array of size "max_num*2" */
void find_neighbours(LDP ld,
                     int i,
                     int max_num,
                     int *indexes,
                     size_t *num_found)
{
    *num_found = 0;

    int up = i;
    while ((up + 1 <= i + max_num) && (up + 1 < ld->nrays) &&
           ld_valid_ray(ld, up + 1) &&
           (ld->cluster[up + 1] == ld->cluster[i]))
    {
        up += 1;
        indexes[(*num_found)++] = up;
    }
    int down = i;
    while ((down >= i - max_num) && (down - 1 >= 0) &&
           ld_valid_ray(ld, down - 1) &&
           (ld->cluster[down - 1] == ld->cluster[i]))
    {
        down -= 1;
        indexes[(*num_found)++] = down;
    }
}

void ld_compute_orientation(LDP ld, int size_neighbourhood, double sigma)
{
    int i;
    for (i = 0; i < ld->nrays; i++)
    {
        if (!ld_valid_ray(ld, i) || (ld->cluster[i] == -1))
        {
            ld->alpha[i] = GSL_NAN;
            ld->cov_alpha[i] = GSL_NAN;
            ld->alpha_valid[i] = 0;
            continue;
        }

        int neighbours[size_neighbourhood * 2];
        size_t num_neighbours;
        find_neighbours(ld, i, size_neighbourhood, neighbours, &num_neighbours);

        if (0 == num_neighbours)
        {
            ld->alpha[i] = GSL_NAN;
            ld->cov_alpha[i] = GSL_NAN;
            ld->alpha_valid[i] = 0;
            continue;
        }

        /*		printf("orientation for i=%d:\n",i); */
        double thetas[num_neighbours];
        double readings[num_neighbours];
        size_t a = 0;
        for (a = 0; a < num_neighbours; a++)
        {
            thetas[a] = ld->theta[neighbours[a]];
            readings[a] = ld->readings[neighbours[a]];
            /* printf(" j = %d theta = %f rho = %f\n", neighbours[a],
             * thetas[a],readings[a]); */
        }

        double alpha = 42, cov0_alpha = 32;
        filter_orientation(ld->theta[i], ld->readings[i], num_neighbours,
                           thetas, readings, alpha, cov0_alpha);

        ROFL_VAR1(std::isnan(GSL_NAN));
        if (std::isnan(alpha))
        {
            ld->alpha[i] = GSL_NAN;
            ld->cov_alpha[i] = GSL_NAN;
            ld->alpha_valid[i] = 0;
        }
        else
        {
            ld->alpha[i] = alpha;
            ld->cov_alpha[i] = cov0_alpha * sigma * sigma;
            ld->alpha_valid[i] = 1;
        }
        /* printf("---------- i = %d alpha = %f sigma=%f cov_alpha = %f\n", i,
         * alpha, ld->cov_alpha[i]);*/
    }
}