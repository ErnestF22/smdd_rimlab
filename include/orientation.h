#include <math.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>

#include <csm/csm.h>
#include <csm/csm_all.h>

#include <egsl/egsl_macros.h>

#include <rofl/common/macros.h>

void find_neighbours(LDP ld,
                     int i,
                     int max_num,
                     int *indexes,
                     size_t *num_found);

void filter_orientation(double theta0,
                        double rho0,
                        size_t n,
                        const double *thetas,
                        const double *rhos,
                        double *alpha,
                        double *cov0_alpha);

/** Requires the "cluster" field to be set */
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
                           thetas, readings, &alpha, &cov0_alpha);

        if (gsl_isnan(alpha))
        {
            ld->alpha[i] = GSL_NAN;
            ld->cov_alpha[i] = GSL_NAN;
            ld->alpha_valid[i] = 0;
        }
        else
        {
            ld->alpha[i] = alpha;
            ld->cov_alpha[i] = cov0_alpha * square(sigma);
            ld->alpha_valid[i] = 1;
        }
        /* printf("---------- i = %d alpha = %f sigma=%f cov_alpha = %f\n", i,
         * alpha, ld->cov_alpha[i]);*/
    }
}

/** A very cool algorithm for finding the orientation */

void filter_orientation(double theta0,
                        double rho0,
                        size_t n,
                        const double *thetas,
                        const double *rhos,
                        double *alpha,
                        double *cov0_alpha)
{
    // ROFL_VAR3(theta0, rho0, n);

    egsl_push();
    /* Y = L x + R epsilon */
    val Y = zeros(n, 1);
    val L = ones(n, 1);
    val R = zeros(n, n + 1);

    size_t i;
    for (i = 0; i < n; i++)
    {
        // ROFL_VAR3(i, thetas[i] - theta0, theta0);

        if (fabs(thetas[i] - theta0) < 1e-6)
        {
            ROFL_VAR3(i, thetas[i], theta0);
            continue;
        }

        // if(!gsl_finite(f1))
        // continue;

        *egsl_atmp(Y, i, 0) = (rhos[i] - rho0) / (thetas[i] - theta0);
        *egsl_atmp(R, i, 0) = -1 / (thetas[i] - theta0);
        *egsl_atmp(R, i, i + 1) = +1 / (thetas[i] - theta0);
    }

    // ROFL_VAR1("Passed first for cycle in filter_orientation()");

    val eRinv = inv(m(R, tr(R)));
    val vcov_f1 = inv(m3(tr(L), eRinv, L));
    val vf1 = m4(vcov_f1, tr(L), eRinv, Y);

    // ROFL_VAR3(eRinv.gslm->data, vcov_f1.gslm->data, vf1.gslm->data);

    double cov_f1 = *egsl_atmp(vcov_f1, 0, 0);
    double f1 = *egsl_atmp(vf1, 0, 0);

    *alpha = theta0 - atan(f1 / rho0);

    // ROFL_VAR2(rho0, atan(f1 / rho0));

    if (cos(*alpha) * cos(theta0) + sin(*alpha) * sin(theta0) > 0)
        *alpha = *alpha + M_PI;

    double dalpha_df1 = rho0 / (square(rho0) + square(f1));
    double dalpha_drho = -f1 / (square(rho0) + square(f1));

    *cov0_alpha = square(dalpha_df1) * cov_f1 + square(dalpha_drho);

    if (!gsl_finite(*alpha))
    {
        // egsl_print("Y", Y);
        // egsl_print("L", L);
        // egsl_print("R", R);
        // egsl_print("eRinv", eRinv);
        // egsl_print("vcov_f1", vcov_f1);

        // printf("   f1 = %f cov =%f \n", f1, cov_f1);
        // printf("   f1/rho = %f \n", f1 / rho0);
        // printf("   atan = %f \n", atan(f1 / rho0));
        // printf("   theta0= %f \n", theta0);
    }

    egsl_pop();
    /*
    //	printf("dalpha_df1 = %f dalpha_drho = %f\n",dalpha_df1,dalpha_drho);
    //	printf("f1 = %f covf1 = %f alpha = %f cov_alpha = %f\n
    ",f1,cov_f1,*alpha,*cov0_alpha);
    //	printf("sotto = %f\n ",(square(rho0) + square(f1)));

    //	printf("   alpha = %f sigma= %f°\n", *alpha,
    rad2deg(0.01*sqrt(*cov0_alpha)));

        printf("l= ");
        gsl_matrix_fprintf(stdout, l, "%f");
        printf("\ny= ");
        gsl_matrix_fprintf(stdout, y, "%f");
        printf("\nr= ");
        gsl_matrix_fprintf(stdout, r, "%f");
        printf("\ninv(r*r)= ");
        gsl_matrix_fprintf(stdout, Rinv, "%f");
        printf("\nf1 = %lf ",f1);
        printf("\ncov_f1 = %lf ",cov_f1);
    */

    // ROFL_VAR1("Exiting filter_orientation()");
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

void computeFIM(Eigen::Matrix3f &fim,
                float sigma,
                float theta_robot,
                const std::vector<double> &alphas,
                const std::vector<double> &phis,
                const std::vector<double> &ranges)
{
    int n = ranges.size(); // = alphas.size();
    float sigmaSq = sigma * sigma;
    ROFL_VAR2(n, sigma);

    // TODO: number of for cycles can most likely be limited

    // betas
    std::vector<float> betas(n, 0.0f);
    for (int i = 0; i < n; ++i)
    {
        float alpha_i = alphas[i];
        float phi_i = phis[i];
        float beta_i = alpha_i - (theta_robot + phi_i);
        betas[i] = beta_i;
    }

    // building matrix
    fim = Eigen::Matrix3f::Zero(); // resetting matrix
    for (int i = 0; i < n; ++i)
    {
        // ROFL_VAR1(i);
        float r = ranges[i];
        Eigen::Matrix3f fim_i(Eigen::Matrix3f::Zero());

        Eigen::Vector2f v_alpha_i;
        float alpha_i = alphas[i];
        float beta_i = betas[i];
        float r_i = ranges[i];
        v_alpha_i << cos(alpha_i), sin(alpha_i);

        float tbi = tan(beta_i);
        float cbi = cos(beta_i);

        Eigen::Matrix2f fim_i_11(Eigen::Matrix2f::Zero());
        fim_i_11 = v_alpha_i * v_alpha_i.transpose() / cos(beta_i);
        Eigen::Vector2f fim_i_12(Eigen::Vector2f::Zero());
        fim_i_12 = r_i * (tbi / cbi) * v_alpha_i;
        // Eigen::RowVector2f fim_i_21 = fim_i_12.transpose();
        float fim_i_22 = r_i * r_i * tbi * tbi;

        fim_i.block(0, 0, 2, 2) =
            fim_i_11; // block() params startRow, startCol, numRows, numCols
        fim_i.block(0, 2, 2, 1) =
            fim_i_12; // block() params startRow, startCol, numRows, numCols
        fim_i.block(2, 0, 1, 2) =
            fim_i_12.transpose(); // block() params startRow, startCol,
                                  // numRows, numCols
        fim_i(2, 2) = fim_i_22;

        // adding fim_i to fim
        fim += fim_i;
    }

    fim *= (1 / sigmaSq);
}