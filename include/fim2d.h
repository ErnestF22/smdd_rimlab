#include <rofl/common/profiler.h>
#include <eigen3/Eigen/Dense>

void cartesian2polar(const std::vector<Eigen::Vector2f> &pts,
                     std::vector<float> &rhos,
                     std::vector<float> &thetas)
{
    int n = pts.size();
    rhos.clear();
    thetas.clear();
    thetas.reserve(n);
    rhos.reserve(n);

    for (int i = 0; i < n; ++i)
    {
        float x = pts[i](0);
        float y = pts[i](1);
        float rho = sqrt(x * x + y * y);
        float theta = atan2(y, x);
        rhos[i] = rho;
        thetas[i] = theta;
    }
}

void computeFIM2D(Eigen::Matrix3f &fim,
                  float sigma,
                  float thetaRobot,
                  const std::vector<float> &alphas,
                  const std::vector<float> &phis,
                  const std::vector<float> &ranges)
{
    int n = ranges.size(); // = alphas.size();
    float sigmaSq = sigma * sigma;
    ROFL_VAR2(n, sigma);

    rofl::ScopedTimer timer("FIM_2D");

    // TODO: number of for cycles can most likely be limited

    // betas
    std::vector<float> betas(n, 0.0f);

    // building matrix
    fim = Eigen::Matrix3f::Zero(); // resetting matrix
    for (int i = 0; i < n; ++i)
    {
        float alphaI = alphas[i];
        float phiI = phis[i];
        float betaI = alphaI - (thetaRobot + phiI);

        if (!std::isfinite(alphaI) || !std::isfinite(betaI) ||
            !std::isfinite(phiI))
        {
            continue;
        }
        betas[i] = betaI;

        // ROFL_VAR1(i);
        float r = ranges[i];
        Eigen::Matrix3f fimI(Eigen::Matrix3f::Zero());

        Eigen::Vector2f vAlphaI;
        float rI = ranges[i];
        if (!std::isfinite(rI))
        {
            continue;
        }
        vAlphaI << cos(alphaI), sin(alphaI);

        float tbi = tan(betaI);
        float cbi = cos(betaI);

        // ROFL_VAR5(i, alphaI, betaI, rI, sigma);

        if (!std::isfinite(betaI) || fabs(cbi) < 1e-6)
        {
            ROFL_VAR2(betaI, cbi);
            continue;
        }

        Eigen::Matrix2f fimI11(Eigen::Matrix2f::Zero());
        fimI11 = vAlphaI * vAlphaI.transpose() / (cbi * cbi);
        Eigen::Vector2f fimI12(Eigen::Vector2f::Zero());
        fimI12 = rI * (tbi / cbi) * vAlphaI;
        // Eigen::RowVector2f fimI_21 = fimI12.transpose();
        float fimI22 = rI * rI * tbi * tbi;

        fimI.block(0, 0, 2, 2) =
            fimI11; // block() params startRow, startCol, numRows, numCols
        fimI.block(0, 2, 2, 1) =
            fimI12; // block() params startRow, startCol, numRows, numCols
        fimI.block(2, 0, 1, 2) =
            fimI12.transpose(); // block() params startRow, startCol,
                                // numRows, numCols
        fimI(2, 2) = fimI22;

        // adding fimI to fim
        fim += fimI;
    }

    fim *= (1 / sigmaSq);
}