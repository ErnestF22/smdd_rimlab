#ifndef DSD_UTILS_CUH_
#define DSD_UTILS_CUH_

#include <rofl/common/macros.h>
#include <cmath>
#include <eigen3/Eigen/Dense>

#include <eigen3/Eigen/Dense>

#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <rofl/common/macros.h>

#include <boost/lexical_cast.hpp>

#include <find_peaks_matrix.h>

#include <thrust/host_vector.h>

// void vonMisesStats3d(std::vector<double> &vmm,
//                      int nSamples,
//                      const VectorVector3 &mus,
//                      const VectorMatrix3 &sigmas,
//                      const std::vector<double> &weights)
// {
//     VectorVector3 musVonMises;
//     std::vector<double> kappas;
//     std::vector<double> weightsVmm;
//     double phiTmp;

//     size_t numGaussians = mus.size();
//     musVonMises.resize(2 * numGaussians);
//     kappas.resize(2 * numGaussians);
//     weightsVmm.resize(2 * numGaussians);

//     int validDiagCtr = 0;
//     for (int i = 0; i < numGaussians; ++i)
//     {
//         ROFL_VAR1(i)
//         Matrix3 eigvalsMat, eigvecs;
//         bool isDiVehiclealid = dsd::diagonalize3d(sigmas[i], eigvalsMat, eigvecs);
//         if (!isDiVehiclealid)
//         {
//             musVonMises[2 * i] = Vector3::Zero();
//             kappas[2 * i] = 0;
//             musVonMises[2 * i + 1] = Vector3::Zero();
//             kappas[2 * i + 1] = kappas[2 * i];
//             weightsVmm[2 * i] = 0;
//             weightsVmm[2 * i + 1] = 0;
//             continue;
//         }
//         auto eigvals = eigvalsMat.diagonal();
//         Eigen::Index minIdx, maxIdx;
//         double lmin = eigvals.minCoeff(&minIdx);
//         double lmax = eigvals.maxCoeff(&maxIdx);
//         if (fabs(lmin) > 1e-5 && fabs(lmax) > 1e-5 && std::isfinite(lmin) &&
//             std::isfinite(lmax))
//         {
//             ROFL_ASSERT_VAR2((lmin / lmax) >= 0, lmin, lmax);
//         }
//         else
//             continue;

//         validDiagCtr++;
//         ROFL_VAR1(validDiagCtr);
//         musVonMises[2 * i] = eigvecs.col(minIdx);
//         if (lmax > 10.0 * lmin)
//         {
//             kappas[2 * i] = 10.0;
//         }
//         else
//         {
//             kappas[2 * i] = lmax / lmin;
//         }

//         musVonMises[2 * i + 1] = -eigvecs.col(minIdx);
//         kappas[2 * i + 1] = kappas[2 * i];
//         weightsVmm[2 * i] = weights[i];
//         weightsVmm[2 * i + 1] = weights[i];

//         // ROFL_VAR3(i, phis[i], kappas[i]);
//     }
//     // ROFL_ASSERT(0)
//     ROFL_VAR1(numGaussians)
//     vonMisesMixture3d(vmm, nSamples, musVonMises, kappas, weightsVmm);
// }

__device__ double vonMisesFisherDevice(const double3 x, const double3 mu, const double k)
{
    // in this context k is the weight
    double c3k = k / (4 * M_PI * sinh(k));
    double scalarProd = mu.x * x.x + mu.y * x.y + mu.z * x.z;
    double e = exp(k * scalarProd);
    return c3k * e;
}

// void vonMisesMixture3d(std::vector<double> &vmm,
//                        int nSamples,
//                        const std::vector<Vector3> &mu,
//                        const std::vector<double> &k,
//                        const std::vector<double> &w)
__global__ void vmm3d_kernel(double *vmm, int nSamples, //!! mu size needs to be 3 times the size of k, w
                             const double3 *mu, const double *k, const int mukwSz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    // int stride = blockDim.x * gridDim.x;

    // int totalSzIn = nrowsIn*ncols; //matrix is considered of size nrows*ncols, with nrows = sumNaturalsUpToN(numPts)

    // int numRows = mukwSz;

    // int numCols = 2 * nSamples * nSamples;
    int sizePhis = 2 * nSamples; // TODO: this can maybe be passed directly
    // int sizeThetas = nSamples;   // TODO: this can maybe be passed directly

    int i = index % mukwSz;
    int j = (index - i) / mukwSz;

    // if (i >= mukwSz)
    //     printf("i >= mukwSz!! %d \t %d\n", i, mukwSz);

    // if (j >= numCols)
    //     printf("j >= numCols!! %d \t %d\n", j, numCols);

    int jPhi = j % sizePhis;
    int jTheta = (j - jPhi) / sizePhis;

    // if (jPhi >= sizePhis)
    //     printf("jPhi >= sizePhis!!\n");
    // if (jTheta >= sizeThetas)
    //     printf("jTheta >= sizeThetas!! %d \t %d\n", jTheta, sizeThetas);

    // double dtheta = M_PI / nSamples;
    double dtheta = M_PI / nSamples;

    // int vmmSz = nSamples * (2 * nSamples);

    // for
    double theta = dtheta * jTheta;
    // int idx = itheta * nSamples + iphi; //idx should be j in this context
    double phi = dtheta * jPhi;

    double3 xThetaPhi;
    double st = sin(theta);
    double ct = cos(theta);
    double cp = cos(phi);
    double sp = sin(phi);
    xThetaPhi.x = st * cp;
    xThetaPhi.y = st * sp;
    xThetaPhi.z = ct;

    double vmf = vonMisesFisherDevice(xThetaPhi, mu[i], k[i]);
    if (std::isfinite(vmf))
    {
        // printf("index %d \t i %d \t j %d \t jTheta %d \t jPhi %d \t theta %f \t phi %f \t vmf %f\n",
        //     index, i, j, jTheta, jPhi, theta, phi, vmf);
        vmm[index] = vmf;
        // printf("vmf %f\n", vmf);
    }
    else
    {
        vmm[index] = 0.0; // TODO: avoid this else -> cudaMemset() before calling kernel?
        // printf("vmf NOT FINITE!\n");
    }

    // if (std::isfinite(vmf))
    // {
    //     vmm[j] += w[i] * vmf;
    // }
}

__global__ void vmm3d_summation_kernel(double *vmm, int nSamples, //!! mu size needs to be 3 times the size of k, w
                                       const double *vmmAll, const double *w, const int mukwSz)
{
    int jndex = blockIdx.x * blockDim.x + threadIdx.x; //=blockIdx.x
    // int stride = blockDim.x * gridDim.x;

    // int totalSzIn = nrowsIn*ncols; //matrix is considered of size nrows*ncols, with nrows = sumNaturalsUpToN(numPts)

    // int numRows = mukwSz;

    // int numCols = 2 * nSamples * nSamples;
    // int sizePhis = 2 * nSamples;
    // int sizeThetas = nSamples;

    int jndexStart = jndex * mukwSz;

    for (int mukw = jndexStart; mukw < jndexStart + mukwSz; ++mukw)
    {
        // int iAll = jndex + numCols * mukw;
        double vmmij = vmmAll[mukw];
        vmm[jndex] += w[mukw - jndexStart] * vmmij;
    }
}

// int numCols = 2 * nSamples * nSamples;
// int sizePhis = 2 * nSamples; // TODO: this can maybe be passed directly
// int sizeThetas = nSamples;   // TODO: this can maybe be passed directly
// for (int idx = 0; idx < totalVmmSz; ++idx)
// {
//     int i = idx % mukwSz;
//     int j = (idx - i) / mukwSz;

//     if (i >= mukwSz)
//         printf("i >= mukwSz!! %d \t %d\n", i, mukwSz);

//     if (j >= numCols)
//         printf("j >= numCols!! %d \t %d\n", j, numCols);

//     int jPhi = j % sizePhis;
//     int jTheta = (j - jPhi) / sizePhis;

//     // printf("i %d \tj %d \tjPhi %d \t jTheta %d\n", i, j, jPhi, jTheta);

//     // double dtheta = M_PI / nSamples;
//     // double dtheta = M_PI / nSamples;

//     // int vmmSz = nSamples * (2 * nSamples);

//     // for
//     // double theta = dtheta * jTheta;
//     // int idx = itheta * nSamples + iphi; //idx should be j in this context
//     // double phi = dtheta * jPhi;

//     // int jPhi = j % sizeThetas;
//     // int jTheta = (j - jPhi) / sizeThetas;

//     // // double dtheta = M_PI / nSamples;
//     // double dtheta = M_PI / nSamples;

//     // int vmmSz = nSamples * (2 * nSamples);
//     // int vmmSz = nSamples * (2 * nSamples);

//     // for
//     // double theta = dtheta * jTheta;
//     // // int idx = itheta * nSamples + iphi; //idx should be j in this context
//     // double phi = dtheta * jPhi;

//     // double3 xThetaPhi;
//     // double st = sin(theta);
//     // double ct = cos(theta);
//     // double cp = cos(phi);
//     // double sp = sin(phi);
//     // xThetaPhi.x = st * cp;
//     // xThetaPhi.y = st * sp;
//     // xThetaPhi.z = ct;

//     double vmmHidx = vmmHost[idx];
//     // ROFL_VAR4(idx, i, j, vmmHidx);
//     if (std::isfinite(vmmHidx))
//         vmm[j] += weightsVmm[i] * vmmHidx;
// }

#endif /*DSD_UTILS_CUH_*/