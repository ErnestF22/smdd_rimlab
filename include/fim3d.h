#ifndef FIM_3D_H_
#define FIM_3D_H_

#include <pcl/common/geometry.h>
#include <pcl/features/normal_3d.h>

#include <rofl/common/macros.h>
#include <rofl/common/profiler.h>
#include <eigen3/Eigen/Dense>

#include <dsd_utils.h>

void cpm(const Eigen::Vector3f &v, Eigen::Matrix3f &m)
{
    m = Eigen::Matrix3f::Zero();
    m << 0.0f, -v(2), v(1), v(2), 0.0f, -v(0), -v(1), v(0), 0.0f;
}

void computeFIM3D(dsd::Matrix6 &fim,
                  pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud,
                  pcl::PointCloud<pcl::Normal>::Ptr &normals,
                  double sigma)
{
    dsd::Vector3 nI, pI;
    dsd::Vector6 grad;
    double cosI, rhoI;
    size_t counter;

    rofl::ScopedTimer timer("FIM_3D");

    ROFL_ASSERT_VAR2(cloud->size() == normals->size(), cloud->size(),
                     normals->size());

    fim.setZero();
    counter = 0;
    for (int i = 0; i < cloud->size(); ++i)
    {
        if (!std::isfinite(normals->at(i).normal_x) ||
            !std::isfinite(normals->at(i).normal_y) ||
            !std::isfinite(normals->at(i).normal_z))
        {
            continue;
        }
        nI << normals->at(i).normal_x, normals->at(i).normal_y,
            normals->at(i).normal_z;
        pI << cloud->at(i).x, cloud->at(i).y, cloud->at(i).z;
        rhoI = pI.norm();
        cosI = nI.dot(pI) / rhoI;

        if (fabs(cosI) > 0.10)
        {
            grad.head<3>() = -pI.cross(nI);
            grad.tail<3>(3) = -nI;
            grad = grad * 1.0 / cosI;
            fim += grad * grad.transpose() / (sigma * sigma);
            counter++;
        }
    }
    if (counter > 0)
    {
        fim *= 1.0 / counter;
    }
}

#endif