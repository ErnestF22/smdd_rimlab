#include <iostream>

#include <rofl/common/macros.h>
#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

#include "csv_utils.h"
#include "dsd_utils.h"

void createCorridor2d(double xsize,
                      double ysize,
                      double xres,
                      int wallNum,
                      dsd::VectorVector2 &cloud);

void createCorridor3d(double xsize,
                      double ysize,
                      double zsize,
                      double xres,
                      double yres,
                      double zres,
                      int wallNum,
                      dsd::VectorVector3 &cloud);

void saveCsv2d(const std::string &filename, dsd::VectorVector2 &cloud);

void saveCsv3d(const std::string &filename, dsd::VectorVector3 &cloud);

int main(int argc, char **argv)
{
    std::string filenameCfg, filenameOut;
    std::string cloudFormat;
    int ndim, nwall;
    double xsize, ysize, zsize, xres, yres, zres;
    dsd::VectorVector2 cloud2;
    dsd::VectorVector3 cloud3;
    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    params.getParam<std::string>("out", filenameOut, std::string("cloud.csv"));
    params.getParam<std::string>("format", cloudFormat, std::string("csv"));
    params.getParam<int>("ndim", ndim, int(2));
    params.getParam<int>("nwall", nwall, int(2));
    params.getParam<double>("xsize", xsize, double(10.0));
    params.getParam<double>("ysize", ysize, double(1.0));
    params.getParam<double>("zsize", zsize, double(2.5));
    params.getParam<double>("xres", xres, double(0.05));
    params.getParam<double>("yres", yres, double(0.05));
    params.getParam<double>("zres", zres, double(0.10));

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    if (ndim == 2)
    {
        createCorridor2d(xsize, ysize, xres, nwall, cloud2);
        if (cloudFormat == "csv")
        {
            saveCsv2d(filenameOut, cloud2);
        }
        else
        {
            ROFL_ERR("invalid file format \"" << cloudFormat << "\"");
        }
    }
    else if (ndim == 3)
    {
        createCorridor3d(xsize, ysize, zsize, xres, yres, zres, nwall, cloud3);
        if (cloudFormat == "csv")
        {
            saveCsv3d(filenameOut, cloud3);
        }
        else
        {
            ROFL_ERR("invalid file format \"" << cloudFormat << "\"");
        }
    }

    return 0;
}

void createCorridor2d(double xsize,
                      double ysize,
                      double xres,
                      int wallNum,
                      dsd::VectorVector2 &cloud)
{
    dsd::Vector2 p;
    double xmin = -0.5 * xsize;
    double ymin = -0.5 * ysize;
    double ymax = 0.5 * ysize;
    int xnum = ceil(xsize / xres);
    cloud.clear();
    for (int ix = 0; ix < xnum; ++ix)
    {
        p << xmin + ix * xres, ymin;
        cloud.push_back(p);
        if (wallNum > 1)
        {
            p << xmin + ix * xres, ymax;
            cloud.push_back(p);
        }
    }
}

void createCorridor3d(double xsize,
                      double ysize,
                      double zsize,
                      double xres,
                      double yres,
                      double zres,
                      int wallNum,
                      dsd::VectorVector3 &cloud)
{
    dsd::Vector3 p;
    double xmin = -0.5 * xsize;
    double ymin = -0.5 * ysize;
    double ymax = 0.5 * ysize;
    double zmin = 0.0;
    double zmax = zsize;
    int xnum = ceil(xsize / xres);
    int ynum = ceil(ysize / yres);
    int znum = ceil(zsize / zres);
    cloud.clear();
    for (int ix = 0; ix < xnum; ++ix)
    {
        double x = xmin + ix * xres;

        for (int iz = 0; ix < xnum; ++ix)
        {
            double z = zmin + iz * zres;
            p << x, ymin, z;
            cloud.push_back(p);
            if (wallNum > 1)
            {
                p << x, ymax, z;
                cloud.push_back(p);
            }
        }

        if (wallNum > 2)
        {
            for (int iy = 0; iy < ynum; ++ix)
            {
                double y = ymin + iy * yres;
                p << x, y, zmin;
                cloud.push_back(p);
                if (wallNum > 3)
                {
                    p << x, y, zmax;
                    cloud.push_back(p);
                }
            }
        }
    }
}

void saveCsv2d(const std::string &filename, dsd::VectorVector2 &cloud)
{
    std::ofstream file(filename);
    file << "x, y\n";
    for (auto p : cloud)
    {
        file << p(0) << ", " << p(1) << "\n";
    }
    file.close();
}

void saveCsv3d(const std::string &filename, dsd::VectorVector3 &cloud)
{
    std::ofstream file(filename);
    file << "x, y, z\n";
    for (auto p : cloud)
    {
        file << p(0) << ", " << p(1) << ", " << p(2) << "\n";
    }
    file.close();
}