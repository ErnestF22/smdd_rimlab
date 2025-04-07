
#include "dsd_utils.h"

#include <fstream>

namespace fs = std::filesystem;

class BinReader {
   public:
    enum class LidarType { VELODYNE, OUSTER, LIVOX };

    const size_t CLOUD_BIN_VELODYNE_SIZE = 22;
    const size_t CLOUD_BIN_OUSTER_SIZE = 30;
    const size_t CLOUD_BIN_LIVOX_SIZE = 18;

    struct CloudBinVelodyne {
        float x;          // 0-3
        float y;          // 4-7
        float z;          // 8-11
        float intensity;  // 12-15
        uint16_t ring;    // 16-17
        float timestamp;  // 18-21
    };

    struct CloudBinOuster {
        float x;                // 0-3
        float y;                // 4-7
        float z;                // 8-11
        float intensity;        // 12-15
        uint64_t timestamp;     // 16-19
        uint16_t reflectivity;  // 20-21
        uint16_t ring;          // 22-23
        uint16_t ambient;       // 24-25
        uint64_t range;         // 26-29
    };

    struct CloudBinLivox {
        float x;          // 0-3
        float y;          // 4-7
        float z;          // 8-11
        float intensity;  // 12-15
        uint8_t tag;      // 16
        uint8_t line;     // 17
    };

    int readCloudBin(const std::string& infile, LidarType type) {
        dsd::Vector3 point;
        points_.clear();
        std::fstream input(infile.c_str(), std::ios::in | std::ios::binary);
        if (!input.good()) {
            ROFL_ERR("Could not read file: " << infile);
            return -1;
        } else {
            ROFL_MSG("reading " << infile);
            input.seekg(0, std::ios::beg);

            // ROFL_VAR3(sizeof(CloudBinVelodyne), sizeof(CloudBinOuster),
            //           sizeof(CloudBinLivox));
            // ROFL_VAR4(sizeof(float), sizeof(uint16_t), sizeof(uint32_t),
            //           sizeof(uint64_t))
            for (int i = 0; input.good() && !input.eof(); i++) {
                if (type == LidarType::VELODYNE) {
                    CloudBinVelodyne record;
                    input.read((char*)&record, CLOUD_BIN_VELODYNE_SIZE);
                    // ROFL_VAR7(i, record.x, record.y, record.z,
                    // record.intensity,
                    //           record.ring, record.timestamp);
                    point(0) = record.x;
                    point(1) = record.y;
                    point(2) = record.z;
                } else if (type == LidarType::OUSTER) {
                    CloudBinOuster record;
                    input.read((char*)&record, CLOUD_BIN_OUSTER_SIZE);
                    point(0) = record.x;
                    point(1) = record.y;
                    point(2) = record.z;
                } else if (type == LidarType::LIVOX) {
                    CloudBinLivox record;
                    input.read((char*)&record, CLOUD_BIN_LIVOX_SIZE);
                    point(0) = record.x;
                    point(1) = record.y;
                    point(2) = record.z;
                }
                // values[3] is the intensity of the beam
                if (std::isfinite(point(0)) && std::isfinite(point(1)) &&
                    std::isfinite(point(2)) && !isVehicleShapePt(point))
                    points_.push_back(point);
                // ROFL_VAR2(i, points_.size());
            }
            input.close();
            return points_.size();
        }
    }

    BinReader() {}

    ~BinReader() {}

    void setVehiclePtsMinMax(const dsd::Vector2 xVeh,
                             dsd::Vector2 yVeh,
                             dsd::Vector2 zVeh) {
        xVeh_ = xVeh;
        yVeh_ = yVeh;
        zVeh_ = zVeh;
    }

    /**
     * x, y, z contain minimum coordinate for each axis at position 0 and,
     * respectively, maximum coordinate at position 1
     */
    bool isVehicleShapePt(const dsd::Vector3& pt) {
        if (pt(0) < xVeh_(0) || pt(0) > xVeh_(1))
            return false;

        if (pt(1) < yVeh_(0) || pt(1) > yVeh_(1))
            return false;

        if (pt(2) < zVeh_(0) || pt(2) > zVeh_(1))
            return false;

        return true;
    }

    const dsd::VectorVector3& getCloud() const { return points_; }

    dsd::VectorVector3 points_;
    dsd::Vector2 xVeh_;
    dsd::Vector2 yVeh_;
    dsd::Vector2 zVeh_;
};