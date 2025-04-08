#include <rofl/common/param_map.h>

#include <eigen3/Eigen/Dense>
#include <limits>

#include <orientation_eigen.h>

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

    // TODO: number of for cycles can most likely test_fisher_matrix_gsl_eigen.cppbe limited

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
        ROFL_VAR1(i);
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

int main(int argc, char *argv[])
{
    /*Fisher’s information matrix is defined by the first derivatives
    of the ray-tracing function r with respect to t and θ*/

    rofl::ParamMap params;

    std::string filenameCfg;
    std::string filePath;
    float tx, ty, theta_robot;

    float sigma;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>("in", filePath,
                                 std::string("sample.csv"));
    params.getParam<float>("sigma", sigma, 1.0f);
    params.getParam<float>("tx", tx, 0.0f);
    params.getParam<float>("ty", ty, 0.0f);
    params.getParam<float>("theta_robot", theta_robot, 0.0f);

    // adapt tildes
    params.adaptTildeInPaths();
    params.getParam<std::string>("in", filePath, std::string("sample.csv"));

    std::cout << "-------\n"
              << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::string logPath = filePath;
    std::cout << "--LogPath: " << logPath << std::endl;

    if (!std::filesystem::exists(logPath))
    {
        std::cout << "Invalid log path!" << std::endl;
        return 0;
    }

    LogReader reader(false);
    reader.openFile(logPath);
    std::cout << "--Filesize: " << reader.getFileSize() << std::endl;
    std::cout << "--Timestamp: " << reader.getTimestamp() << std::endl;

    Scan scan_main;
    float angleIncMain = 0.0;
    while (true)
    {
        LogFrame frame;
        Scan scan;
        bool readCommRecv = reader.readFrame(FrameType::COMM_RECV, frame);
        if (readCommRecv)
        {
            frame.type = FrameType::COMM_RECV;
            Deserializer des(frame.data);
            uint32_t remoteIp = des.read32u();
            uint16_t remotePort = des.read16u();
            des.skip(4); // uint32_t recvBytes = des.read32u();

            uint32_t magicNum = des.read32u();
            uint16_t pckType = des.read16u();
            uint16_t version = des.read16u();

            if (magicNum != PCK_MAGIC_NUMBER)
            {
                throw("Wrong magic number!");
            }
            if (version > 2)
            {
                throw("Wrong version!");
            }

            if (pckType == 2 && version == 2)
            {
                // Vehicle pose
                float VehiclePoseX = des.read32f();
                float VehiclePoseY = des.read32f();
                float VehiclePoseT = des.read32f();
                float VehiclePoseCov_m00 = des.read32f();
                float VehiclePoseCov_m10 = des.read32f();
                float VehiclePoseCov_m20 = des.read32f();
                float VehiclePoseCov_m11 = des.read32f();
                float VehiclePoseCov_m21 = des.read32f();
                float VehiclePoseCov_m22 = des.read32f();

                // Vehicle velocity
                float VehicleVx = des.read32f();
                float VehicleVy = des.read32f();
                float VehicleOmega = des.read32f();
                float VehicleSpdCov_m00 = des.read32f();
                float VehicleSpdCov_m10 = des.read32f();
                float VehicleSpdCov_m20 = des.read32f();
                float VehicleSpdCov_m11 = des.read32f();
                float VehicleSpdCov_m21 = des.read32f();
                float VehicleSpdCov_m22 = des.read32f();

                // Scan info
                uint8_t scanner = des.read8u();
                uint32_t scanId = des.read32u();
                uint64_t timestamp = des.read64u();
                uint32_t locModes = des.read32u();
                float laser2VehicleX = des.read32f();
                tx = -laser2VehicleX;
                float laser2VehicleY = des.read32f();
                ty = -laser2VehicleY;
                float laser2VehicleT = des.read32f(); // suppose T stands for theta
                theta_robot = -laser2VehicleT;
                float laserHeight = des.read32f();
                float frequency = des.read32f();
                float angleMin = des.read32f();
                float angleInc = des.read32f();
                angleIncMain = angleInc;

                uint16_t numBeams = des.read16u();
                float rangeScale = des.read32f();

                // Beams
                std::vector<float> ranges;
                std::vector<point> points;
                ranges.resize(numBeams);
                points.resize(numBeams);
                for (size_t i = 0; i < numBeams; ++i)
                {
                    auto &r = ranges[i];
                    auto &p = points[i];
                    uint16_t rangeRead = des.read16u();
                    if (rangeRead == uint16_t(65535))
                    {
                        r = std::numeric_limits<float>::infinity(); // invalid
                                                                    // range
                        p.x = p.y = std::numeric_limits<float>::infinity();
                    }
                    else
                    {
                        r = float(rangeRead) * rangeScale;
                        float theta = float(i) * angleInc + angleMin;
                        p.x = r * cosf(theta);
                        p.y = r * sinf(theta);
                        p.phi = theta;
                    }
                }

                scan.VehicleOmega = VehicleOmega;
                scan.VehicleVx = VehicleVx;
                scan.VehicleVy = VehicleVy;
                scan.points = points;
                scan.ranges = ranges;
                scan.id = scanId;
            }
            scan_main = scan;
        }

        bool readCommSent = reader.readFrame(FrameType::COMM_SENT, frame);
        if (readCommSent)
        {
            Deserializer des(frame.data);
            des.seek(14);
            uint32_t magic = des.read32u();
            if (magic == PCK_MAGIC_NUMBER)
            {
                des.seek(18);
                uint16_t pckType = des.read16u();
                if (pckType == 8)
                {
                    des.seek(24);
                    uint8_t errLength = des.read8u();
                    des.skip(errLength);
                    scan.VehicleX_gt = des.read32f();
                    scan.VehicleY_gt = des.read32f();
                    scan.VehicleT_gt = des.read32f();
                }
            }
        }
        if (!readCommRecv && !readCommSent)
        {
            break;
        }
        std::cout << "Scan id: " << scan.id << " "
                  << "scanPose_gt: " << scan.VehicleX_gt << " " << scan.VehicleY_gt
                  << " " << scan.VehicleT_gt << "\t" << "scanV: " << scan.VehicleVx
                  << " " << scan.VehicleVy << " " << scan.VehicleOmega << " "
                  << std::endl;
        std::cout << "---------------------------------------------"
                  << std::endl;
    } // end of log reader

    std::vector<double> ranges(
        scan_main.ranges.begin(),
        scan_main.ranges
            .end()); // scan_main.ranges is actually a vector of floats

    int n = ranges.size();
    std::cout << "scan_main ranges size " << n << std::endl;

    // phis
    std::vector<double> phis(n, 0.0f);
    for (int i = 0; i < n; ++i)
    {
        double phi_i = scan_main.points[i].phi; // this would actually be phi
        phis[i] = phi_i;
    }

    // alphas
    // for (int i = 0; i < n; ++i)
    // {
    //     // TODO: fix here -> alpha is actually the direction of the normal to
    //     the surface at the sensed point float alpha_i =
    //     scan_main.points[i].phi; alphas[i] = alpha_i;
    // }
    // FILTER ORIENTATION FROM CENSI'S CSM
    // double *alpha, *cov0_alpha;
    double rho_robot = sqrt(tx * tx + ty * ty);
    std::vector<double> alphas(n, std::numeric_limits<double>::quiet_NaN()),
        cov0_alphas(n, 0.0f);
    int whsize = 5;
    // double theta0 = theta_robot;
    for (int i = 0; i < n; i++)
    {
        std::vector<double> ranges_part;
        std::vector<double> phis_part; //(phis.begin(), phis.begin() + 100);
        int jmin = std::max<int>(0, i - whsize);
        int jmax = std::min<int>(n - 1, i + whsize);
        for (int j = jmin; j <= jmax; ++j)
        {
            if (j != i && isValid(ranges[i]))
            {
                ROFL_VAR2(i, j);
                ranges_part.push_back(ranges[j]);
                phis_part.push_back(phis[j]);
            }
        }
        double theta0 = phis[i];
        double rho0 = ranges[i];
        double alpha = 42, cov0_alpha = 32;
        if (!isValid(rho0))
        {
            continue;
        }
        int wsize = ranges_part.size();
        filter_orientation(theta0, rho0, wsize, phis_part.data(),
                           ranges_part.data(), alpha, cov0_alpha);
        std::cout << "alpha " << alpha << std::endl;
        std::cout << "cov0_alpha " << cov0_alpha << std::endl;
        std::cout << "angleInc " << angleIncMain << std::endl;
        alphas[i] = alpha;
        cov0_alphas[i] = cov0_alpha;
    }

    Eigen::Matrix3f fimMat(Eigen::Matrix3f::Zero());
    computeFIM(fimMat, sigma, theta_robot, alphas, phis, ranges);

    return 0;
}
