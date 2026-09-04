#include <Eigen/Dense> // Eigen::Matrix3d, Eigen::Vector3d
#include <ceres/rotation.h> // ceres::RotationMatrixToAngleAxis

#include <algorithm>  // std::replace
#include <cmath>      // std::cos, std::sin
#include <cstdio>     // std::printf
#include <fstream>    // std::ifstream
#include <sstream>    // std::istringstream
#include <string>     // std::string, std::getline
#include <vector>     // std::vector

// Direction3D 用三个分量保存一个三维方向。
struct Direction3D
{
    double x;
    double y;
    double z;
};

// CalibrationObservation 保存CSV中的一行标定观测，不包含表头。
struct CalibrationObservation
{
    int id;
    double x_m;
    double y_m;
    double z_m;
    double horizontal_deg;
    double vertical_deg;
};

constexpr double kPi = 3.14159265358979323846;       //这里定义的kPi是 很精确没办法修改的量

// 把设备输出的水平角、垂直角转换成设备坐标系中的单位方向。
// 参数单位均为度；返回方向使用 x向前、y向左、z向上的约定。右手定则
Direction3D anglesToBearing(double horizontal_deg, double vertical_deg)
{
    const double theta_rad = horizontal_deg * kPi / 180.0;        //水平角转弧度
    const double phi_rad = vertical_deg * kPi / 180.0;            //垂直角转弧度
    const double cos_phi = std::cos(phi_rad);

    return Direction3D{
        cos_phi * std::cos(theta_rad),            //x方向
        cos_phi * std::sin(theta_rad),            //y方向
        std::sin(phi_rad),                         //z方向
    };
}

// 解析一行简单的逗号分隔CSV；成功时填充 observation 并返回 true。
/*在 csv_line 的全部字符中  把逗号 ',' 替换为空格 ' '*/
bool parseObservation(std::string csv_line, CalibrationObservation& observation)
{
    /* 遍历一下，csv中每一行的逗号 ',' 替换为空格 ' ' */
    std::replace(csv_line.begin(), csv_line.end(), ',', ' ');   //replace直接修改原字符串 `csv_line`，不会生成新字符串副本

    std::istringstream row_stream(csv_line); //istringstream 是从字符串中读取数据，和cin 从键盘读取数据一样，跳过空格

    //相当 从 row_stream 读取一个数，放到 observation
    return static_cast<bool>(
        row_stream >> observation.id
                   >> observation.x_m
                   >> observation.y_m
                   >> observation.z_m
                   >> observation.horizontal_deg
                   >> observation.vertical_deg);
}

// 第2～4步：根据t0计算地图方向，根据设备角度计算设备方向，最后通过Wahba/SVD求R0。
// observations：CSV中的全部原始观测。
// translation_initial：设备原点在地图坐标系中的平移初值t0，单位为米。
// 返回值：设备坐标系到LiDAR地图坐标系的初始旋转矩阵R0。
Eigen::Matrix3d estimateInitialRotation(
    const std::vector<CalibrationObservation>& observations,
    const Eigen::Vector3d& translation_initial)
{
    Eigen::Matrix3d direction_correlation = Eigen::Matrix3d::Zero();

     
    for (const CalibrationObservation& observation : observations)
    {
        // 第2步：由地图点P_i^M和设备位置初值t0计算地图单位方向d_i^M。
        const Eigen::Vector3d point_map(
            observation.x_m,
            observation.y_m,
            observation.z_m);
        const Eigen::Vector3d bearing_map =
            (point_map - translation_initial).normalized();

        // 第3步：由水平角和垂直角计算设备坐标系单位方向d_i^S。
        const Direction3D bearing_sensor_components = anglesToBearing(
            observation.horizontal_deg,
            observation.vertical_deg);
        const Eigen::Vector3d bearing_sensor(
            bearing_sensor_components.x,
            bearing_sensor_components.y,
            bearing_sensor_components.z);

        // 第4步：构造了H  累加H = sum_i(d_i^M * (d_i^S)^T)。其中很多部分省略了
        direction_correlation += bearing_map * bearing_sensor.transpose();
    }
    // 遍历每一行观测，计算地图方向和设备方向，并累加到构造 H矩阵中。

    
    // 对H进行SVD分解：H = U * Sigma * V^T。
    /*第一个参数是待分解的矩阵，第二个参数是分解选项，这里表示计算完整的U和V矩阵*/
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
        direction_correlation,
        Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Eigen::Matrix3d matrix_u = svd.matrixU();
    const Eigen::Matrix3d matrix_v = svd.matrixV();

    // 保证det(R0)=+1，避免得到镜像矩阵。
    Eigen::Matrix3d determinant_correction = Eigen::Matrix3d::Identity();
    if ((matrix_u * matrix_v.transpose()).determinant() < 0.0)
    {
        determinant_correction(2, 2) = -1.0;
    }

    return matrix_u * determinant_correction * matrix_v.transpose();
}

/*
    计算一条观测的三维方向残差 r_i。
    observation：CSV中的一行地图点和水平/垂直角观测。
    rotation_map_from_sensor：候选旋转矩阵R，负责把设备坐标系方向转到地图坐标系。
    translation_map_from_sensor：候选平移向量t，表示设备原点在地图坐标系中的位置。
    返回值：r_i = d_i^S x d_i^S_predicted；两条单位射线重合时返回零向量。
*/
Eigen::Vector3d calculateBearingResidual(
    const CalibrationObservation& observation,
    const Eigen::Matrix3d& rotation_map_from_sensor,
    const Eigen::Vector3d& translation_map_from_sensor)
{
    // CSV水平角、垂直角对应的实测设备单位方向d_i^S。
    const Direction3D measured_direction_components = anglesToBearing(
        observation.horizontal_deg,
        observation.vertical_deg);
    const Eigen::Vector3d measured_bearing_sensor(
        measured_direction_components.x,
        measured_direction_components.y,
        measured_direction_components.z);

    // CSV中的LiDAR地图点P_i^M。
    const Eigen::Vector3d point_map(
        observation.x_m,
        observation.y_m,
        observation.z_m);

    // R负责“设备到地图”，因此用R^T把(P_i^M-t)转换回设备坐标系。
    const Eigen::Vector3d point_sensor =
        rotation_map_from_sensor.transpose() *
        (point_map - translation_map_from_sensor);

    // 去掉未知距离lambda_i，得到候选外参预测的设备单位方向d_i^S_predicted。
    const Eigen::Vector3d predicted_bearing_sensor =
        point_sensor.normalized();

    // 测量方向和预测方向的叉积作为残差，重合时返回零向量。
    return measured_bearing_sensor.cross(predicted_bearing_sensor);
}

// 对全部观测累加总平方残差 J(R, t) = sum_i ||r_i(R, t)||^2。
// observations：CSV中的全部观测；R和t：同一组待评价的候选外参。
// 返回值：标量J。J越小，说明这组外参与全部方向观测越一致。
double calculateTotalSquaredResidual(
    const std::vector<CalibrationObservation>& observations,
    const Eigen::Matrix3d& rotation_map_from_sensor,
    const Eigen::Vector3d& translation_map_from_sensor)
{
    double total_squared_residual = 0.0;

    for (const CalibrationObservation& observation : observations)
    {
        const Eigen::Vector3d residual = calculateBearingResidual(
            observation,
            rotation_map_from_sensor,
            translation_map_from_sensor);

        // squaredNorm() = r_x^2 + r_y^2 + r_z^2，也就是||r_i||^2。
        total_squared_residual += residual.squaredNorm();
    }

    return total_squared_residual;
}


int main(int argc, char* argv[])
{
    // argv[0] 是程序名；另外需要CSV路径和人工测得的tx0、ty0、tz0。
    if (argc != 5)
    {
        std::printf("Usage: %s <calibration_csv> <tx0_m> <ty0_m> <tz0_m>\n", argv[0]);
        return 1;
    }

    // argv[0]：程序名  argv[1]：CSV路径；argv[2..4]：平移初值的三个分量，单位为米。
    const std::string csv_path = argv[1];
    double tx_initial_m = 0.0;
    double ty_initial_m = 0.0;
    double tz_initial_m = 0.0;

    if (std::sscanf(argv[2], "%lf", &tx_initial_m) != 1 ||
        std::sscanf(argv[3], "%lf", &ty_initial_m) != 1 ||
        std::sscanf(argv[4], "%lf", &tz_initial_m) != 1)
    {
        std::printf("Invalid translation initial value. tx0, ty0, tz0 must be numbers.\n");
        return 1;
    }

    // 构造平移初值向量 t0。程序读入的 设备坐标初值 传给 translation_initial
    const Eigen::Vector3d translation_initial(
        tx_initial_m,
        ty_initial_m,
        tz_initial_m);

    std::ifstream input_file(csv_path);
    if (!input_file.is_open())        //判断文件是否成功打开，如果没有成功打开，就打印错误信息并返回错误码 1
    {
        /* csv_path.c_str()：把 C++ std::string 转成 C 语言 const char*，因为 printf 只接受 C 风格字符串。*/
        std::printf("Failed to open CSV: %s\n", csv_path.c_str());
        return 1;
    }

    // 读取表头行（跳过列名）
    std::string header;
    if (!std::getline(input_file, header))   // 一次性读取一整行，直到换行符为止
    {
        std::printf("CSV is empty: %s\n", csv_path.c_str());
        return 1;
    }

    std::vector<CalibrationObservation> observations;     //csv 读取的数据都放在 observations 中
    std::string line;
    int line_number = 1;

    while (std::getline(input_file, line))
    {
        ++line_number;
        if (line.empty())   //遇到空行就跳过
        {
            continue;
        }

        CalibrationObservation observation{};

        /*一边执行函数parseObservation判断，一边填充 observation ，填充 地图坐标系xyz 与 设备水平/垂直角度  */
        if (!parseObservation(line, observation))    //某一行 CSV 解析失败，就打印错误信息并返回错误码 1
        {
            std::printf("Invalid CSV row at line %d: %s\n", line_number, line.c_str());
            return 1;
        }
        observations.push_back(observation);
    }

    //填充完 observations 后，检查是否为空，如果为空就打印错误信息并返回错误码 1
    if (observations.empty())
    {
        std::printf("CSV contains no calibration observations: %s\n", csv_path.c_str());
        return 1;
    }
    if (observations.size() < 3)
    {
        std::printf("At least 3 calibration observations are required.\n");
        return 1;
    }

    // 第2～4步：根据全部点—方向对应关系和 初始偏差t0 ，通过SVD求初始旋转矩阵R0。
    const Eigen::Matrix3d rotation_initial = estimateInitialRotation(
        observations,
        translation_initial);

    /*
        第5步：把SVD得到的3x3初始旋转矩阵R0转换为Ceres使用的旋转向量phi0。
        rotation_initial.data() 指向R0的9个double元素；rotation_vector_initial.data()
        提供3个double元素，接收 phi0 = theta0 * a0，单位为弧度。
    */ 
    Eigen::Vector3d rotation_vector_initial;   //定义旋转向量 phi0
    //把第一个参数 旋转矩阵R0 转换为旋转向量 phi0，存储在 rotation_vector_initial 中
    ceres::RotationMatrixToAngleAxis(
        rotation_initial.data(),
        rotation_vector_initial.data());

    // 第6步的初值：前三个元素是旋转向量phi0（弧度），后三个元素是平移t0（米）。
    // 下一阶段Ceres会直接修改这个数组，因此这里不能写成const。
    double extrinsic_parameters[6] = {
        rotation_vector_initial.x(), // [0] = phi_x0
        rotation_vector_initial.y(), // [1] = phi_y0
        rotation_vector_initial.z(), // [2] = phi_z0
        translation_initial.x(),     // [3] = tx0
        translation_initial.y(),     // [4] = ty0
        translation_initial.z()      // [5] = tz0
    };

    // 使用当前初值x0，对CSV第一行观测计算一次残差，暂时不进行优化。
    const Eigen::Vector3d first_initial_residual =
        calculateBearingResidual(
            observations.front(),
            rotation_initial,
            translation_initial);

    const double initial_total_squared_residual =
        calculateTotalSquaredResidual(
            observations,
            rotation_initial,
            translation_initial);

    const double rotation_angle_initial_rad = rotation_vector_initial.norm();  //计算旋转角度

    // 下面的逻辑 是 如果旋转角度大于一个很小的阈值，就计算旋转轴 a0 = phi0 / theta0，
    // 否则旋转轴为零向量，相当于任意方向旋转了
    Eigen::Vector3d rotation_axis_initial = Eigen::Vector3d::Zero();        
    if (rotation_angle_initial_rad > 1e-12)
    {
        rotation_axis_initial = rotation_vector_initial / rotation_angle_initial_rad;
    }

    

    std::printf("Calibration observations: %zu\n", observations.size());
    std::printf("Initial translation t0: [%.6f, %.6f, %.6f] m\n",
                translation_initial.x(),
                translation_initial.y(),
                translation_initial.z());
    std::printf("Initial rotation matrix R0_map_from_sensor:\n");
    std::printf("[%.9f %.9f %.9f]\n",
                rotation_initial(0, 0),
                rotation_initial(0, 1),
                rotation_initial(0, 2));
    std::printf("[%.9f %.9f %.9f]\n",
                rotation_initial(1, 0),
                rotation_initial(1, 1),
                rotation_initial(1, 2));
    std::printf("[%.9f %.9f %.9f]\n",
                rotation_initial(2, 0),
                rotation_initial(2, 1),
                rotation_initial(2, 2));
    std::printf("Initial rotation vector phi0: [%.9f, %.9f, %.9f] rad\n",
                rotation_vector_initial.x(),
                rotation_vector_initial.y(),
                rotation_vector_initial.z());
    std::printf("Initial rotation angle theta0: %.9f rad (%.6f deg)\n",
                rotation_angle_initial_rad,
                rotation_angle_initial_rad * 180.0 / kPi);
    std::printf("Initial rotation axis a0: [%.9f, %.9f, %.9f]\n",
                rotation_axis_initial.x(),
                rotation_axis_initial.y(),
                rotation_axis_initial.z());
    std::printf("Initial parameters [phi_x, phi_y, phi_z, tx, ty, tz]:\n");
    std::printf("[%.9f, %.9f, %.9f, %.6f, %.6f, %.6f]\n",
                extrinsic_parameters[0],
                extrinsic_parameters[1],
                extrinsic_parameters[2],
                extrinsic_parameters[3],
                extrinsic_parameters[4],
                extrinsic_parameters[5]);
    std::printf("First initial bearing residual r1: [%.9f, %.9f, %.9f]\n",
                first_initial_residual.x(),
                first_initial_residual.y(),
                first_initial_residual.z());
    std::printf("Initial total squared residual J: %.12f\n",
                initial_total_squared_residual);

    return 0;
}





