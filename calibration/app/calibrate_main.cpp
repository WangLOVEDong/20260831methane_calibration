#include <Eigen/Dense> // Eigen::Matrix3d, Eigen::Vector3d

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
    Direction3D bearing_sensor;    //射线方向
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

/*输入参数：
    observation ： CSV中一行标定观测数据。包括 ：地图坐标系xyz 与 设备水平/垂直角度
    rotation_map_from_sensor: 将设备坐标系向量旋转到地图坐标系的候选旋转矩阵 R。
    translation_map_from_sensor: 设备原点在地图坐标系中的候选位置 t，单位为米。
*/ 
Eigen::Vector3d calculateBearingResidual(
    const CalibrationObservation& observation,
    const Eigen::Matrix3d& rotation_map_from_sensor,
    const Eigen::Vector3d& translation_map_from_sensor)
{
    // 遥测设备角度换算出的实测射线方向 d^S。
    const Eigen::Vector3d measured_bearing_sensor(
        observation.bearing_sensor.x,
        observation.bearing_sensor.y,
        observation.bearing_sensor.z);

    // CSV给出的LiDAR地图点 P^M。
    const Eigen::Vector3d point_map(
        observation.x_m,
        observation.y_m,
        observation.z_m);

    // R负责“设备到地图”，所以这里用R的转置，把(P^M-t)变回设备坐标系。
    const Eigen::Vector3d point_sensor =
        rotation_map_from_sensor.transpose() *
        (point_map - translation_map_from_sensor);

    /*去掉未知距离，只保留候选外参预测出的单位射线方向。
      .normalized() 返回**原向量除以自身二范数 (模长)**之后得到的**新的单位向量**。     
      此时得到的 predicted_bearing_sensor 就是估计 d^S 射线向量 */ 
    const Eigen::Vector3d predicted_bearing_sensor = point_sensor.normalized();    
    
    /* 然后实际观测的射线向量和 上面估计的射线向量求叉乘，得到残差  如果两个向量完全重合，叉乘结果为零  */
    return measured_bearing_sensor.cross(predicted_bearing_sensor);
}


int main(int argc, char* argv[])
{
    // argv[0] 是程序名，因此 argc == 2 表示用户额外传入了一个CSV路径。
    if (argc != 2)         // 如果用户没有传入参数，argc 就不等于 2，就会返回错误码 1
    {
        std::printf("Usage: %s <calibration_csv>\n", argv[0]);
        return 1;
    }
    // argv[0]：你运行的**可执行程序名字**，永远自带
    // argv[1]：你运行的**可执行程序的第一个参数**，就是你在命令行中输入的第一个参数
    /*  此时
        argc = 2
        argv[0] = "./calib_parser"
        argv[1] = "/home/data/calibration.csv"
    */
    const std::string csv_path = argv[1];    //csv路径 例如：/home/data/calibration.csv

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
        // anglesToBearing 根据水平角度和垂直角度计算射线方向
        observation.bearing_sensor = anglesToBearing(observation.horizontal_deg, observation.vertical_deg);
        observations.push_back(observation);
    }

    //填充完 observations 后，检查是否为空，如果为空就打印错误信息并返回错误码 1
    if (observations.empty())
    {
        std::printf("CSV contains no calibration observations: %s\n", csv_path.c_str());
        return 1;
    }

    //打印第一行的  统计数据
    const CalibrationObservation first = observations.front();
    std::printf("Parsed calibration rows: %zu\n", observations.size());
    std::printf("First map point: [%.6f, %.6f, %.6f] m\n",
                first.x_m, first.y_m, first.z_m);
    std::printf("First angles: horizontal=%.2f deg, vertical=%.2f deg\n",
                first.horizontal_deg, first.vertical_deg);
    std::printf("First bearing d^S: [%.6f, %.6f, %.6f]\n",
                first.bearing_sensor.x,
                first.bearing_sensor.y,
                first.bearing_sensor.z);

    /*初始化搜索：
        给出 旋转角度偏差（0，0，-20到20）只考虑 z轴旋转，其他两个轴不考虑。
        给出 平移偏差（-0.2到0.2）只考虑 xyz平移。
        根据给出的旋转角度求出旋转矩阵R
    */

    // 先用一个最简单的候选值验证残差代码：R0为单位旋转，t0为零平移。
    // 下一步再把这里替换成上面计划的角度和平移搜索。
    const Eigen::Matrix3d rotation_initial = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d translation_initial = Eigen::Vector3d::Zero();

    const Eigen::Vector3d first_residual = calculateBearingResidual(
        first,
        rotation_initial,
        translation_initial);

    std::printf("First residual: [%.6f, %.6f, %.6f]\n",
                first_residual.x(),
                first_residual.y(),
                first_residual.z());
    std::printf("First residual norm: %.6f\n", first_residual.norm());

    return 0;
}
