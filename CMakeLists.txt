# 要求本机安装的 CMake 版本不低于 3.16。
cmake_minimum_required(VERSION 3.16)

# 定义项目名称 methane_calibration，并声明本项目使用 C++ 语言。
# methane_calibration 是项目名；CXX 在 CMake 中表示 C++。
project(methane_calibration LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 使用 C++17 标准编译本项目。
set(CMAKE_CXX_STANDARD 17)

# ON 表示必须支持 C++17；不允许编译器自动退回到更旧的 C++ 标准。
set(CMAKE_CXX_STANDARD_REQUIRED ON)




# 查找本机安装的 Eigen3 矩阵库。
# REQUIRED 表示如果没有找到 Eigen3，CMake 配置阶段立即报错。
find_package(Eigen3 REQUIRED)

# 把 src/main.cpp 编译成名为 methane_app 的可执行文件。
# methane_app 是前面创建项目时使用的基础练习程序。
# add_executable(methane_app src/main.cpp)

# 把标定程序的源文件编译成名为 methane_calibrate 的可执行文件。
add_executable(methane_calibrate
    calibration/app/calibrate_main.cpp
)

# 让 methane_calibrate 可以使用 Eigen3 的头文件和矩阵运算功能。
# 当前只学习到Wahba/SVD初始化阶段，暂时不链接Ceres。
target_link_libraries(methane_calibrate PRIVATE Eigen3::Eigen)
