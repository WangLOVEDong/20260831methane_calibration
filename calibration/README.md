# LiDAR 地图甲烷遥测设备外参标定

本目录实现一个独立的 **外参标定模块**。它根据 LiDAR 地图中选取的三维点，以及甲烷遥测设备对应的水平角、垂直角，求出该设备相对 LiDAR 地图坐标系的外参：

```text
旋转矩阵 R_map_from_sensor: 设备坐标系 -> LiDAR 地图坐标系
平移向量 t_map_from_sensor: 设备原点在 LiDAR 地图坐标系中的位置
```

每台甲烷设备必须分别运行一次标定。例如设备 A 使用自己的 CSV 和平移初值，设备 B 也使用自己的 CSV 和平移初值。两台设备不能混用同一份观测数据。

> 当前状态：已完成模拟数据上的 SVD 初值、Ceres/LM 优化和 YAML 导出。真实双设备标定、PCD 读取和三维甲烷融合仍未实现。

## 1. 目录说明

```text
calibration/
├── app/
│   └── calibrate_main.cpp          # 标定程序源码
├── data/
│   ├── simulated_device_a.csv      # 无噪声模拟观测
│   └── simulated_device_a_ground_truth.yaml
│                                    # 仅用于模拟验证的真值，真实标定不使用
├── output/                         # 程序运行时自动创建，保存最新 YAML 外参
└── README.md                       # 本使用教程
```

## 2. 坐标系与观测模型

设备坐标系约定为：

```text
x 轴：设备前方
y 轴：设备左侧
z 轴：设备上方
水平角：逆时针为正
垂直角：向上为正
```

一条水平角 `theta` 和垂直角 `phi` 对应设备坐标系中的单位方向：

```text
d_sensor = [cos(phi) cos(theta), cos(phi) sin(theta), sin(phi)]^T
```

甲烷遥测设备不提供距离。对于第 `i` 条观测，LiDAR 地图点 `P_i_map` 位于设备射线上：

```text
P_i_map = t_map_from_sensor + lambda_i * R_map_from_sensor * d_i_sensor
```

其中：

- `P_i_map`：CSV 给出的 LiDAR 地图点，单位米。
- `d_i_sensor`：由该行水平角、垂直角换算得到的设备单位方向。
- `R_map_from_sensor`：待求的 `3 x 3` 旋转矩阵。
- `t_map_from_sensor`：待求的 `3 x 1` 平移向量，表示设备原点的位置。
- `lambda_i`：该点沿射线的未知距离。程序通过单位化方向消除它，不需要设备测距。

## 3. 标定算法流程

程序的输入是人工提供的平移初值 `t0` 和 CSV 观测。计算步骤如下：

```text
1. 由 t0 和每个地图点计算地图单位方向 d_i_map。
2. 由水平角、垂直角计算设备单位方向 d_i_sensor。
3. 构造 H = sum(d_i_map * d_i_sensor^T)，通过 Wahba/SVD 求初始旋转 R0。
4. 将 R0 转为旋转向量 phi0。
5. 以 [phi_x, phi_y, phi_z, tx, ty, tz] 为六维待优化参数。
6. 对每条观测使用方向叉乘残差：
   r_i = d_i_sensor x normalize(R(phi)^T * (P_i_map - t))。
7. Ceres 使用中心差分计算数值雅可比，并使用 LM 算法迭代优化。
8. 输出最终 R*、t*、最终总平方残差 J，并写入 YAML。
```

残差长度接近两条方向夹角的弧度值。全部观测的总平方残差为：

```text
J = sum(||r_i||^2)
```

`J` 越小，说明这组外参与标定观测越一致；但真实数据中不应期待 `J = 0`，因为存在角度噪声、选点误差和 LiDAR 点误差。

## 4. 依赖与首次配置

程序使用 C++17、Eigen3 和 Ceres Solver。Ubuntu 上首次安装依赖时可执行：

```bash
sudo apt update
sudo apt install cmake g++ libeigen3-dev libceres-dev
```

在项目根目录执行 CMake 配置。项目根目录指包含 `CMakeLists.txt` 的目录：

```bash
cmake -S . -B build
```

参数含义：

- `-S .`：源代码目录是当前目录。
- `-B build`：将 CMake 生成的构建文件放进 `build/`。

配置成功后，编译标定程序：

```bash
cmake --build build --target methane_calibrate
```

生成的可执行文件是：

```text
build/methane_calibrate
```

## 5. CSV 输入格式

CSV 第一行必须是表头。每一行是一组“地图点 - 设备方向”对应关系：

```csv
id,x_m,y_m,z_m,horizontal_deg,vertical_deg
1,5.377098,-1.065153,0.551369,-25.0,-10.0
2,6.198768,0.176950,0.639457,-10.0,-10.0
```

列含义：

- `id`：观测编号。
- `x_m, y_m, z_m`：该标定目标在 LiDAR 地图坐标系中的坐标，单位米。
- `horizontal_deg`：甲烷设备的水平偏转角，单位度。
- `vertical_deg`：甲烷设备的垂直偏转角，单位度。

至少需要 3 条观测才能运行；建议采集 10 到 20 条。标定点应覆盖不同的水平角、垂直角和空间位置，避免所有点几乎共线、共面或位于相近方向。

## 6. 运行标定

程序的完整用法是：

```bash
./build/methane_calibrate <标定CSV> <tx0_m> <ty0_m> <tz0_m> [输出YAML路径]
```

参数说明：

- `<标定CSV>`：输入 CSV 路径。
- `<tx0_m> <ty0_m> <tz0_m>`：设备原点在 LiDAR 地图坐标系的平移初值 `t0`，单位米。
- `[输出YAML路径]`：可选参数。未填写时，程序根据 CSV 文件名自动生成输出路径。

使用无噪声模拟数据的示例：

```bash
./build/methane_calibrate calibration/data/simulated_device_a.csv 1.50 -1.10 1.10
```

使用带角度噪声的 CSV 时，为避免终端对特殊字符产生歧义，建议用双引号包住路径：

```bash
./build/methane_calibrate "calibration/data/噪声0.1°版本.csv" 1.10 -1.30 1.00
```

也可以显式指定外参输出文件，例如真实设备 A：

```bash
./build/methane_calibrate calibration/data/device_a.csv 1.20 -0.80 0.90 calibration/output/device_a_extrinsic.yaml
```

## 7. YAML 输出与覆盖规则

若不指定 `[输出YAML路径]`，程序使用下列规则：

```text
calibration/data/device_a.csv
    -> calibration/output/device_a_extrinsic.yaml
```

同一个 CSV 再次运行时，会覆盖同名 YAML。这是设计行为：`output/` 中保留的是该数据文件当前最新的标定结果，不会因重复运行积累大量过期结果。

程序自动创建 `calibration/output/`。YAML 包含：

```yaml
source_csv: "calibration/data/device_a.csv"
observation_count: 15
final_total_squared_residual: 0.000000000000

coordinate_convention:
  sensor_x: forward
  sensor_y: left
  sensor_z: up

rotation_vector_rad:
  - phi_x
  - phi_y
  - phi_z

rotation_map_from_sensor:
  rows: 3
  cols: 3
  data:
    - R00
    - R01
    # 其余元素按行优先顺序继续保存

translation_map_from_sensor_m:
  - tx
  - ty
  - tz
```

后续 `fusion` 模块将读取其中的 `rotation_map_from_sensor` 和 `translation_map_from_sensor_m`，把设备射线转换到 LiDAR 地图坐标系。

## 8. 两台设备的推荐操作

设备 A 与设备 B 的机械安装位置和姿态不同，因此必须独立标定：

```text
设备 A 的 CSV + 设备 A 的 t0
    -> device_a_extrinsic.yaml

设备 B 的 CSV + 设备 B 的 t0
    -> device_b_extrinsic.yaml
```

不要把两台设备的观测混到一个 CSV 中，也不要让设备 B 覆盖设备 A 的 YAML 文件。

## 9. 模拟数据验证

`simulated_device_a.csv` 是无噪声模拟数据。`simulated_device_a_ground_truth.yaml` 只用于验证程序是否恢复了模拟真值，它不是实际标定的输入。

项目还包含不同角度噪声版本的 CSV。比较它们时，可以观察：

- 优化前后 `J` 是否显著下降。
- 噪声增大后，最终 `J` 是否整体增大。
- 仅在模拟环境中，才可用 ground truth 计算 `R*`、`t*` 与真值的误差。

真实数据没有 ground truth 时，应额外准备未参与标定的验证点，检查预测方向与实际方向的夹角误差。

## 10. 常见问题

### CMake 报错：找不到 Eigen 或 Ceres

确认已安装第 4 节的依赖，然后重新执行：

```bash
cmake -S . -B build
```

### 移动或重命名项目目录后，CMake 报路径错误

旧的 `build/` 中保存了旧目录路径。将旧构建目录改名或删除后，再重新配置：

```bash
mv build build_before_move
cmake -S . -B build
```

### 最终残差较大

依次检查：平移初值 `t0` 是否接近设备原点、CSV 坐标系是否与程序约定一致、水平角和垂直角符号是否正确、标定点是否分布充分、以及是否混入了两台设备的数据。

## 11. 下一模块

下一步将在项目的 `fusion/` 目录实现三维融合。其输入是：

```text
PCD 点云 / LiDAR 地图
设备 A 的外参 YAML
设备 B 的外参 YAML
两台甲烷设备随时间变化的角度与浓度数据
```

融合模块的职责是读取外参、将甲烷设备的射线变换到地图坐标系，并将浓度信息关联到三维点云或体素；它不应重复执行本目录的外参标定。
