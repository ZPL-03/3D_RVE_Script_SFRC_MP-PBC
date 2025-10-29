# 3D RVE Model for Short Fiber Reinforced Composites (SFRC)

## English

### Overview

This script provides a comprehensive solution for automated parametric modeling of Representative Volume Elements (RVE) for Short Fiber Reinforced Composites (SFRC). It implements the classical **Random Sequential Adsorption (RSA)** algorithm with advanced optimizations including spatial grid acceleration and intelligent mirroring strategies for periodic boundary conditions.

The script covers the complete workflow from geometry creation, material definition, and mesh generation to the application of periodic boundary conditions (PBC) using the ABAQUS Micromechanics Plugin.

### Key Features

#### 🚀 **Efficient Fiber Distribution Generation**
- **Spatial Grid Acceleration**: Reduces interference detection complexity from O(N²) to approximately O(N)
- **Analytical Interference Detection**: Precise 27-box method ensures minimum fiber spacing
- **Performance Statistics**: Detailed reports on spatial grid search efficiency
- **Controllable Volume Fraction**: Target fiber volume fraction with CSV export capability

#### 🔄 **Periodic Geometry Modeling**
- **Smart Mirroring Strategy**: Creates mirror copies only for fibers near boundaries in the required directions
- **Boolean Operations**: Precise fiber-matrix cutting and merging
- **Unified Boundary Trimming**: Ensures geometric periodicity of the final model

#### 🎯 **Material Property Definition**
- **Orthotropic Fiber Material**: Complete specification with E1, E2, E3, ν12, ν13, ν23, G12, G13, G23
- **Isotropic Matrix Material**: Young's modulus and Poisson's ratio
- **Intelligent Cell Classification**: Geometry-based automatic material assignment
- **Per-Fiber Local Coordinate Systems**: Accurate material orientation for each fiber

#### 🌐 **Mesh Generation**
- **Parametric Mesh Control**: Automated 3D solid element generation
- **Robust Fallback Strategy**: Attempts C3D10 (quadratic tetrahedral) first, falls back to C3D4 (linear tetrahedral) if needed
- **Quality Assurance**: Comprehensive error handling and reporting

#### 📐 **Periodic Boundary Conditions**
- **Micromechanics Plugin Integration**: Automatic PBC application for non-periodic meshes
- **Complete Setup**: Creates reference points, constraint equations, analysis steps, and loads
- **Homogenization Ready**: Prepares model for subsequent material property homogenization

### Algorithm Details

#### RSA with Spatial Grid Acceleration

The script implements an enhanced RSA algorithm:

1. **Random Generation**: Completely random fiber positions and orientations
2. **Sequential Placement**: Attempts to place fibers one by one
3. **Interference Check**: Validates against existing fibers using spatial grid
4. **Accept/Reject**: Accepts if no interference, otherwise generates new position
5. **Grid Optimization**: O(N) complexity for most operations instead of O(N²)

#### Smart Mirroring Strategy

For fibers near boundaries:
- Analyzes which boundaries the fiber is close to
- Creates mirror copies only in necessary directions
- Significantly reduces geometry count and improves performance
- Maintains full periodicity with minimal computational overhead

### Requirements

#### Software Requirements
- **ABAQUS**: Version 2023 or higher
- **Micromechanics Plugin**: Must be installed and enabled in ABAQUS/CAE
- **Python**: 2.7 (built-in with ABAQUS)

#### Plugin Installation
The Micromechanics Plugin must be installed in ABAQUS. The script expects it at:
```python
d:/Abaqus2023/Plugins
```
Modify the path in the script if your installation differs.

### Installation

1. Clone this repository:
```bash
git clone https://github.com/ZPL-03/3D_RVE_Model_SFRC_MP-PBC.git
cd 3D_RVE_Model_SFRC_MP-PBC
```

2. Ensure the Micromechanics Plugin is properly installed in your ABAQUS installation.

3. Open ABAQUS/CAE and load the script via File → Run Script.

### Usage

#### Basic Usage

1. **Configure Parameters**: Open `3D_RVE_Model_SFRC_MP-PBC.py` and modify the parameters in the configuration section:

```python
# RVE Geometry (mm)
LENGTH = 0.5
WIDTH = 0.5
HEIGHT = 0.5

# Fiber Geometry (mm)
FIBRE_LENGTH = 0.25
FIBRE_DIAMETER = 0.025

# Fiber Distribution
MIN_DISTANCE_FACTOR = 1.5
TARGET_VF = 0.05  # Target volume fraction (5%)
MAX_ITERATIONS = 100000

# Material Properties (MPa)
FIBER_E1 = 230000.0  # Longitudinal Young's modulus
MATRIX_E = 3500.0    # Matrix Young's modulus
```

2. **Run the Script**: In ABAQUS/CAE:
   - File → Run Script → Select `3D_RVE_Model_SFRC_MP-PBC.py`
   - Or via command line:
     ```bash
     abaqus cae noGUI=3D_RVE_Model_SFRC_MP-PBC.py
     ```

3. **Monitor Progress**: The script provides detailed progress updates:
   - Fiber placement statistics
   - Geometry creation steps
   - Material assignment progress
   - Mesh generation status
   - PBC application confirmation

4. **Review Results**: After completion:
   - Model name: `RVE-VfXXX` (e.g., `RVE-Vf005` for 5% volume fraction)
   - Job name: `Job-1` (ready for submission)
   - CSV file: Fiber center coordinates (if enabled)

#### Advanced Configuration

**Mesh Control:**
```python
GLOBAL_SEED_SIZE = 0.01      # Global seed size (mm)
DEVIATION_FACTOR = 0.05      # Curvature deviation
MIN_SIZE_FACTOR = 0.05       # Minimum element size factor
```

**Boundary Threshold:**
```python
BOUNDARY_THRESHOLD_FACTOR = 0.5  # Multiple of fiber length
```

**CSV Export:**
```python
EXPORT_COORDINATES = True     # Enable/disable CSV export
CSV_FILENAME = None          # Auto-generate or specify name
```

**Analysis Job:**
```python
JOB_NAME = 'Job-1'
NUM_CPUS = 8                 # Number of CPU cores
```

### Parameter Guidelines

#### Volume Fraction Recommendations

| Volume Fraction | Difficulty | Recommended Settings |
|----------------|-----------|---------------------|
| 0-5% | Low | Default settings work well |
| 5-10% | Medium | Increase MAX_ITERATIONS to 200000 |
| >10% | High | Consider using advanced algorithms |

#### Fiber Spacing

`MIN_DISTANCE_FACTOR` controls minimum spacing between fibers:
- **1.0**: Fibers can touch (minimum)
- **1.5**: Recommended for most cases
- **2.0**: Conservative spacing, reduces achievable volume fraction

### Output Files

1. **ABAQUS CAE Model**: `RVE-VfXXX.cae`
2. **Fiber Coordinates**: `fiber_centers_VfXXX.csv` (if enabled)
   - Format: `fiber_id, center_x, center_y, center_z, angle_z, angle_x`
3. **Analysis Job**: `Job-1` ready for submission

### Workflow Diagram

```
┌─────────────────────────┐
│  Parameter Setup        │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  Clean Old Models       │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  RSA Fiber Generation   │
│  (with Spatial Grid)    │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  Periodic Mirroring     │
│  (Smart Strategy)       │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  Boolean Operations     │
│  (Fiber-Matrix)         │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  Material Definition    │
│  & Cell Classification  │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  Mesh Generation        │
│  (C3D10 → C3D4 fallback)│
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  Apply PBC              │
│  (Micromechanics Plugin)│
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  Model Ready            │
└─────────────────────────┘
```

### Troubleshooting

#### Common Issues

1. **"Micromechanics Plugin not found"**
   - Ensure the plugin is installed
   - Verify the plugin path in the script (line 76)
   - Check plugin is enabled in ABAQUS/CAE

2. **"Mesh generation failed"**
   - The script automatically tries fallback strategy
   - If both C3D10 and C3D4 fail, reduce TARGET_VF or increase MIN_DISTANCE_FACTOR
   - Check geometry quality manually

3. **"Cannot achieve target volume fraction"**
   - Increase MAX_ITERATIONS
   - Reduce MIN_DISTANCE_FACTOR (carefully)
   - Consider larger RVE size

4. **"Script runs very slowly"**
   - Reduce MAX_ITERATIONS
   - Increase BOUNDARY_THRESHOLD_FACTOR
   - Use lower TARGET_VF for testing

### Performance Tips

- **Start Small**: Test with low volume fractions (1-2%) first
- **Grid Size**: Spatial grid is automatically optimized
- **Iterations**: MAX_ITERATIONS = 100000 is sufficient for Vf ≤ 5%
- **Mesh Density**: Balance between accuracy and computational cost

### Examples

#### Example 1: Basic 5% Volume Fraction Model
```python
LENGTH = WIDTH = HEIGHT = 0.5
FIBRE_LENGTH = 0.25
FIBRE_DIAMETER = 0.025
TARGET_VF = 0.05
```
Expected fibers: ~21 fibers

#### Example 2: High Fiber Content (10%)
```python
TARGET_VF = 0.10
MAX_ITERATIONS = 200000
MIN_DISTANCE_FACTOR = 1.3
```
Expected fibers: ~42 fibers

#### Example 3: Conservative Spacing
```python
MIN_DISTANCE_FACTOR = 2.0
TARGET_VF = 0.03
```
Ensures maximum fiber separation

### Validation

The model has been validated for:
- ✅ Volume fractions: 1% - 10%
- ✅ Various fiber aspect ratios (L/D = 5-20)
- ✅ Periodic boundary conditions
- ✅ Homogenization analysis compatibility

### Citation

If you use this code in your research, please cite:

```bibtex
@software{3D_RVE_SFRC_2025,
  author = {ZPL-03},
  title = {3D RVE Model for Short Fiber Reinforced Composites},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/ZPL-03/3D_RVE_Model_SFRC_MP-PBC}
}
```

### Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Contact

- GitHub: [@ZPL-03](https://github.com/ZPL-03)
- Issues: [GitHub Issues](https://github.com/ZPL-03/3D_RVE_Model_SFRC_MP-PBC/issues)

### Acknowledgments

- ABAQUS Micromechanics Plugin developers
- RSA algorithm research community

---

## 中文

### 概述

本脚本为短纤维增强复合材料(SFRC)的代表性体积单元(RVE)提供了全自动参数化建模解决方案。实现了经典的**随机序贯添加(RSA)算法**，并结合空间网格加速和智能镜像策略优化周期性边界条件处理。

脚本覆盖从几何创建、材料定义、网格划分到周期性边界条件(PBC)施加的完整流程，通过ABAQUS Micromechanics Plugin实现自动化。

### 主要特点

#### 🚀 **高效纤维分布生成**
- **空间网格加速**：将干涉检测复杂度从O(N²)降至约O(N)
- **解析干涉检测**：精确的27盒法确保纤维间最小间距
- **性能统计**：详细的空间网格搜索效率报告
- **可控体积分数**：目标纤维体积分数，可导出CSV数据

#### 🔄 **周期性几何建模**
- **智能镜像策略**：仅为靠近边界的纤维在必要方向创建镜像
- **布尔运算**：纤维与基体的精确切割与合并
- **统一边界裁剪**：确保最终模型的几何周期性

#### 🎯 **材料属性定义**
- **正交各向异性纤维材料**：完整的E1, E2, E3, ν12, ν13, ν23, G12, G13, G23规范
- **各向同性基体材料**：杨氏模量和泊松比
- **智能单元分类**：基于几何的自动材料赋予
- **每根纤维局部坐标系**：精确的纤维材料方向

#### 🌐 **网格划分**
- **参数化网格控制**：自动三维实体单元生成
- **鲁棒回退策略**：优先尝试C3D10(二阶四面体)，失败时回退至C3D4(一阶四面体)
- **质量保证**：全面的错误处理和报告

#### 📐 **周期性边界条件**
- **Micromechanics Plugin集成**：自动为非周期性网格施加PBC
- **完整设置**：创建参考点、约束方程、分析步和载荷
- **均匀化就绪**：为后续材料性能均匀化做好准备

### 算法详情

#### 空间网格加速的RSA算法

脚本实现了增强的RSA算法：

1. **随机生成**：完全随机的纤维位置和方向
2. **序贯放置**：逐个尝试放置纤维
3. **干涉检查**：使用空间网格验证与现有纤维的干涉
4. **接受/拒绝**：无干涉则接受，否则生成新位置
5. **网格优化**：大多数操作的复杂度为O(N)而非O(N²)

#### 智能镜像策略

对于靠近边界的纤维：
- 分析纤维靠近哪些边界
- 仅在必要方向创建镜像副本
- 显著减少几何体数量，提升性能
- 以最小计算开销维持完整周期性

### 环境要求

#### 软件要求
- **ABAQUS**：2023版本或更高
- **Micromechanics Plugin**：必须在ABAQUS/CAE中安装并启用
- **Python**：2.7（ABAQUS内置）

#### 插件安装
必须在ABAQUS中安装Micromechanics Plugin。脚本期望其位于：
```python
d:/Abaqus2023/Plugins
```
如安装路径不同，请修改脚本中的路径。

### 安装

1. 克隆仓库：
```bash
git clone https://github.com/ZPL-03/3D_RVE_Model_SFRC_MP-PBC.git
cd 3D_RVE_Model_SFRC_MP-PBC
```

2. 确保Micromechanics Plugin已正确安装在您的ABAQUS中。

3. 打开ABAQUS/CAE，通过 文件 → 运行脚本 加载脚本。

### 使用方法

#### 基本使用

1. **配置参数**：打开 `3D_RVE_Model_SFRC_MP-PBC.py` 并修改配置区的参数：

```python
# RVE几何参数 (mm)
LENGTH = 0.5
WIDTH = 0.5
HEIGHT = 0.5

# 纤维几何参数 (mm)
FIBRE_LENGTH = 0.25
FIBRE_DIAMETER = 0.025

# 纤维分布参数
MIN_DISTANCE_FACTOR = 1.5
TARGET_VF = 0.05  # 目标体积分数 (5%)
MAX_ITERATIONS = 100000

# 材料属性 (MPa)
FIBER_E1 = 230000.0  # 纵向杨氏模量
MATRIX_E = 3500.0    # 基体杨氏模量
```

2. **运行脚本**：在ABAQUS/CAE中：
   - 文件 → 运行脚本 → 选择 `3D_RVE_Model_SFRC_MP-PBC.py`
   - 或通过命令行：
     ```bash
     abaqus cae noGUI=3D_RVE_Model_SFRC_MP-PBC.py
     ```

3. **监控进度**：脚本提供详细的进度更新：
   - 纤维放置统计
   - 几何创建步骤
   - 材料赋予进度
   - 网格生成状态
   - PBC施加确认

4. **查看结果**：完成后：
   - 模型名称：`RVE-VfXXX`（如 `RVE-Vf005` 表示5%体积分数）
   - 作业名称：`Job-1`（可直接提交）
   - CSV文件：纤维中心坐标（如已启用）

#### 高级配置

**网格控制：**
```python
GLOBAL_SEED_SIZE = 0.01      # 全局种子尺寸 (mm)
DEVIATION_FACTOR = 0.05      # 曲率偏差
MIN_SIZE_FACTOR = 0.05       # 最小单元尺寸因子
```

**边界阈值：**
```python
BOUNDARY_THRESHOLD_FACTOR = 0.5  # 纤维长度的倍数
```

**CSV导出：**
```python
EXPORT_COORDINATES = True     # 启用/禁用CSV导出
CSV_FILENAME = None          # 自动生成或指定名称
```

**分析作业：**
```python
JOB_NAME = 'Job-1'
NUM_CPUS = 8                 # CPU核心数
```

### 参数指南

#### 体积分数建议

| 体积分数 | 难度 | 推荐设置 |
|---------|------|---------|
| 0-5% | 低 | 默认设置即可 |
| 5-10% | 中 | 增加MAX_ITERATIONS至200000 |
| >10% | 高 | 考虑使用高级算法 |

#### 纤维间距

`MIN_DISTANCE_FACTOR` 控制纤维间最小间距：
- **1.0**：纤维可接触（最小值）
- **1.5**：大多数情况推荐值
- **2.0**：保守间距，降低可达到的体积分数

### 输出文件

1. **ABAQUS CAE模型**：`RVE-VfXXX.cae`
2. **纤维坐标**：`fiber_centers_VfXXX.csv`（如已启用）
   - 格式：`fiber_id, center_x, center_y, center_z, angle_z, angle_x`
3. **分析作业**：`Job-1` 可直接提交

### 工作流程图

```
┌─────────────────────────┐
│  参数设置               │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  清理旧模型             │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  RSA纤维生成            │
│  (空间网格加速)         │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  周期性镜像             │
│  (智能策略)             │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  布尔运算               │
│  (纤维-基体)            │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  材料定义               │
│  & 单元分类             │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  网格划分               │
│  (C3D10 → C3D4回退)     │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  施加PBC                │
│  (Micromechanics Plugin)│
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  模型就绪               │
└─────────────────────────┘
```

### 问题排查

#### 常见问题

1. **"未找到Micromechanics Plugin"**
   - 确保插件已安装
   - 验证脚本中的插件路径（第76行）
   - 检查插件在ABAQUS/CAE中已启用

2. **"网格生成失败"**
   - 脚本自动尝试回退策略
   - 如C3D10和C3D4均失败，降低TARGET_VF或增加MIN_DISTANCE_FACTOR
   - 手动检查几何质量

3. **"无法达到目标体积分数"**
   - 增加MAX_ITERATIONS
   - 减小MIN_DISTANCE_FACTOR（谨慎操作）
   - 考虑增大RVE尺寸

4. **"脚本运行非常慢"**
   - 减少MAX_ITERATIONS
   - 增加BOUNDARY_THRESHOLD_FACTOR
   - 测试时使用较低的TARGET_VF

### 性能提示

- **从小开始**：先用低体积分数(1-2%)测试
- **网格尺寸**：空间网格自动优化
- **迭代次数**：MAX_ITERATIONS = 100000 对于Vf ≤ 5%已足够
- **网格密度**：在精度和计算成本间平衡

### 示例

#### 示例1：基本5%体积分数模型
```python
LENGTH = WIDTH = HEIGHT = 0.5
FIBRE_LENGTH = 0.25
FIBRE_DIAMETER = 0.025
TARGET_VF = 0.05
```
预期纤维数：~21根

#### 示例2：高纤维含量(10%)
```python
TARGET_VF = 0.10
MAX_ITERATIONS = 200000
MIN_DISTANCE_FACTOR = 1.3
```
预期纤维数：~42根

#### 示例3：保守间距
```python
MIN_DISTANCE_FACTOR = 2.0
TARGET_VF = 0.03
```
确保最大纤维分离度

### 验证

该模型已针对以下情况验证：
- ✅ 体积分数：1% - 10%
- ✅ 各种纤维长径比(L/D = 5-20)
- ✅ 周期性边界条件
- ✅ 均匀化分析兼容性

### 引用

如果您在研究中使用此代码，请引用：

```bibtex
@software{3D_RVE_SFRC_2025,
  author = {ZPL-03},
  title = {短纤维增强复合材料三维RVE模型},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/ZPL-03/3D_RVE_Model_SFRC_MP-PBC}
}
```

### 贡献

欢迎贡献！请随时提交Pull Request。

### 许可证

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件。

### 联系方式

- GitHub: [@ZPL-03](https://github.com/ZPL-03)
- Issues: [GitHub Issues](https://github.com/ZPL-03/3D_RVE_Model_SFRC_MP-PBC/issues)

### 致谢

- ABAQUS Micromechanics Plugin开发团队
- RSA算法研究社区

---

## Star History

如果此项目对您有帮助，请给个⭐️！

If this project helps you, please give it a ⭐️!
