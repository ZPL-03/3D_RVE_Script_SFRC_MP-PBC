# -*- coding: utf-8 -*-
"""
短纤维增强复合材料RVE参数化建模脚本-混合策略 (集成Micromechanics Plugin)
================================================================

适用场景: 各种纤维体积分数的纤维填充
算法特点: 混合策略（快速RSA + 模拟退火优化）+ 空间网格加速 + 智能镜像策略

功能概述:
---------
本程序实现了短纤维增强复合材料代表性体积单元(RVE)的全自动参数化建模，通过集成多种先进算法，从纤维分布、几何处理到边界条件施加，全程追求计算性能与建模精度。
其核心优势在于采用了一套高度优化的混合纤维生成策略，专门为解决高体积分数下的填充效率与成功率问题而设计。

核心算法:
---------
**分层混合策略（Tiered Hybrid Strategy）**

第一层级 - 快速随机放置（Fast Random Placement）:
  • 对每根待放置的纤维，首先尝试N次快速廉价的随机放置
  • 在RVE空间尚不拥挤时非常高效
  • 算法简单，计算成本低

第二层级 - 模拟退火优化（Simulated Annealing）:
  • 当快速放置失败后，启动计算成本更高但更强大的模拟退火算法
  • 通过温度控制的随机搜索，在局部空间寻找合适位置
  • 高温时接受较差解以探索空间，随温度降低逐渐收敛
  • 自适应扰动：温度越高，位置和角度扰动越大
  • 内置重新加热机制，避免陷入局部最优
  • 能量函数考虑：干涉惩罚 + 边界惩罚（鼓励内部放置）

主要功能模块:
-------------
1. 高效随机纤维分布生成
   - 采用“快速随机填充”与“模拟退火搜索”相结合的混合生成策略，兼顾前期速度与后期填充成功率
   - 引入空间网格数据结构(Spatial Grid)，用于邻近纤维的快速查询，大幅降低干涉检测的计算复杂度
   - 基于解析法的全方位干涉检测（空间线段最短距离），精确高效地确保纤维间最小间距
   - 采用分层生成策略，优先填充RVE内部区域，以提高放置效率
   - 性能统计报告，可输出关于空间网格搜索效率的详细信息
   - 可控的纤维体积分数，纤维中心坐标数据可导出至CSV文件

2. 周期性几何建模
   - 采用智能镜像策略，仅在纤维靠近的边界方向上创建必要的镜像副本，大幅减少几何体数量
   - 通过布尔运算实现纤维与基体的精确切割与合并
   - 对合并后的纤维集合进行统一的边界裁剪，确保最终模型的几何周期性

3. 材料属性定义
   - 为纤维定义正交各向异性材料，为基体定义各向同性材料
   - 采用基于几何距离的智能算法对布尔运算后的体单元(Cell)进行分类，准确识别出属于每根原始纤维的几何域
   - 为每根纤维自动创建独立的局部坐标系，以精确指派其材料方向
   - 自动创建材料、截面并将其赋予对应的几何集合

4. 网格划分
   - 参数化的网格控制，自动生成三维实体单元网格
   - 优先尝试生成二阶四面体单元(C3D10)，当划分失败时自动回退并尝试生成一阶四面体单元(C3D4)

5. 周期性边界条件 (PBC)
   - 创建分析作业（Job）
   - 通过脚本调用 ABAQUS Micromechanics Plugin 施加周期性边界条件，能够处理随机纤维分布导致的非周期性网格
   - 自动创建参考点、约束方程、分析步和载荷，为后续的均匀化分析做好准备

6. 模型管理
   - 在脚本开始时自动清理旧模型，确保建模环境的纯净
   - 根据目标体积分数(Vf)对最终生成的模型进行智能命名

版本: v1.0
日期: 2025.10.09
适用软件: ABAQUS 2023及以上版本 (需确保Micromechanics Plugin已安装并启用)
Python版本: 2.7 (ABAQUS内置)
"""

from abaqus import *
from abaqusConstants import *
from caeModules import *
from mesh import ElemType
import random
import math
import time
import os
import sys

# 导入Micromechanics Plugin所需的模块
sys.path.insert(0, r'd:/Abaqus2023/Plugins')
import microMechanics
from microMechanics.mmpBackend import Interface
from microMechanics.mmpBackend.mmpInterface.mmpRVEConstants import *


# =================================================================
#                 全局参数配置区
# =================================================================

# RVE几何参数 (单位: mm)
LENGTH = 0.5  # RVE在X方向的长度
WIDTH = 0.5  # RVE在Y方向的宽度
HEIGHT = 0.5  # RVE在Z方向的高度

# 纤维几何参数 (单位: mm)
FIBRE_LENGTH = 0.25  # 纤维长度
FIBRE_DIAMETER = 0.025  # 纤维直径
FIBRE_RADIUS = FIBRE_DIAMETER / 2.0  # 纤维半径

# 纤维分布参数
MIN_DISTANCE_FACTOR = 1.5  # 纤维间最小距离因子 (直径的倍数)
TARGET_VF = 0.08  # 目标体积分数
MAX_ITERATIONS = 5000  # 纤维放置的最大尝试次数

# 边界判断参数（统一阈值）
BOUNDARY_THRESHOLD_FACTOR = 0.5  # 边界阈值因子（纤维长度的倍数）
BOUNDARY_THRESHOLD = BOUNDARY_THRESHOLD_FACTOR * FIBRE_LENGTH  # 边界判断阈值(mm)

# 模拟退火参数
SA_INITIAL_TEMP = 1.0  # 初始温度
SA_COOLING_RATE = 0.95  # 冷却速率
SA_MIN_TEMP = 0.001  # 最小温度
SA_MAX_ATTEMPTS_PER_FIBER = 200  # 每根纤维的最大尝试次数
SA_POSITION_PERTURBATION = 0.1  # 位置扰动基准范围(mm) - 随温度自适应
SA_ANGLE_PERTURBATION = 1.5  # 角度扰动基准范围(弧度) - 随温度自适应

# 分层生成参数
TARGET_INNER_RATIO = 0.7  # 期望内部纤维的比例概率

# 混合生成策略参数 (分层策略)
FAST_PLACEMENT_ATTEMPTS = 666  # 在调用SA之前，先尝试一定次数的快速随机放置
SA_REHEAT_THRESHOLD = 150  # 连续失败这么多次后，进行一次再加热

# 纤维材料属性（正交各向异性）
FIBER_E1 = 230000.0  # 纤维方向的杨氏模量 (纵向)
FIBER_E2 = 15000.0  # 垂直纤维方向的杨氏模量 (径向)
FIBER_E3 = 15000.0  # 垂直纤维方向的杨氏模量 (径向)
FIBER_NU12 = 0.2  # 主泊松比
FIBER_NU13 = 0.2  # 主泊松比
FIBER_NU23 = 0.25  # 横截面内的泊松比
FIBER_G12 = 15000.0  # 剪切模量 (纵向-径向)
FIBER_G13 = 15000.0  # 剪切模量 (纵向-径向)
FIBER_G23 = 7000.0  # 剪切模量 (径向-径向)

# 基体材料属性（各向同性）
MATRIX_E = 3500.0  # 杨氏模量
MATRIX_NU = 0.35  # 泊松比

# 网格参数
GLOBAL_SEED_SIZE = 0.01  # 全局网格种子尺寸
DEVIATION_FACTOR = 0.05  # 曲面的网格偏差因子
MIN_SIZE_FACTOR = 0.05  # 最小单元尺寸因子

# CSV导出参数
EXPORT_COORDINATES = True  # 是否导出纤维中心坐标
CSV_FILENAME = None  # 导出的CSV文件名 (None表示自动生成)

# 分析作业参数
JOB_NAME = 'Job-1'  # 分析作业的名称
NUM_CPUS = 8        # 分析时使用的CPU核心数


# =================================================================
#                 空间网格结构（周期性边界支持 + 智能优化）
# =================================================================
class SpatialGridPeriodic:
    """周期性空间网格 - 用于快速邻近查询

    采用网格划分方法加速纤维干涉检测，支持周期性边界条件。
    每个网格单元记录其影响范围内的纤维索引，查询时只需检查相关网格单元，大幅减少计算量。
    智能周期性检查：只在必要时检查镜像位置。
    性能统计：记录查询次数和早期拒绝次数。
    """

    def __init__(self, length, width, height, cell_size):
        """初始化空间网格

        参数:
            length, width, height: RVE尺寸
            cell_size: 网格单元尺寸
        """
        self.length = length
        self.width = width
        self.height = height
        self.cell_size = cell_size

        # 计算各方向网格数量
        self.nx = int(math.ceil(length / cell_size))
        self.ny = int(math.ceil(width / cell_size))
        self.nz = int(math.ceil(height / cell_size))

        # 网格字典：键为网格索引，值为纤维索引列表
        self.grid = {}
        # 纤维列表：存储(纤维数据, 端点)元组
        self.fibers = []

        # 性能统计
        self.query_count = 0
        self.early_reject_count = 0

    def _get_cell_index(self, x, y, z):
        """获取点所在的网格索引

        使用模运算实现周期性边界
        """
        i = int(x / self.cell_size) % self.nx
        j = int(y / self.cell_size) % self.ny
        k = int(z / self.cell_size) % self.nz
        return (i, j, k)

    def add_fiber(self, fiber_data, endpoints):
        """添加纤维到空间网格

        参数:
            fiber_data: [centx, centy, centz, angle_z, angle_x]
            endpoints: ((x1,y1,z1), (x2,y2,z2))
        """
        fiber_idx = len(self.fibers)
        self.fibers.append((fiber_data, endpoints))

        centx, centy, centz = fiber_data[0], fiber_data[1], fiber_data[2]

        # 计算纤维影响半径
        influence_radius = FIBRE_LENGTH / 2.0 + MIN_DISTANCE_FACTOR * FIBRE_DIAMETER

        # 获取影响范围内的所有网格单元
        affected_cells = self._get_affected_cells_periodic(
            (centx, centy, centz), influence_radius
        )

        # 在所有受影响的网格单元中注册该纤维
        for cell in affected_cells:
            if cell not in self.grid:
                self.grid[cell] = []
            self.grid[cell].append(fiber_idx)

    def remove_fiber(self, fiber_idx):
        """从空间网格移除纤维

        参数:
            fiber_idx: 纤维索引
        """
        if fiber_idx >= len(self.fibers):
            return

        # 从所有网格单元中移除该纤维索引
        for cell_key in self.grid.keys():
            if fiber_idx in self.grid[cell_key]:
                self.grid[cell_key].remove(fiber_idx)

    def _get_affected_cells_periodic(self, center, radius):
        """获取影响范围内的网格单元（考虑周期性）

        参数:
            center: 中心点坐标 (x, y, z)
            radius: 影响半径

        返回:
            set: 受影响的网格单元索引集合
        """
        cx, cy, cz = center
        cells = set()

        # 计算影响范围（网格单元数）
        cell_range = int(math.ceil(radius / self.cell_size)) + 1
        ci, cj, ck = self._get_cell_index(cx, cy, cz)

        # 遍历周围的网格单元
        for di in range(-cell_range, cell_range + 1):
            for dj in range(-cell_range, cell_range + 1):
                for dk in range(-cell_range, cell_range + 1):
                    # 使用模运算实现周期性
                    i = (ci + di) % self.nx
                    j = (cj + dj) % self.ny
                    k = (ck + dk) % self.nz
                    cells.add((i, j, k))

        return cells

    def get_nearby_fibers_periodic(self, center, search_radius):
        """获取附近的纤维索引（智能周期性）

        只在纤维靠近边界时才检查镜像位置，大幅减少不必要的计算

        参数:
            center: 查询中心点 (x, y, z)
            search_radius: 搜索半径

        返回:
            set: 附近纤维的索引集合
        """
        self.query_count += 1
        nearby = set()

        centx, centy, centz = center

        # 智能判断是否需要考虑周期性边界
        need_periodic_x = (centx < search_radius) or (centx > self.length - search_radius)
        need_periodic_y = (centy < search_radius) or (centy > self.width - search_radius)
        need_periodic_z = (centz < search_radius) or (centz > self.height - search_radius)

        # 如果不靠近任何边界，只检查中心位置
        if not (need_periodic_x or need_periodic_y or need_periodic_z):
            self.early_reject_count += 1
            affected_cells = self._get_affected_cells_periodic(center, search_radius)
            for cell in affected_cells:
                if cell in self.grid:
                    nearby.update(self.grid[cell])
            return nearby

        # 靠近边界时，根据位置决定检查哪些镜像
        offsets_to_check = [(0, 0, 0)]

        for dx in ([-self.length, 0, self.length] if need_periodic_x else [0]):
            for dy in ([-self.width, 0, self.width] if need_periodic_y else [0]):
                for dz in ([-self.height, 0, self.height] if need_periodic_z else [0]):
                    if dx != 0 or dy != 0 or dz != 0:
                        offsets_to_check.append((dx, dy, dz))

        # 对每个需要检查的位置进行查询
        for dx, dy, dz in offsets_to_check:
            mirror_center = (
                center[0] + dx,
                center[1] + dy,
                center[2] + dz
            )

            affected_cells = self._get_affected_cells_periodic(mirror_center, search_radius)

            for cell in affected_cells:
                if cell in self.grid:
                    nearby.update(self.grid[cell])

        return nearby

    def get_statistics(self):
        """获取性能统计信息"""
        early_reject_rate = (self.early_reject_count * 100.0 / self.query_count
                             if self.query_count > 0 else 0)
        avg_fibers_per_cell = (len(self.fibers) * 1.0 / len(self.grid)
                               if self.grid else 0)
        return {
            'total_queries': self.query_count,
            'early_rejects': self.early_reject_count,
            'early_reject_rate': early_reject_rate,
            'total_cells': len(self.grid),
            'avg_fibers_per_cell': avg_fibers_per_cell
        }


# =================================================================
#                 纤维干涉检测（完整27盒子）
# =================================================================
def shortest_distance_between_segments(p1, q1, p2, q2):
    """计算空间中两条线段 p1-q1 和 p2-q2 之间的最短距离（解析法）

    采用向量代数方法，比离散采样法更精确、更高效。
    算法来源: David Eberly (Geometric Tools, LLC), "Distance Between Two Lines in 3D"

    参数:
        p1, q1: 第一条线段的端点 ((x,y,z), (x,y,z))
        p2, q2: 第二条线段的端点 ((x,y,z), (x,y,z))

    返回:
        float: 最短距离
    """
    # 向量表示
    u = (q1[0] - p1[0], q1[1] - p1[1], q1[2] - p1[2])
    v = (q2[0] - p2[0], q2[1] - p2[1], q2[2] - p2[2])
    w = (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2])

    # 点积
    a = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]  # dot(u, u)
    b = u[0] * v[0] + u[1] * v[1] + u[2] * v[2]  # dot(u, v)
    c = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]  # dot(v, v)
    d = u[0] * w[0] + u[1] * w[1] + u[2] * w[2]  # dot(u, w)
    e = v[0] * w[0] + v[1] * w[1] + v[2] * w[2]  # dot(v, w)

    D = a * c - b * b
    sD, tD = D, D

    # 计算无限长直线上的最近点参数 sN/sD, tN/tD
    # s = sN/sD, t = tN/tD
    if D < 1e-7:  # 平行线
        sN = 0.0
        sD = 1.0
        tN = e
        tD = c
    else:
        sN = (b * e - c * d)
        tN = (a * e - b * d)
        if sN < 0.0:
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD:
            sN = sD
            tN = e + b
            tD = c

    if tN < 0.0:
        tN = 0.0
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = a
        else:
            sN = -d
        sD = a
    elif tN > tD:
        tN = tD
        if (-d + b) < 0.0:
            sN = 0.0
        elif (-d + b) > a:
            sN = a
        else:
            sN = -d + b
        sD = a

    sc = 0.0 if abs(sN) < 1e-7 else sN / sD
    tc = 0.0 if abs(tN) < 1e-7 else tN / tD

    # 计算距离向量
    dP_x = w[0] + (sc * u[0]) - (tc * v[0])
    dP_y = w[1] + (sc * u[1]) - (tc * v[1])
    dP_z = w[2] + (sc * u[2]) - (tc * v[2])

    return math.sqrt(dP_x ** 2 + dP_y ** 2 + dP_z ** 2)


def check_interference_complete(spatial_grid, endpoints):
    """
    使用空间网格加速的27盒子干涉检测 ("空间线段最短距离"解析法)

    **检测策略**：
    1. 使用空间网格快速查询新纤维附近的纤维子集
    2. 对这个子集中的每根纤维，执行完整的27盒子干涉检测
       - 现有纤维：保持在中心盒子的原始位置
       - 新纤维：在27个位置检测 (1个中心位置 + 26个镜像位置)
       - 计算新纤维在每个位置与现有纤维之间的最短几何距离 (模拟退火算法需要用到)

    参数:
        spatial_grid: 空间网格对象
        endpoints: 新纤维端点 ((x1,y1,z1), (x2,y2,z2))

    返回:
        (是否无干涉, 最小距离): tuple
            是否无干涉: bool
            最小距离: float
    """
    c, d = endpoints[0], endpoints[1]
    center_new = (
        (c[0] + d[0]) / 2.0,
        (c[1] + d[1]) / 2.0,
        (c[2] + d[2]) / 2.0
    )

    min_allowed_distance = MIN_DISTANCE_FACTOR * FIBRE_DIAMETER
    fiber_half_length = FIBRE_LENGTH / 2.0
    search_radius = fiber_half_length + min_allowed_distance

    # 使用空间网格快速获取附近的纤维
    nearby_indices = spatial_grid.get_nearby_fibers_periodic(center_new, search_radius)

    if not nearby_indices:
        return True, float('inf')

    # 定义完整的27个周期性位置（包括中心位置(0,0,0)和26个镜像位置）
    offsets = []
    for dx in [-LENGTH, 0, LENGTH]:
        for dy in [-WIDTH, 0, WIDTH]:
            for dz in [-HEIGHT, 0, HEIGHT]:
                offsets.append((dx, dy, dz))

    min_dist_found = float('inf')

    # 对附近的每个纤维进行完整的27盒子检测
    for fiber_idx in nearby_indices:
        _, existing_endpoints = spatial_grid.fibers[fiber_idx]
        a, b = existing_endpoints[0], existing_endpoints[1]

        # 对新纤维的所有27个周期性位置进行检测
        for dx, dy, dz in offsets:
            c_mir = (c[0] + dx, c[1] + dy, c[2] + dz)
            d_mir = (d[0] + dx, d[1] + dy, d[2] + dz)

            dist = shortest_distance_between_segments(a, b, c_mir, d_mir)

            if dist < min_dist_found:
                min_dist_found = dist

            if dist < min_allowed_distance:
                return False, dist

    return True, min_dist_found


# =================================================================
#                 模拟退火算法实现
# =================================================================
class SimulatedAnnealingOptimizer:
    """模拟退火算法 - 用于纤维放置策略

    通过温度控制的随机搜索方法，寻找合适的纤维放置位置。
    高温时接受较差的解以探索更大空间，随温度降低逐渐收敛。
    """

    def __init__(self, spatial_grid):
        """初始化模拟退火器

        参数:
            spatial_grid: 空间网格对象
        """
        self.spatial_grid = spatial_grid
        self.temperature = SA_INITIAL_TEMP

    def is_in_inner_zone(self, centx, centy, centz):
        """判断纤维中心是否在内部区域（远离边界）

        参数:
            centx, centy, centz: 纤维中心坐标

        返回:
            bool: True表示在内部区域
        """
        return (centx > BOUNDARY_THRESHOLD and centx < LENGTH - BOUNDARY_THRESHOLD and
                centy > BOUNDARY_THRESHOLD and centy < WIDTH - BOUNDARY_THRESHOLD and
                centz > BOUNDARY_THRESHOLD and centz < HEIGHT - BOUNDARY_THRESHOLD)

    def generate_candidate_fiber(self, prefer_inner=True):
        """生成候选纤维

        参数:
            prefer_inner: 是否优先在内部区域生成

        返回:
            list: [centx, centy, centz, angle_z, angle_x]
        """
        if prefer_inner and random.random() < TARGET_INNER_RATIO:
            # 优先在内部区域生成
            centx = random.uniform(BOUNDARY_THRESHOLD, LENGTH - BOUNDARY_THRESHOLD)
            centy = random.uniform(BOUNDARY_THRESHOLD, WIDTH - BOUNDARY_THRESHOLD)
            centz = random.uniform(BOUNDARY_THRESHOLD, HEIGHT - BOUNDARY_THRESHOLD)
        else:
            # 在整个RVE中生成
            centx = random.uniform(0, LENGTH)
            centy = random.uniform(0, WIDTH)
            centz = random.uniform(0, HEIGHT)

        angle_z = random.uniform(0, 2 * math.pi)
        angle_x = random.uniform(0, 2 * math.pi)

        return [centx, centy, centz, angle_z, angle_x]

    def perturb_fiber(self, fiber_data):
        """扰动纤维位置和方向 (自适应)

        扰动范围会根据当前温度进行缩放，实现
        高温时大范围探索，低温时精细调整。

        参数:
            fiber_data: 原始纤维数据

        返回:
            list: 扰动后的纤维数据
        """
        new_fiber = fiber_data[:]

        # 自适应因子：温度越高，扰动越大
        adaptive_factor = self.temperature / SA_INITIAL_TEMP

        # 位置扰动
        pos_perturb = SA_POSITION_PERTURBATION * adaptive_factor
        new_fiber[0] += random.uniform(-pos_perturb, pos_perturb)
        new_fiber[1] += random.uniform(-pos_perturb, pos_perturb)
        new_fiber[2] += random.uniform(-pos_perturb, pos_perturb)

        # 确保在RVE内 (使用模运算实现周期性环绕)
        new_fiber[0] = new_fiber[0] % LENGTH
        new_fiber[1] = new_fiber[1] % WIDTH
        new_fiber[2] = new_fiber[2] % HEIGHT

        # 角度扰动
        ang_perturb = SA_ANGLE_PERTURBATION * adaptive_factor
        new_fiber[3] += random.uniform(-ang_perturb, ang_perturb)
        new_fiber[4] += random.uniform(-ang_perturb, ang_perturb)

        # 角度归一化
        new_fiber[3] = new_fiber[3] % (2 * math.pi)
        new_fiber[4] = new_fiber[4] % (2 * math.pi)

        return new_fiber

    def compute_endpoints(self, fiber_data):
        """计算纤维端点坐标

        参数:
            fiber_data: [centx, centy, centz, angle_z, angle_x]

        返回:
            tuple: ((x1,y1,z1), (x2,y2,z2))
        """
        centx, centy, centz, angle_z, angle_x = fiber_data

        x1 = centx + 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.cos(angle_x)
        y1 = centy + 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.sin(angle_x)
        z1 = centz + 0.5 * FIBRE_LENGTH * math.cos(angle_z)

        x2 = centx - 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.cos(angle_x)
        y2 = centy - 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.sin(angle_x)
        z2 = centz - 0.5 * FIBRE_LENGTH * math.cos(angle_z)

        endpoints = ((x1, y1, z1), (x2, y2, z2))

        return endpoints

    def compute_energy(self, fiber_data, min_dist):
        """计算能量函数

        能量越低表示状态越好。包含两部分：
        1. 干涉惩罚：距离过近时能量增加
        2. 边界惩罚：鼓励内部放置

        参数:
            fiber_data: 纤维数据
            endpoints: 纤维端点
            min_dist: 与最近纤维的距离

        返回:
            float: 能量值
        """
        # 干涉惩罚
        min_allowed_dist = MIN_DISTANCE_FACTOR * FIBRE_DIAMETER
        if min_dist < min_allowed_dist:
            interference_penalty = 1000.0 * (1.0 - min_dist / min_allowed_dist)
        else:
            interference_penalty = 0.0

        # 边界惩罚（鼓励内部放置）
        centx, centy, centz = fiber_data[0], fiber_data[1], fiber_data[2]
        if self.is_in_inner_zone(centx, centy, centz):
            boundary_penalty = 0.0
        else:
            boundary_penalty = 0.1

        total_energy = interference_penalty + boundary_penalty

        return total_energy

    def accept_move(self, old_energy, new_energy):
        """Metropolis准则判断是否接受新状态

        参数:
            old_energy: 旧状态能量
            new_energy: 新状态能量

        返回:
            bool: 是否接受新状态
        """
        if new_energy < old_energy:
            return True

        delta_e = new_energy - old_energy
        # 避免温度过低时 exp() 溢出
        if self.temperature < 1e-6:
            return False
        probability = math.exp(-delta_e / self.temperature)

        return random.random() < probability

    def optimize_fiber_placement(self, fiber_data, max_attempts=SA_MAX_ATTEMPTS_PER_FIBER):
        """使用模拟退火算法放置单根纤维

        参数:
            fiber_data: 初始纤维数据
            max_attempts: 最大尝试次数

        返回:
            tuple: (最佳纤维数据, 最佳端点, 是否有效)
        """
        current_fiber = fiber_data[:]
        current_endpoints = self.compute_endpoints(current_fiber)
        no_interference, current_min_dist = check_interference_complete(
            self.spatial_grid, current_endpoints
        )

        # 提前退出，如果初始安置已经有效
        if no_interference:
            return current_fiber, current_endpoints, True

        current_energy = self.compute_energy(current_fiber, current_min_dist)

        best_fiber = current_fiber[:]
        best_endpoints = current_endpoints
        best_energy = float('inf')  # 从高能量开始

        for attempt in range(max_attempts):
            # 生成新候选
            new_fiber = self.perturb_fiber(current_fiber)
            new_endpoints = self.compute_endpoints(new_fiber)

            # 检查干涉
            no_interference_new, new_min_dist = check_interference_complete(
                self.spatial_grid, new_endpoints
            )
            new_energy = self.compute_energy(new_fiber, new_min_dist)

            # 判断是否接受
            if self.accept_move(current_energy, new_energy):
                current_fiber = new_fiber
                current_energy = new_energy

                # 更新最佳解
                if no_interference_new:
                    return new_fiber, new_endpoints, True
                elif new_energy < best_energy:
                    best_energy = new_energy
                    best_fiber = new_fiber[:]
                    best_endpoints = new_endpoints

        # 如果循环结束仍未找到有效解，则返回失败
        return best_fiber, best_endpoints, False

    def cool_down(self):
        """降温"""
        self.temperature *= SA_COOLING_RATE
        if self.temperature < SA_MIN_TEMP:
            self.temperature = SA_MIN_TEMP

    def reheat(self):
        """重新加热（用于跳出局部最优）"""
        self.temperature = min(SA_INITIAL_TEMP, self.temperature * 10.0)


# =================================================================
#                 智能纤维生成器（混合策略）
# =================================================================
def generate_fibers_with_optimization(spatial_grid, target_vf, verbose=True):
    """使用分层策略和模拟退火生成纤维

    采用"快速随机尝试 + 精细模拟退火"的分层策略：
    1. Tier 1 (Fast Placement): 对每一根待放置的纤维，首先尝试N次快速、廉价的随机放置。这在RVE空间尚不拥挤时非常高效。
    2. Tier 2 (Simulated Annealing): 如果N次快速尝试全部失败，说明空间已经变得拥挤，此时再启动计算成本更高但更强大的模拟退火算法进行局部精细搜索。

    这种策略确保了计算资源被用在最需要的地方，大幅提升了在高体积分数下的填充效率和成功率。

    参数:
        spatial_grid: 空间网格对象
        target_vf: 目标体积分数
        verbose: 是否打印详细信息

    返回:
        list: 纤维列表，每个元素为 [centx, centy, centz, angle_z, angle_x]
    """
    fibre_list = []
    volum_base = LENGTH * WIDTH * HEIGHT
    volum_fiber = FIBRE_LENGTH * math.pi * FIBRE_RADIUS ** 2

    optimizer = SimulatedAnnealingOptimizer(spatial_grid)

    iteration = 0
    consecutive_failures = 0

    start_time = time.time()
    last_report_time = start_time
    last_reported_fiber_count = -1

    if verbose:
        print("\n" + "=" * 70)
        print("FIBER GENERATION WITH TIERED STRATEGY")
        print("=" * 70)
        print("Strategy: Tier 1 (Fast Random) -> Tier 2 (Simulated Annealing)")
        print("Fast Placement Attempts per Fiber: %d" % FAST_PLACEMENT_ATTEMPTS)
        print("SA Reheat Threshold (consecutive failures): %d" % SA_REHEAT_THRESHOLD)
        print("Target Vf: %.2f%%" % (target_vf * 100))
        print("Boundary threshold: %.4f mm (%.1fx fiber length)" %
              (BOUNDARY_THRESHOLD, BOUNDARY_THRESHOLD_FACTOR))
        print("Starting generation...\n")

    while iteration < MAX_ITERATIONS:
        iteration += 1

        current_vf = len(fibre_list) * volum_fiber / volum_base
        if current_vf >= target_vf:
            if verbose:
                print("\nTarget volume fraction reached!")
            break

        placed_successfully = False
        strategy = "N/A"

        # 第1层级：快速随机放置 ---
        strategy = "Fast Random"
        for _ in range(FAST_PLACEMENT_ATTEMPTS):
            candidate_fiber = optimizer.generate_candidate_fiber(prefer_inner=True)
            candidate_endpoints = optimizer.compute_endpoints(candidate_fiber)
            is_valid, _ = check_interference_complete(spatial_grid, candidate_endpoints)

            if is_valid:
                fibre_list.append(candidate_fiber)
                spatial_grid.add_fiber(candidate_fiber, candidate_endpoints)
                placed_successfully = True
                break

        # 第2层级：模拟退火（如果第1层失败）
        if not placed_successfully:
            strategy = "Simulated Annealing"

            # 检查是否需要重新加热
            if consecutive_failures > 0 and consecutive_failures % SA_REHEAT_THRESHOLD == 0:
                optimizer.reheat()
                if verbose:
                    print("\n  [Reheating Triggered] T=%.4f (at %d consecutive failures)" % (
                        optimizer.temperature, consecutive_failures))

            # 动态调整内部生成偏好
            progress = current_vf / target_vf
            prefer_inner = progress < (1.0 - TARGET_INNER_RATIO)

            initial_candidate = optimizer.generate_candidate_fiber(prefer_inner=prefer_inner)
            optimized_fiber, optimized_endpoints, is_valid = optimizer.optimize_fiber_placement(initial_candidate)

            if is_valid:
                fibre_list.append(optimized_fiber)
                spatial_grid.add_fiber(optimized_fiber, optimized_endpoints)
                placed_successfully = True

        # 放置后逻辑
        if placed_successfully:
            consecutive_failures = 0
            # 仅当通过"模拟退火"策略成功放置纤维时，才进行降温。
            if strategy == "Simulated Annealing":
                optimizer.cool_down()
        else:
            consecutive_failures += 1

        # 详细且不刷屏的日志记录逻辑
        if verbose:
            current_fiber_count = len(fibre_list)
            time_since_last_report = time.time() - last_report_time
            # 如果成功放置了10的倍数根纤维，或距离上次报告已超过10秒，则打印日志
            if (current_fiber_count > last_reported_fiber_count and current_fiber_count % 10 == 0) or \
                    (time_since_last_report > 10.0):

                elapsed = time.time() - start_time
                report_str = "  Iter %-6d | Fibers: %-4d | Vf: %7.4f%% | Strategy: %-19s" % \
                             (iteration, current_fiber_count, current_vf * 100, strategy)

                # 仅在模拟退火策略时显示温度
                if strategy == "Simulated Annealing":
                    report_str += " | T: %.4f" % optimizer.temperature
                else:
                    report_str += " | T: %-8s" % "N/A"

                if not placed_successfully:
                    report_str += " | [Working... %d fails]" % consecutive_failures

                report_str += " | Time: %.1fs" % elapsed
                print(report_str)
                last_report_time = time.time()
                last_reported_fiber_count = current_fiber_count

        if consecutive_failures > (SA_REHEAT_THRESHOLD * 4):  # 如果长时间重加热后仍无进展，则停止
            if verbose:
                print("\n  Warning: Too many consecutive failures even after reheating, stopping...")
            break

    elapsed = time.time() - start_time
    final_vf = len(fibre_list) * volum_fiber / volum_base

    if verbose:
        print("\n" + "=" * 70)
        print("FIBER GENERATION COMPLETED")
        print("=" * 70)
        print("Total iterations: %d" % iteration)
        print("Fibers generated: %d" % len(fibre_list))
        print("Final Vf: %.4f%%" % (final_vf * 100))
        success_rate = (len(fibre_list) * 100.0 / iteration) if iteration > 0 else 0
        print("Success rate: %.2f%%" % success_rate)
        print("Time elapsed: %.2f seconds" % elapsed)

        # 性能统计报告
        grid_stats = spatial_grid.get_statistics()

        print("\n" + "-" * 70)
        print("PERFORMANCE STATISTICS")
        print("-" * 70)
        print("Spatial Grid:")
        print("  Total queries: %d" % grid_stats['total_queries'])
        print("  Early reject rate: %.2f%% (%d early rejects)" % (
            grid_stats['early_reject_rate'],
            grid_stats['early_rejects']
        ))
        print("  Total cells used: %d" % grid_stats['total_cells'])
        print("  Avg fibers per cell: %.2f" % grid_stats['avg_fibers_per_cell'])
        print("-" * 70)
        print("=" * 70 + "\n")

    return fibre_list


# =================================================================
#                 其他辅助函数
# =================================================================
def pointToLineSegmentDistance(point, line_end1, line_end2):
    """计算点到线段的最短距离

    采用向量投影方法计算点到线段的距离。

    参数:
        point: 点坐标 (px, py, pz)
        line_end1: 线段端点1 (x1, y1, z1)
        line_end2: 线段端点2 (x2, y2, z2)

    返回:
        float: 最短距离
    """
    px, py, pz = point
    x1, y1, z1 = line_end1
    x2, y2, z2 = line_end2

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1

    length_sq = dx * dx + dy * dy + dz * dz

    if length_sq == 0:
        return math.sqrt((px - x1) ** 2 + (py - y1) ** 2 + (pz - z1) ** 2)

    t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy + (pz - z1) * dz) / length_sq))

    proj_x = x1 + t * dx
    proj_y = y1 + t * dy
    proj_z = z1 + t * dz

    return math.sqrt((px - proj_x) ** 2 + (py - proj_y) ** 2 + (pz - proj_z) ** 2)


def determineMirrorOffsets(fiber_data):
    """(智能策略) 判断纤维需要在哪些方向上创建周期性镜像

    仅在纤维靠近的边界方向上创建镜像，而不是创建全部26个。
    这能显著减少后续建模步骤中的几何体数量，提升性能。

    参数:
        fiber_data: [centx, centy, centz, angle_z, angle_x]

    返回:
        list: 必需的镜像偏移列表 [(dx, dy, dz), ...] 或 []
    """
    centx, centy, centz, _, _ = fiber_data

    # 独立判断每个轴是否靠近边界
    near_x = centx < BOUNDARY_THRESHOLD or centx > LENGTH - BOUNDARY_THRESHOLD
    near_y = centy < BOUNDARY_THRESHOLD or centy > WIDTH - BOUNDARY_THRESHOLD
    near_z = centz < BOUNDARY_THRESHOLD or centz > HEIGHT - BOUNDARY_THRESHOLD

    if not (near_x or near_y or near_z):
        return []

    # 根据需要，构建每个方向的偏移列表
    x_offsets = [-LENGTH, 0, LENGTH] if near_x else [0]
    y_offsets = [-WIDTH, 0, WIDTH] if near_y else [0]
    z_offsets = [-HEIGHT, 0, HEIGHT] if near_z else [0]

    mirror_offsets = []
    # 组合所有必需的偏移
    for dx in x_offsets:
        for dy in y_offsets:
            for dz in z_offsets:
                # 跳过(0,0,0)这个原始位置
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                mirror_offsets.append((dx, dy, dz))

    return mirror_offsets


def buildAllFiberData_ShortFiber(fibre_list):
    """(智能策略) 构建完整纤维数据（含必要的周期性镜像和原始ID）用于Cell分类

    对于靠近边界的纤维，仅在其靠近的方向上生成镜像数据。
    每个数据点都包含原始纤维的ID，确保Cell分类时能正确识别。

    参数:
        fibre_list: 原始纤维列表

    返回:
        list: 分类必要的纤维数据列表 [(original_id, center, end1, end2), ...]
    """
    all_fiber_data = []

    for i, fiber in enumerate(fibre_list):
        centx, centy, centz, angle_z, angle_x = fiber

        # 计算纤维端点
        x1 = centx + 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.cos(angle_x)
        y1 = centy + 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.sin(angle_x)
        z1 = centz + 0.5 * FIBRE_LENGTH * math.cos(angle_z)

        x2 = centx - 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.cos(angle_x)
        y2 = centy - 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.sin(angle_x)
        z2 = centz - 0.5 * FIBRE_LENGTH * math.cos(angle_z)

        # 独立判断每个轴是否靠近边界
        near_x = centx < BOUNDARY_THRESHOLD or centx > LENGTH - BOUNDARY_THRESHOLD
        near_y = centy < BOUNDARY_THRESHOLD or centy > WIDTH - BOUNDARY_THRESHOLD
        near_z = centz < BOUNDARY_THRESHOLD or centz > HEIGHT - BOUNDARY_THRESHOLD

        # 根据需要，构建每个方向的位置列表（含原始位置0）
        x_positions = [-LENGTH, 0, LENGTH] if near_x else [0]
        y_positions = [-WIDTH, 0, WIDTH] if near_y else [0]
        z_positions = [-HEIGHT, 0, HEIGHT] if near_z else [0]

        positions_to_create = []
        for dx in x_positions:
            for dy in y_positions:
                for dz in z_positions:
                    positions_to_create.append((dx, dy, dz))

        # 为所有必需的位置（原始+镜像）创建纤维数据
        for dx, dy, dz in positions_to_create:
            mirror_center = (centx + dx, centy + dy, centz + dz)
            mirror_end1 = (x1 + dx, y1 + dy, z1 + dz)
            mirror_end2 = (x2 + dx, y2 + dy, z2 + dz)
            all_fiber_data.append((i, mirror_center, mirror_end1, mirror_end2))

    return all_fiber_data


def exportFiberCentersToCSV(fiber_list, filename, rveSize, fiberRadius, fiberLength, target_Vf):
    """将纤维中心坐标导出为CSV文件

    参数:
        fiber_list: 纤维列表
        filename: 输出文件名
        rveSize: RVE尺寸 [length, width, height]
        fiberRadius: 纤维半径
        fiberLength: 纤维长度
        target_Vf: 目标体积分数

    返回:
        str: 文件路径，失败返回None
    """
    try:
        work_dir = os.getcwd()
        filepath = os.path.join(work_dir, filename)

        with open(filepath, 'w') as f:
            # 写入文件头信息
            f.write("# Short Fiber RVE - Fiber Center Coordinates\n")
            f.write("# Generated: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
            f.write("# RVE Size: %.6f x %.6f x %.6f mm\n" % (rveSize[0], rveSize[1], rveSize[2]))
            f.write("# Fiber Radius: %.6f mm\n" % fiberRadius)
            f.write("# Fiber Length: %.6f mm\n" % fiberLength)
            f.write("# Target Vf: %.4f%%\n" % (target_Vf * 100))
            f.write("# Total Fibers: %d\n" % len(fiber_list))
            f.write("#\n")
            f.write("Fiber_ID,Center_X,Center_Y,Center_Z,Angle_Z_rad,Angle_X_rad\n")

            # 写入纤维数据
            for i, fiber_data in enumerate(fiber_list, start=1):
                centx, centy, centz, angle_z, angle_x = fiber_data
                f.write("%d,%.8f,%.8f,%.8f,%.8f,%.8f\n" %
                        (i, centx, centy, centz, angle_z, angle_x))

        print("\n" + "=" * 60)
        print("CSV Export: SUCCESS")
        print("File: %s" % filepath)
        print("=" * 60 + "\n")
        return filepath

    except Exception as e:
        print("\nWARNING: CSV export failed: %s\n" % str(e))
        return None


def getCellCenterFromVertices(cell):
    """通过顶点计算Cell中心（备用方法）

    当无法直接获取Cell质心时，通过顶点坐标平均值估算。

    参数:
        cell: ABAQUS Cell对象

    返回:
        tuple: (center_x, center_y, center_z) 或 None
    """
    try:
        vertices = cell.getVertices()
        if not vertices or len(vertices) == 0:
            return None

        x_coords, y_coords, z_coords = [], [], []

        for vertex in vertices:
            # 检查 vertex 是否有效以及 pointOn 属性是否存在
            if hasattr(vertex, 'pointOn') and vertex.pointOn:
                coord = vertex.pointOn[0]
                if coord:
                    x_coords.append(coord[0])
                    y_coords.append(coord[1])
                    z_coords.append(coord[2])

        if len(x_coords) > 0 and len(y_coords) > 0 and len(z_coords) > 0:
            center_x = sum(x_coords) / len(x_coords)
            center_y = sum(y_coords) / len(y_coords)
            center_z = sum(z_coords) / len(z_coords)
            return (center_x, center_y, center_z)
        else:
            return None

    except Exception as e:
        print("      Warning: Error in getCellCenterFromVertices: %s" % str(e))
        return None


def classifyCellsImproved_ShortFiber(all_cells, fibre_list, rveVolume, original_fiber_count):
    """
    基于几何位置的智能cell分类算法

    核心思路:
    1. 最大体积cell必为主基体
    2. 构建含原始纤维ID的周期性镜像纤维数据 (智能策略，仅包含必要的镜像)
    3. 对所有剩余cell逐个进行几何判断:
       - 计算cell质心到所有纤维(含镜像)轴线的最短距离
       - 记录下最近的纤维的原始ID
       - 如果最短距离 < 纤维半径, 则将该cell归类到对应的纤维ID下
       - 否则, 将其归类为基体碎片

    参数:
        all_cells: 部件中的所有体单元
        fibre_list: 纤维数据列表 [[centx, centy, centz, angle_z, angle_x], ...]
        rveVolume: RVE总体积 (用于体积分数验证)
        original_fiber_count: 原始RVE内的纤维数量 (用于最终体积验证)

    返回:
        fiber_cells_by_id: 按纤维ID分类的cell字典 {fiber_id: [cell1, cell2], ...}
        matrix_cells: 基体单元列表(含主基体和碎片)
        actual_vf: 实际体积分数
        all_fiber_cells: 所有纤维单元的平铺列表
    """
    print("  Starting intelligent cell classification...")
    print("  Total cells to classify: %d" % len(all_cells))

    if not all_cells:
        print("  WARNING: No cells found in the part to classify.")
        return {}, [], 0.0, []

    # 步骤1: 按体积排序, 最大的必为主基体
    sorted_cells = sorted(all_cells, key=lambda c: c.getSize(), reverse=True)
    matrix_cells = [sorted_cells[0]]
    potential_cells = sorted_cells[1:]
    print("  Largest cell identified as main matrix (Volume: %.6e)" % sorted_cells[0].getSize())

    # 步骤2: 构建包含周期性镜像和原始ID的纤维数据列表
    all_fiber_data_with_ids = buildAllFiberData_ShortFiber(fibre_list)
    print("  Built fiber data list with %d entries (including periodic mirrors)" % len(all_fiber_data_with_ids))

    # 步骤3: 对所有potential_cells进行几何判断
    fiber_cells_by_id = {}
    matrix_fragments = []

    for idx, cell in enumerate(potential_cells):
        # 获取cell质心 (3级回退策略获取单元中心点)
        cell_center = None
        # 策略1: 尝试 getCentroid()
        try:
            centroid = cell.getCentroid()
            if centroid and len(centroid) >= 3:
                cell_center = (centroid[0], centroid[1], centroid[2])
        except:
            pass  # 静默失败，进入下一策略

        # 策略2: 尝试通过顶点平均计算
        if cell_center is None:
            cell_center = getCellCenterFromVertices(cell)

        # 策略3: 最终保障，使用 pointOn 属性
        if cell_center is None:
            try:
                cell_center = cell.pointOn[0]
            except:
                pass

        if cell_center is None:
            raise Exception("FATAL ERROR: Cannot determine center for potential cell %d" % (idx + 1))

        # 计算cell质心到所有纤维轴线的最短距离, 并记录最近的纤维ID
        min_dist = float('inf')
        closest_fiber_original_id = -1
        for original_id, _, end1, end2 in all_fiber_data_with_ids:
            dist = pointToLineSegmentDistance(cell_center, end1, end2)
            if dist < min_dist:
                min_dist = dist
                closest_fiber_original_id = original_id

        # 判断归属: 距离 < 纤维半径 → 纤维, 否则 → 基体碎片
        if min_dist < FIBRE_RADIUS:
            if closest_fiber_original_id not in fiber_cells_by_id:
                fiber_cells_by_id[closest_fiber_original_id] = []
            fiber_cells_by_id[closest_fiber_original_id].append(cell)
        else:
            matrix_fragments.append(cell)

    matrix_cells.extend(matrix_fragments)

    # 构建所有纤维cell的平铺列表
    all_fiber_cells = []
    for fiber_id in fiber_cells_by_id:
        all_fiber_cells.extend(fiber_cells_by_id[fiber_id])

    # 步骤4: 输出分类结果
    print("\n  Classification Results:")
    print("  Fiber cells identified and grouped for %d fibers" % len(fiber_cells_by_id))
    print("  Matrix cells total: %d (1 main + %d fragments)" % (len(matrix_cells), len(matrix_fragments)))

    # 步骤5: 体积验证
    theoretical_fiber_volume = original_fiber_count * FIBRE_LENGTH * math.pi * FIBRE_RADIUS ** 2
    theoretical_vf = theoretical_fiber_volume / rveVolume if rveVolume > 0 else 0

    fiber_total_volume = sum([c.getSize() for c in all_fiber_cells])
    actual_vf = fiber_total_volume / rveVolume if rveVolume > 0 else 0

    deviation = abs(actual_vf - theoretical_vf) / theoretical_vf * 100 if theoretical_vf > 0 else 0

    print("\n  Volume Fraction Validation:")
    print("  Theoretical Vf: %.4f%%" % (theoretical_vf * 100))
    print("  Actual Vf (from geometry): %.4f%%" % (actual_vf * 100))
    print("  Deviation: %.2f%%" % deviation)

    return fiber_cells_by_id, matrix_cells, actual_vf, all_fiber_cells


def createCellSequence(part, cells):
    """创建Cell序列

    将Cell列表转换为ABAQUS所需的Cell序列对象。

    参数:
        part: Part对象
        cells: Cell列表

    返回:
        CellArray: ABAQUS Cell序列
    """
    if not cells:
        return part.cells[0:0]

    indices = []
    for cell in cells:
        for i, part_cell in enumerate(part.cells):
            if cell == part_cell:
                indices.append(i)
                break

    if indices:
        return part.cells[indices[0]:indices[0] + 1] if len(indices) == 1 else reduce(
            lambda a, b: a + b, [part.cells[i:i + 1] for i in indices]
        )
    else:
        return part.cells[0:0]


# =================================================================
#                 主建模流程
# =================================================================
script_start_time = time.time()

print("\n" + "=" * 70)
print("Short Fiber RVE - Parametric Modeling Script")
print("Using Tiered Hybrid Strategy and Micromechanics Plugin for PBC")
print("=" * 70 + "\n")

volum_base = LENGTH * WIDTH * HEIGHT
volum_fiber = FIBRE_LENGTH * math.pi * FIBRE_RADIUS ** 2

print("Configuration:")
print("  RVE: %.3f x %.3f x %.3f mm" % (LENGTH, WIDTH, HEIGHT))
print("  Fiber: L=%.3f mm, D=%.4f mm" % (FIBRE_LENGTH, FIBRE_DIAMETER))
print("  Target Vf: %.2f%%" % (TARGET_VF * 100))
print("  Boundary threshold: %.4f mm (%.1fx fiber length)" %
      (BOUNDARY_THRESHOLD, BOUNDARY_THRESHOLD_FACTOR))
print("")

# =================================================================
# 步骤1: 模型清理与创建
# =================================================================
print("Step 1: Cleaning up existing models and creating new one...")
step_start_time = time.time()

target_model_name = 'RVE-Vf-%d' % int(round(TARGET_VF * 100))
# 为初始创建创建一个唯一的临时名称，以避免在清理过程中发生冲突
temp_model_name_for_creation = 'TEMP_MODEL_FOR_CREATION_' + str(int(time.time()))

# 首先创建临时模型
mdb.Model(name=temp_model_name_for_creation, modelType=STANDARD_EXPLICIT)

# 获取要删除的所有其他模型名称的列表
models_to_delete = [m for m in mdb.models.keys() if m != temp_model_name_for_creation]

if models_to_delete:
    print("  Deleting old models:")
    for m in models_to_delete:
        print("    - %s" % m)
        del mdb.models[m]
    print("  Old models deleted.")
else:
    print("  No old models to delete.")

# 将新的空模型重命名为最终目标名称
mdb.models.changeKey(fromName=temp_model_name_for_creation, toName=target_model_name)
myModel = mdb.models[target_model_name]
print("  New model created and named: '%s'" % target_model_name)
print("  Model management completed in %.2f seconds\n" % (time.time() - step_start_time))


# 步骤2: 创建基体
print("Step 2: Creating base matrix...")
step_start_time = time.time()
myPart = myModel.Part(name="Part-base", dimensionality=THREE_D, type=DEFORMABLE_BODY)
mySketch = myModel.ConstrainedSketch(name="sketch-1", sheetSize=2.0)
mySketch.rectangle(point1=(0, 0), point2=(LENGTH, WIDTH))
myPart.BaseSolidExtrude(sketch=mySketch, depth=HEIGHT)
print("  Base matrix created in %.2f seconds\n" % (time.time() - step_start_time))

# 步骤3: 使用混合策略生成纤维
print("Step 3: Generating fibers with tiered strategy...")
step_start_time = time.time()

# 基于实际搜索半径设置初始化空间网格
influence_radius = FIBRE_LENGTH / 2.0 + MIN_DISTANCE_FACTOR * FIBRE_DIAMETER
grid_cell_size = influence_radius * 0.8  # 稍小于搜索半径
spatial_grid = SpatialGridPeriodic(LENGTH, WIDTH, HEIGHT, grid_cell_size)
print("  Spatial grid initialized: %d x %d x %d cells (cell_size=%.4f mm)" %
      (spatial_grid.nx, spatial_grid.ny, spatial_grid.nz, grid_cell_size))

# 生成纤维
fibre = generate_fibers_with_optimization(spatial_grid, TARGET_VF, verbose=True)

num_fibers = len(fibre)
theoretical_vf = num_fibers * volum_fiber / volum_base

# CSV导出
if EXPORT_COORDINATES and len(fibre) > 0:
    if CSV_FILENAME is None:
        CSV_FILENAME = "FiberCenters_Vf%d_%s.csv" % (
            int(TARGET_VF * 100), time.strftime("%Y%m%d_%H%M%S")
        )
    exportFiberCentersToCSV(fibre, CSV_FILENAME,
                            [LENGTH, WIDTH, HEIGHT],
                            FIBRE_RADIUS, FIBRE_LENGTH, TARGET_VF)

# 步骤4: 创建纤维模板
print("Step 4: Creating fiber template...")
step_start_time = time.time()
myPart3 = myModel.Part(name="Part-fibre-solid", dimensionality=THREE_D, type=DEFORMABLE_BODY)
mySketch3 = myModel.ConstrainedSketch(name="sketch-3", sheetSize=2.0)
mySketch3.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(FIBRE_RADIUS, 0.0))
myPart3.BaseSolidExtrude(sketch=mySketch3, depth=FIBRE_LENGTH)
print("  Template created in %.2f seconds\n" % (time.time() - step_start_time))

# =================================================================
# 步骤5: 在装配中创建纤维实例（含周期性镜像）
# =================================================================
print("Step 5: Creating fiber instances with periodic mirrors in assembly...")
step_start_time = time.time()

myAssembly = myModel.rootAssembly
instance_counter = 0

for num in range(len(fibre)):
    fiber_data = fibre[num]
    x = fiber_data[0]
    y = fiber_data[1]
    z = fiber_data[2]
    angle_z = fiber_data[3]
    angle_x = fiber_data[4]

    # (智能策略) 判断该纤维需要在哪些镜像位置创建副本
    mirror_offsets = determineMirrorOffsets(fiber_data)

    # 将原始位置(0,0,0)和所有镜像位置组合
    positions_to_create = [(0, 0, 0)]
    positions_to_create.extend(mirror_offsets)

    for offset_idx, (dx, dy, dz) in enumerate(positions_to_create):
        instance_counter += 1
        instance_name = 'Part-fibre-solid-%d' % instance_counter

        myAssembly.Instance(name=instance_name, part=myPart3, dependent=ON)

        myAssembly.rotate(instanceList=(instance_name,),
                          axisPoint=(0, 0, 0), axisDirection=(0, 1, 0),
                          angle=angle_z * 360 / (2 * math.pi))

        myAssembly.rotate(instanceList=(instance_name,),
                          axisPoint=(0, 0, 0), axisDirection=(0, 0, 1),
                          angle=angle_x * 360 / (2 * math.pi))

        # 平移到正确位置（原始位置或镜像位置）
        final_x = x + dx - 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.cos(angle_x)
        final_y = y + dy - 0.5 * FIBRE_LENGTH * math.sin(angle_z) * math.sin(angle_x)
        final_z = z + dz - 0.5 * FIBRE_LENGTH * math.cos(angle_z)

        myAssembly.translate(instanceList=(instance_name,),
                             vector=(final_x, final_y, final_z))

assembly_time = time.time() - step_start_time
print("  Fiber instances created: %d total instances (including mirrors) in %.2f seconds" %
      (instance_counter, assembly_time))
print("    Original fibers: %d" % len(fibre))
print("    Mirror instances: %d\n" % (instance_counter - len(fibre) if instance_counter > len(fibre) else 0))

# =================================================================
# 步骤6: 将所有纤维合并为单个部件
# =================================================================
print("Step 6: Merging all fiber instances...")
step_start_time = time.time()

if instance_counter > 0:
    instances = []
    for ins in myAssembly.instances.values():
        instances.append(ins)

    myAssembly.InstanceFromBooleanMerge(name='Part-fibre-all',
                                        instances=tuple(instances),
                                        keepIntersections=ON,
                                        originalInstances=DELETE,
                                        domain=GEOMETRY)
    print(
        "  All fibers (including mirrors) merged into single part in %.2f seconds\n" % (time.time() - step_start_time))
else:
    print("  No fibers to merge, creating an empty 'Part-fibre-all'.")
    myModel.Part(name='Part-fibre-all', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    myAssembly.Instance(name='Part-fibre-all-1', part=myModel.parts['Part-fibre-all'], dependent=ON)

# 步骤7: 创建基准面
print("Step 7: Creating datum planes for trimming...")
step_start_time = time.time()
myModel.parts['Part-fibre-all'].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XYPLANE)
myModel.parts['Part-fibre-all'].DatumPlaneByPrincipalPlane(offset=HEIGHT, principalPlane=XYPLANE)
myModel.parts['Part-fibre-all'].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=YZPLANE)
myModel.parts['Part-fibre-all'].DatumPlaneByPrincipalPlane(offset=LENGTH, principalPlane=YZPLANE)
myModel.parts['Part-fibre-all'].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XZPLANE)
myModel.parts['Part-fibre-all'].DatumPlaneByPrincipalPlane(offset=WIDTH, principalPlane=XZPLANE)
myModel.parts['Part-fibre-all'].DatumAxisByPrincipalAxis(principalAxis=XAXIS)
myModel.parts['Part-fibre-all'].DatumAxisByPrincipalAxis(principalAxis=YAXIS)
print("  Datum planes created in %.2f seconds\n" % (time.time() - step_start_time))

# 步骤8: 裁剪到边界
print("Step 8: Trimming merged fibers to RVE boundaries...")
step_start_time = time.time()

# Z方向裁剪
myModel.ConstrainedSketch(
    gridSpacing=0.16, name='__profile__', sheetSize=6.52,
    transform=myModel.parts['Part-fibre-all'].MakeSketchTransform(
        sketchPlane=myModel.parts['Part-fibre-all'].datums[3],
        sketchPlaneSide=SIDE1,
        sketchUpEdge=myModel.parts['Part-fibre-all'].datums[8],
        sketchOrientation=RIGHT,
        origin=(LENGTH / 2, WIDTH / 2, HEIGHT)
    )
)
myModel.parts['Part-fibre-all'].projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=myModel.sketches['__profile__']
)
myModel.sketches['__profile__'].rectangle(point1=(-2, 2), point2=(2, -2))
myModel.parts['Part-fibre-all'].CutExtrude(
    flipExtrudeDirection=ON,
    sketch=myModel.sketches['__profile__'],
    sketchOrientation=RIGHT,
    sketchPlane=myModel.parts['Part-fibre-all'].datums[3],
    sketchPlaneSide=SIDE1,
    sketchUpEdge=myModel.parts['Part-fibre-all'].datums[8]
)
del myModel.sketches['__profile__']

myModel.ConstrainedSketch(
    gridSpacing=0.16, name='__profile__', sheetSize=6.52,
    transform=myModel.parts['Part-fibre-all'].MakeSketchTransform(
        sketchPlane=myModel.parts['Part-fibre-all'].datums[2],
        sketchPlaneSide=SIDE1,
        sketchUpEdge=myModel.parts['Part-fibre-all'].datums[8],
        sketchOrientation=RIGHT,
        origin=(LENGTH / 2, WIDTH / 2, 0.0)
    )
)
myModel.parts['Part-fibre-all'].projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=myModel.sketches['__profile__']
)
myModel.sketches['__profile__'].rectangle(point1=(-2, 2), point2=(2, -2))
myModel.parts['Part-fibre-all'].CutExtrude(
    flipExtrudeDirection=OFF,
    sketch=myModel.sketches['__profile__'],
    sketchOrientation=RIGHT,
    sketchPlane=myModel.parts['Part-fibre-all'].datums[2],
    sketchPlaneSide=SIDE1,
    sketchUpEdge=myModel.parts['Part-fibre-all'].datums[8]
)
del myModel.sketches['__profile__']

# X方向裁剪
myModel.ConstrainedSketch(
    gridSpacing=0.14, name='__profile__', sheetSize=5.94,
    transform=myModel.parts['Part-fibre-all'].MakeSketchTransform(
        sketchPlane=myModel.parts['Part-fibre-all'].datums[4],
        sketchPlaneSide=SIDE1,
        sketchUpEdge=myModel.parts['Part-fibre-all'].datums[9],
        sketchOrientation=RIGHT,
        origin=(0.0, WIDTH / 2, HEIGHT / 2)
    )
)
myModel.parts['Part-fibre-all'].projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=myModel.sketches['__profile__']
)
myModel.sketches['__profile__'].rectangle(point1=(-2, 2), point2=(2, -2))
myModel.parts['Part-fibre-all'].CutExtrude(
    flipExtrudeDirection=OFF,
    sketch=myModel.sketches['__profile__'],
    sketchOrientation=RIGHT,
    sketchPlane=myModel.parts['Part-fibre-all'].datums[4],
    sketchPlaneSide=SIDE1,
    sketchUpEdge=myModel.parts['Part-fibre-all'].datums[9]
)
del myModel.sketches['__profile__']

myModel.ConstrainedSketch(
    gridSpacing=0.13, name='__profile__', sheetSize=5.57,
    transform=myModel.parts['Part-fibre-all'].MakeSketchTransform(
        sketchPlane=myModel.parts['Part-fibre-all'].datums[5],
        sketchPlaneSide=SIDE1,
        sketchUpEdge=myModel.parts['Part-fibre-all'].datums[9],
        sketchOrientation=RIGHT,
        origin=(LENGTH, WIDTH / 2, HEIGHT / 2)
    )
)
myModel.parts['Part-fibre-all'].projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=myModel.sketches['__profile__']
)
myModel.sketches['__profile__'].rectangle(point1=(-2, 2), point2=(2, -2))
myModel.parts['Part-fibre-all'].CutExtrude(
    flipExtrudeDirection=ON,
    sketch=myModel.sketches['__profile__'],
    sketchOrientation=RIGHT,
    sketchPlane=myModel.parts['Part-fibre-all'].datums[5],
    sketchPlaneSide=SIDE1,
    sketchUpEdge=myModel.parts['Part-fibre-all'].datums[9]
)
del myModel.sketches['__profile__']

# Y方向裁剪
myModel.ConstrainedSketch(
    gridSpacing=0.11, name='__profile__', sheetSize=4.79,
    transform=myModel.parts['Part-fibre-all'].MakeSketchTransform(
        sketchPlane=myModel.parts['Part-fibre-all'].datums[7],
        sketchPlaneSide=SIDE1,
        sketchUpEdge=myModel.parts['Part-fibre-all'].datums[8],
        sketchOrientation=RIGHT,
        origin=(LENGTH / 2, WIDTH, HEIGHT / 2)
    )
)
myModel.parts['Part-fibre-all'].projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=myModel.sketches['__profile__']
)
myModel.sketches['__profile__'].rectangle(point1=(-2, 2), point2=(2, -2))
myModel.parts['Part-fibre-all'].CutExtrude(
    flipExtrudeDirection=ON,
    sketch=myModel.sketches['__profile__'],
    sketchOrientation=RIGHT,
    sketchPlane=myModel.parts['Part-fibre-all'].datums[7],
    sketchPlaneSide=SIDE1,
    sketchUpEdge=myModel.parts['Part-fibre-all'].datums[8]
)
del myModel.sketches['__profile__']

myModel.ConstrainedSketch(
    gridSpacing=0.1, name='__profile__', sheetSize=4.26,
    transform=myModel.parts['Part-fibre-all'].MakeSketchTransform(
        sketchPlane=myModel.parts['Part-fibre-all'].datums[6],
        sketchPlaneSide=SIDE1,
        sketchUpEdge=myModel.parts['Part-fibre-all'].datums[8],
        sketchOrientation=RIGHT,
        origin=(LENGTH / 2, 0.0, HEIGHT / 2)
    )
)
myModel.parts['Part-fibre-all'].projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=myModel.sketches['__profile__']
)
myModel.sketches['__profile__'].rectangle(point1=(-2, 2), point2=(2, -2))
myModel.parts['Part-fibre-all'].CutExtrude(
    flipExtrudeDirection=OFF,
    sketch=myModel.sketches['__profile__'],
    sketchOrientation=RIGHT,
    sketchPlane=myModel.parts['Part-fibre-all'].datums[6],
    sketchPlaneSide=SIDE1,
    sketchUpEdge=myModel.parts['Part-fibre-all'].datums[8]
)
del myModel.sketches['__profile__']

print("  Trimming completed in %.2f seconds\n" % (time.time() - step_start_time))

# =================================================================
# 步骤9: 布尔运算 - 切割基体并合并
# =================================================================
print("Step 9: Performing Boolean operations (creates Composite Part and Instance)...")
step_start_time = time.time()

myModel.rootAssembly.DatumCsysByDefault(CARTESIAN)
myModel.rootAssembly.Instance(dependent=ON, name='Part-base-1',
                              part=myModel.parts['Part-base'])

myModel.rootAssembly.InstanceFromBooleanCut(
    cuttingInstances=(myModel.rootAssembly.instances['Part-fibre-all-1'],),
    instanceToBeCut=myModel.rootAssembly.instances['Part-base-1'],
    name='Matrix',
    originalInstances=DELETE)

myModel.rootAssembly.Instance(dependent=ON, name='Part-fibre-all-1',
                              part=myModel.parts['Part-fibre-all'])

myModel.rootAssembly.InstanceFromBooleanMerge(
    domain=GEOMETRY,
    instances=(myModel.rootAssembly.instances['Matrix-1'],
               myModel.rootAssembly.instances['Part-fibre-all-1']),
    keepIntersections=ON,
    name='Composite',
    originalInstances=DELETE)

print("  Boolean operations completed in %.2f seconds - Composite Part with periodic geometry created\n" % (
        time.time() - step_start_time))

# =================================================================
# 步骤10: 清理临时部件
# =================================================================
print("Step 10: Cleaning up temporary parts...")
step_start_time = time.time()

del myModel.parts['Part-base']
del myModel.parts['Part-fibre-solid']

print("  Temporary parts deleted in %.2f seconds\n" % (time.time() - step_start_time))

# 步骤11-12: 材料、截面与分类
print("Step 11-12: Creating materials and sections...")
step_start_time = time.time()
fiberMaterial = myModel.Material(name='Material-Fiber')
fiberMaterial.Elastic(
    type=ENGINEERING_CONSTANTS,
    table=((FIBER_E1, FIBER_E2, FIBER_E3,
            FIBER_NU12, FIBER_NU13, FIBER_NU23,
            FIBER_G12, FIBER_G13, FIBER_G23),)
)

matrixMaterial = myModel.Material(name='Material-Matrix')
matrixMaterial.Elastic(table=((MATRIX_E, MATRIX_NU),))

myModel.HomogeneousSolidSection(material='Material-Fiber',
                                name='Section-Fiber', thickness=None)
myModel.HomogeneousSolidSection(material='Material-Matrix',
                                name='Section-Matrix', thickness=None)
print("  Materials and sections created in %.2f seconds\n" % (time.time() - step_start_time))

print("Step 12: Classifying cells and assigning sections/orientations...")
step_start_time = time.time()
p_composite = myModel.parts['Composite']
all_cells = p_composite.cells
rveVolume = LENGTH * WIDTH * HEIGHT

# 使用分类函数，按纤维ID对单元进行分组
fiber_cells_by_id, matrix_cells, actual_vf, all_fiber_cells = classifyCellsImproved_ShortFiber(
    all_cells, fibre, rveVolume, num_fibers
)

print("\n  Assigning sections and precise material orientations for each fiber...")

# 为每一根纤维(或其碎片)创建独立的Set,并赋予截面和方向
for fiber_id, cells_of_this_fiber in fiber_cells_by_id.items():
    if not cells_of_this_fiber:
        continue

    # 1. 为当前纤维创建单元集合
    set_name = 'Set-Fiber-%d' % (fiber_id + 1)
    cell_sequence = createCellSequence(p_composite, cells_of_this_fiber)
    region = p_composite.Set(name=set_name, cells=cell_sequence)

    # 2. 为当前纤维集合赋予截面属性
    p_composite.SectionAssignment(region=region,
                                  sectionName='Section-Fiber',
                                  offset=0.0,
                                  offsetType=MIDDLE_SURFACE,
                                  offsetField='')

    # 3. 为当前纤维创建并指派精确的材料方向
    centx, centy, centz, angle_z, angle_x = fibre[fiber_id]

    # 计算纤维轴向矢量 (材料的1轴)
    v1_x = math.sin(angle_z) * math.cos(angle_x)
    v1_y = math.sin(angle_z) * math.sin(angle_x)
    v1_z = math.cos(angle_z)

    # 计算一个与轴向垂直的矢量 (材料的2轴)，通过与全局Z轴叉乘得到。如果纤维接近Z轴，则与Y轴叉乘
    global_axis = (0.0, 0.0, 1.0) if abs(v1_z) < 0.999 else (0.0, 1.0, 0.0)

    # 叉乘: axis_2 = axis_1 x global_axis
    v2_x = v1_y * global_axis[2] - v1_z * global_axis[1]
    v2_y = v1_z * global_axis[0] - v1_x * global_axis[2]
    v2_z = v1_x * global_axis[1] - v1_y * global_axis[0]

    # 在纤维中心创建局部坐标系
    csys_name = 'Csys-Fiber-%d' % (fiber_id + 1)
    origin_pt = (centx, centy, centz)
    pt_on_axis1 = (centx + v1_x, centy + v1_y, centz + v1_z)
    pt_on_plane12 = (centx + v2_x, centy + v2_y, centz + v2_z)

    try:
        # 步骤A: 创建基准坐标系特征，并捕获返回的Feature对象
        datum_feature = p_composite.DatumCsysByThreePoints(name=csys_name,
                                                           coordSysType=CARTESIAN,
                                                           origin=origin_pt,
                                                           point1=pt_on_axis1,
                                                           point2=pt_on_plane12)

        # 步骤B: 使用Feature对象的.id属性从datums库中获取正确的DatumCsys对象
        local_csys_object = p_composite.datums[datum_feature.id]

        # 步骤C: 将正确的DatumCsys对象传递给MaterialOrientation函数
        p_composite.MaterialOrientation(region=region,
                                        orientationType=SYSTEM,
                                        axis=AXIS_1,
                                        localCsys=local_csys_object,
                                        stackDirection=STACK_1)
    except Exception as e:
        print("  WARNING: Could not create/assign orientation for Fiber-%d. Error: %s" % (fiber_id + 1, str(e)))

print("  Processed %d fibers for section and orientation assignment." % len(fiber_cells_by_id))

# 为所有基体单元创建集合并赋予截面
if len(matrix_cells) > 0:
    matrix_sequence = createCellSequence(p_composite, matrix_cells)
    p_composite.Set(name='Set-Matrix', cells=matrix_sequence)
    p_composite.SectionAssignment(region=p_composite.sets['Set-Matrix'],
                                  sectionName='Section-Matrix',
                                  offset=0.0,
                                  offsetType=MIDDLE_SURFACE,
                                  offsetField='')
    print("  Matrix section assigned to %d cells." % len(matrix_cells))
else:
    print("  WARNING: No matrix cells were identified!")

print("  Cell classification and property assignment completed in %.2f seconds\n" % (time.time() - step_start_time))

# 步骤13: 网格生成
print("Step 13: Generating mesh...")
step_start_time = time.time()
p_composite.seedPart(
    size=GLOBAL_SEED_SIZE,
    deviationFactor=DEVIATION_FACTOR,
    minSizeFactor=MIN_SIZE_FACTOR,
    constraint=FREE
)

all_classified_cells = createCellSequence(p_composite, all_fiber_cells + matrix_cells)
p_composite.setMeshControls(regions=all_classified_cells, elemShape=TET, technique=FREE)

elemType_tet = ElemType(elemCode=C3D10, elemLibrary=STANDARD)
p_composite.setElementType(regions=(all_classified_cells,), elemTypes=(elemType_tet,))

try:
    print("  Attempting to mesh with C3D10 (quadratic tet) elements...")
    p_composite.generateMesh()
    print("  SUCCESS: Mesh generated with C3D10.")
except Exception as e:
    print("\n  WARNING: C3D10 meshing failed. Error: %s" % str(e))
    print("  FALLBACK: Attempting to mesh with C3D4 (linear tet) elements...")
    elemType_linear = ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    p_composite.setElementType(regions=(all_classified_cells,), elemTypes=(elemType_linear,))
    try:
        p_composite.generateMesh()
        print("  SUCCESS: Mesh generated with C3D4.")
    except Exception as e2:
        print("\n  FATAL: Meshing failed with both C3D10 and C3D4. The geometry may be too complex.")
        raise e2

total_elements = len(p_composite.elements)
print("  Mesh generation completed in %.2f seconds" % (time.time() - step_start_time))
print("  Total elements generated: %d\n" % total_elements)

if 'Composite-1' in myModel.rootAssembly.instances.keys():
    del myModel.rootAssembly.instances['Composite-1']

myModel.rootAssembly.Instance(name='Composite-1', part=p_composite, dependent=ON)

# =================================================================
# 步骤14: 创建分析作业并使用 Micromechanics Plugin 施加PBC
# =================================================================
print("Step 14: Creating Analysis Job and Applying PBCs with Micromechanics Plugin...")
step_start_time = time.time()

try:
    # 插件要求Job必须预先存在，因此先创建Job
    print("  Creating analysis job '%s' with %d CPUs as a prerequisite for the plugin..." % (JOB_NAME, NUM_CPUS))
    mdb.Job(name=JOB_NAME, model=myModel.name, description='RVE Homogenization Analysis',
            type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
            scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=NUM_CPUS,
            numDomains=NUM_CPUS, numGPUs=0)
    print("  Job '%s' created successfully." % JOB_NAME)

    # 调用插件函数，参数与GUI设置对应
    print("  Calling Micromechanics Plugin to apply PBC and configure analysis...")
    Interface.Loading.MechanicalModelMaker(
        modelName=myModel.name,
        jobName=JOB_NAME,
        constraintType='PERIODIC',
        drivenField='STRAIN',
        doNotSubmit=True,  # 对应 "Do not Submit Job for analysis" 已勾选，即只施加边界条件不求解
        homogenizeProperties=(True, False, False)  # 设置弹性属性，均质复选框：[弹性性能，热膨胀，密度]
    )
    print("  Micromechanics Plugin executed successfully.")
    print("  PBC, reference points, steps, and loads have been configured in job '%s'." % JOB_NAME)

except Exception as e:
    print("\n  FATAL ERROR: Failed to apply PBC using Micromechanics Plugin.")
    print("  Please ensure the plugin is installed and enabled in Abaqus/CAE.")
    print("  Error details: %s" % str(e))
    raise e

print("\n  Job creation and PBC application completed in %.2f seconds\n" % (time.time() - step_start_time))

# --- TIMER END: 脚本执行总时间 ---
total_script_time = time.time() - script_start_time

# 最终总结
print("=" * 70)
print("RVE MODEL COMPLETED")
print("=" * 70)
print("\nSummary:")
print("  Model: %s" % target_model_name)
print("  RVE: %.3f x %.3f x %.3f mm" % (LENGTH, WIDTH, HEIGHT))
print("  Fibers: %d" % num_fibers)
print("  Theoretical Vf: %.4f%%" % (theoretical_vf * 100))
print("  Actual Vf (from geometry): %.4f%%" % (actual_vf * 100))
print("  Elements: %d" % total_elements)
print("  Total script execution time: %.2f seconds" % total_script_time)
print("\nKey Features:")
print("  - Hybrid generation (Fast Random + Simulated Annealing)")
print("  - Analytical interference detection (vector algebra)")
print("  - Layered filling strategy for high Vf")
print("  - Complete 27-box periodicity for interference checks")
print("  - Smart periodic mirroring (performance-aware)")
print("  - Spatial Grid for fast interference checks")
print("  - Robust PBCs for non-periodic meshes via Micromechanics Plugin")
print("  - Performance statistics reporting")
print("\nNext Steps:")
print("  1. The model, analysis steps, loads, and job ('%s') are ready." % JOB_NAME)
print("  2. Verify mesh quality and boundary conditions in Visualization/Load modules.")
print("  3. You can now directly submit '%s' for analysis." % JOB_NAME)
print("=" * 70 + "\n")