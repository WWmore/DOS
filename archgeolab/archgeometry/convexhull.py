#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 25 18:04:59 2025

@author: wanghui
"""


import numpy as np
from scipy.spatial import ConvexHull

def build_tri_mesh_from_centroid_and_convex_hull(points):
    """
    基于平面点集的重心和凸包，构建三角形网格
    参数:
        points: 平面点集，形状 (N, 2) 或 (N, 3)（平面点）
    返回:
        vertices: 网格顶点（原始点集 + 重心，可选）
        faces: 三角形面索引，形状 (M, 3)（M为凸包顶点数，每个面对应一个凸包边+重心）
    """
    # 1. 处理输入点集，转为2D（若为3D平面点）
    if points.shape[1] == 3:
        points_2d = points[:, :2]  # 假设在XY平面，可根据实际平面调整
    else:
        points_2d = points.copy()

    # 2. 计算点集重心
    centroid = np.mean(points_2d, axis=0)  # (2,) 或 (3,)
    # 将重心加入顶点集（作为公共顶点）
    vertices = np.vstack([points_2d, centroid])
    centroid_idx = len(vertices) - 1  # 重心在顶点集中的索引

    # 3. 计算点集的凸包
    if len(points_2d) < 3:
        raise ValueError("点集数量不足，至少需要3个点构建凸包")
    hull = ConvexHull(points_2d)
    convex_hull_vertices_idx = hull.vertices  # 凸包顶点在原始点集中的索引

    # 4. 构建三角形面：重心 + 凸包相邻顶点对
    faces = []
    num_hull_vertices = len(convex_hull_vertices_idx)
    for i in range(num_hull_vertices):
        # 凸包当前顶点索引
        v1 = convex_hull_vertices_idx[i]
        # 凸包下一个顶点索引（闭合循环）
        v2 = convex_hull_vertices_idx[(i+1) % num_hull_vertices]
        # 三角形面：(重心, v1, v2)
        faces.append([centroid_idx, v1, v2])
    faces = np.array(faces, dtype=int)

    if vertices.shape[1]==2 and points.shape[1]==3:
        new_col = np.full((vertices.shape[0], 1), np.mean(points[:,2])) 
        vertices = np.hstack((vertices,new_col))

#     """
#     对顶点集进行重映射（去重+唯一索引），并生成对应顺序面列
#     修复：移除错误的 view 语法，使用兼容所有NumPy版本的实现
#     参数:
#         vertices: 原始顶点集（可能含重复顶点），形状 (N, 2)/(N, 3)
#         faces: 原始面列（对应原始顶点索引），形状 (M, 3)
#     返回:
#         remapped_vertices: 重映射后的唯一顶点集，形状 (K, 2)/(K, 3)（K ≤ N）
#         remapped_faces: 重映射后的顺序面列（对应新顶点索引），形状 (M, 3)
#         vertex_map: 旧顶点索引 → 新顶点索引的映射关系，形状 (N,)
#     """
#     # 1. 顶点去重并构建映射关系（修复核心：无需view传入形状，直接用reshape+unique）
#     # 先将顶点转换为 float64 类型（统一类型，避免去重误差）
#     vertices_float = vertices.astype(np.float64)
#     # 重塑为 (N, -1) 确保是二维数组（兼容2D/3D顶点），再直接用 np.unique 按行去重
#     # 关键：np.unique(axis=0) 直接支持二维数组按行去重，无需手动展平
#     unique_vertices, inv_indices = np.unique(
#         vertices_float,
#         axis=0,
#         return_inverse=True  # 返回原始顶点在唯一顶点中的索引（旧→新映射）
#     )

#     # 恢复唯一顶点的原始数据类型（可选，保持与输入一致）
#     remapped_vertices = unique_vertices.astype(vertices.dtype)
#     # 构建旧索引→新索引的映射（直接使用 inv_indices，无需额外处理）
#     vertex_map = inv_indices

#     # 2. 基于顶点映射，更新面列索引（保持原有面顺序，逻辑不变）
#     remapped_faces = np.zeros_like(faces)
#     for i in range(faces.shape[0]):
#         for j in range(faces.shape[1]):
#             old_idx = faces[i, j]
#             remapped_faces[i, j] = vertex_map[old_idx]

#     return remapped_vertices, remapped_faces


# def remove_redundant_vertices(vertices, faces):
    """
    剔除未被面使用的冗余顶点，生成拓扑三角网的Vlist和Flist
    参数:
        vertices: 原始顶点集（可能含冗余顶点），形状 (N, 2)/(N, 3)
        faces: 原始面列（对应原始顶点索引），形状 (M, 3)
    返回:
        Vlist: 无冗余有用顶点集（仅保留被面使用的顶点），形状 (K, 2)/(K, 3)（K ≤ N）
        Flist: 更新后的面列（索引对应Vlist），形状 (M, 3)
        useful_idx: 原始顶点集中的有用顶点索引（用于追溯），形状 (K,)
        idx_map: 原始有用索引 → Vlist新索引的映射（形状 (N,)，冗余顶点对应-1）
    """
    # 1. 提取面列中所有出现过的顶点索引（有用索引，去重）
    # faces.flatten() 将面列展平为一维数组，获取所有用到的顶点索引
    useful_idx = np.unique(faces.flatten())  # 去重，得到有用顶点的原始索引，形状 (K,)

    # 2. 构建「原始顶点索引 → Vlist新索引」的映射（冗余顶点标记为-1）
    idx_map = np.full(vertices.shape[0], -1, dtype=int)  # 初始化全为-1（冗余标记）
    for new_idx, old_idx in enumerate(useful_idx):
        idx_map[old_idx] = new_idx  # 有用顶点：原始索引→新索引（0~K-1连续）

    # 3. 生成无冗余顶点集Vlist（按有用索引筛选原始顶点）
    Vlist = vertices[useful_idx]  # 直接索引筛选，保留有用顶点，形状 (K, 2)/(K, 3)

    # 4. 生成更新后面列Flist（将原始面索引替换为Vlist的新索引）
    Flist = np.zeros_like(faces)
    for i in range(faces.shape[0]):
        for j in range(faces.shape[1]):
            old_idx = faces[i, j]
            Flist[i, j] = idx_map[old_idx]  # 替换为Vlist对应的新索引

    return Vlist, Flist

if __name__ == "__main__":
    # 1. 生成测试平面点集
    np.random.seed(42)
    test_points = np.random.rand(20, 2) * 10  # 20个平面点，范围[0,10)

    # 2. 构建三角形网格
    vertices, faces = build_tri_mesh_from_centroid_and_convex_hull(test_points)
    
    print(vertices.shape, faces.shape)
    print(vertices, faces)
