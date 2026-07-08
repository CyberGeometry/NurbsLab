#pragma once

#include "../concepts.h"
#include "../vec3.h"
#include "../nurbs_curve.h"
#include <vector>
#include <optional>

namespace nurbslab::algorithms {

/// @brief 曲线求交结果
template <FloatingPoint T>
struct IntersectionResult {
    T u;            // 曲线1的参数值
    T v;            // 曲线2的参数值
    Vec3<T> point;  // 交点坐标
};

/// @brief AABB包围盒
template <FloatingPoint T>
struct AABB {
    Vec3<T> min_bound;
    Vec3<T> max_bound;
    
    bool overlaps(const AABB& other) const;
    Vec3<T> center() const;
};

/// @brief 计算曲线的包围盒
template <FloatingPoint T = double>
AABB<T> compute_aabb(const NurbsCurve<T>& curve);

/// @brief 计算曲线一段参数区间的包围盒（采样近似）
template <FloatingPoint T = double>
AABB<T> compute_aabb_segment(const NurbsCurve<T>& curve, T u_start, T u_end, int samples = 10);

/// @brief 两条NURBS曲线的求交算法
/// 使用 AABB 包围盒树递归细分 + 牛顿迭代精化
/// 
/// @tparam T 浮点数类型
/// @param curve1 第一条曲线
/// @param curve2 第二条曲线
/// @param tolerance 交点容差
/// @param max_subdivisions 最大递归细分次数
/// @return 所有交点
template <FloatingPoint T = double>
std::vector<IntersectionResult<T>> curve_curve_intersection(
    const NurbsCurve<T>& curve1,
    const NurbsCurve<T>& curve2,
    T tolerance = 1e-6,
    int max_subdivisions = 20);

} // namespace nurbslab::algorithms