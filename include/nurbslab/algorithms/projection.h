#pragma once

#include "../concepts.h"
#include "../vec3.h"
#include "../nurbs_curve.h"
#include <optional>
#include <vector>

namespace nurbslab::algorithms {

/// @brief 投影结果
template <FloatingPoint T>
struct ProjectionResult {
    T u;                        // 投影点对应的参数值
    Vec3<T> point;              // 曲线上的投影点
    T distance;                 // 点到曲线的距离
    int iterations;             // 迭代次数
    bool converged;             // 是否收敛
};

/// @brief 单点投影：将三维空间点投影到NURBS曲线上
/// 使用牛顿-拉夫森迭代法求解 Point Inversion 问题
/// 
/// 数学原理：
/// 寻找参数 u*，使得 ||C(u*) - P|| 最小
/// 等价于求解：C'(u) · (C(u) - P) = 0
/// 
/// @tparam T 浮点数类型
/// @param curve NURBS曲线
/// @param point 待投影的三维点
/// @param u_init 参数初始猜测值
/// @param tolerance 收敛容差（空间距离）
/// @param max_iterations 最大迭代次数
/// @return 投影结果，若不收敛返回 nullopt
template <FloatingPoint T = double>
std::optional<ProjectionResult<T>> project_point_on_curve(
    const NurbsCurve<T>& curve,
    const Vec3<T>& point,
    T u_init,
    T tolerance = 1e-6,
    int max_iterations = 50);

/// @brief 全局投影：采样多个初始值，寻找全局最近点
/// 适用于不知道初始参数猜测的情况
/// 
/// @tparam T 浮点数类型
/// @param curve NURBS曲线
/// @param point 待投影的三维点
/// @param num_samples 初始采样数量
/// @param tolerance 收敛容差
/// @return 最优投影结果
template <FloatingPoint T = double>
std::optional<ProjectionResult<T>> project_point_global(
    const NurbsCurve<T>& curve,
    const Vec3<T>& point,
    int num_samples = 20,
    T tolerance = 1e-6);

} // namespace nurbslab::algorithms