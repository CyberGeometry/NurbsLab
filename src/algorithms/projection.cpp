#include "nurbslab/algorithms/projection.h"
#include <cmath>
#include <limits>
#include <algorithm>

namespace nurbslab::algorithms {

template <FloatingPoint T>
std::optional<ProjectionResult<T>> project_point_on_curve(
    const NurbsCurve<T>& curve,
    const Vec3<T>& point,
    T u_init,
    T tolerance,
    int max_iterations)
{
    T u = u_init;
    T u_min = curve.domain_min();
    T u_max = curve.domain_max();
    
    // 确保初始值在参数域内
    u = std::clamp(u, u_min, u_max);
    
    constexpr T eps = std::numeric_limits<T>::epsilon();
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        // 求值和导数
        auto [C_u, C_prime] = curve.evaluate_with_derivative(u);
        
        // 残差向量
        Vec3<T> diff = C_u - point;
        
        // 空间距离检查（收敛准则1）
        T dist = diff.norm();
        if (dist < tolerance) {
            return ProjectionResult<T>{u, C_u, dist, iter + 1, true};
        }
        
        // 目标函数：f(u) = C'(u) · (C(u) - P)
        T f = dot(C_prime, diff);
        
        // 检查正交性收敛（收敛准则2）
        T c_prime_norm = C_prime.norm();
        if (std::abs(f) < tolerance * std::max(c_prime_norm, T{1})) {
            return ProjectionResult<T>{u, C_u, dist, iter + 1, true};
        }
        
        // 计算二阶导数
        Vec3<T> C_double_prime = curve.evaluate_second_derivative(u);
        
        // f'(u) = C''(u) · (C(u) - P) + C'(u) · C'(u)
        T f_prime = dot(C_double_prime, diff) + dot(C_prime, C_prime);
        
        // 数值稳定性检查：防止除以零
        if (std::abs(f_prime) < eps * 100) {
            // 雅可比矩阵奇异（可能处于拐点或尖点）
            // 尝试使用梯度下降法作为后备
            T step_size = T{0.01} * (u_max - u_min);
            T f_current = f;
            
            // 沿负梯度方向小步搜索
            T u_try = u - std::copysign(step_size, f);
            u_try = std::clamp(u_try, u_min, u_max);
            
            auto [C_try, C_prime_try] = curve.evaluate_with_derivative(u_try);
            Vec3<T> diff_try = C_try - point;
            T f_try = dot(C_prime_try, diff_try);
            
            if (std::abs(f_try) < std::abs(f_current)) {
                u = u_try;
                continue;
            } else {
                // 无法继续迭代
                return std::nullopt;
            }
        }
        
        // 牛顿迭代步
        T delta_u = -f / f_prime;
        
        // 步长限制：防止参数跳跃过大
        T max_step = (u_max - u_min) * T{0.5};
        if (std::abs(delta_u) > max_step) {
            delta_u = std::copysign(max_step, delta_u);
        }
        
        u += delta_u;
        
        // 参数域钳位
        u = std::clamp(u, u_min, u_max);
        
        // 参数收敛检查（收敛准则3）
        if (std::abs(delta_u) < tolerance * T{0.01}) {
            auto [C_final, _] = curve.evaluate_with_derivative(u);
            T final_dist = (C_final - point).norm();
            return ProjectionResult<T>{u, C_final, final_dist, iter + 1, true};
        }
    }
    
    // 未收敛
    return std::nullopt;
}

template <FloatingPoint T>
std::optional<ProjectionResult<T>> project_point_global(
    const NurbsCurve<T>& curve,
    const Vec3<T>& point,
    int num_samples,
    T tolerance)
{
    T u_min = curve.domain_min();
    T u_max = curve.domain_max();
    T step = (u_max - u_min) / static_cast<T>(num_samples - 1);
    
    std::optional<ProjectionResult<T>> best_result;
    T best_distance = std::numeric_limits<T>::max();
    
    // 均匀采样初始点
    for (int i = 0; i < num_samples; ++i) {
        T u_init = u_min + i * step;
        
        auto result = project_point_on_curve(curve, point, u_init, tolerance);
        
        if (result && result->converged && result->distance < best_distance) {
            best_distance = result->distance;
            best_result = result;
        }
    }
    
    return best_result;
}

// 显式实例化
template std::optional<ProjectionResult<double>> project_point_on_curve(
    const NurbsCurve<double>&, const Vec3<double>&, double, double, int);

template std::optional<ProjectionResult<float>> project_point_on_curve(
    const NurbsCurve<float>&, const Vec3<float>&, float, float, int);

template std::optional<ProjectionResult<double>> project_point_global(
    const NurbsCurve<double>&, const Vec3<double>&, int, double);

template std::optional<ProjectionResult<float>> project_point_global(
    const NurbsCurve<float>&, const Vec3<float>&, int, float);

} // namespace nurbslab::algorithms