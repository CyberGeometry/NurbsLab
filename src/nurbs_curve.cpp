#include "nurbslab/nurbs_curve.h"
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <string>

namespace nurbslab {

template <FloatingPoint T>
NurbsCurve<T>::NurbsCurve(int degree,
                          std::span<const PointType> control_points,
                          std::span<const T> weights,
                          KnotVector<T> knots)
    : degree_(degree)
    , control_points_(control_points.begin(), control_points.end())
    , weights_(weights.begin(), weights.end())
    , knots_(std::move(knots))
{
    validate();
}

template <FloatingPoint T>
NurbsCurve<T>::NurbsCurve(int degree,
                          std::span<const PointType> control_points,
                          KnotVector<T> knots)
    : degree_(degree)
    , control_points_(control_points.begin(), control_points.end())
    , weights_(control_points.size(), T{1}) // 所有权重为1
    , knots_(std::move(knots))
{
    validate();
}

template <FloatingPoint T>
void NurbsCurve<T>::validate() const {
    if (degree_ < 0) {
        throw std::invalid_argument("NurbsCurve: degree must be >= 0");
    }
    
    int n = static_cast<int>(control_points_.size());
    if (n < degree_ + 1) {
        throw std::invalid_argument("NurbsCurve: not enough control points for given degree");
    }
    
    if (static_cast<int>(weights_.size()) != n) {
        throw std::invalid_argument("NurbsCurve: weights size must match control points size");
    }
    
    int expected_knots = n + degree_ + 1;
    if (static_cast<int>(knots_.size()) != expected_knots) {
        throw std::invalid_argument(
            "NurbsCurve: knot vector size must be n + p + 1 = " + 
            std::to_string(expected_knots));
    }
    
    // 检查权重为正
    for (size_t i = 0; i < weights_.size(); ++i) {
        if (weights_[i] <= T{0}) {
            throw std::invalid_argument("NurbsCurve: weights must be positive");
        }
    }
}

// ==================== Cox-de Boor 基函数计算 ====================

template <FloatingPoint T>
std::vector<T> NurbsCurve<T>::basis_functions(int span, T u) const {
    int p = degree_;
    std::vector<T> N(p + 1, T{0});
    std::vector<T> left(p + 1);
    std::vector<T> right(p + 1);
    
    N[0] = T{1};
    
    for (int j = 1; j <= p; ++j) {
        left[j] = u - knots_[span + 1 - j];
        right[j] = knots_[span + j] - u;
        
        T saved = T{0};
        for (int r = 0; r < j; ++r) {
            T denom = right[r + 1] + left[j - r];
            if (std::abs(denom) < std::numeric_limits<T>::epsilon()) {
                // 处理退化情况：节点重合导致分母为零
                N[r] = saved;
                saved = T{0};
                continue;
            }
            T temp = N[r] / denom;
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
    
    return N;
}

// ==================== De Boor 求值算法 ====================

template <FloatingPoint T>
typename NurbsCurve<T>::PointType NurbsCurve<T>::evaluate(T u) const {
    // 参数范围钳位
    u = std::clamp(u, domain_min(), domain_max());
    
    int span = knots_.find_span(u);
    auto N = basis_functions(span, u);
    
    // 有理曲线：在齐次坐标下计算
    PointType numerator{T{0}, T{0}, T{0}};
    T denominator = T{0};
    
    for (int i = 0; i <= degree_; ++i) {
        int idx = span - degree_ + i;
        T w = weights_[idx] * N[i];
        numerator = numerator + control_points_[idx] * w;
        denominator += w;
    }
    
    // 防止除以零（理论上不会出现，因为权重为正且基函数非负）
    if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
        return control_points_[span]; // 退化处理
    }
    
    return numerator / denominator;
}

// ==================== 一阶导数 ====================

template <FloatingPoint T>
typename NurbsCurve<T>::PointType NurbsCurve<T>::evaluate_derivative(T u) const {
    u = std::clamp(u, domain_min(), domain_max());
    
    // 对于有理曲线，使用商法则：
    // C(u) = A(u) / w(u)
    // C'(u) = (A'(u) * w(u) - A(u) * w'(u)) / w(u)^2
    
    int span = knots_.find_span(u);
    
    // 计算齐次坐标下的控制点：Pw_i = (w_i * P_i, w_i)
    // 对非有理曲线简化处理
    if (!is_rational()) {
        // 非有理曲线：直接对控制点差分求导
        // C'(u) = sum_{i=0}^{p-1} N'_{i,p-1}(u) * (P_{i+1} - P_i) * p / (knots_{i+p+1} - knots_{i+1})
        
        // 简化实现：使用中心差分（生产环境应使用精确导数公式）
        T h = T{1e-8};
        T u_min = domain_min();
        T u_max = domain_max();
        
        T u_plus = std::min(u + h, u_max);
        T u_minus = std::max(u - h, u_min);
        
        PointType p_plus = evaluate(u_plus);
        PointType p_minus = evaluate(u_minus);
        
        return (p_plus - p_minus) / (u_plus - u_minus);
    }
    
    // 有理曲线：使用中心差分（简化实现）
    T h = T{1e-8};
    T u_min = domain_min();
    T u_max = domain_max();
    
    T u_plus = std::min(u + h, u_max);
    T u_minus = std::max(u - h, u_min);
    
    PointType p_plus = evaluate(u_plus);
    PointType p_minus = evaluate(u_minus);
    
    return (p_plus - p_minus) / (u_plus - u_minus);
}

// ==================== 二阶导数 ====================

template <FloatingPoint T>
typename NurbsCurve<T>::PointType NurbsCurve<T>::evaluate_second_derivative(T u) const {
    u = std::clamp(u, domain_min(), domain_max());
    
    // 使用中心差分计算二阶导数
    T h = T{1e-5}; // 二阶导数需要稍大的步长
    T u_min = domain_min();
    T u_max = domain_max();
    
    T u_plus = std::min(u + h, u_max);
    T u_minus = std::max(u - h, u_min);
    
    PointType p_plus = evaluate(u_plus);
    PointType p_center = evaluate(u);
    PointType p_minus = evaluate(u_minus);
    
    // 二阶中心差分：f''(x) ≈ (f(x+h) - 2f(x) + f(x-h)) / h^2
    T h_actual = (u_plus - u_minus) / T{2};
    return (p_plus - p_center * T{2} + p_minus) / (h_actual * h_actual);
}

// ==================== 同时求值和导数 ====================

template <FloatingPoint T>
std::pair<typename NurbsCurve<T>::PointType, typename NurbsCurve<T>::PointType>
NurbsCurve<T>::evaluate_with_derivative(T u) const {
    return {evaluate(u), evaluate_derivative(u)};
}

// ==================== 辅助函数 ====================

template <FloatingPoint T>
bool NurbsCurve<T>::is_rational() const {
    constexpr T eps = std::numeric_limits<T>::epsilon() * 100;
    for (const auto& w : weights_) {
        if (std::abs(w - T{1}) > eps) return true;
    }
    return false;
}

// 显式实例化
template class NurbsCurve<double>;
template class NurbsCurve<float>;

} // namespace nurbslab