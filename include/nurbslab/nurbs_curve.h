#pragma once

#include "concepts.h"
#include "vec3.h"
#include "knot_vector.h"
#include <vector>
#include <span>
#include <memory>

namespace nurbslab {

/// @brief NURBS 曲线类
/// 非均匀有理B样条曲线，工业级几何内核的核心数据结构
template <FloatingPoint T = double>
class NurbsCurve {
public:
    using PointType = Vec3<T>;

    /// @brief 默认构造
    NurbsCurve() = default;

    /// @brief 完整构造 NURBS 曲线
    /// @param degree 曲线次数（通常 2-5）
    /// @param control_points 控制点序列
    /// @param weights 权重序列（与控制点一一对应）
    /// @param knots 节点向量
    /// @throws std::invalid_argument 参数不匹配时抛出
    NurbsCurve(int degree,
               std::span<const PointType> control_points,
               std::span<const T> weights,
               KnotVector<T> knots);

    /// @brief 简化构造：B样条曲线（所有权重为1）
    NurbsCurve(int degree,
               std::span<const PointType> control_points,
               KnotVector<T> knots);

    // ==================== 核心求值算法 ====================

    /// @brief De Boor 算法求值
    /// @param u 参数值，必须在 [knots_[degree], knots_[m-degree-1]] 范围内
    /// @return 曲线上的点
    PointType evaluate(T u) const;

    /// @brief 一阶导数（切向量）
    /// @param u 参数值
    /// @return 一阶导数向量
    PointType evaluate_derivative(T u) const;

    /// @brief 二阶导数（用于牛顿迭代等）
    PointType evaluate_second_derivative(T u) const;

    /// @brief 同时求值和一阶导数（性能优化，避免重复计算）
    std::pair<PointType, PointType> evaluate_with_derivative(T u) const;

    // ==================== 访问器 ====================

    int degree() const { return degree_; }
    int num_control_points() const { return static_cast<int>(control_points_.size()); }
    const KnotVector<T>& knots() const { return knots_; }
    std::span<const PointType> control_points() const { return control_points_; }
    std::span<const T> weights() const { return weights_; }

    /// 参数域
    T domain_min() const { return knots_[degree_]; }
    T domain_max() const { return knots_[knots_.size() - degree_ - 1]; }

    /// 是否为有理曲线（存在非1权重）
    bool is_rational() const;

private:
    int degree_{0};
    std::vector<PointType> control_points_;
    std::vector<T> weights_;
    KnotVector<T> knots_;

    /// @brief 计算B样条基函数值（Cox-de Boor递推）
    /// @param span 节点区间索引
    /// @param u 参数值
    /// @return 基函数值数组 N_{i,p}(u)，i 从 span-p 到 span
    std::vector<T> basis_functions(int span, T u) const;

    /// @brief 验证构造参数的合法性
    void validate() const;
};

// 显式实例化声明
extern template class NurbsCurve<double>;
extern template class NurbsCurve<float>;

} // namespace nurbslab