#pragma once

#include "concepts.h"
#include <vector>
#include <span>
#include <stdexcept>

namespace nurbslab {

/// @brief NURBS 节点向量类
/// 节点向量是非递减的实数序列，定义了B样条基函数的支撑区间
template <FloatingPoint T = double>
class KnotVector {
public:
    /// @brief 默认构造
    KnotVector() = default;

    /// @brief 从已有节点序列构造
    /// @param knots 非递减节点序列
    /// @throws std::invalid_argument 若节点序列不满足非递减条件
    explicit KnotVector(std::span<const T> knots);

    /// @brief 从初始化列表构造
    KnotVector(std::initializer_list<T> knots);

    /// @brief 生成均匀节点向量
    /// @param degree 次数
    /// @param num_control_points 控制点数量
    /// @return 均匀节点向量
    static KnotVector uniform(int degree, int num_control_points);

    /// @brief 生成准均匀节点向量（端点重复度 = degree+1）
    /// @param degree 次数
    /// @param num_control_points 控制点数量
    /// @return 准均匀节点向量
    static KnotVector clamped(int degree, int num_control_points);

    // 访问器
    T operator[](size_t i) const { return knots_[i]; }
    size_t size() const { return knots_.size(); }
    T min() const { return knots_.front(); }
    T max() const { return knots_.back(); }
    std::span<const T> data() const { return knots_; }

    /// @brief 查找节点区间索引，使得 knots[i] <= u < knots[i+1]
    /// @param u 参数值
    /// @return 节点区间索引
    int find_span(T u) const;

    /// @brief 计算节点 u 的重数（重复次数）
    int multiplicity(T u) const;

    /// @brief 插入节点（用于细化）
    /// @param u 待插入节点值
    void insert(T u);

    /// @brief 验证节点向量的合法性
    bool is_valid() const;

    /// @brief 判断是否为 clamped 节点向量
    bool is_clamped(int degree) const;

private:
    std::vector<T> knots_;
};

// 显式实例化声明
extern template class KnotVector<double>;
extern template class KnotVector<float>;

} // namespace nurbslab