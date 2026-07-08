#include "nurbslab/algorithms/intersection.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace nurbslab::algorithms {

// ==================== AABB 实现 ====================

template <FloatingPoint T>
bool AABB<T>::overlaps(const AABB& other) const {
    // 分离轴定理：任一轴不重叠则不相交
    if (max_bound.x < other.min_bound.x || min_bound.x > other.max_bound.x) return false;
    if (max_bound.y < other.min_bound.y || min_bound.y > other.max_bound.y) return false;
    if (max_bound.z < other.min_bound.z || min_bound.z > other.max_bound.z) return false;
    return true;
}

template <FloatingPoint T>
Vec3<T> AABB<T>::center() const {
    return (min_bound + max_bound) / T{2};
}

template <FloatingPoint T>
AABB<T> compute_aabb(const NurbsCurve<T>& curve) {
    auto cpts = curve.control_points();
    
    Vec3<T> min_bound{
        std::numeric_limits<T>::max(),
        std::numeric_limits<T>::max(),
        std::numeric_limits<T>::max()
    };
    Vec3<T> max_bound{
        std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::lowest()
    };
    
    for (const auto& pt : cpts) {
        min_bound.x = std::min(min_bound.x, pt.x);
        min_bound.y = std::min(min_bound.y, pt.y);
        min_bound.z = std::min(min_bound.z, pt.z);
        max_bound.x = std::max(max_bound.x, pt.x);
        max_bound.y = std::max(max_bound.y, pt.y);
        max_bound.z = std::max(max_bound.z, pt.z);
    }
    
    return {min_bound, max_bound};
}

template <FloatingPoint T>
AABB<T> compute_aabb_segment(const NurbsCurve<T>& curve, T u_start, T u_end, int samples) {
    Vec3<T> min_bound{
        std::numeric_limits<T>::max(),
        std::numeric_limits<T>::max(),
        std::numeric_limits<T>::max()
    };
    Vec3<T> max_bound{
        std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::lowest()
    };
    
    T step = (u_end - u_start) / static_cast<T>(samples);
    for (int i = 0; i <= samples; ++i) {
        T u = u_start + i * step;
        auto pt = curve.evaluate(u);
        
        min_bound.x = std::min(min_bound.x, pt.x);
        min_bound.y = std::min(min_bound.y, pt.y);
        min_bound.z = std::min(min_bound.z, pt.z);
        max_bound.x = std::max(max_bound.x, pt.x);
        max_bound.y = std::max(max_bound.y, pt.y);
        max_bound.z = std::max(max_bound.z, pt.z);
    }
    
    return {min_bound, max_bound};
}

// ==================== 递归求交核心 ====================

namespace detail {

template <FloatingPoint T>
struct ParamInterval {
    T start;
    T end;
    T mid() const { return (start + end) / T{2}; }
    T length() const { return end - start; }
};

/// @brief 递归细分求交
template <FloatingPoint T>
void subdivide_intersect(
    const NurbsCurve<T>& curve1,
    const NurbsCurve<T>& curve2,
    ParamInterval<T> interval1,
    ParamInterval<T> interval2,
    T tolerance,
    int depth,
    int max_depth,
    std::vector<IntersectionResult<T>>& results)
{
    // 终止条件：细分到足够小
    if (depth >= max_depth || 
        (interval1.length() < tolerance && interval2.length() < tolerance)) {
        
        // 牛顿迭代精化
        T u = interval1.mid();
        T v = interval2.mid();
        
        // 简化牛顿迭代：固定步数
        for (int iter = 0; iter < 20; ++iter) {
            auto [P1, dP1] = curve1.evaluate_with_derivative(u);
            auto [P2, dP2] = curve2.evaluate_with_derivative(v);
            
            Vec3<T> diff = P1 - P2;
            T dist = diff.norm();
            
            if (dist < tolerance) {
                // 检查是否重复
                bool is_duplicate = false;
                for (const auto& r : results) {
                    if (std::abs(r.u - u) < tolerance * 10 && 
                        std::abs(r.v - v) < tolerance * 10) {
                        is_duplicate = true;
                        break;
                    }
                }
                if (!is_duplicate) {
                    Vec3<T> mid_point = (P1 + P2) / T{2};
                    results.push_back({u, v, mid_point});
                }
                return;
            }
            
            // 构建 2x2 雅可比矩阵
            // J = [dP1/du, -dP2/dv]
            // 求解 J * [du, dv]^T = -diff
            
            T a = dot(dP1, dP1);
            T b = -dot(dP1, dP2);
            T c = dot(dP2, dP2);
            T e = -dot(dP1, diff);
            T f = dot(dP2, diff);
            
            T det = a * c - b * b;
            if (std::abs(det) < std::numeric_limits<T>::epsilon() * 100) {
                break; // 奇异，停止迭代
            }
            
            T delta_u = (c * e - b * f) / det;
            T delta_v = (a * f - b * e) / det;
            
            u += delta_u;
            v += delta_v;
            
            // 钳位到参数域
            u = std::clamp(u, curve1.domain_min(), curve1.domain_max());
            v = std::clamp(v, curve2.domain_min(), curve2.domain_max());
        }
        
        return;
    }
    
    // 计算包围盒
    AABB<T> box1 = compute_aabb_segment(curve1, interval1.start, interval1.end);
    AABB<T> box2 = compute_aabb_segment(curve2, interval2.start, interval2.end);
    
    // 包围盒不相交则剪枝
    if (!box1.overlaps(box2)) {
        return;
    }
    
    // 选择较长区间进行细分
    if (interval1.length() >= interval2.length()) {
        ParamInterval<T> left1{interval1.start, interval1.mid()};
        ParamInterval<T> right1{interval1.mid(), interval1.end};
        
        subdivide_intersect(curve1, curve2, left1, interval2, tolerance, depth + 1, max_depth, results);
        subdivide_intersect(curve1, curve2, right1, interval2, tolerance, depth + 1, max_depth, results);
    } else {
        ParamInterval<T> left2{interval2.start, interval2.mid()};
        ParamInterval<T> right2{interval2.mid(), interval2.end};
        
        subdivide_intersect(curve1, curve2, interval1, left2, tolerance, depth + 1, max_depth, results);
        subdivide_intersect(curve1, curve2, interval1, right2, tolerance, depth + 1, max_depth, results);
    }
}

} // namespace detail

template <FloatingPoint T>
std::vector<IntersectionResult<T>> curve_curve_intersection(
    const NurbsCurve<T>& curve1,
    const NurbsCurve<T>& curve2,
    T tolerance,
    int max_subdivisions)
{
    std::vector<IntersectionResult<T>> results;
    
    // 全局包围盒快速检测
    AABB<T> global_box1 = compute_aabb(curve1);
    AABB<T> global_box2 = compute_aabb(curve2);
    
    if (!global_box1.overlaps(global_box2)) {
        return results; // 无交点
    }
    
    detail::ParamInterval<T> interval1{curve1.domain_min(), curve1.domain_max()};
    detail::ParamInterval<T> interval2{curve2.domain_min(), curve2.domain_max()};
    
    detail::subdivide_intersect(
        curve1, curve2, interval1, interval2, 
        tolerance, 0, max_subdivisions, results);
    
    return results;
}

// 显式实例化
template AABB<double> compute_aabb(const NurbsCurve<double>&);
template AABB<float> compute_aabb(const NurbsCurve<float>&);
template AABB<double> compute_aabb_segment(const NurbsCurve<double>&, double, double, int);
template AABB<float> compute_aabb_segment(const NurbsCurve<float>&, float, float, int);

template std::vector<IntersectionResult<double>> curve_curve_intersection(
    const NurbsCurve<double>&, const NurbsCurve<double>&, double, int);
template std::vector<IntersectionResult<float>> curve_curve_intersection(
    const NurbsCurve<float>&, const NurbsCurve<float>&, float, int);

} // namespace nurbslab::algorithms