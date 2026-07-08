#pragma once

#include "../concepts.h"
#include "../vec3.h"
#include "../knot_vector.h"
#include <vector>
#include <stdexcept>
#include <cmath>

namespace nurbslab::algorithms {

/// @brief NURBS 求值算法集合
template <FloatingPoint T = double>
class Evaluate {
public:
    /// @brief Cox-de Boor 递推公式计算 B 样条基函数
    /// @param i 基函数索引
    /// @param p 次数
    /// @param u 参数值
    /// @param knots 节点向量
    static T basis_function(int i, int p, T u, const KnotVector<T>& knots) {
        if (p == 0) {
            // 特殊处理：当 u 等于最后一个节点时，包含右端点
            if (i == static_cast<int>(knots.size()) - 2 && 
                std::abs(u - knots[i + 1]) < 1e-15) {
                return T{1};
            }
            return (u >= knots[i] && u < knots[i + 1]) ? T{1} : T{0};
        }

        T left = T{0}, right = T{0};
        T denom1 = knots[i + p] - knots[i];
        T denom2 = knots[i + p + 1] - knots[i + 1];

        if (std::abs(denom1) > 1e-15) {
            left = (u - knots[i]) / denom1 * basis_function(i, p - 1, u, knots);
        }
        if (std::abs(denom2) > 1e-15) {
            right = (knots[i + p + 1] - u) / denom2 * basis_function(i + 1, p - 1, u, knots);
        }
        return left + right;
    }

    /// @brief De Boor 算法求 NURBS 曲线上的点
    /// @param degree 曲线次数
    /// @param ctrl_pts 控制点
    /// @param weights 权重
    /// @param knots 节点向量
    /// @param u 参数值
    static Vec3<T> nurbs_curve_point(int degree, 
                                     const std::vector<Vec3<T>>& ctrl_pts, 
                                     const std::vector<T>& weights, 
                                     const KnotVector<T>& knots, 
                                     T u) 
    {
        int span = knots.find_span(u);
        std::vector<Vec3<T>> tmp(degree + 1);
        std::vector<T> tmp_w(degree + 1);

        // 提取相关控制点和权重
        for (int i = 0; i <= degree; ++i) {
            tmp[i] = ctrl_pts[span - degree + i];
            tmp_w[i] = weights[span - degree + i];
        }

        // De Boor 递推
        for (int r = 1; r <= degree; ++r) {
            for (int i = degree; i >= r; --i) {
                int idx = span - degree + i;
                T alpha = T{0};
                T denom = knots[idx + degree - r + 1] - knots[idx];
                if (std::abs(denom) > 1e-15) {
                    alpha = (u - knots[idx]) / denom;
                }
                
                // 齐次坐标下的线性插值
                tmp[i] = tmp[i - 1] * (T{1} - alpha) * tmp_w[i - 1] + tmp[i] * alpha * tmp_w[i];
                tmp_w[i] = tmp_w[i - 1] * (T{1} - alpha) + tmp_w[i] * alpha;
            }
        }

        // 除以权重得到真实三维坐标
        if (std::abs(tmp_w[degree]) < 1e-15) return {0, 0, 0};
        return tmp[degree] / tmp_w[degree];
    }
};

} // namespace nurbslab