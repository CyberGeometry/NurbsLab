#pragma once

#include "concepts.h"
#include <cmath>
#include <iostream>

namespace nurbslab {

/// @brief 三维向量/点，支持基本运算
template <FloatingPoint T = double>
struct Vec3 {
    T x{0}, y{0}, z{0};

    constexpr Vec3() = default;
    constexpr Vec3(T x, T y, T z) : x(x), y(y), z(z) {}

    /// 向量加法
    constexpr Vec3 operator+(const Vec3& rhs) const {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }

    /// 向量减法
    constexpr Vec3 operator-(const Vec3& rhs) const {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }

    /// 标量乘法
    constexpr Vec3 operator*(T scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    /// 标量除法
    constexpr Vec3 operator/(T scalar) const {
        return {x / scalar, y / scalar, z / scalar};
    }

    /// 点积
    constexpr T dot(const Vec3& rhs) const {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    /// 叉积
    constexpr Vec3 cross(const Vec3& rhs) const {
        return {
            y * rhs.z - z * rhs.y,
            z * rhs.x - x * rhs.z,
            x * rhs.y - y * rhs.x
        };
    }

    /// 向量长度
    T norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    /// 向量长度的平方（避免开方，性能优化）
    constexpr T norm_squared() const {
        return x * x + y * y + z * z;
    }

    /// 单位化
    Vec3 normalized() const {
        T len = norm();
        if (len < 1e-15) return {0, 0, 0};
        return *this / len;
    }

    /// 距离
    T distance_to(const Vec3& other) const {
        return (*this - other).norm();
    }

    /// 输出
    friend std::ostream& operator<<(std::ostream& os, const Vec3& v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }

    /// 相等判断（带容差）
    bool equals(const Vec3& other, T tol = 1e-9) const {
        return std::abs(x - other.x) < tol &&
               std::abs(y - other.y) < tol &&
               std::abs(z - other.z) < tol;
    }
};

/// 类型别名
using Vec3d = Vec3<double>;
using Vec3f = Vec3<float>;

/// 全局点积函数
template <FloatingPoint T>
constexpr T dot(const Vec3<T>& a, const Vec3<T>& b) {
    return a.dot(b);
}

/// 全局叉积函数
template <FloatingPoint T>
constexpr Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) {
    return a.cross(b);
}

} // namespace nurbslab