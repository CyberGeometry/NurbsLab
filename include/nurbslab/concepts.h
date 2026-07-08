#pragma once

#include <concepts>
#include <type_traits>

namespace nurbslab {

/// @brief 浮点数概念，约束模板只接受 float/double/long double
template <typename T>
concept FloatingPoint = std::is_floating_point_v<T>;

/// @brief 可求值几何对象概念
template <typename Geo, typename T>
concept Evaluatable = requires(Geo g, T t) {
    { g.evaluate(t) } -> std::same_as<typename Geo::PointType>;
};

/// @brief 可求导几何对象概念
template <typename Geo, typename T>
concept Differentiable = requires(Geo g, T t) {
    { g.evaluate(t) } -> std::same_as<typename Geo::PointType>;
    { g.evaluate_derivative(t) } -> std::same_as<typename Geo::PointType>;
};

} // namespace nurbslab