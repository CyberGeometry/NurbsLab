#include "nurbslab/knot_vector.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

namespace nurbslab {

template <FloatingPoint T>
KnotVector<T>::KnotVector(std::span<const T> knots) 
    : knots_(knots.begin(), knots.end()) 
{
    if (!is_valid()) {
        throw std::invalid_argument("KnotVector: 节点序列必须非递减");
    }
}

template <FloatingPoint T>
KnotVector<T>::KnotVector(std::initializer_list<T> knots)
    : knots_(knots)
{
    if (!is_valid()) {
        throw std::invalid_argument("KnotVector: 节点序列必须非递减");
    }
}

template <FloatingPoint T>
KnotVector<T> KnotVector<T>::uniform(int degree, int num_control_points) {
    int n = num_control_points;
    int m = n + degree + 1; // 节点总数
    std::vector<T> knots(m);
    
    for (int i = 0; i < m; ++i) {
        knots[i] = static_cast<T>(i);
    }
    
    return KnotVector(knots);
}

template <FloatingPoint T>
KnotVector<T> KnotVector<T>::clamped(int degree, int num_control_points) {
    int n = num_control_points;
    int m = n + degree + 1;
    std::vector<T> knots(m);
    
    // 前 degree+1 个节点为 0
    for (int i = 0; i <= degree; ++i) {
        knots[i] = T{0};
    }
    
    // 中间节点均匀分布
    int num_internal = m - 2 * (degree + 1);
    for (int i = 1; i <= num_internal; ++i) {
        knots[degree + i] = static_cast<T>(i) / static_cast<T>(num_internal + 1);
    }
    
    // 后 degree+1 个节点为 1
    for (int i = 0; i <= degree; ++i) {
        knots[m - 1 - i] = T{1};
    }
    
    return KnotVector(knots);
}

template <FloatingPoint T>
int KnotVector<T>::find_span(T u) const {
    int n = static_cast<int>(knots_.size()) - 2;
    
    // 特殊情况：u 等于最后一个节点
    if (u >= knots_[n]) return n - 1;
    if (u <= knots_[0]) return 0;
    
    // 二分查找
    int low = 0;
    int high = n;
    int mid = (low + high) / 2;
    
    while (u < knots_[mid] || u >= knots_[mid + 1]) {
        if (u < knots_[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    
    return mid;
}

template <FloatingPoint T>
int KnotVector<T>::multiplicity(T u) const {
    constexpr T eps = std::numeric_limits<T>::epsilon() * 100;
    int mult = 0;
    for (const auto& k : knots_) {
        if (std::abs(k - u) < eps) {
            ++mult;
        }
    }
    return mult;
}

template <FloatingPoint T>
void KnotVector<T>::insert(T u) {
    auto it = std::lower_bound(knots_.begin(), knots_.end(), u);
    knots_.insert(it, u);
}

template <FloatingPoint T>
bool KnotVector<T>::is_valid() const {
    if (knots_.size() < 2) return false;
    
    // 检查非递减
    for (size_t i = 1; i < knots_.size(); ++i) {
        if (knots_[i] < knots_[i - 1]) {
            return false;
        }
    }
    return true;
}

template <FloatingPoint T>
bool KnotVector<T>::is_clamped(int degree) const {
    constexpr T eps = std::numeric_limits<T>::epsilon() * 100;
    
    if (static_cast<int>(knots_.size()) < 2 * (degree + 1)) return false;
    
    // 检查前 degree+1 个节点是否相同
    for (int i = 1; i <= degree; ++i) {
        if (std::abs(knots_[i] - knots_[0]) > eps) return false;
    }
    
    // 检查后 degree+1 个节点是否相同
    T last = knots_.back();
    for (int i = 1; i <= degree; ++i) {
        if (std::abs(knots_[knots_.size() - 1 - i] - last) > eps) return false;
    }
    
    return true;
}

// 显式实例化
template class KnotVector<double>;
template class KnotVector<float>;

} // namespace nurbslab