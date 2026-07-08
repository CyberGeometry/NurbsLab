#include "nurbslab/nurbs_curve.h"
#include <iostream>
#include <cmath>
#include <numbers>

int main() {
    // 使用 9 个控制点逼近完整圆 (4段二次圆弧)
    std::vector<nurbslab::Vec3d> ctrl_pts;
    std::vector<double> weights;
    
    double r = 2.0;
    double w_corner = 1.0;
    double w_mid = std::sqrt(2.0) / 2.0; // 中间点的有理权重

    for (int i = 0; i <= 8; ++i) {
        double angle = i * std::numbers::pi / 4.0;
        ctrl_pts.push_back({r * std::cos(angle), r * std::sin(angle), 0});
        weights.push_back((i % 2 == 0) ? w_corner : w_mid);
    }

    // 准均匀节点向量
    auto knots = nurbslab::KnotVector<double>::clamped(2, 9);
    nurbslab::NurbsCurve<double> circle(2, ctrl_pts, weights, knots);

    // 在 u=0.25 处求值，理论上应该是 (0, 2, 0)
    auto pt = circle.evaluate(0.25);
    std::cout << "Point at u=0.25: " << pt << "\n";
    std::cout << "Distance to origin: " << pt.norm() << "\n";

    return 0;
}