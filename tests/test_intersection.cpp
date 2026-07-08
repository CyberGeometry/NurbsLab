#include <gtest/gtest.h>

#include "nurbslab/algorithms/intersection.h"
#include "nurbslab/nurbs_curve.h"
#include "nurbslab/knot_vector.h"
#include "nurbslab/vec3.h"

#include <vector>
#include <cmath>

using namespace nurbslab;
using namespace nurbslab::algorithms;

// ==========================================
// 辅助函数：创建一个简单的半圆弧 NURBS 曲线
// 用于测试求值、投影和求交
// ==========================================
NurbsCurve<double> CreateSemiCircleCurve2() {
    // 控制点 (表示一个半圆)
    std::vector<Vec3d> points = {
        {1.0, 0.0, 0.0},  // P0
        {1.0, 1.0, 0.0},  // P1
        {-1.0, 1.0, 0.0}, // P2
        {-1.0, 0.0, 0.0}  // P3
    };

    // 权重 (用于表示圆弧)
    std::vector<double> weights = {1.0, 0.70710678, 0.70710678, 1.0};

    // 节点向量 (Clamped, Degree 2)
    // 对于 4 个控制点，degree 2，节点数应为 4+2+1 = 7
    // 标准半圆节点: {0,0,0, 0.5, 1,1,1} (size 7)
    KnotVector<double> knots({0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0});

    return NurbsCurve<double>(2, points, weights, knots);
}


TEST(IntersectionTest, LineCurveIntersection) {
    // 曲线 1: 半圆弧
    auto curve1 = CreateSemiCircleCurve2();

    // 曲线 2: 一条直线 y = 0.5 (穿过半圆)
    std::vector<Vec3d> line_points = {
        {-2.0, 0.5, 0.0},
        { 2.0, 0.5, 0.0}
    };
    std::vector<double> line_weights = {1.0, 1.0};
    KnotVector<double> line_knots({0.0, 0.0, 1.0, 1.0});
    
    auto curve2 = NurbsCurve<double>(1, line_points, line_weights, line_knots);

    // 执行求交
    auto intersections = curve_curve_intersection(curve1, curve2, 1e-6, 20);

    // 直线 y=0.5 应该与半圆有两个交点
    ASSERT_EQ(intersections.size(), 2);

    bool found_p1 = false;
    bool found_p2 = false;

    for (const auto& res : intersections) {
        // 检查 Y 坐标
        EXPECT_NEAR(res.point.y, 0.5, 1e-4);
        
        // 检查 X 坐标 (正或负)
        if (std::abs(res.point.x - 0.866025) < 1e-4) found_p1 = true;
        if (std::abs(res.point.x + 0.866025) < 1e-4) found_p2 = true;
    }

    EXPECT_TRUE(found_p1);
    EXPECT_TRUE(found_p2);
}