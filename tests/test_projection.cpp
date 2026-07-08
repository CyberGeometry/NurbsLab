#include <gtest/gtest.h>

#include "nurbslab/algorithms/projection.h"
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
NurbsCurve<double> CreateSemiCircleCurve3() {
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


TEST(ProjectionTest, ProjectPointGlobal) {
    auto curve = CreateSemiCircleCurve3();

    // 待投影点：位于圆心上方 (0, 2, 0)
    // 理论投影点应该是圆弧顶点 (0, 1, 0)，参数 u = 0.5
    Vec3d point_to_project{0.0, 2.0, 0.0};

    auto result_opt = project_point_global(curve, point_to_project, 20, 1e-6);

    ASSERT_TRUE(result_opt.has_value());
    
    auto result = result_opt.value();
    
    // 检查参数 u 是否接近 0.5
    EXPECT_NEAR(result.u, 0.5, 1e-4);
    
    // 检查投影点坐标是否接近 (0, 1, 0)
    EXPECT_NEAR(result.point.x, 0.0, 1e-5);
    EXPECT_NEAR(result.point.y, 1.0, 1e-5);
    EXPECT_NEAR(result.point.z, 0.0, 1e-5);
    
    // 检查距离 (点到圆弧顶点的距离应为 1.0)
    EXPECT_NEAR(result.distance, 1.0, 1e-5);
}

TEST(ProjectionTest, ProjectPointNewton) {
    auto curve = CreateSemiCircleCurve3();
    
    // 待投影点：(0.5, 1.5, 0)
    // 初始猜测 u = 0.4 (靠近顶点)
    Vec3d point{0.5, 1.5, 0.0};

    auto result_opt = project_point_on_curve(curve, point, 0.4, 1e-6, 50);

    ASSERT_TRUE(result_opt.has_value());
    EXPECT_TRUE(result_opt->converged);
}