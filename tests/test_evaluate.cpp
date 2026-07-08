#include <gtest/gtest.h>


#include "nurbslab/algorithms/evaluate.h"
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
NurbsCurve<double> CreateSemiCircleCurve1() {
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

// ==========================================
// 1. 求值测试 (Evaluation)
// ==========================================
TEST(CurveEvaluationTest, EvaluatePointOnSemiCircle) {
    auto curve = CreateSemiCircleCurve1();

    // 测试起点 (u=0)
    Vec3d start = curve.evaluate(0.0);

    EXPECT_NEAR(start.x, 1.0, 1e-6);
    EXPECT_NEAR(start.y, 0.0, 1e-6);
    EXPECT_NEAR(start.z, 0.0, 1e-6);

    // 测试终点 (u=1)
    Vec3d end = curve.evaluate(1.0);
    EXPECT_NEAR(end.x, -1.0, 1e-6);
    EXPECT_NEAR(end.y, 0.0, 1e-6);
    EXPECT_NEAR(end.z, 0.0, 1e-6);

    // 测试顶点 (u=0.5, 应该是 0, 1, 0)
    Vec3d mid = curve.evaluate(0.5);
    EXPECT_NEAR(mid.x, 0.0, 1e-6);
    EXPECT_NEAR(mid.y, 1.0, 1e-6);
    EXPECT_NEAR(mid.z, 0.0, 1e-6);
}

TEST(CurveEvaluationTest, EvaluateDerivative) {
    auto curve = CreateSemiCircleCurve1();
    
    // 在顶点处 (0,1,0)，切线应该是水平的 (-x 方向)
    Vec3d deriv = curve.evaluate_derivative(0.5);
    
    // 检查 Y 分量是否接近 0
    EXPECT_NEAR(deriv.y, 0.0, 1e-5);
    // 检查 X 分量是否为负
    EXPECT_LT(deriv.x, 0.0);
}
