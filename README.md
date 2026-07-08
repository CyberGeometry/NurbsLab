# NurbsLab

**NurbsLab** 是一个现代化的 C++20 几何算法库，专注于 **NURBS (非均匀有理B样条)** 曲线的底层数学计算与几何操作。

本库采用现代 C++ 范式（Concepts, `std::span`），提供类型安全、高性能且数学严谨的几何内核，适用于 CAD/CAM、机器人轨迹规划、计算机图形学及物理仿真等领域。

---

##  核心特性

- **现代 C++20 架构**：使用 `Concepts` 约束模板参数，利用 `std::span` 实现零拷贝视图，拒绝内存泄漏。
- **工业级数学内核**：基于 Cox-de Boor 递推与 De Boor 算法，支持任意阶次 NURBS 曲线的精确求值。
- **丰富的几何算法**：内置曲线求导、点到曲线投影（牛顿迭代）、曲线求交等高级算法。
- **高精度与类型安全**：支持 `float` / `double` / `long double` 模板切换，内置容差与合法性校验。
- **轻量与易用**：Header-friendly 设计，CMake 一键集成，开箱即用。

---

##  项目结构

```text
NurbsLab/
├── include/nurbslab/      # 核心头文件 (Concepts, Vec3, KnotVector, NurbsCurve)
├── src/                   # 算法实现 (求值, 投影, 求交)
├── tests/                 # 单元测试 (Catch2 / GTest)
├── examples/              # 示例代码 (如：圆逼近)
├── docs/                  # 数学公式推导与文档
├── CMakeLists.txt         # 构建配置
└── README.md
```

---

##  快速开始

### 1. 环境要求
- C++20 兼容编译器 (GCC 10+, Clang 10+, MSVC 2019+)
- CMake >= 3.20

### 2. 构建与测试
```bash
git clone https://github.com/yourusername/NurbsLab.git
cd NurbsLab
mkdir build && cd build
cmake .. -DNURBSLAB_BUILD_TESTS=ON
make -j$(nproc)
ctest --output-on-failure
```

### 3. 代码示例

```cpp
#include "nurbslab/nurbs_curve.h"
#include <iostream>

int main() {
    // 1. 定义控制点与权重
    std::vector<nurbslab::Vec3d> ctrl_pts = {
        {0, 0, 0}, {1, 2, 0}, {3, 2, 0}, {4, 0, 0}
    };
    std::vector<double> weights = {1.0, 1.0, 1.0, 1.0};

    // 2. 生成准均匀节点向量 (Clamped)
    auto knots = nurbslab::KnotVector<double>::clamped(3, 4);

    // 3. 构造三次 NURBS 曲线
    nurbslab::NurbsCurve<double> curve(3, ctrl_pts, weights, knots);

    // 4. 求值与求导
    auto [point, tangent] = curve.evaluate_with_derivative(0.5);
    
    std::cout << "Point at u=0.5: " << point << "\n";
    std::cout << "Tangent: " << tangent << "\n";

    return 0;
}
```

---

##  支持的算法

| 模块 | 功能 | 状态 |
| :--- | :--- | :--- |
| **Evaluate** | De Boor 求值、一阶/二阶导数 |  已完成 |
| **Projection** | 空间点到 NURBS 曲线的最近点投影 |  已完成 |
| **Intersection** | 曲线-曲线求交 (空间/平面) |  开发中 |
| **Knot Vector** | 节点插入、细化、合法性校验 |  已完成 |

---

##  贡献指南

欢迎提交 Issue 和 Pull Request！在提交代码前，请确保：
1. 新增算法必须附带单元测试。
2. 遵循现有的代码注释规范（Doxygen 风格）。
3. 数学公式请在 `docs/math_notes.md` 中给出推导过程。

---

##  License

本项目基于 **MIT License** 开源。详见 [LICENSE](LICENSE) 文件。
