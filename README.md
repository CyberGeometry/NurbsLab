#  NurbsLab

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
├── CMakeLists.txt         # 构建配置
└── README.md
