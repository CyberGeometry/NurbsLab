# NurbsLab 数学笔记

## 1. Cox-de Boor 递推公式
NURBS 的基函数 $N_{i,p}(u)$ 由以下递推关系定义：

$$
N_{i,0}(u) = 
\begin{cases} 
1 & \text{if } u_i \le u < u_{i+1} \\ 
0 & \text{otherwise} 
\end{cases}
$$

$$
N_{i,p}(u) = \frac{u - u_i}{u_{i+p} - u_i} N_{i,p-1}(u) + \frac{u_{i+p+1} - u}{u_{i+p+1} - u_{i+1}} N_{i+1,p-1}(u)
$$

> **注意**：当分母为 0 时，该项定义为 0。

## 2. De Boor 算法
给定参数 $u$，首先找到所在的节点区间 $[u_i, u_{i+1})$。
对于非有理 B 样条，递推公式为：
$$ P_i^{(r)} = (1 - \alpha_i^{(r)}) P_{i-1}^{(r-1)} + \alpha_i^{(r)} P_i^{(r-1)} $$
其中 $\alpha_i^{(r)} = \frac{u - u_i}{u_{i+p-r+1} - u_i}$。

对于 **NURBS**，需要在齐次坐标下进行线性插值，最后一步进行透视除法（除以权重 $w$）。