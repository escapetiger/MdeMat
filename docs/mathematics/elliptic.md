# Poisson 的 LDG 与 PDG

连续问题采用正椭圆符号：

$$
-\Delta u=f\quad\text{in }\Omega,\qquad q=\nabla u. \tag{P1}
$$

面向数值内核的新入口将问题数据与离散空间显式分开：

```matlab
problem = mdemat.solver.elliptic.poissonProblem(@(x) -2+zeros(1,size(x,2)), ...
    boundaryType="dirichlet", boundaryValue=@(x) x(1,:).^2);
[solution, report] = mdemat.solver.elliptic.solvePoisson(problem, space, method="ldg");
values = space.eval(referenceNodes, solution.U);
gradients = space.eval(referenceNodes, solution.Q);
```

`poissonProblem` 只持有物理源项和边界数据，`space` 是调用方构造的 FEM 离散空间；`method="pdg"` 选择消元系统。`solvePoisson` 直接组合公共装配与线性代数内核，不构造 `Scheme/State/ConfigParser/Visualizer`。返回的 `solution` 是全部系数的纯记录，`report` 含独立残差、物理均值、周期相容性积分、乘子、规模与阶段时间，不持有可变空间/工作区句柄。

`bc=[]` 选择周期配对；非空 `bc(x)` 选择 Dirichlet $u=g_D$。`exact.U/Q` 仅用于分析和绘图，可以省略。$U\in\mathbb R^{N_x}$，$Q\in\mathbb R^{N_x\times d}$。`Q(:,r)` 是第 $r$ 个正梯度分量。参考梯度回调应返回 `d × Npoints`，而非其转置。

令 `AdjointAssembly.assembleMatrix(I)` 返回的块为 $D_{ab}^{(r)}$。具体定义是质量归一化的“面通量减体积分”，不是直接采用连续微分矩阵：

$$
D_{12}^{(r)}=G_{F,12}^{(r)}-G_{V,12}^{(r)},\quad
D_{21}^{(r)}=G_{F,21}^{(r)}-G_{V,21}^{(r)},\quad
D_{11}^{(r)}=G_{F,11}^{(r)}. \tag{P2}
$$

`xPenaltyType` 是 $2\times2$ 的通量选择表；默认 `{'boundary','right';'left',''}` 使用成对左右单边通量，并在 Dirichlet 边界提供必要罚项。`central`、`all/allw`、`boundary/boundaryw` 的系数由 `computeInternalFluxPenalty`、`computeDirichletBoundaryFluxPenalty` 给出；改变这些选项会改变 $D_{11}$ 和伴随块，不能把所有配置都称为同一个对称内罚方法。

`D11` 的内部自身块为负、邻居块为正，是负跳跃生成元。因此正椭圆方程必须对它取负；扩散方程直接使用其负耗散符号。LDG 的实际线性系统为

$$
\begin{aligned}
-\sum_rD_{11}^{(r)}U-\sum_rD_{12}^{(r)}Q_r
 &=F+\sum_rb_{11}^{(r)},\\
-D_{21}^{(r)}U+Q_r&=b_{21}^{(r)}.
\end{aligned} \tag{P3}
$$

$F=\mathrm{SourceAssembly}(f)$，$b_{ab}^{(r)}=\mathrm{AdjointAssembly.vector}(I,g_D)_{ab}^{(r)}$；周期情形这些边界载荷为零。完整未知量是 `[U;Q1;...;Qd]`，总长度 $(d+1)N_x$。`LdgScheme.step` 必须保留 `reshape(W,Nx,d+1)` 的全部 `2:end` 列。

消去 $Q_r=D_{21}^{(r)}U+b_{21}^{(r)}$ 得

$$
K U=F+b_E,\quad
K=\sum_r(-D_{11}^{(r)}-D_{12}^{(r)}D_{21}^{(r)}),\quad
b_E=\sum_r(D_{12}^{(r)}b_{21}^{(r)}+b_{11}^{(r)}). \tag{P4}
$$

`EllipticAssembly.assembleMatrix(I)` 返回的是 **$G=-K$**，即扩散生成元；`PdgScheme` 对它取负得到 $K$，并将 `assembleVector` 返回的 $b_E$ 加到右端。这个符号也解释了扩散求解器为什么直接使用 $G$。Dirichlet 的 `assembleVector` 依赖之前 `assembleMatrix` 缓存的 $D_{12}$。

周期问题要求相容性 $\int_\Omega f\,dx=0$，并额外指定零均值。设常数函数系数为 $c\in\mathbb R^{N_x}$，物理积分权为

$$
\ell_{K,a}=|\det J_K|\sum_qw_q\phi_a(\xi_q),\qquad
\int_\Omega u_h\,dx=\ell^TU. \tag{P5}
$$

质量归一化方程的增广形式是

$$
\begin{bmatrix}K&c\\\ell^T&0\end{bmatrix}
\begin{bmatrix}U\\\lambda\end{bmatrix}
=\begin{bmatrix}F\\0\end{bmatrix}. \tag{P6}
$$

LDG 将 $c,\ell$ 在梯度对应位置补零。非均匀网格上 $c$ 和 $\ell$ 不相同，不能使用 `c.'` 代替物理均值行。对于不相容载荷，乘子会吸收常数分量；得到可解增广系统并不意味着原周期 PDE 已满足相容性。

| 公式 | 代码入口 |
| --- | --- |
| P1、边界类别 | `PoissonScheme.configure` |
| P2、通量选择 | `AdjointAssembly.assembleMatrix`、`computePenalty` 及其体/面装配方法 |
| P3、未知量形状 | `LdgScheme.setLhsTerm/setRhsTerm/step` |
| P4、消元、符号 | `EllipticAssembly.assembleMatrix/assembleVector`、`PdgScheme.setLhsTerm/setRhsTerm` |
| P5–P6、周期增广 | `LdgScheme.setLhsTerm`、`PdgScheme.setLhsTerm` |
| P2–P6、独立数值入口 | `mdemat.solver.elliptic.poissonProblem/solvePoisson` |

回归入口为 `tests/regression/solvers/TestSolverContracts.m`：独立多分量解包、二维制造解 $u=x^2+2y^2+3x-4y+1/2$（$f=-6$、$q=(2x+3,4y-4)$）、LDG/PDG 一致性与非均匀网格周期零均值。制造解同时检查梯度符号与两个分量，避免仅比较两个可能共享错误的实现。

`TestPoissonInterface` 对新的纯入口分别验证一维/二维 Dirichlet 制造解、非均匀周期四次制造解、全部梯度分量、离散不相容载荷拒绝，以及三层网格上独立高阶物理积分测得的正弦解收敛。新入口默认拒绝超过 `compatibilityTolerance` 的周期离散 forcing 积分，避免把乘子吸收的常数缺陷当作原 PDE 的解；旧 Scheme 的兼容行为仍在 P6 后明确说明。`residualTolerance` 与相容性阈值是求解验收配置，参考解不属于求解输入。

`TestDgAssemblyIdentities` 还直接代入解析多项式验证体/面弱形式，并用弱梯度核内的末单元最高 Legendre 模态验证边界罚项的耗散符号。只有制造解连续时跳跃为零，无法暴露反号罚项。内部单元集合为空时显式跳过装配，避免将空选择误当作整个网格重复添加通量。
