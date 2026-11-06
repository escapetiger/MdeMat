# 数学模型与代码定位

这里记录实际实现的离散方程。公式编号是定位锚点，方法名是实现入口；更改离散符号、缩放、边界或时间推进时应同时更新对应公式。

| 内容 | 文档 | 实现 |
| --- | --- | --- |
| Poisson、LDG/PDG 消元、周期约束 | [elliptic.md](elliptic.md) | `mdemat.solver.elliptic` |
| 迎风与半拉格朗日输运、扩散 | [transport-diffusion.md](transport-diffusion.md) | `mdemat.solver.hyperbolic`、`mdemat.solver.parabolic` |
| 辐射宏微分解、速度耦合、epsilon 缩放 | [radiation.md](radiation.md) | `mdemat.application.radiation` |
| Hermite 等离子体、Poisson、电场三角递推 | [plasma.md](plasma.md) | `mdemat.application.plasma` |
| 公共时间积分及代数约束 | [time-integration.md](time-integration.md) | `mdemat.common.ode` |
| 物理采样、分量范数与收敛阶 | [assessment.md](assessment.md) | `common.fem` 采样、`tools.analyze` 归约、应用层分布重构 |
| 低秩加权度量、动态秩及输运消费者 | [low-rank.md](low-rank.md) | `mdemat.common.lowrank` |
| 有限差分的物理符号与边界载荷 | [finite-difference.md](finite-difference.md) | `mdemat.common.fdm` |
| 独立无限介质输运参考 | [reference-transport.md](reference-transport.md) | `benchmarks/ganapol` |

公共记号：空间维数为 $d$，单元数为 $N_e$，局部自由度为 $N_l$，全局自由度为 $N_x=N_eN_l$。不连续空间为

$$
V_h=\{u_h:u_h|_K\in\mathcal Q_p(K)\text{ 或 }\mathcal P_p(K)\},\qquad
u_h|_K=\sum_{a=1}^{N_l}U_{K,a}\phi_a\circ F_K^{-1}. \tag{N1}
$$

`BH1OrthotopeElement.nodal/modal` 的 `order` 表示每方向基函数数；可表示的最高单方向次数为 `order-1`。模态 `pattern='P'/'Q'` 决定多项式截断。`FiniteElementSpace` 按单元连续存储局部系数，多个场按列存储。`U(:)` 依次排列整个第一场、整个第二场，不能按节点交错解包。

| 数据 | 数组形状 |
| --- | --- |
| 函数回调的坐标 `x` | `d × Npoints` |
| 多分量函数回调结果 | `Ncomponents × Npoints` |
| 空间系数 | `Nx × Ncomponents` |
| `XDisc.eval(xRef, coefficients)` | `(Ne*NpointsPerElement) × Ncomponents` |
| 速度基函数值 | `NvelocityDofs × NvelocityPoints` |
| 重构分布函数 | `NspacePoints × NvelocityPoints` |

参考积分节点为 $(\xi_q,w_q)$，$B_{aq}=\phi_a(\xi_q)$。空间局部质量为

$$
\widehat M=\overline B\operatorname{diag}(w)B^T,\qquad
M_K=|\det J_K|\widehat M. \tag{N2}
$$

装配中的 `Approximator.project` 已施加参考质量矩阵逆；`AdjointAssembly`、`UpwindAssembly`、`SourceAssembly` 等返回系数空间中的质量归一化算子/载荷。因此时间推进通常使用单位质量，辐射及等离子体中的 `M` 另表示物理缩放。不要再次对这些算子应用空间质量逆。

`FiniteElementAssembly.createMatrix/createVector` 负责局部到全局索引；体积分在 `assembleConstantVolumeMatrix/assembleVolumeMatrix` 中，面通量在 `assembleConstantFluxMatrix/assembleFluxMatrix` 中。`MultiplierAssembly.assembleMatrix(a)` 表示系数空间的乘法投影 $M_K^{-1}\int_K\overline{\phi_i}a\phi_j\,dx$，通常不同于节点系数逐点相乘。`WeightedAdjointAssembly` 是带 `weight/dWeight` 的加权派生装配；当前 AP Hermite 求解器使用普通 `AdjointAssembly` 加显式平衡电场乘法项。

Galerkin 测试函数采用复共轭：体积项、面通量、半拉格朗日转移及边界载荷均含 $\overline{\phi_i}$，试函数仍为 $\phi_j$。因此统一换基 $\phi'_i=i\phi_i$ 后系数 $U'=-iU$ 表示同一物理函数，质量归一化矩阵不变，载荷变为 $-i$ 倍；`TestComplexGalerkinGauge` 对 nodal/modal、回调/常数路径、加权算子和半拉格朗日边界分别验证此协变关系。物理积分 $\ell^TU=\int u_h\,dx$ 仍是对系数的线性函数，$\ell$ 不含测试函数共轭，周期 gauge 的下行使用非共轭转置。

`WeightedAdjointAssembly` 的第一行表示 $A^*=-c\partial_x(a\,\cdot)/a$，第二行表示 $A=c\,a\partial_x((\cdot)/a)$，因此权重恒为一时第一行与普通 `AdjointAssembly` 的符号相反。`AdjointAssembly` 的 tensor callback 返回 `d^2 × points`（列主序矩阵分量），标量样本解释为各向同性张量；常张量 callback 的所有块和 Dirichlet 载荷均与数值张量路径相等。变张量的 `G21` 定义为对 $D\,u$ 的保守弱导数，含面上 $Dn$ 系数，不能直接称为一般变系数的 $D\nabla u$；专用扩散入口使用显式常张量模型。

这些文档说明离散定义和已有回归所覆盖的行为；类名中的 AP、滤波器名称或一个算例的成功，不构成所有参数下稳定性、守恒、正性或渐近保持性的证明。
