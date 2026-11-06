# DG–Hermite 等离子体求解器

当前 `DghsScheme` 与 `ApDghsScheme` 是一维空间、一维速度的 Hermite 展开，Poisson 路径使用周期约束。基类存在 `bc` 字段不代表这里已经实现 Dirichlet 电势/动理学流入边界。

## 标准 Hermite 表示

以温度 $T_0>0$ 的 Maxwellian–Hermite 基 $\psi_k(v)$ 展开

$$
f_h(x,v)=\sum_{k=0}^{N_h-1}D_k(x)\psi_k(v),\qquad
v\psi_k=\sqrt{T_0}\big(\sqrt{k+1}\psi_{k+1}+\sqrt{k}\psi_{k-1}\big),\quad
\partial_v\psi_k=-\sqrt{(k+1)/T_0}\psi_{k+1}. \tag{V1}
$$

`L2OrthotopeElement.hermite` 构造速度基，`HermiteState.Dofs.D` 为 `Nx × Nh`，MATLAB 列 `k+1` 对应数学模态 $k$。`distribution` 将空间系数先评估，再使用速度基重构。末端采用截断 $D_{N_h}=0$。

模型的尺度约定是

$$
\varepsilon\partial_tf+v\partial_xf+E\partial_vf
=\tau_0^{-1}\partial_v(vf+T_0\partial_vf),\qquad
E=-\partial_xP,\quad-\partial_{xx}P=\rho-\rho_i,\quad\rho=D_0. \tag{V2}
$$

令 `AdjointAssembly` 的一维块为 $D_{11},D_{12},D_{21}$，并记 $T_a=\mathrm{MultiplierAssembly}(a)$。标准线性子步为

$$
\varepsilon\dot D_k
=-\sqrt{kT_0}D_{21}D_{k-1}
-\sqrt{(k+1)T_0}D_{12}D_{k+1}
-\frac{k}{\tau_0}T_1D_k. \tag{V3}
$$

第一/末模态略去不存在的邻居。`setTransportTerm` 存储 `Dual12=-D12`、`Dual21=D21`、`Dual11=-D11`；`addFreeTransport/addCollision` 实现 V3，`M{1}=epsilon*I`。这里没有把可选 Hermite 滤波项计入连续碰撞算子。

电势及电场由

$$
\begin{bmatrix}
-D_{11}&D_{12}&c\\D_{21}&I&0\\\ell^T&0&0
\end{bmatrix}
\begin{bmatrix}P\\E\\\lambda\end{bmatrix}
=\begin{bmatrix}D_0-\rho_i\\0\\0\end{bmatrix} \tag{V4}
$$

给出；$c$ 是常数函数的系数，$\ell$ 是物理积分权，定义见 [P5](elliptic.md)。由第二行得到 $E=-D_{21}P$，消元的正椭圆矩阵为 $-D_{11}-D_{12}D_{21}$。罚项存在时 $D_{11}$ 的符号不能反转。

## 电场子步的严格三角解

冻结本子步的电场后，剩余系统是

$$
\dot D_0=0,\qquad
\dot D_k=\frac{\sqrt{k}}{\varepsilon\sqrt{T_0}}T_ED_{k-1},\quad k=1,\ldots,N_h-1. \tag{V5}
$$

设 $J_{k,k-1}=\sqrt{k}$，则全部模态堆叠向量 $d=D(:)$ 的生成元是

$$
H=J\otimes\frac{T_E}{\varepsilon\sqrt{T_0}}. \tag{V6}
$$

三个实际子步与完整矩阵方程一一对应：

$$
\begin{aligned}
\text{BE: }&(I-hH)d^+=d,\\
\text{中点: }&(I-\tfrac h2H)d^+=(I+\tfrac h2H)d,\\
\text{SDIRK2: }&(I-\gamma hH)y_1=d,\\
&(I-\gamma hH)y_2=d+(1-\gamma)hHy_1,\\
&d^+=d+h\big((1-\gamma)Hy_1+\gamma Hy_2\big),\quad\gamma=1-1/\sqrt2.
\end{aligned} \tag{V7}
$$

矩阵下三角，因而可以按模态由低到高代入求解，无需构造完整 $H$。每个隐式阶段必须在 `for k=2:nh` 中逐列更新；对所有列同时赋值会读取尚未更新的高阶依赖，漏掉 $h^2,h^3,\ldots$ 耦合。BE 和其他子步必须使用相同的 `epsilon` 因子。

`DghsScheme.step` 组合线性半步、电场整步、线性半步，随后可选滤波、重解 Poisson，并记录完整 `[D,P,E]` 历史。二阶精度需要子步本身满足相应精度，参见 [T6](time-integration.md)。

## 平衡加权的 AP 表示

`ApHermiteState` 使用

$$
f_h=\sqrt{\rho_\infty}\sum_kD_k\psi_k,\quad
d_\infty=(\sqrt{\rho_\infty},0,\ldots)^T,\quad
E_{\mathrm{pert}}=\sqrt{T_0}\sqrt{\rho_\infty}\,Q,\quad
P_{\mathrm{pert}}=T_0\Psi/\sqrt{\rho_\infty}. \tag{V8}
$$

这一变量替换与原模型匹配要求平衡数据满足 $E_\infty=-\phi_\infty'$、$T_0\rho_\infty'=E_\infty\rho_\infty$。`Config.phiInf/rhoInf/sqrtRhoInf/EInf` 需要来自同一平衡；随意组合回调不自动产生平衡保持性质。

对 A/B 两种通量配置分别装配 $D^A,D^B$，令 $T_\infty=T_{E_\infty}$。代码的加权块为

$$
X_{12}^{A/B}=-\sqrt{T_0}D_{12}^{A/B}-T_\infty/(2\sqrt{T_0}),\quad
X_{21}^{A/B}=\sqrt{T_0}D_{21}^{A/B}-T_\infty/(2\sqrt{T_0}),\quad
X_{11}^{A}=\sqrt{T_0}D_{11}^{A},\quad X_{11}^{B}=-\sqrt{T_0}D_{11}^{B}. \tag{V9}
$$

`addFreeTransport` 的上邻块为 $\sqrt{k+1}X_{12}^A$，下邻块为 $-\sqrt{k}X_{21}^A$。碰撞对角块为 $-k\varepsilon T_1/\tau$，质量为 $\varepsilon I$，所以实际模态衰减率是 $k/\tau$。旧算例采用 `tau=epsilon*tau0` 与 V2 的衰减率匹配，比较两个方案时必须同时检查这两个参数。

设 $T_\rho=T_{\rho_\infty}$，AP Poisson 系统为

$$
\begin{bmatrix}
X_{11}^B&-X_{12}^B&c\\X_{21}^B&T_\rho&0\\\ell_w^T&0&0
\end{bmatrix}
\begin{bmatrix}\Psi\\Q\\\lambda\end{bmatrix}
=\begin{bmatrix}D_0-\sqrt{\rho_\infty}\\0\\0\end{bmatrix},\quad
\ell_w^T\Psi=\int_\Omega\frac{\Psi_h}{\sqrt{\rho_\infty}}dx. \tag{V10}
$$

`poissonGauge` 用 `MultiplierAssembly(1/sqrtRhoInf).'*integrationWeights()` 构造 $\ell_w$，而非把两个函数的系数直接点乘。耦合时间推进使用 [T4–T5](time-integration.md)：未知量 `[D(:);Psi;Q;lambda]`；$P$ 只在第一高阶模态对应的 Q 块放入 $T_\rho$，$A$ 在密度块放入 $-I$，$B$ 就是 V10 的约束矩阵，$C=[-\sqrt{\rho_\infty};0;0]$。对应方法名保留原拼写 `addPoissonCoulping`。

分裂的非线性剩余电场子步为仿射三角系统

$$
\dot d=H_Q(d-d_\infty),\qquad
H_Q=J\otimes\frac{T_{\sqrt{\rho_\infty}Q}}{\varepsilon}. \tag{V11}
$$

`stepNonlinearBe/stepNonlinearSdirk2` 使用同一乘法投影；V7 将每次斜率 `H*y` 替换为 `H_Q*(y-dInf)` 即得到仿射公式。耦合路径的 `assembleL/assembleS/assembleF/updatePoissonCoupling` 则构造分块非线性矩阵与源项；部分平衡源项使用节点系数乘法。它们需要独立的时间收敛及平衡残差检验，不能仅凭分裂子步测试宣称两条路径完全等价。

| 公式/组件 | 代码入口 |
| --- | --- |
| V1–V2、重构与物理矩 | `HermiteState.distribution/updateHistory`、`VlasovScheme.setup` |
| V3 | `DghsScheme.setMassTerm/setTransportTerm/addFreeTransport/addCollision` |
| V4 | `DghsScheme.setPoissonLhs/setPoissonRhs/stepPoisson` |
| V5–V7 | `DghsScheme.stepNonlinearBe/stepNonlinearImp/stepNonlinearSdirk2` |
| V8 | `ApHermiteState.distribution/equilibrium/updateHistory` |
| V9 | `ApDghsScheme.setTransportTerm/addFreeTransport/addCollision` |
| V10、约束耦合 | `ApDghsScheme.poissonGauge/setPoissonLhs/addPoissonCoulping/stepLinearized/stepCoupled` |
| V11 | `ApDghsScheme.stepNonlinearBe/stepNonlinearSdirk2` |

`tests/regression/solvers/TestHermiteForce.m` 以五个 Hermite 模态、非单位 epsilon 和非对角空间力矩阵，对每种电场阶段与独立 Kronecker 全局方程求解结果比较；还检查 AP 仿射平衡固定点、加权物理 gauge 以及标准 Vlasov 与正椭圆 Poisson 的符号一致性。它验证上述代数实现，不证明连续模型正性、滤波后能量守恒或任意 epsilon 下 AP 收敛；这些仍需指定极限尺度、网格/步长及误差指标的数值研究。

## 应用问题记录、物理诊断与独立验证

`mdemat.application.plasma.problem(xSpace,vSpace,initial,tFinal,dt,method='ap'|'dghs',config=...)`
记录空间、Hermite 表示、初值及尺度；`initial(X,nh)` 返回 `nh × N` 模态样本。
`solve(problem,observe=...,callback=...)` 调用原 DGHS/APDGHS 的初始化、Poisson、线性/碰撞、电场和滤波步骤；默认不构造图形、分析器、计时器或输出文件。
AP 入口要求原有的节点型空间平衡乘积形式；标准 DGHS 可使用节点型或模态型空间。
速度基在入口处验证为按 $k=0,\ldots,N_h-1$ 排列的归一化 Hermite 基与中心位于零的 Maxwellian 之积。
Hermite 自变量必须为 $v/\sqrt{T_0}$；显式 `config.T0` 必须是正、有限、实数标量，并与速度表示温度一致。
普通 Legendre 空间、未归一化/错缩放 Hermite 或不同温度组合会在问题构造时拒绝，而不是延迟到推进时访问不存在的 `.Rhs`。
问题记录顶层使用封闭 schema，未知、缺失字段或不匹配的模型/版本均拒绝；`config` 内部为兼容 v4 保持开放扩展元数据语义。
速度模型仍在实线上，`nh` 是谱截断；`config.vBBox` 指定绘图与显式物理 `F` 评估的有限速度区间，
不改变 Hermite 演化模型或矩诊断的积分域。`F` 的评估使用该截断区间上的正 Gauss–Legendre
积分权与物理空间权，区别于全实线或 Maxwellian 加权的系数范数，定义见 [采样评估](assessment.md)。

物理诊断在空间积分节点先重构函数，禁止把一般空间系数直接相乘。
令 $s=\sqrt{\rho_\infty}$，标准 DGHS 取 $s=1$，则

$$
\mathcal M=\int sD_0\,dx,\qquad
\mathcal P=\sqrt{T_0}\int sD_1\,dx,\qquad
\mathcal K=\frac{T_0}{2}\int s(D_0+\sqrt2D_2)\,dx. \tag{V12}
$$

`updateHistory` 用物理 Jacobian 和体积积分权评估 V12；截断中缺少的高阶模态取零。
标准 DGHS 的 `potentialEnergy` 为 $\int E^2/2$；AP 历史保留缩放量 $\int E_{pert}^2/(2T_0)$。
两者的 `totalEnergy` 均为各自已存储动能和场诊断之和。因此 AP 的这个字段在非单位温度下不能直接当作未缩放物理总能量；零背景场时，未缩放场能需将其 `potentialEnergy` 乘 $T_0$。
非均匀外部平衡还需计入背景势作用，不能直接套用均匀无碰撞能量守恒。

| 独立应用参考 | 验证入口和范围 |
| --- | --- |
| 非单位温度均匀 Maxwellian $M_{T_0}$ | `TestApplicationPlasma.homogeneousMaxwellianIsStationary`；两种方法、$T_0=0.7,1.5$，直接比较 Gaussian 公式、物理质量和动能 |
| 均匀二阶 Hermite 碰撞模态 | `collisionRelaxationHasIndependentSecondOrderTimeReference`；$D_2=0.1e^{-2t/\tau}$，$\varepsilon=1,\tau_0=\tau=0.3$，时间步减半要求误差至少缩小三倍 |
| 平移 Gaussian 的速度谱截断 | `shiftedGaussianConvergesWithHermiteTruncation`；$N_h=3,5,8$，速度样本区间 $[-4,4]$，直接比较 $M_1(v-0.4)$，每次截断误差至少缩小十倍 |
| 连续 Poisson 场及非线性质量 | `poissonCouplingMatchesContinuousFieldAndPreservesMass`；$\rho=1+0.02\cos x,\rho_i=1$ 的参考 $P=0.02\cos x,E=0.02\sin x$，随后进行真实场耦合推进 |
| 非均匀连续平衡的空间一致性 | `nonuniformContinuousEquilibriumConvergesUnderSpatialRefinement`；$\rho_\infty=e^{0.1\cos x}$，固定时间步、P1 空间单元数加倍，细网格误差界 `2e-3`，要求误差缩小至少 `2.5` 倍；以最新集中证据确认 |
| 三种 Hermite 滤波 | `hermiteFiltersDampOnlyUnresolvedModes`；`hou-li/tc/quasi-tc` 检查高阶模态衰减和零阶质量，未把它们当作物理碰撞或能量守恒 |

平移 Gaussian 的参考系数来自生成函数，而非内部离散矩阵：

$$
M_{T_0}(v-u)=M_{T_0}(v)
\sum_{k=0}^{\infty}\frac{(u/\sqrt{T_0})^k}{k!}
\operatorname{He}_k(v/\sqrt{T_0}),\qquad
D_k=\frac{(u/\sqrt{T_0})^k}{\sqrt{k!}}. \tag{V13}
$$

均匀无碰撞场景无空间或时间演化误差，因而可单独观察 Hermite 截断。
非均匀场景取 $T_0=\varepsilon=1$、$\tau=0.4$，
$\phi_\infty=-a\cos x$、$E_\infty=-a\sin x$、$\rho_i=\rho_\infty+a\cos x$，
从连续关系得到静态 $f_\infty=\rho_\infty M_1$，没有用内部矩阵制造体源。
其时间步固定为 `0.002`、终止时间 `0.02`，空间单元数为 `8/16`；这是非均匀平衡一致性检验，不是小 $\varepsilon$ 渐近定理。
执行证据见 [应用数值验证](../evidence/deep-numerical-verification.json)。这些测试验证了所列有限尺度和周期模型，不证明一般非均匀平衡、耦合/分裂路径完全等价、完整相空间正性或任意 $\varepsilon\to0$ 下的 AP 定理。
