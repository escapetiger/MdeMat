# 辐射输运的宏微 DG

`MmDgScheme` 表达的模型族为

$$
\varepsilon^{s_t}\partial_tf+v\cdot\nabla_x f
=\varepsilon^{s_s}\sigma_s(E\rho-f)
-\varepsilon^{s_a}\sigma_af+\varepsilon^{s_a}q,
\quad\rho=\int_Vf\,dv,\quad\int_VE\,dv=1. \tag{R1}
$$

`timeScale/scatteringScale/absorptionScale` 分别是 $s_t,s_s,s_a$；当前源项使用与吸收相同的尺度。`scattering/absorption=[]` 省略对应项，`source=[]` 省略源项。纯离散速度的散射使用 `Config.E` 与速度积分权，因此必须选择与离散速度测度一致的归一化。空间及速度降维后的测度不能混用。

## 表示与投影

设宏基数 $m$、微速度自由度 $n$，缩放指数为 $p_i$。实际表示为

$$
f_h(x,v)=\sum_{i=1}^m\varepsilon^{p_i}U_i(x)\phi_i(v)
+\sum_{j=1}^n\varepsilon^{p_{m+j}}G_j(x)\psi_j(v). \tag{R2}
$$

`Dofs.U` 为 `Nx × m`，`Dofs.G` 为 `Nx × n`，`SumSpace.LhsWeight/RhsWeight` 存储 $\varepsilon^{p_i}$。`setDecomposition('')` 取零指数，`partial` 取 `[0,1,...,1]`，`full` 采用 `Levels` 给出的模态层级。`distribution` 与 `kineticReconstruct` 必须通过 `evalLhs/evalRhs/eval` 应用相同权重；`macroFit/microFit` 的逆映射必须除回这些权重。

混合表示的宏系数由宏积分投影取得，微残差为原函数在微节点的值减去宏重构；这是数值投影的定义。不同速度积分规则不自动保证离散微残差的每个宏矩均为机器零。

设 `EU=VLhsBasisValues.'`，`EG=macroBasis.eval(vMicro).'`；$W_U=\operatorname{diag}(w_U)$，$W_G=\operatorname{diag}(w_G/E(v_G))$，$V_U^{(r)},V_G^{(r)}$ 是第 $r$ 个速度分量的对角矩阵。令 $S$ 只保留最高宏层级。`Pattern.freeTransport(:,:,r)` 的混合分块精确定义为

$$
\begin{aligned}
T_{UU}^{(r)}&=E_U^TW_UV_U^{(r)}E_U, &
T_{UG}^{(r)}&=S E_G^TW_GV_G^{(r)},\\
T_{GU}^{(r)}&=V_G^{(r)}E_G-E_GT_{UU}^{(r)}, &
T_{GG}^{(r)}&=V_G^{(r)}-E_GS E_G^TW_GV_G^{(r)}.
\end{aligned} \tag{R3}
$$

纯宏模式只使用 $T_{UU}$，纯微模式使用 $T_{GG}^{(r)}=V_G^{(r)}$。`setFreeTransport` 将绝对值小于 `1e-13` 的项设零；该阈值会改变极小耦合项。

## 行缩放与离散块

按 `[U(:);G(:)]` 排列 $W\in\mathbb R^{N_x(m+n)}$。缩放前第 $i$ 行的幂次为

$$
e_{i,t}=s_t+p_i,\quad e_{i,j}=p_j\ (T_{ij}\ne0),\quad
e_{i,s}=s_s+p_i,\quad e_{i,a}=s_a+p_i. \tag{R4}
$$

零耦合和关闭项使用 $+\infty$；有宏密度时其散射幂次也设为 $+\infty$。每行减去有限最小值 $r_i=\min e_{i,*}$，得到 `Pattern.scaling = epsilon.^(e-r)`，再将无效项置零。这是方程行预缩放，不是额外物理项。

$$
M\dot W=L W+S(t),\quad
M_{ii}=\varepsilon^{e_{i,t}-r_i}I_{N_x}. \tag{R5}
$$

`addFreeTransport` 对宏–宏上三角、宏–微使用 $D_{12}^{(r)}$，宏–宏下三角和微–宏使用 $D_{21}^{(r)}$，系数均为 $-\varepsilon^{p_j-r_i}T_{ij}^{(r)}$。微–微块通过 `UpwindAssembly.assembleMatrix(a,v_j)` 构造。非周期且启用首个罚项时，还在宏对角块增加由 `D11`、模态奇偶性和 $T$ 决定的跳跃修正；`addFreeTransport` 中四个具名分区就是对应定位点，不能用一个统一的中心通量公式覆盖它们。

令 $C_s=\mathrm{MultiplierAssembly}(\sigma_s)$、$C_a=\mathrm{MultiplierAssembly}(\sigma_a)$。有宏密度时散射对第一宏密度模态为零，对其余模态加入 $-\varepsilon^{e_{i,s}-r_i}C_s$。纯微情形另加恢复密度的秩一速度耦合：

$$
(L_s)_{ij}=\varepsilon^{e_{i,s}-r_i}C_s
\big(Ew_j-\delta_{ij}\big),\qquad
(L_a)_{ii}=-\varepsilon^{e_{i,a}-r_i}C_a. \tag{R6}
$$

`addSource` 在当前时间调用 `macroFit/microFit`，再乘吸收列对应的行缩放权。`addBc` 以 `macroTraceEvaluate/microTraceEvaluate` 构造流入数据：宏项使用伴随通量，微项使用逐速度流入方向。周期时不添加这些边界载荷。`step` 向公共积分器传 `M/L/S`，解包后依次执行可选滤波、可选正性限制和物理重构。

| 公式/组件 | 代码入口 |
| --- | --- |
| R1、参数及边界类别 | `RadiationScheme.setup/configure`、`MmDgScheme.setup` |
| R2、拟合、重构 | `MacroMicroState.setLevels/setDecomposition/setScaling`、`macroImpl/microImpl/distribution/kineticReconstruct` |
| R3 | `MmDgScheme.setFreeTransport` |
| R4 | `MmDgScheme.setScaling` |
| R5 | `MmDgScheme.setMassTerm/setLinearTerm/step` |
| R6 | `MmDgScheme.addCollision` |
| 体源、流入 | `MmDgScheme.addSource/addBc`、`MacroMicroState.*TraceEvaluate` |
| 诊断 | `MacroMicroState.density/updateHistory` |

`History.mass` 使用空间物理积分及 `density`。`History.freeEnergy` 当前记录未缩放 $U,G$ 的空间系数平方和的一半；在缩放或非正交微基下，它不是通用的物理 $\int f^2/(2E)\,dx\,dv$，不能据此断言物理熵耗散。限制器对负单元平均值进行截断时会改变质量；滤波与限制器均不提供普适的守恒/正性保证。

`TestSolverContracts` 覆盖状态复制不共享可变权重、纯微缩放重构、拟合/重构互逆和密度不依赖过期重构缓存。下面的应用验证补充了指定两速度尺度下的扩散极限参考；缩放代数及这一组有限参数试验均不构成所有速度表示、边界和参数下的 AP 定理。

## 纯计算入口与验证范围

`mdemat.application.radiation.problem(xSpace,vSpace,initial,tFinal,dt,config=...)`
返回显式问题记录；`vSpace` 为宏/微 `SumSpace`，可使用纯宏、纯微或组合表示。
`initial(X,V,t)` 对成对的物理位置与速度返回 `1 × N` 样本，有限标量表示各向同性常数。
初值、源项及流入回调中的 `V` 都是 `VLhsNodes/VRhsNodes` 转换后的物理 Cartesian 速度。
`state.distribution(xRef,v)` 与 `state.kineticReconstruct(v=...)` 则直接评估速度基，
`v` 必须使用该基的 chart 坐标，可直接取对应 `Element.Volume.Nodes`：圆的 chart 是角度
`theta`；默认球面 `cossph` chart 是 `[-cos(vartheta);phi]`，不是三行 Cartesian 向量。
一维 symmetry 区间及两离散速度的 chart 与物理速度分量相同。
纯问题入口要求 `L2SphereElement` 角表示和原有 product-form 宏基；混合表示共享物理域与
chart 维数。二维 symmetry 的宏 `modal(2,...)` 与微 `nodal(3,...)` 实际同属三维球面，
因此按实际几何兼容，而不按 factory 输入维数或 reduction 标签机械拒绝。
通用 `MacroMicroState` 的重构能力保留，不在其内部施加应用入口的额外限制。
`solve(problem,observe=...,callback=...)` 创建独立的状态和可变缩放权，复用上述全部数值组件。
默认不构造计时器、可视化器、分析器，不需要精确解、参考文件或输出目录。
显式 observer 的 `F` 分量对应选定角 chart 节点的物理样值，在空间方向仍存有限元系数；
`U/G` 是对应宏/微系数或矩序列。旧 `run/converge` 的适配器把物理回调转换成这些对应量，
统一向通用采样观察器提供 `fields × physical-points`，不把 raw RBF 系数当分布样值。
`internal.execution.compute` 依次调用 `initialize`、`preStep/step/postStep`；
报告和输出由显式 `run/converge` 或调用者回调承担。
问题记录的顶层字段为封闭 schema：未知、缺失字段或版本/模型不匹配均拒绝。
`config` 内部保留 v4 既有的开放扩展元数据语义；它不是顶层拼写错误的容错入口。

| 验证及独立参考 | 代码入口和预声明界 |
| --- | --- |
| 均匀 $f=1$，纯宏/纯微/组合及圆、球角表示 | `TestApplicationRadiation.uniformEquilibriumUsesEveryVelocityRepresentation`；重构与物理质量误差界分别为 `2e-10`、`3e-10` |
| 均匀散射、吸收、常体源 | `homogeneousScatteringAbsorptionAndSourceMatchExactSolution`；直接比较下式，误差界 `8e-7` |
| 两速度扩散尺度、固定空间网格与时间步 | `diffusionScalingHasAnIndependentTelegraphReference`；$\varepsilon=0.1,0.01$，$\Delta t=0.005$，密度与 R7 误差界 `3e-3` |
| 纯微周期自由输运及时间依赖流入 | `pureMicroFreeStreamingAndTimeDependentInflow`；连续参考 $1+0.1\cos(x-vt)$、$1+0.2(x-vt)$，P1 空间界 `2e-3`，可精确表示的线性流入界 `2e-10` |
| 真正二维、三维空间的均匀平衡 | `uniformEquilibriumInTwoAndThreeSpatialDimensions`；每轴两个单元、Q1 空间、圆/S2 组合角表示，重构与质量界 `3e-9` |
| 滤波及密度限制 | `angularFiltersDampHigherModesAndPreserveMass`、`densityLimiterUsesPhysicalSamplesAndPreservesAdmissibleMeans`；检查显式模态衰减、物理质量及体积节点密度 |

执行记录见 [应用数值验证](../evidence/deep-numerical-verification.json)；新增空间/流入场景以最新集中测试证据为准。
均匀两速度参考为
$f_\pm(t)=q/\sigma_a+(1-q/\sigma_a)e^{-\sigma_at}
\pm0.2e^{-(\sigma_s+\sigma_a)t}$，测试使用 $\sigma_s=0.7,\sigma_a=0.2,q=0.3$。

对 $v=\pm1$、$E=1/2$、$\sigma_s=1$、$\sigma_a=q=0$，令
$r=(f_++f_-)/2$、$j=(f_+-f_-)/2$，则

$$
\begin{aligned}
r_t+j_x/\varepsilon&=0, &j_t+r_x/\varepsilon&=-j/\varepsilon^2,\\
\lambda&=-\frac{2k^2}{1+\sqrt{1-4\varepsilon^2k^2}}, &r&=1+a e^{\lambda t}\cos(kx),\\
j&=-\varepsilon\lambda a e^{\lambda t}\sin(kx)/k, &f_\pm&=r\pm j.
\end{aligned} \tag{R7}
$$

该实慢模态参考要求 $2\varepsilon k\le1$；测试取 $k=1,a=0.1,t_f=0.05$。
物理密度为 $\rho=2r$，热极限为 $2(1+a e^{-k^2t}\cos(kx))$。
这是指定离散速度、周期边界、尺度和误差界下的验证，不替代连续角空间、一般非均匀散射或其它边界的极限研究。

限制器只针对第一宏密度模态在体积积分节点的值。节点系数顺序与积分节点顺序可以不同，代码先用 $B^TC$ 评估，再用设计矩阵逆映射。
Zhang–Shu 分支对非负单元平均值采用
$\theta=\min(1,\bar u/(\bar u-u_{\min}))$，保持单元平均值并限制这些采样点；负平均值单元被置零，会改变质量。
这不保证整个单元或完整相空间分布非负。`cutoff` 的模态分支只截断平均系数，节点分支截断采样值，均不能宣称普适的守恒限制器。

`rot` 滤波每接受一步使用 $\sigma_k=[1+(s_N^{-1}-1)\eta_k^4]^{-1}$；固定强度下改变步数会改变总衰减。
`exp` 使用 $\sigma_k=\exp(\log(s_N)\Delta t\eta_k^p)$，强度按单位时间解释。
这里的滤波验证不作为时间精度或物理熵耗散的证明。


## 可控的圆周近似研究

`toolbox/examples/radiation/ex00_func_approx.m` 保留四个目标函数、原 modal/nodal/hybrid
扫描和固定 profile 参数；`ex05_generalization.m` 保留七个输运角剖面、五组表示和全部原
DOF 预算。两个入口返回 chart 网格、目标样本、重构值和误差，可显式指定
`ReferencePoints`、`MakePlots=false`、`OutputDirectory=""` 来禁止图形与文件。
前者另支持 `ModalOrders/NodalDofs/HybridConfigurations`，后者支持 `Dofs`。
非空输出目录显式保存结果；覆盖已有结果必须指定 `Overwrite=true`。

拟合与评估都使用角 chart。宏投影满足
$M c=\overline B W f$，节点残差通过原 design solve 求系数；不把 RBF 系数当成节点值，
也不改变原角节点顺序、`pi/6` 偏移或 RBF 形状参数。不同采样网格使用同一个目标函数。
`ex05` 的窄束采用圆周截断 Gaussian 的解析归一化：

$$
g_\sigma(\theta)=
\frac{\exp[-\operatorname{wrap}(\theta-\theta_0)^2/(2\sigma^2)]}
{\sigma\sqrt{2\pi}\,\operatorname{erf}[\pi/(\sqrt2\sigma)]},
\qquad\operatorname{wrap}(\alpha)\in[-\pi,\pi),\quad\theta_0=\pi/6. \tag{R8}
$$

归一化不随 reference/modal/nodal 的采样点数改变。原有
$\sigma=0.20,0.10,0.05$ 均保留；粗角积分仍可能欠分辨窄束，不能把离散归一化冒充连续质量正确。
`TestAngularStudies` 调用这两个实际入口：独立 Fourier 多项式投影及原节点插值界为
`2e-12`；`ex05` 平滑输运剖面 DOF `9 -> 17` 要求相对 $L^2$ 误差低于 `2e-3` 且改善至少
20 倍，其连续平均值独立取
$1+\cos(\pi x_0)\cos(\pi y_0)J_0(\pi\sqrt2)$。
