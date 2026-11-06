# 输运与扩散

空间、质量归一化及数组约定见 [公共记号](README.md)。三个求解器共享 `SpatialState.Dofs.U : Nx × 1` 与 `MolScheme` 时间线。

## Euler 迎风 DG

`UwdgScheme` 对应保守输运问题

$$
\partial_tu+\nabla\cdot(a u)=0,\qquad u|_{\Gamma_-}=g_D,\quad
\Gamma_-=\{x\in\partial\Omega:a\cdot n<0\}. \tag{A1}
$$

当前 `advection` 是数值常向量/逐单元常向量。取内侧 $u^-$、外侧 $u^+$，默认 `bias=a` 的面通量为

$$
\widehat{a u}\cdot n=(a\cdot n)\{u\}
+\tfrac12|a\cdot n|(u^- -u^+),\qquad
(\partial_tu_h,\phi)_K=(au_h,\nabla\phi)_K
-\langle\widehat{a u}\cdot n,\phi\rangle_{\partial K}. \tag{A2}
$$

周期面连接相对边；Dirichlet 仅在流入方向引入外侧 $g_D(x,t)$。`UpwindAssembly.assembleMatrix(field,bias)` 可把输运系数和选边方向分开，辐射微观块会使用此能力。实际系数系统是

$$
\dot U=L U+S(t),\quad
L=-\mathrm{UpwindAssembly.matrix}(a),\quad
S=-\mathrm{UpwindAssembly.vector}(a,a,g_D(\cdot,t)). \tag{A3}
$$

| 公式 | 代码入口 |
| --- | --- |
| A1、初值投影 | `AdvectionScheme.setup/configure`、`UwdgScheme.initialize` |
| A2 | `UpwindAssembly.assembleVolumeMatrix/assembleFluxMatrix/assembleFluxVector` |
| A3 | `UwdgScheme.setLinearTerm/setSourceTerm/computeSourceTerm/step` |

## 半拉格朗日 DG

`SldgScheme` 使用特征足点 $X(t_n;t_{n+1},x)$，对固定输运速度满足 $u(t_{n+1},x)=u(t_n,X(t_n;t_{n+1},x))$。`SemiLagrangianAssembly` 对平移后的试验/测试单元交集积分，得到历史转移矩阵 $T_{n+1,n}\in\mathbb R^{N_x\times N_x}$ 和流入边界载荷 $b_{n+1,n}\in\mathbb R^{N_x}$：

$$
M_{n+1}U^{n+1}=T_{n+1,n}U^n+b_{n+1,n}. \tag{A4}
$$

这里的 $M,T$ 是质量回调协议，不应理解成普通的时变 $M(t)\dot U=0$。`computeM(tBegin,tEnd)` 返回 `[boundaryVector, transferMatrix]`；同一时刻调用返回左端质量。周期边界向量为空，非周期边界在时间区间内由特征流入提供。`SldgScheme.step` 只向公共时间推进器传 `M=obj.M`，多步/多阶段方法按自己的阶段和历史时间调用该协议。

转移块在交集上积分 $\overline{\phi_i(x)}\phi_j(X)$；完全出域的足点集合不会产生内部块。对跨越多个边界平面的特征，边界数据取距离目标时间最近的有效交点，即最后进入计算域的时刻；不能按面循环顺序覆盖交点。`assembleVector` 也可在尚未调用 `assembleMatrix` 时独立裁剪。

| 公式 | 代码入口 |
| --- | --- |
| A4、两时刻协议 | `SldgScheme.setMassTerm/computeM/step` |
| 单元交集、周期环绕、流入积分 | `SemiLagrangianAssembly.assembleMatrix/assembleVector` |
| CFL 或指定步长 | `SldgScheme.preStep`、`DynamicTimeline.setTimeStep` |

## LDG 扩散

`diffusion.LdgScheme` 对应

$$
\partial_tu=\nabla\cdot(D\nabla u),\qquad
u|_{\partial\Omega}=g_D(x,t)\text{ 或周期边界}. \tag{D1}
$$

$D$ 是常数/逐单元常数扩散张量；`EllipticAssembly` 的形状按 `d × d × Ne` 解读。使用 [P2–P4](elliptic.md) 的消元块，质量归一化系统为

$$
\dot U=G_DU+b_D(t),\quad
G_D=\sum_r(D_{12}^{(r)}D_{21}^{(r)}+D_{11}^{(r)}),\quad
b_D=\mathrm{EllipticAssembly.vector}(D,g_D). \tag{D2}
$$

`LdgScheme.setLinearTerm` 设置 `L=G_D`，`computeSourceTerm` 直接使用边界载荷，`step` 调用 `TDisc.step(L=L,S=S)`。不要套用 Poisson 的负号；扩散的 $G_D$ 在正确稳定配置下应是耗散生成元。

格式的空间收敛阶、离散质量守恒和能量耗散依赖网格、通量及时间推进。周期离散常数模态/总质量、扩散生成元谱与制造解收敛是独立验收指标；显式推进还需要相应稳定步长。这里没有宣称任意通量选项和任意步长均稳定或保持正性。

## 独立计算入口

`mdemat.solver.hyperbolic.advectionProblem/solveAdvection` 与 `mdemat.solver.parabolic.diffusionProblem/solveDiffusion` 将物理数据、显式空间和调用方选择的 ODE 工作区分开。`solveAdvection(...,method="uwdg"/"sldg")` 分别使用 A3/A4 的原有公共装配；`solveDiffusion` 使用 D2。边界回调接受 `(x,t)`，周期模型不接收额外 Dirichlet 数据。它们返回系数纯记录和质量/时间报告，默认没有观测或文件行为。所有单步/多步公共 ODE 类均作为显式工作区传入；数值元数据不隐含参考解。

`TestEvolutionInterfaces` 验证周期迎风三层空间收敛/质量、SLDG 常数传递、扩散正弦精确衰减/质量、双向单单元长步长的纯流入，以及不等边界 $[0,2]\times[-1,3]$、混合符号速度 $(2,-4)$ 下二维特征的独立进入时刻 $t_{\rm entry}=1-\min(x/2,(3-y)/4)$。`TestOdeScientific` 进一步验证 SLDG 指数流入数据的双时刻质量回调，初值和边界均来自独立的行波精确解。
