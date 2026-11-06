# 公共时间积分与代数约束

接口入口为 `mdemat.common.ode.OdeIntegrator`。Euler 框架下常质量问题为

$$
M\dot U=L(U,t)U+F(U,t)+S(U,t),\qquad U\in\mathbb R^N. \tag{T1}
$$

`M=[]` 表示单位矩阵，其他空项表示零；`L/F/S` 可为常值或按实际签名接收 `U`、`(U,t)`、`(U,tBegin,tEnd)` 的回调。`evalOp` 只分派参数，不额外计算空间离散或修改物理量。

单时刻 `M(t)` 对应守恒积方程

$$
\frac{d}{dt}\bigl(M(t)U(t)\bigr)=R(U,t). \tag{T1a}
$$

每个阶段先按自己的时刻求 $M(t_i)U_i$，再以该阶段状态和时刻评估右端；$R=0$ 保持的是 $M(t)U(t)$。若模型原本为 $M(t)\dot U=R$，调用方必须把 $\dot M(t)U$ 加入 T1a 的右端，不能直接复用而忽略乘积求导。质量接口另支持半拉格朗日两时刻 `[b,T]=M(tBegin,tEnd)`：将历史向量变成 `T*Uold+b`，同一时刻调用提供左端质量。这是历史转移协议，由 `evalMLhs/evalMRhs` 分派。

常质量 DIRK 的数学定位是

$$
M U_i=M U_n+h\sum_{j<i}a_{ij}R_j+h a_{ii}R(U_i,t_i),\quad
M U_{n+1}=M U_n+h\sum_i b_i R_i. \tag{T2}
$$

`DirkIntegrator.A/b/c` 存储表，`stage/setLhs/setRhs` 构造阶段系统。BE 是单阶段全隐式；隐式中点、SDIRK2/3/4、ESDIRK3 具有各自表。BE/DIRK 的 `F` 与依赖当前状态的 `L/S/P` 在当前阶段隐式求解，随后 DIRK 的阶段斜率包含 `LU+F+S+P*lambda` 全部项。`OdeIntegrator.implicitResidual/solveImplicitStage` 复用 JFNK，并独立复算代数残差；`SolveReports` 记录系数空间欧氏残差及阈值，`NonlinearStatistics` 保留每次 Newton 的残差历史、停止状态，以及内部 Krylov 的 flag/相对残差/迭代数/残差向量，超过阈值抛出 `mdemat:ode:NonlinearStageFailed`，不会返回未收敛的旧阶段。`NonlinearTolerance/NonlinearMaxIterations` 是显式工作区设置。不同物理缩放需要调用方选择适当容差，代数残差不是统一物理误差界。

IMEX 将 $R_I=LU$、$R_E=F+S$ 分开：

$$
M U_i=M U_n+h\sum_{j<i}a^I_{ij}R_{I,j}
+h\sum_{j<i}a^E_{ij}R_{E,j}+h a^I_{ii}L_iU_i. \tag{T3}
$$

`ImexrkIntegrator.AI/AE/bI/bE/c` 与 `Ars111/222/443Integrator` 对应 T3；其 `F/S` 显式处理，`L` 是隐式线性算子。ARS 的显式初始阶段计入 `c/AE/bE`，隐式表 `AI/bI` 少一个阶段；内部第 `i` 个隐式阶段对应 `c(i+1)`，不能用 `c(i)` 评估非自治系数。显式 RK 家族包括 FE、中点/Heun、SSPRK2/3、经典 RK4，其全部右端项显式处理。

BDF 使用历史质量/转移后的数值量：

$$
M_{n+1}U_{n+1}=-\sum_{j=1}^{k}\alpha_j\mathcal T_{n+1,n+1-j}(U_{n+1-j})
+h\beta\,[L(U_{n+1})U_{n+1}+F(U_{n+1},t_{n+1})+S(U_{n+1},t_{n+1})]. \tag{T7}
$$

单时刻质量下 $\mathcal T$ 是对应历史时刻的 `M(tHistory)*UHistory`；双时刻回调下是 `transfer*UHistory+boundary`。`BdfIntegrator.setRhs/stage` 与 `Bdf2/3Integrator.setCoefficients` 实现 T7。BDF2 使用 BE 启动；BDF3 默认 `Startup="sdirk3"` 用匹配三阶的前两步初值，随后采用三步公式。`Startup="legacy"` 保留旧 BE/BDF2 启动，该固定步长策略通常只给全局二阶，不能因后续三步公式三阶就宣称整体三阶。公开 `F` 为当前状态的隐式项，连同状态相关 `L/S/P` 做完整非线性残差求解；BDF3 的两种启动也传递同一个 `F`。`TestBdfNonlinear` 以独立非线性衰减、变质量制造解、增广约束及显式重建后的轨迹验证二/三阶和历史。质量回调每次按真实时间评估，固定步长不代表质量矩阵不变。

带约束的等离子体阶段使用微分–代数系统

$$
M\dot U=LU+P\Lambda+F+S,\qquad A U+B\Lambda=C, \tag{T4}
$$

对应阶段矩阵

$$
\begin{bmatrix}M-h a_{ii}L&-h a_{ii}P\\A&B\end{bmatrix}
\begin{bmatrix}U_i\\\Lambda_i\end{bmatrix}
=\begin{bmatrix}R_i\\C_i\end{bmatrix}. \tag{T5}
$$

其中 $M,L\in\mathbb R^{N\times N}$，$P\in\mathbb R^{N\times N_c}$，$A\in\mathbb R^{N_c\times N}$，$B\in\mathbb R^{N_c\times N_c}$，$C\in\mathbb R^{N_c}$。`BeIntegrator`、`DirkIntegrator`、`ImexrkIntegrator` 的增广路径负责 T5，`AugmentSolver` 用于末步约束恢复。`LinearSolver` 只解装配后的代数方程，物理约束不属于其职责。

省略 `M` 时单位质量只补在 $N$ 个微分变量上，不能把乘子也当作带质量的状态。显式 RK 由当前阶段的 $AU+B\Lambda=C$ 恢复乘子，要求该恢复子问题可解；这里的回归覆盖可逆 $B$ 的 index-1 约束，不宣称一般高指数 DAE 的阶数或稳定性。

`LieTrotterIntegrator` 顺序组合子流；`StrangIntegrator` 为两算子提供

$$
\Phi_h^{\mathrm{Strang}}=\Phi_{h/2}^{(1)}\circ\Phi_h^{(2)}\circ\Phi_{h/2}^{(1)}. \tag{T6}
$$

二阶结论要求足够光滑的子问题及相应精度的子步；若中间使用一阶非对称 BE 近似，不能只因外层叫 Strang 就宣称整体二阶。`SeparableOdeIntegrator` 管理因子及共享时间线，各应用的 `step` 负责实际组合。

历史 `U0` 是 `NStages × NSteps` 的 cell：行表示阶段，列表示时间历史。`update` 推进历史，`add` 写入当前系数，`advance` 推进时间线。`MolScheme.precompute` 只在 `[Timeline.Now, Timeline.Next]` 的非零区间准备观测插值，取 `U0{1,historyIndex}`，并按所有字段的完整打包偏移提取启用观测的字段。`TestSolverContracts.observationHistoryUsesTimeIntervalAndPackedOffsets` 覆盖这些契约。

标准验收应分别检验已知标量/矩阵解的时间收敛、刚性极限、约束残差和变步长历史；应用回归不能替代这些公共积分器验收。Hermite 电场子步另用独立完整三角线性系统校验，见 [plasma.md](plasma.md)。

独立执行入口 `mdemat.common.ode.integrateLinear` 接受初始系数、工作区及显式 `mass/linear/source`；可选存储轨迹与 `maxSteps` 中断。完成步边界的 checkpoint 由 `snapshotIntegrator/restoreIntegrator` 按固定方法白名单保存/重建，包含阶段×历史 `U0`、质量数据、BDF/RK 系数、全部时间与步长队列、非线性设置及递归 split 因子；重建求解器缓存，物理回调由调用方重新提供。`Phase="completed-step"` 不接受阶段内部或半个 split 的恢复。结果与报告仅持有数值记录。

`TestOdeScientific` 对所有具体时间方法分别做非自治强制收敛，另验证非线性隐式方程、IMEX 加性非线性强制、变化质量守恒积、非交换算子的 Lie/Strang 阶数、增广约束和 SLDG 流入；`TestEvolutionInterfaces/TestIntegratorPersistence` 比较连续执行与中断恢复的完整轨迹、BDF历史和递归 split 状态。这些验收直接比较独立精确解或矩阵指数，不以实现自己的输出作为参考。

`TestOdeScientific.allUnsplitMethodsRespectSemiLagrangianStageAndHistoryTransfers` 还对全部十七种非 split 方法验证双时刻质量回调：仿射行波 $u(t,x)=1+2(x-at)$ 在正/负速度下同时穿越单元和 Dirichlet 边界，所有阶段和多步历史均须组合正确的转移与流入载荷。这里没有把 split 子流包装成通用半拉格朗日质量模型。

## 指数积分的整函数

公共 `expseq(Z,n)` 返回 $\phi_0$ 到 $\phi_n$：

$$
\phi_k(Z)=\sum_{j=0}^{\infty}\frac{Z^j}{(j+k)!},\qquad
\phi_k(Z)=\frac1{(k-1)!}\int_0^1e^{(1-s)Z}s^{k-1}\,ds\quad(k\ge1). \tag{T8}
$$

零矩阵有 $\phi_k(0)=I/k!$，奇异或幂零矩阵也有唯一整函数值；递推式 $Z\phi_{k+1}=\phi_k-I/k!$ 在核空间上不唯一，不能用奇异线性求解替代定义。实现对小范数自适应求和，其余情况计算上三角增广矩阵指数的首块行：对角首块为 $Z$，相邻上对角块为 $I$，其余块为零。该首行正是 $[e^Z,\phi_1(Z),\ldots,\phi_n(Z)]$。输出为稠密矩阵，增广维数为 $(n+1)\dim Z$；大规模问题应另选矩阵函数作用算法，当前接口没有宣称稀疏存储复杂度。旧的求逆递推私有函数已由这个无逆算法替代。

`TestPhiFunctions` 以零/大范数幂零矩阵的有限解析级数、复标量积分、奇异及非正规矩阵的独立变常数积分验证 T8；同时检查 single/sparse 输入、标量负无穷刚性极限，以及非整数/非有限阶数和未定义矩阵值的拒绝。

矩阵无关差分 Jacobian 在小残差区间可能遇到浮点舍入，GMRES 可报告停滞或 `tooSmallTolerance` 警告；这类内部诊断保留在上述统计中，不被屏蔽，也不会自动放宽验收阈值。阶段最终是否接受仍由重新计算的完整非线性方程残差决定，内部线性迭代 flag 不能单独证明阶段成功或失败。

`krylovStats.rawIterations` 保留 MATLAB 原始计数。GMRES 的 `linearSteps` 为返回迭代点的 $(outer-1)restart+inner$（零步情况为零），其它方法保留原标量迭代数；这与完整残差轨迹覆盖的总工作量应分别解读，失败时返回的可能是较早的最佳迭代点。计数约定见 [MathWorks GMRES 文档](https://www.mathworks.com/help/matlab/ref/gmres.html)，`TestNewtonDiagnostics` 以重启 GMRES 的独立成功求解/残差轨迹验证。
