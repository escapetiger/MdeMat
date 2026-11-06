# 工具层优化设计与实施计划

这是已落实到源码的工具层重构说明及验收合同，替代首轮“主要保留类、增加适配”的方案。
范围覆盖当前 main 的工具文件、局部/private 函数、公开属性/选项及实际调用链；只关注 v4，
不重新引入旧版本或 legacy。

**当前总体验收：PASS。** 新接口及调用方已落地，原生完整研究工作流、生产覆盖、
[测试](evidence/tools-required-tests.json)、[符号测试](evidence/tools-symbolic-tests.json)、
[分析告警复核](evidence/tools-analyzer-review.json)、[打包](evidence/tools-package-report.json)、
[全新安装](evidence/tools-installed-package.json) 和 [项目核验](evidence/tools-native-project.json)
均已完成。任务及指标以 [tasks.csv](tasks.csv)、[status.csv](status.csv) 为准；历史探针和旧产物
保留其各自的源码身份，不计作当前验收。

[逐文件方案](tools-file-review.csv) 逐项给出保留、替换、合并、迁移或移除，以及目标路径、
兼容影响和验收。[首轮证据](evidence/tools-layer-audit.json) 保留原始观察；
[深化证据](evidence/tools-design-review.json) 记录本轮新增发现、原生探针和源码身份核对。
数量统一维护在 [status.csv](status.csv)，实施状态与依赖统一在 [tasks.csv](tasks.csv)。

本轮目录归纳统一使用 `mdemat.tools.analyze`，其中包含 `Timer`、`weightedError`、
`convergenceRates`、`summarizeConvergence`、`benchmark` 和 `memorySnapshot`。
原 `analysis` 包已更名，`timing.Timer` 已迁入；空的 `except` 和 `timing` 目录已移除。
现有诊断错误 ID 和数值文件格式保持稳定，历史证据保留当时的原始名称。

## 已实现的职责划分

工具层应保留通用计算、值结构和图形/计时资源管理。数值方法、State、空间采样、实验循环
和重启语义由对应拥有者负责。取消没有独立不变量或资源身份的类，直接使用 MATLAB 值和函数。

| 原设计 | 当前结构 | 理由 |
|---|---|---|
| GraphPlotter / PolytopalMeshPlotter / SeparableGraphPlotter | data.mesh 相应类的外置 plot 方法 | 绘图依赖的数据结构已在 data 定义；不再在工具层重复建一个类型体系 |
| Dataset / Database | Type+Data 的普通 struct；struct-of-struct 集合；纯 selectComponent 函数 | 不需要独立对象身份，内置字段操作已能表达存取；选择与网格重排分离 |
| Strategy 及其四个子类 | plotStyle 值记录 + 无持久状态的 plotField 函数 | 配置和算法不需要 handle，也无需策略继承树 |
| FigureManager / TiledPlotter / Visualizer | Visualizer 为唯一图形资源 owner；布局/样式是值，算法是小函数 | 统一图窗/轴/图元生命周期，取消重复登记和状态同步 |
| ConfigParser | 不继承 handle 的 Config 值类 + 纯 readConfig 函数 | 一个私有 key→(Value,Validator) 结构维护成对不变量，文件解析独立 |
| Metric / Analyzer | 纯数值归约、收敛汇总函数 + 上层观察记录和采样函数 | 取消无执行合同的基类、FEM factory、隐藏实验控制及绘图状态 |
| timing.Record / Timer | Timer 为唯一计时 handle；记录为普通 struct 快照 | 同一时刻只有一个活动 token，无需每条记录独立可变对象 |
| except.assert / verify / backtrack | 原生 error/warning/assert/MException + 显式完整 ID | 原生设施已覆盖能力；取消目录/堆栈猜身份和二次格式化 |

这不是“类越少越好”的一般规则。Config 有成对存储和受控更新的不变量，因此保留小值类；
Timer、Visualizer 拥有跨调用资源，因此保留 handle。简单样本、样式、布局和分析结果则是值。
不添加 ConfigEntry、Schema、Resolver、Adapter、MetricProvider、TimerScope 或 renderer 类树。

重构前 Scheme 已继承 handle 并组合 ConfigParser，并没有继承 ConfigParser。现在它组合
Config 值，向运行代码提供只读快照；此次替换的是组合成员及配置生命周期。

## 原始问题与修复依据

首轮 F01–F17 的容器漏检、hash hex 错误、Timer 重名、复数损坏接受、坐标/分量/时间错配、
图形异常恢复、极端误差溢出和职责耦合仍适用，证据见首轮 JSON 与逐文件表。
原始新增发现如下。“复现”记录重构前的独立反例；修复是否通过验收由当前回归与任务证据说明，
不能仅凭源码变更或该历史表判断。

| 编号 | 证据与影响 | 处理拥有者 |
|---|---|---|
| F18 | **复现**：Scheme.loadFromStruct 第二键失败后，第一键已改变；文件→覆盖→configure 也没有整体事务 | Scheme 候选构造/校验后提交；Config.update 原子返回新值 |
| F19 | **复现**：配置改成 fe、tFinal=2，实际 TDisc 仍为 BeIntegrator、Final=1 | 数值配置与运行对象一起冻结；禁止实例运行后静默更改模型参数 |
| F20 | **复现**：OPTION paths 与 MATLABPATH 共享字段，触发原始 UNIQUE 错误 | readConfig 的资源描述放第二输出，OPTION 值保持独立 |
| F21 | **复现**：选 [2,4] 后 sequence 只存一列，报告仍称两组件，标成 U2，draw 越界 | 分量选择、聚合方式与结果身份明确分开，取消 Analyzer 的隐藏状态 |
| F22 | **复现绝对 Metric 路径；相对路径静态同构**：TYPO 归约未拒绝，返回原始样本而非范数 | tools 的严格归约枚举；raw samples 是单独结果，不是未知 norm fallback |
| F23 | **复现**：合法无边 SeparableGraph 的降采样把 0×2 改为 []，即使 ShowEdges=false 也报错 | data 绘图保持空连通性 shape，避免为显示重建数据对象 |
| F24 | **复现**：物理字段 x 覆盖 Dataset.Data.x 坐标；v 同样存在命名冲突 | 每字段 view 内部存几何，不与物理字段同级混装 |
| F25 | **复现**：改变纯观察 components 改变 checkpoint 数值身份；合法 nested exact 被拒绝 | execution/模型明确描述数值身份和回调路径；通用 IO 不推测角色 |
| F26 | **复现**：严格递减的相邻大尺度 log 值舍入相同，收敛率 Inf 却标 Defined | 对接近正数用 log1p 的稳定 log 比；病态程度与是否有定义分开 |
| F27 | **复现**：Vlasov 接受 exact 函数，configure 却只注册 struct，参考解静默丢失 | plasma 自己定义物理参考到被测字段的适配，工具不猜 D/F 的含义 |
| F28 | **命名转换复现、未跑完整保存工作流**：0.11 与 0.14 均成为 t0p1 | state 用 accepted-step 序号；figure 用观察帧序号；精确时间写元数据 |
| F29 | **设计/数学合同**：旧 sequence 的 L1/max 与全分量 L1/LInf 不同；值 struct 内嵌 handle 仍共享 | 明确范数及浅值语义，不能机械替类后改变数学或声称递归深拷贝 |

其他静态风险必须随对应改动处理：捕获整个 Scheme 的 validator、构造/initialize/resume 所需键
不同、同列数却不同速度基的补零、采样返回形状中的 squeeze、坏 IO 描述的预算及大 MAT 提交。
此前的全量通过证据只对应旧快照；这些边界已纳入新回归并通过本轮集成验收，适用范围以
当前测试和正式证据为准。

## 数据结构的绘图接口

当前目录和完整类名：

```text
toolbox/+mdemat/+data/+mesh/
  @Graph/Graph.m
  @Graph/plot.m
  @SeparableGraph/SeparableGraph.m
  @SeparableGraph/plot.m
  @PolytopalMesh/PolytopalMesh.m
  @PolytopalMesh/plot.m
  Triangulation.m
```

```matlab
h = plot(graph, axes=ax, showVertices=true);
h = plot(separableGraph, axes=ax, maxVertices=1000);
h = plot(triangulation, axes=ax, showBoundary=true);
```

Triangulation 继承 PolytopalMesh.plot，不复制实现；SeparableGraph 的普通 Graph fallback 由
正常方法分派替代。返回新建图元句柄，绘图不修改数据对象，也不向其添加 Figure/Axes/Style
或缓存属性。Graph 现有顶点为 vertices×dimension，不为显示全仓改存储方向。
新 axes 选项显式指定对象；缺省 gca 仅保留交互便利，副作用在 help 中说明。

外置方法要求 classdef 本体也移动进 @Class；旧 Graph.m 等不能与新路径并存。
构造器和 getter/setter 仍在 classdef。本轮临时命名空间原型已验证外置方法、继承和显式轴；
项目类迁移、产物重建及独立进程安装验证均已完成。[MathWorks R2025b 外置方法规则](https://www.mathworks.com/help/releases/R2025b/matlab/matlab_oop/methods-in-separate-files.html)。

边绘制只保留一个实现，PolytopalMesh 可复用其已有 Graph 的 plot；确需共享的数值 helper
置于少量 data.mesh.internal 函数，不未经验证假定 sibling private 对 @Class 可见。
端点保留以预算可行为前提：最小顶点数为各因子 min(2,n) 的乘积，预算不足稳定拒绝，
不能同时承诺二维首尾点完整且 maxVertices=3。抽样后仍为诱导子图，不补造原图没有的边。
保留边界、标签、图例、投影、着色、采样预算和空/退化输入政策。当前 Triangulation 不走的
一般多面体分支只承诺稳定拒绝或明确支持，不为未来通用性新增多面体算法。

## 样本集合与显示

### 用值结构代替 Dataset / Database

```matlab
view = struct('Values', values, 'ComponentIds', [2,4]);
view.Axes = {x,y};
fields.u = view;
datasets.RUN = struct('Type', "numeric", 'Data', fields);
group = mdemat.tools.visual.selectComponent(datasets, "u", 4);
```

Type+Data 是普通 scalar struct，集合是普通 scalar struct-of-struct。Type 为
fixed/exact/numeric 来源和默认样式描述；它不再驱动工具去采样 FEM、参考解或 DOF。
数据由上层完整提交。物理变量 x/v/Values 都是 Data 的字段，其坐标在各自 view 内。
不为 fieldnames/isfield/numel 再提供一套 getter/setter/计数缓存。

selectComponent 是纯选择：按字段和显式 ComponentIds 找行，返回保留数据集 key 的普通 group，
附可选缺失清单。默认跳过不含目标字段/分量的数据集，显式 required 集合缺失则报错。
不解析 u2，不排序坐标，不 reshape 成网格，不插值，不创建图窗，也不返回另一层对象。
空集合用 isempty(fieldnames(datasets))；构造中防止 cell 参数展开成意外 struct 数组。
集合 ID 明确为合法 MATLAB 标识符；任意显示文字放布局 label，不用 makeValidName 静默合并。

| 容器 | 选择理由 |
|---|---|
| scalar struct / struct-of-struct | 当前小型、已命名的集合优先；没有额外类身份 |
| struct 数组 | 真有任意字符串 ID 或显式记录序列需求时选用，不同时维护第二份索引真相 |
| table | 用于收敛摘要/导出；不用于不同坐标和异形大数组的热路径容器 |
| dictionary | 任意 key、大规模动态查询有明确需求时评估；它本身是值语义，不能误称 handle |
| 自定义 Dataset/Database 类 | 已移除；其存取与集合操作由普通值字段和选择函数表达 |

选择 struct 首先减少别名与职责复杂度，不声称必然快于类/dictionary。MATLAB 可共享未修改
数值数组，修改和索引成本需测；struct 内放 handle 仍浅共享。[官方值/句柄语义](https://www.mathworks.com/help/releases/R2025b/matlab/matlab_oop/comparing-handle-and-value-classes.html)。

### 每字段 view 合同

* Values 为 C×N；ComponentIds 为 1×C 唯一正整数，只在该物理字段内有意义。
  [2,4] 意味着两行分别对应原分量 2、4，不能再次取第 4 行。
* Coordinates 为 d×N，或 Axes 为 1×d 的坐标轴 cell；恰好一种，禁止两份互相矛盾的坐标。
  Axes 为严格单调的 tensor 轴，N 为各轴长度乘积；按 ndgrid 的列主序逐分量展平。
  不要求为 Axes 再生成完整 Coordinates。
* Point view 保持列顺序和重复点。DG 元素/侧迹及必要 Connectivity 由上层明确提供，
  不自动 unique、平均、全局网格化或跨单元连线。重复 tensor 迹采用 Point view。
* 空 C×0、单点轴、非均匀轴、complex/sparse 的接受或显示政策逐项明确；不以 squeeze 猜维数。
* Times/Samples 在上层先配对；低秩/高维对象按请求点采样。工具不导入它们的类型或展开全场。

Coordinates→Axes 的转换由 common.fem 的 plotSamplingPlan/plotFieldSamples 负责：依据真实 Grid
拓扑和参考节点建立 Permutation，验证完整 tensor 对应后同时置换 Values，再输出 Axes。
不把任意点云 unique 一遍就宣称结构网格；不完整或有 DG 重复迹时保留 Point view。
Point view 默认可作显式 points 显示；1D 连线需要明确顺序/连通性，2D 场面需要 Connectivity。
3D/slice1d 的体插值要求合法 tensor，或由上层直接提供已采样切片；不支持时稳定拒绝，
不偷偷插值或丢失原有结构网格 image/slice 路径。

### 显示算法与资源分开

```matlab
style = mdemat.tools.visual.plotStyle();
h = mdemat.tools.visual.plotField(ax, datasets, "u", 4, mode="2d", style=style, datasetIds="RUN");
```

style 为值记录，plotField 是有显式图形副作用、无持久状态的函数。1D、2D、3D slice、物理
对角线采用局部/私有函数，颜色边界只有一个实现，不再继承 Strategy。
普通工具独立调用时，调用者拥有返回图元；已有颜色方案和有效样式选项均有映射。
原生研究绘图通过 curveStyle(style,colorIndex,patternIndex,kind=...) 取得属性值，
再用 namedargs2cell 传入 plot，与 plotField 共用配色/线型/marker 选择。
布局绑定含固定有序 DatasetIds：1d/slice1d 共用一轴，2d/3d 每个数据集对应一个独立轴，
轴数必须与绑定数量一致。样本缺失只清理该绑定的自有 artists 并显示无数据，不压缩轴或改变
其他数据集的颜色/标签身份；新增数据集通过显式重配布局加入。共享色限只据当前有效样本计算。
独立 plotField 调用也使用明确 DatasetIds 顺序，不用每帧临时 fieldnames 位置绑定身份。

Visualizer 仅保留私有图窗/布局/axes/自建 artists 和生命周期；接收普通 layout specs、style
及 render(datasets,time=...,title=...)。不再持有 Density、Exacts、Timeline、RefXNodes、
空间 factory、State、插值历史或独立 Database。采样选择属于上层，绘图绑定明确字段/分量 ID。
布局和定位算法放小函数，避免合并成巨型类。保留跨格布局、自定义绘图 callback、自动定位、
全局标题与共享色限；TileSpacing/Padding 是普通布局值。自定义 callback 必须返回其
新建/更新的 artists（句柄向量或 data.plot 的分层句柄 struct）以声明所有权。
内置 plotField/createLayout 对未提交的新资源立即登记并在失败时撤销；临时状态用 onCleanup 恢复，必要时使用
自有暂存容器，不以借用轴的所有 Children 推断所有权。自定义 callback 仅运行在 session
专属自有轴中（借用图窗内可新建自有轴），不直接运行在借用轴上；失败清理该专属轴并允许
重建。返回句柄必须属于所给轴，不能认领他人对象。任意 callback 在其他位置的副作用不属于
自动回滚承诺；新增图元后 throw、返回错误归属句柄都要注入验收。

正常和异常退出恢复临时 NextPlot/根单位，InRender 用 onCleanup 复位；只删除自建图元。
新核不清空借用轴的他人内容；旧交互 plot 的 hold-off replace 语义作为明确接口选择保留。
配置/图窗/布局变化先验证再创建，并检查 axes 归属。close 尊重回调，veto 保留登记并报告；
显式 force 仅作用于自有图窗。finalize 终止 session，重复清理可重试；reset 先解决未关闭
自有资源，再重新创建布局。借用图窗不关闭。避免位掩码组合产生“已结束但还初始化”的状态。

非均匀轴不能仅用 imagesc 的两端点；先明确节点/像素中心/单元值，再选 image/surface/patch，
防止默认 pcolor 丢末行末列。保留正确的 3D 排列、物理对角线/有符号弧长；log→linear 必须
显式恢复。NaN/Inf/常数/复数色限有明确政策。重复渲染图元数量有界，是否复用图元由量测决定。

## Config 数据与文件解析

tools.environment.Config 是值类，私有 Entries 是 scalar struct，每键只保存
Value、Validator。没有独立 Entry 类，不维护可分开写入的两张表，不重载 subsref/subsasgn。

```matlab
cfg = mdemat.tools.environment.Config();
cfg = cfg.define("dt", [], @isValidStep);
cfg = cfg.set("dt", 0.01);
cfg2 = cfg.update(struct('dt',0.02));
[value, validator] = cfg.get("dt");
cfg.validate(["dt"]);
values = cfg.toStruct();                    % 仅值，不含 validator
[patch, info] = mdemat.tools.environment.readConfig(file);
cfg = cfg.update(patch);
```

最小方法为 define/get/has/keys/set/update/validate/toStruct。define 维护完整 pair；重复定义
须显式 override，set/update 仅换值。修改必须赋回；不自动修改别名中的 Config。
不为来源、事件、全局缓存、阶段、Required 再建元数据框架。validator 为只读、只调用一次的
logical 标量 predicate（每次显式键校验调用一次）；[] 表示无 validator；不按 nargout 猜
无返回 mustBe* 协议。显式再次 validate 是一次新校验，不承诺一次实例生命周期只调用一次。

组装期允许 deferred 默认，set/update 验写入项，validate(keys) 验当前操作所需项。
[] 可以是合法空值、未选方案或待提供数据，不能统一当 missing。不能在 Scheme 构造时
无条件全键 validate：resume 可不需要 ic，而重新 initialize 需要。需求由操作拥有者给定。
跨字段关系用拥有者的纯函数校验完整候选，不让 validator 闭包捕获 Scheme/可变配置。
值语义不递归复制 Value 内的 Timer 等 handle，也不回滚 validator 的任意外部副作用。

优先级固定父默认→子 override→文件→显式 patch→纯规范化。先合并最终候选，再检查最终
有效值和跨字段约束，避免文件的中间值与最终覆盖顺序偶然决定合法性；语法错误仍立即拒绝。
Config.update 自身原子；Scheme 的 file+patch+规范化也须在本地候选上完成，成功才提交。

readConfig 只返回 OPTION 的值；info 提供来源/行号和 Paths，避免 paths 键冲突。
它负责 UTF-8/BOM、文件行号、源相对资源、小型允许语法与重复键规则，不 addpath、不 eval、
不创建 Config、不注册 validator、不实例化 Scheme。新读取入口拒绝重复 OPTION；旧文件
若确有 last-wins 需求应显式迁移，不把继承 override 与文件重复混为一谈。
连幂及一元符号/负指数按 MATLAB 子集定义，不承诺任意 MATLAB 表达式。

新 Config 拒绝未声明键；自定义键通过 define 明确增加。现有 schemeName、vBasisType、
useSplitting、customParam 和 configure 派生字段不能直接丢弃或被锁死：T200 逐项区分实验
元数据、数值输入、派生描述、观察/输出设置与运行资源，按拥有者迁移/注册。未知拼写不能
被静默忽略。纯解析仍可以接受任意合法 OPTION 名。
每个实际键核对类型、标量/向量 shape、单位、空值、有限性和操作前置条件；不能仅把旧宽松
validator 原样搬进新类，也不一刀切禁止具有数学含义的 Inf 边界。

Scheme 的 config 输入接受 scalar struct patch 或 scalar Config。struct 只提供已注册键的值；
显式 Config 可声明扩展键，例如：

```matlab
cfg = mdemat.tools.environment.Config();
cfg = cfg.define("tFinal", 1, []);
cfg = cfg.define("customParam", 42, @(v) isscalar(v) && isfinite(v));
scheme = mdemat.solver.hyperbolic.UwdgScheme(config=cfg);
```

导入已知键时只采用 Value，并重新使用 Scheme 拥有者的 validator；外来 validator 不能替换
现有数值规则。新扩展键保留完整 pair，不能占用派生/运行资源保留名。外来 Config 的键视为
显式输入，不把它当包含未解析默认值的模板；partial/resume 仅提供本操作需要覆盖的键即可。

Scheme 组合 Config；数值计算读取一次生成的只读普通 struct 快照，避免每次访问都重建视图。
Scheme 构造完成后数值配置冻结；运行对象依据该快照创建。重新设置方法/终止时间须创建新的 Scheme，
不能只更新 Config 字段。显式 resume 的允许变化单独定义/验收，不引入 Config 通知框架。
冻结点为构造结束，冻结全部数值绑定，不维护逐键依赖图。deferred 默认仍可留在 partial 实例，
它可用于不需该键的 resume/读取；initialize 若需要尚不合法的 ic 等则稳定拒绝，用户补齐
创建前 Config 后新建实例。这样不会在已有积分器上半途填入新的数值模型。
观察/输出资源独立管理，不因换显示组件而重建积分器。过渡用 Config/Validators 直接写入和
configure 原地写法必须全部映射后移除；数据结构正确不等于派生数值对象已经同步。
派生缓存和可变运行状态归 workspace，不回写冻结 Config。冻结的是数值绑定，不是假定其中
任意 callback/handle 都被递归冻结；实际模型值/修订仍由拥有者管理。

## 分析：采样、归约、实验和显示分离

Metric 及两个 FEM Metric 类已退出；其采样能力迁为 common.fem 普通函数，
不保留另一个空基类。Analyzer 也已退出工具层：它的采样/参考/组件选择/level 循环迁至
execution 或研究入口，观察数据为普通 struct 数组。工具只保留数值归约、收敛汇总与格式化。

```text
上层选择 field / component IDs / reference / physical representation
  → common.fem 采样计划 + 物理值/权重
  → tools.analyze.weightedError
  → 上层明确归约并收集 observations
  → tools.analyze.summarizeConvergence
  → 普通结果 / 显式 visual 绘图 / IO 值记录
```

samplingPlan/nestedSamplingPlan 是普通函数，结果只含明确布局的节点、物理坐标/权重和
必要的 NestLI/对应关系；不携带 Space/Integrator handle 或缓存 factory。保留非均匀 Jacobian、
实际坐标对齐、同网格、已有多级整数比例参考比较；不能把全部研究缩成只支持一次二分。
一般非嵌套网格用显式转移/共同点策略或稳定拒绝，不自动投影改变问题。
新采样选项明确每方向点数，标量表示各方向同一数量；不继续让一个 Density 在显示中 floor、
在分析中 ceil 而含义不同。仓外旧调用若采用总密度，迁移时须显式换算为每方向点数，
再核对原显示/积分规则；不把相同标量在新接口下解释为原来的总点数。
显示内点采样与范数积分的规则由上层分别选择，不因共享采样计划就强制相同；用于严格范数
的权重须非负，不能把有负权的积分近似无说明地当作范数。

令 e=U-Uref、非负物理权重 w 包含 Jacobian，q 仅遍历正权重点：

\[
 a_{1,c}=\sum_qw_q|e_{cq}|,\quad
 a_{2,c}=\sqrt{\sum_qw_q|e_{cq}|^2},\quad
 a_{\infty,c}=\max_q|e_{cq}|.
\]

weightedError 提供 ComponentL1/ComponentL2/ComponentLInf，以及已有的
L1=Σc a1,c、L2=||a2||2、LInf=maxc a∞,c；不按体积额外归一化。
旧 sequence 为 Ep=||ap||2：L1 与 LInf 均不同于上面的全分量归约，必须在结果中明确标为
component 2-norm after spatial norm。数值 [2,4]+sequence 的身份就是所选集合，不叫 U2。
原始 samples 独立于严格 norm 枚举；未知名字不能返回原始数组冒充一个误差数。

模式补零只有在同一嵌套基、归一化、尺度和模态 ID 都已证明一致时成立。一般 kinetic 比较
须由 application 重构到共同物理点及测度；ApHermite 的 sqrtRhoInf、宏微分解和速度 Gram
不能由工具猜测。系数平方和只有在对应正交归一测度下才等于相应物理 L2。

旧 Richardson Metric 实际给 interlevel difference，没有除 r^p-1。保留旧差值的数学结果，
准确命名；真正细网格误差估计 ||uh-uH||/(r^p-1) 须显式给 p、r 和渐近假设，不自动改变旧结果。
每字段可选独立 exact/interlevel/unassessed 策略，不因任一字段有 exact 就假设全部都有。
观察记录包含 Scale、Field、ComponentIds、SpatialNorm、ComponentNorm、Comparison、Reference
和 Error；汇总按这些显式身份对齐，不以拼接名称猜分量，不隐式补齐尚未计算的 level。

显式数值选择保留有序身份：[2,4]、[4,2] 和 all 是不同序列。研究中需要逐 level 增加速度
模态数时，`ComponentIds='all'` 只用于 `ComponentNorm='sequence'`，表示“该 level 的全部
分量”；同时每条记录必须提供 `ResolvedComponentIds=1:N`。汇总按 all 的研究身份合并，
另用逐 level 的 ResolvedComponentIds 保存真实支持，避免 3/6/12 个模态被误拆成多条单点曲线。
all 只属于分析观察身份；显示 view 的 ComponentIds 仍须是各行实际的数值 ID。
这只表达选择政策，不证明不同基、归一化或物理重构之间可比；后者仍由模型/研究拥有者校验。

收敛表的 NComponents 记录各 level 的分量数，Properties.UserData 保存数学身份和完整的
ResolvedComponentIds；同表的空间范数由各自误差/阶数列表示。struct table key 只是有界、唯一的
显示标签，可能为合法标识符而转换、截短或加后缀，不能从 key 反解字段、组件或参考来源，
也不能把 key 当作持久化的身份元数据。相同身份下不同范数列必须具有相同尺度和实际支持。

\[
 p_i=\frac{\log(E_i/E_{i+1})}{\log(h_i/h_{i+1})}.
\]

实现对相近正数用 log1p((a-b)/b)，相差大时用 log(a)-log(b)，同时避免比值溢出与 log 消减。
零误差未定义，负阶数保留，接近尺度的病态与浮点结果有限性分别报告，不用任意容差造阶数。
weightedError 避免可避免的中间差值/加权溢出；真实不可表示的 LInf 保持明确非有限状态。
参考阈值、零权重、RelativeDefined 和比值有限性分别定义，保留 finite 输入合同。

### 收敛图保留为显式数组函数

旧 Analyzer.draw/drawImpl 的 log-log 曲线和指定参考阶线由 plotConvergence 接替，
多字段布局仍是普通 layout specs 或原生 tiledlayout，不再由分析对象持有图窗。

```matlab
series = mdemat.tools.analyze.summarizeConvergence(observations);
fig = figure; ax = axes(fig);
[h, fit] = mdemat.tools.visual.plotConvergence(ax, ...
    series(1).Scales, series(1).Errors, referenceOrder=2);
% 同一数学身份的拟合：fitted=true, fitRange=[firstIndex,lastIndex]
```

已知返回的收敛表含 L2_Error 时，也可显式传入 `table.Scale` 和 `table.L2_Error`；
绘图函数本身不识别 table、Summary、模型或 State。多范数共轴由调用者逐条传数组，
使用固定 colorIndex 和 label；不同物理字段、分量集合、归约和参考来源不得混成一条序列。
原表格的首列现在是数值 Scale，组件集合及范数身份显式命名，属于源 API 迁移。

拟合仅使用含首尾的 fitRange 中正的有限误差。设锚点为其中首个有效样本 a，
采用 x=log(h/ha)、y=log(E/Ea) 的中心化最小二乘求经验阶 p；相邻数值的 log 比用
log1p 避免消减。指定参考阶 q 的曲线严格为 Ea*(h/ha)^q；默认使用首个正误差归一，
与原有效参考线一致。显式 referenceOrder 优先于拟合参考线；不把经验拟合当作精度证明。
零/NaN 不进入 log 视图，负值/Inf 拒绝；参考线的有效不同尺度不足或不能表示时稳定报告。
返回图元和 Order/LogAmplitude/FitIndices/ReferenceOrder/Anchor 普通值，供独立核对。

## 其余工具与持久化

* base 保留小函数，不合成 Utils 类；全元素/文本列表/空值/子类与 exact class 合同统一。
* except 三文件退役，调用点改原生设施和完整语义 ID。保留有效既有 ID，变化逐项映射；
  先查每个 warning 是否可继续。已格式化消息以 error(id,'%s',message) 传递，保留原 cause。
* Timer 仅私有持有 ActiveIndex/Token，记录是 Name/Duration/IsActive 值快照。保留同名累计与
  互斥，Duration 默认已停止区间，墙钟/零值/活动报告明确。cleanup 只处理自己启动的区间；
  不为 scope 再造类，不允许快照改坏 owner。报告与格式化尽量分开。
* hash 修逐字节两位 hex；所测 uint8 raw digest 正确不代表全部类型已验。LinearSolver 真正
  复用依据是 isequaln，停止强制构建无消费者的诊断 hash；按需身份不阻塞无 JVM 的常规求解。
* 根/版本/环境为短函数；修 nested checkout、类型/失败上下文；不建 EnvironmentManager。
  benchmark 保留显式重复/预热与调用者副作用；memorySnapshot 只称原生进程快照，未知平台
  明确 Unavailable，不冒充峰值。已安装产品不等于依赖或许可证可用。

IO 保留受限 primitive grammar、精确类型/空 shape、UTF-16 文本、little-endian/base64、
显式 callback、默认拒覆盖和同目录独占提交。Config/validator、State handle、图形和 table
不会因改成 struct 外层就自动可保存；RunRecord 只接收显式值及分离报告，不替实验判准确性。

读取：文件限额→MAT whos 元数据/限额 manifest 或 JSON 解析→版本/完整描述/预算校验
→按引用加载/解码→重建值→上层语义。complex/sparse 为标量标志；missing 是逐元素 mask，
保留合法 JSON 空数组。实虚部独立精确字节数；单叶/累计字节、shape/容器/节点/深度及
稀疏 CSC 列指针均计预算。解析本身仍分配，不宣称任意输入零分配或完整内存隔离。
MAT 不支持的大变量先明确拒绝，临时文件变量清单核验后才提交；不在普通 CI 分配超大数组。

State schema/回调由 State 拥有者负责，checkpoint schema/完成步/ODE 历史由 execution 负责。
旧 readStateRecord/validateStateRecord 退出 tools。State.load 只认 State；checkpoint 的
“只提取 State”用途改为 execution 显式提取再交相应 fromRecord，保留用途而不让 data 识别
执行格式。不把上行依赖简单搬到 data helper，也不留永久 tools→data 代理。

checkpoint 数值身份由模型拥有者明确选取，观察组件/exact/标题不影响数值一致性；真实模型
参数变化仍拒绝。必要的嵌套模型 callback 用明确路径重绑，参考解不进入数值身份；不保存
validator、闭包源码，也不以 func2str 证明模型一致。旧 v1 按旧布局读取后，在 execution
投影到同一数值身份再比较；不直接把新描述与旧的含观察字段完整结构 isequaln，也不改盘上 v1。
execution 自动命名的 checkpoint/state 输出用 accepted-step、图像用观察帧序号命名；
State.save(file) 仍使用调用者显式文件名，data 不获知积分步。精确时间在元数据中，
不用显示舍入字符串作为唯一身份。

既有合法 v1 数值/State/空间/ODE/采样研究记录继续可读；内部新 view 不自动改
mdemat.sampledStudy 的 Points/Values、轴序或其他 wire keys。修坏输入接受不需要改正常 codec。
确需扩展 checkpoint 的 callback/描述能力时单独升内层 schema，并实现旧 reader/迁移；
不只改版本 metadata，不默默改 v1 语义。源 API 重构与文件兼容分别验收。

## 实施、验证与效率

沿用 T200–T290；源码迁移与交付验收已完成，准确状态和证据只在 tasks.csv 维护。
下表保留每项验收责任，供后续修改核对；历史探针执行不计作修复回归。

| 任务 | 交付关口 |
|---|---|
| T200 | 每个类/方法/属性/选项到新函数或值字段的映射；实际配置键分类、错误 ID、文件 fixture、独立反例；冻结 source API 变更与文件兼容策略 |
| T210 | 谓词修复；全仓原生异常迁移；失败/警告与稳定 ID；清除 backtrack |
| T220 | Timer 单 owner/快照、hash 正确且可选诊断、环境边界与 cleanup |
| T230 | value Config/readConfig；Scheme 候选事务、按操作验证、静态 predicate、冻结真实运行配置 |
| T240 | IO 预算/精确 codec/MAT 提交；State/checkpoint 职责和身份/回调/命名；旧文件适配 |
| T250 | 取消 Metric/Analyzer；FEM/应用采样、纯数值归约与汇总；保留各范数/参考研究语义 |
| T260 | data 的 plot 方法；普通样本集合与选择；取消 Strategy/多重资源 owner；几何/时间/样式正确性 |
| T270 | 所有生产、示例、测试、help、build/包入口迁移；旧类退出默认生产包；工具无上层类型/语义依赖 |
| T280 | 受控性能/内存/图元基准；仅实施有收益的局部优化，没收益则保留简单实现 |
| T290 | 全数值家族/配置/重启/原生测试、analyzer/覆盖、构建和全新安装；来源与证据闭合 |

T230/T240/T250/T260 在基础合同就绪后可并行，每项都自带局部回归；不把测试拖到最后。
首先设计新接口和独立 oracle，再迁调用方，最后移除重复类。已无不可替代能力的类不因
“至少保留一个周期”而成为永久架构例外；仓内调用一次迁移，源 API 变化有清晰迁移表。
若明确出现仓外兼容需求，提供单独、有边界的迁移工具，不纳入 tools 正式依赖图。
源 API 采用 2.0 主版本迁移；文件 schema 与源 API 版本分别处理，合法 v1 文件继续按原契约读取。
版本 metadata、重建产物和全新安装已核对相符，当前完整验收为 PASS。

后续修改须保留以下验收覆盖：

* 全工具文件覆盖、所有公开方法/选项能力映射；删除结构身份断言时保留原有有效行为要求。
* struct 副本隔离/显式赋回，内嵌 handle 浅共享；Config 第二键失败、file+override/跨字段失败；
  构造/initialize/resume 差异、默认/扩展/派生键；配置与实际 TDisc 一致，不测字段就算成功。
* 外置方法路径/继承、0×2 无边图、数据不被绘图修改、边端点/连接关系、投影/边界/图例。
* x/v/含数字字段、[2]/[2,4]/[4,2]、分量聚合身份、空集合/空 shape、非均匀内部位置/末行列、
  DG 双侧迹、u(t)=t 变步长配对；保留 3D 仿射切片/矩形对角线独立答案。
* frame/axis/图窗替换、异常重试、veto/reset、借用资源、log→linear；自定义 callback 和跨格布局。
* 分离支持函数的 sequence 范数、错误 norm 拒绝、物理基/Gram/重构、非恒定 sqrtRhoInf、
  非均匀 Jacobian、多级共同坐标；相近大尺度、极值误差、零 reference/权重与有限性。
* v1 golden MAT/JSON、primitive/Unicode/missing、坏 flags/实虚字节/shape、低预算、无大分配
  大小判断、保存/提交故障旧目标保全；观察改变仍可恢复、模型改变拒绝、嵌套 callback 和命名唯一。
* 所有原 numerical family、FDM/低秩、完整研究工作流、无绘图求解、新进程外部 cwd 安装；
  @Class 被打包且无旧类并存。不能只跑新纯函数测试而放弃应用集成。

T270 同步真实 MATLAB 项目成员和打包清单，以及来源/声明审计的路径与职责合并映射；
不能只移动文件却让旧 namespace 映射把 ledger 刷回不存在路径，或把有意的类删除记成能力遗漏。

统一风格沿用 [standards.md](standards.md)：英文 help、UTF-8/LF/四空格；类/公开属性与值记录
字段 UpperCamelCase，函数/方法/变量/name-value lowerCamelCase，常量 UPPER_SNAKE_CASE。
物理字段和已有 wire keys 保留身份。不用全仓排版掩盖语义变化，不为了风格引入更多 wrapper。

效率只在相同答案下量测：当前容器与普通值集合的构造/分量选择/更新、FEM 采样与归约、
图元批绘/更新、JSON/MAT 字节量与叶子数、hash/metadata/Timer 开销。区分 dense/sparse/complex、
大数组与小叶子，报告实际原生进程、输入、预热/重复、随机流和测量范围；内存区分快照/采样
最大值/真实峰值，普通单测无硬件时间阈值。没有证据就不增加全局缓存或强称 struct 提速。

源码中的值结构、外置方法、纯分析函数和单资源 owner 已实施，仓内调用已同步迁移。
本轮完整验收与交付关口已闭合，接受状态、源码及产物身份记录在对应正式证据中。
首轮机器证据和重构前探针仍保留，用于解释修复依据，不是新 API 的测试入口。

历史探针 [runToolsDesignReview.m](evidence/runToolsDesignReview.m) 对应其记录的旧源码身份；
在当前源码直接运行时，已删除类及接口的结果不再代表原始反例。当前使用 runTests 和构建关口，
不能把历史探针计数当作修复后的通过数量。效率结论只采纳受控原生量测；没有实测支持时，
保持简单实现，不声称少用类就必然提速。
