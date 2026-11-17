# 上边界压强：从“定死”改为“由内部决定”——导师公式的解读与讨论

**Status: PLAN** — 2026-10-27 的讨论记录，尚无代码改动、没有新的测量。导师手写公式对应的闭合在代码里已经存在（`ISO_OUTER_P_NEUMANN=2`）；已有的测量见 [outer_bc_pressure_temperature_swap_recap.md](../studies/boundaries/outer_bc_pressure_temperature_swap_recap.md)，但那次测量早于 hydrostatic reconstruction（2026-10-13）和 single-T 上边界（2026-09-22），需要在当前 release 数值方法下重跑（§7）。

## 1. 公式在说什么

导师写的两行（`b` = boundary，即上边界 ghost 中心；`n` = 最后一个物理格子 N，也就是顶格；$\Delta x$ = 两者中心的距离；$T$ 写成 $p = T\rho$，即 $T$ 已吸收了 $(1+x)k_B/m_H$）：

$$
p_b = p(T_b,\rho_b) = T_b\rho_b, \qquad p_n - p_b = \tfrac12(\rho_n+\rho_b)\,g\,\Delta x
$$

$$
\frac{p_n}{\Delta x} - \tfrac12\rho_n g \;=\; \frac{T_b\rho_b}{\Delta x} + \tfrac12\rho_b g
$$

第一行是两个条件：

1. **EOS**：ghost 的压强由 ghost 的温度和密度决定。
2. **离散静力平衡（梯形）**：从顶格 N 到 ghost，压强按 $dp/dx = -\rho g$ 下降一个格距。

第二行是把 (1) 代入 (2) 再除以 $\Delta x$，目的是解出唯一未知数 $\rho_b$：

$$
\rho_b = \frac{p_n - \tfrac12\rho_n g\Delta x}{T_b + \tfrac12 g\Delta x}, \qquad p_b = T_b\rho_b .
$$

也就是说：**给定 $T_b$，ghost 的 $\rho_b$ 和 $p_b$ 都由顶格的 live 状态 $(p_n,\rho_n)$ 加静力平衡推出来**，外部只提供一个温度。你的理解是对的：密度和压强都不再定死，唯一定死的是温度。

两点细节：

- 在 Saha EOS 下 $T_b$ 里的 $(1+x)$ 依赖 $\rho_b$，所以需要不动点迭代。`hse_rung`（`scenarios/model_column.cpp`）已经这样做了。
- release 的 hydrostatic reconstruction 平衡的是**对数形式**：$\ln p_b - \ln p_n = -\tfrac{g\Delta x}{2}\left(\rho_n/p_n + \rho_b/p_b\right)$。导师的梯形式是它的一阶展开，两者相差 $O((\Delta x/H)^2)$。顶格 $\Delta x \approx 276$ m，$H \approx 1300$ km，所以数值上可以忽略。但 ghost 必须用对数形式构造，静止柱才是重构的精确不动点，否则上边界会重新出现一个残差。`hse_rung` 在 `kHseRecon` 下用的正是对数形式。顺带一提，$\rho_b/p_b = 1/T_b$ 只依赖 $T_b$（和 $x_b$），所以 $T_b$ 固定时 ghost 那一项几乎是常数。

速度公式里没写。代码里 ghost 速度取顶格速度，并按 $0.1\,c_s$ 做 Mach cap。

**它在代码里已经存在**：`ISO_OUTER_P_NEUMANN=2` 正是这个闭合，即“与 release 相同的 EOS-closed 静力阶梯，但锚在 live 顶格上，而不是固定的 `kOuterPRef`”。release 现在的 hydro ghost 温度固定为 22 kK（`ISO_HYDRO_T_DECOUPLE=0`），这正是公式里的 $T_b$。所以导师的方案 = release preset + `ISO_OUTER_P_NEUMANN=2`，不需要写新代码。

## 2. 现在的 release 边界条件

| | 上边界（2153 km，TR 底部，subsonic **outflow**） | 下边界（1600 km，中色球截断，subsonic **inflow**） |
|---|---|---|
| ghost 压强 | **固定** `kOuterPRef`（IC 时从顶格单侧静力外推一次，≈ 0.0103 Pa），第二个 ghost 再上一个静力阶 | **固定** 储层静力阶梯 `kInnerPRes0/1`（IC 时捕获） |
| ghost 温度 | **固定** 22 kK（hydro ghost 和 conduction 的 Dirichlet 点是同一个） | **固定** `kInnerTRef`（conduction 默认 Dirichlet） |
| ghost 密度 | EOS$(p_g, T_g)$，派生量 | EOS，派生量 |
| ghost 速度 | 顶格速度外推，按 $0.1\,c_s$ 截断 | $V = 0$ |

所以你记得没错：**上边界的温度和压强都是定死的**。

## 3. 是不是超定了？

（本节在讨论中修正过一次：最初的版本把条件贴上了“hydro”和“thermal”的标签，说温度“只算热传导条件”。这个说法不对，因为 $T$ 同时出现在流体和热传导里。下面是修正后的说法。）

### 3.1 纯 Euler：$T$ 本身就是一个合法的流体条件

把线性化的 1-D Euler 扰动写成三个特征振幅：$w_\pm = \delta p \pm Z\,\delta u$（声波，速度 $u\pm c$）和 $w_s = \delta p - c^2\delta\rho$（熵波，速度 $u$），其中 $Z=\rho c$。由 $\delta T/T = \delta p/p - \delta\rho/\rho$ 和 $c^2=\gamma p/\rho$ 可得

$$
\frac{\delta T}{T} \;=\; \frac{\gamma-1}{\gamma}\,\frac{w_+ + w_-}{2p} \;+\; \frac{w_s}{\gamma p}.
$$

在 subsonic outflow 面上，$w_+$ 和 $w_s$ 往外走，只有 $w_-$ 往里走。令 $\delta T = 0$，得到 $w_- = -w_+ - 2w_s/(\gamma-1)$。只要 $\gamma \ne 1$，这就把入射波定下来了。所以**对纯 Euler 来说，只固定 $T$ 是一个合法的单一外部条件**。温度确实“既是流体也是热”，这个直觉是对的。代价是它会反射：往外走的熵扰动会被转成往回走的声波。

这个式子也说明了 $T$ 的局限：系数是 $(\gamma-1)/\gamma$。在等温极限 $\gamma\to1$ 下，声波不带温度扰动，固定 $T$ 就完全控制不了入射声波。

### 3.2 加上热传导以后，$T$ 只能用一次

热传导让能量方程变成关于 $T$ 的二阶方程，于是每个边界多需要一个条件。按不完全抛物系统的计数方法（Strikwerda 1977）：抛物部分每端 1 个；剩下的双曲部分 $(\rho,u)$ 的特征速度是 $u\pm c_T$（$c_T$ 是**等温**声速），subsonic 时每端各进来 1 条。所以 subsonic outflow 面上需要 **2 个独立条件**（纯 Euler 是 1 个）。*用 $c_T$ 作剩余双曲部分速度的这一步是我自己的推导，值得和导师核对。*

“2 个独立条件”并不说哪个属于流体、哪个属于热传导，因为能量方程只有一个，$T$ 两边都管。它说的是：固定 $T$ 只占一个名额，不能占两次，另一个名额必须由与 $T$ 独立的量来填（$p$、$\rho$、入射 Riemann 不变量……）。从 3.1 的式子也能直观看出原因：热传导把短波推向等温（$\gamma\to1$ 的行为），恰好是 $T$ 看不见声波的那个极限。

### 3.3 定常计数：$T$ 和 $p$ 是两个独立参数

更直接的论证来自定常解。定常方程是 $\rho u = J$、$d(\rho u^2+p)/ds = -\rho g$，以及关于 $T$ 的二阶能量方程。积分常数有 4 个：$J$、压强水平、能量通量、$T$ 的另一个积分常数。所以需要 4 个边界数据。

| | 上边界给 | 下边界给 | 个数 | 结果 |
|---|---|---|---|---|
| release | $T_\text{top}$, $p_\text{top}$ | $T_\text{base}$, $p_\text{base}$ | 4 | $J$ 被确定 |
| 导师方案 | $T_\text{top}$（外推的 $p$ 在定常时自动满足，没有信息量） | $T_\text{base}$, $p_\text{base}$ | 3 | **差 1 个：$J$ 不被边界确定** |

这和我们自己的解析表一致：$J = J(T_\text{top}, p_\text{top}, T_c)$ 里 $T_\text{top}$ 和 $p_\text{top}$ 是两个独立参数，给了 $T$ 并不等于给了 $p$。

### 3.4 release 是不是超定

不是。上边界给 $p$ 和 $T$，正好是 2 个独立条件；定常计数也正好是 4 个。`outer_temperature_bc_single_t_recap.md` 里“formally imposes two conditions”只按纯 Euler 计数，没有把热传导算进去。

另外，离散上 ghost-cell + Riemann solver 是**弱施加**：ghost 给的是整个状态，solver 只取特征线允许进来的信息。在 outflow 面上，ghost 的 $T$（即 ghost 密度）主要落在往外走的接触间断上，只通过声阻抗进入 $p^*$。如果流动反转成 inflow，熵波变成入射，ghost 的 $T$ 就直接决定进来的熵。

## 4. 导师的方案在计数上是什么

上边界有两个条件：

1. $T_b$ 固定（Dirichlet）；
2. $p_b$ 由顶格静力外推。这是一个条件（Neumann 型，$dp/ds = -\rho g$），**不是“什么都没给”**。在离散上它等价于“入射声波振幅为零”：hydrostatic reconstruction 之后 ghost 和顶格之间没有压强跳跃，速度也相同，Riemann 问题不产生入射波。

所以按 3.2 的计数它是**够数的、适定的**，属于标准的 non-reflecting outflow 条件。问题出在它选出的是哪个解。Non-reflecting 条件在定常意义下没有信息量（3.3 表中只有 3 个），于是平均压强水平没有外部的恢复力，会慢慢漂。这是 CFD 里的已知现象：Rudy & Strikwerda (1980, JCP 36, 55) 和 Poinsot & Lele (1992, JCP 101, 104) 都因此在 non-reflecting 条件里加了一个弛豫项 $K(p - p_\infty)$，用来把平均压强拉回一个给定值。

这时压强水平只剩底部储层来锚定：

$$
p_\text{top} \approx p_\text{base}\exp\!\left(-\int_\text{base}^\text{top}\frac{g}{(1+x)k_BT/m_H}\,ds\right).
$$

$T$ 在这里确实起作用：它决定了标高，也就是 $p_\text{top}/p_\text{base}$ 这个**比值**。但压强的**水平**来自 $p_\text{base}$。所以“由动态过程决定的压强”实际上是**由 1600 km 储层压强 + 柱内温度结构 + 演化历史决定的**。外部压强参数并没有消失，只是从 TR 顶部（AWSoM-R 给）挪到了中色球截断处（IC 给），而后者同样没人能告诉我们。

## 5. 已有的测量（早于当前数值方法）

`outer_bc_pressure_temperature_swap_recap.md` 用 `ISO_OUTER_T_HYDRO_WALL=1` + `ISO_OUTER_P_NEUMANN=2` 跑过 4000 s（N=500/R4，22 kK）。和 release（固定 back-pressure）比：

| | release | 只固定 T（导师方案） |
|---|---:|---:|
| `V_top` [m/s] | +9.50 | +0.49（**约 19 倍抑制**） |
| 顶面 `f_total` [kg m⁻² s⁻¹] | +2.70e-10 | +1.47e-11 |
| 底面 `f_total` | +2.71e-10（同号，一致上流） | **−2.42e-10**（质量从底部流出） |
| `p_top` 漂移 | ~0 | +5.8 %，仍在上升 |
| 柱质量 | 99.994 % | 98.97 % |
| 100 kK 墙 | 色球保持 6.5 kK | `p_top` +128 %，TR 下移到 ~1900 km，柱质量 −22.7 % |

这些现象和 §3.3 / §4 一致：

- `p_top` 上升：柱被 conduction 加热、标高变大，而唯一的压强锚在底部，所以顶部压强按上面的积分上升。100 kK 时柱子加热得更多，漂移也更大。
- 底部流出：加热膨胀需要释放，而唯一“握住压强”的面是底部储层，于是质量从底部排出，而不是从顶部蒸发。
- 蒸发被抑制：蒸发本质上是“TR 上方压强给定，加热的色球必须把质量推出去”。上方不再给压强，这个机制就消失了。

**一个重要的保留**：这组测量做于 2026-09 月，当时 hydro ghost T 还是 decoupled，重构还是普通 `ln p` 而不是 hydrostatic reconstruction。导师的动机很可能正是“limiter 修好以后，上边界的静力梯度不会再被误当成振荡去 limit，那就可以重新试一下只定温度”。§4 的论证不依赖 limiter，所以我预计结论不变，但这一点**没有测过**，应该重跑（§7）。

## 6. 和 AWSoM-R 耦合的含义

想让 AWSoM-R 不给压强，理由是成立的：AWSoM-R TR 顶部的压强来自它自己的解析静态 TR 闭合（Sokolov 2021 Eqs 19–22），而这个闭合正是我们要替代的，拿它当输入有循环之嫌。

但物理上信息流是这样的：在 loop 里，冕区压强由冕区能量平衡决定（RTV $p\propto T^3/L$），色球是响应的储层。冕区“要质量”这个信号，正是沿 $u-c$ 特征线、以压强的形式向下传的。如果交界面只传 $T$，我们的柱子就感受不到冕区压强。并且我们的上边界位置是固定的（2153 km），TR 不能像真实大气那样上下移动来适应不同的冕区压强，所以 $p_\text{TR}$ 会被锁死在 C7 层结给出的值附近。而 AWSoM-R 那边 TR 顶压强的实际范围大约是 $10^{-4}$–$10^{-1}$ Pa。

可以考虑的几个方向：

1. **保持 release**：AWSoM-R 给 $(T, p)$，我们回 $J$（即现有 $J(T,p)$ 查表的设计）。
2. **特征耦合**：AWSoM-R 只给入射的 Riemann 不变量（线性化为 $\delta u - \delta p/Z$），$p_b$ 由内外两侧共同决定，既不定死也不纯外推。这在概念上最接近导师“压强由动态决定”的意图，同时保留了那个与 $T$ 独立的外部条件。注意：下边界的 characteristic BC 已经测过并被拒（对蒸发影响 ≤ 0.4 %）。在低 Mach 数下它大概率退化为固定 back-pressure，上边界需要它自己的证据。
3. **导师方案**（只给 $T$）：需要先确认重跑后蒸发不再被抑制，并且要接受 `kInnerPRes` 成为新的自由压强参数。
4. **用热通量代替压强**（未验证的想法）：在定常计数里，$\{T_\text{top}, q_\text{top}, T_\text{base}, p_\text{base}\}$ 同样是 4 个数据，并且无辐射时 $J = q_\text{top}/(h_\text{top}-h_c)$，于是 $p_\text{top}$ 变成输出，AWSoM-R 只需要给 $T$ 和 $q$。但在含时问题里，同一个边界上同时给 $T$ 和 $\partial T/\partial s$ 是热方程的 Cauchy 数据，不能直接施加。可行的做法大概是方向 2/§4 的弛豫：缓慢调整 ghost 压强，使顶部热通量趋向 AWSoM-R 给的 $q$。另外，AWSoM-R 的 $q$ 是否比它的 TR 压强更独立于静态 TR 闭合，还要对照代码确认。

## 7. 建议的下一步（便宜）

在当前 release（hydrostatic reconstruction + single-T）下重跑一次导师方案，和现有 release 4000 s 基线对比。release preset 已是 `ISO_HYDRO_T_DECOUPLE=0`，所以只需要多一个环境变量：

```bash
CHROMO_T_END=4000 ISO_OUTER_P_NEUMANN=2 scripts/run_chromo_realtime.sh \
    outputs/model_column/outer_bc_swap/pneumann2_hserecon_4000s.txt \
    full no-ionization model_column - 20.0 no-cooling
```

判据：`V_top`、顶/底面 `f_total` 是否同号且相等、`p_top` 漂移、柱质量。如果仍然是约 19 倍抑制加底部流出，§4 的计数解释就成立，可以带着这张表和导师讨论特征耦合（§6 方向 2）。如果抑制消失了，说明旧的失败来自非平衡的重构，那就需要更新 swap recap。

## 8. 重跑结果（2026-10-27，当前 release 数值：hydrostatic reconstruction + MC3 修正）

Mac mini，`outputs/model_column/outer_bc_swap_hse_recon/`。设置与 10-13 的 T_wall sweep 重跑（`twall_sweep_hse_recon`，`ISO_HYDRO_T_DECOUPLE=1`）完全相同，只改外边界。导师方案 = `ISO_OUTER_T_HYDRO_WALL=1 ISO_OUTER_P_NEUMANN=2`；对照组 = `ISO_OUTER_T_HYDRO_WALL=1`（只固定 ghost 温度，保留 `kOuterPRef`）。

| | release | 只固定 T（对照） | 导师方案 |
|---|---:|---:|---:|
| 22 kK，4000 s：`V_top`（最后 500 s 平均）[m/s] | 9.639 | 9.649 | **0.734**（速度低 13 倍；顶面质量通量 3500–4000 s 平均低 12.5 倍） |
| 22 kK：顶面 / 底面 `f_total` | +2.70e-10 / +2.71e-10 | 同 release | +2.17e-11 / **−2.23e-10** |
| 22 kK：柱质量，`p_top` 漂移 | 99.995 %，0 | 同 release | 99.00 %，**+5.6 %** |
| 100 kK，1480 s：`V_top`（前 500 s 平均）[m/s] | 1271 | 1271 | **59**（速度低 21 倍；顶面质量通量 980–1480 s 平均低 11.5 倍） |
| 100 kK：顶面 / 底面 `f_total` | +7.50e-09 / +7.25e-09 | +7.49e-09 / +7.85e-09 | +5.9e-10 / **−6.6e-09** |
| 100 kK：柱质量，`p_top` 漂移 | 99.85 %，0 | 99.85 %，0 | **84.3 %，+88 %** |

结论：limiter 修好以后，导师方案的失败方式和旧测量相同（22 kK 的抑制从 19 倍变为 13 倍，量级不变）：蒸发被压低一个量级，质量从底部流出，`p_top` 持续上升。所以失败来自 §3–4 的计数问题，而不是旧的重构。100 kK 的两组 run 为省时间在 1500 s 截断。

100 kK 时导师方案下，TR 底（7.5 kK 等温线）从 2140 km 在 200 s 内沉到 2003 km，1480 s 时到 1956 km；固定压强时同期只下沉 20 km（2120 km）。

## 9. 附：热墙下 TR 的网格锁定脉冲（同日测量）

- 现象：≥ 45 kK 时 TR 速度出现约为局地速度 2 % 的脉冲；每个脉冲都对应 7 kK 等温线跨进下一个格子（100 kK、N500 时周期约 40 s）。7.0–7.5 kK 层宽度远小于格子（100 kK 时约 61 m，格子 276 m）。
- 面导热：TR 底冷的一侧（6.7 kK，约 60 % 电离）的导热以**中性氢**为主（κ_n ≈ 0.6，κ_e ≈ 0.03 W/m/K），热的一侧完全电离（κ_n ≈ 0）。所以底部那个面上：调和平均（默认）0.17，算术平均 0.36，而 `ISO_COND_KIRCHHOFF=1`（只对电子部分做 Kirchhoff，中性部分仍用调和平均）反而**降到** 0.064。它让脉冲变大（rms 31 → 114 m/s，平均顶面通量只变 0.5 %），这与“调和平均卡住热流”的假设方向一致，而不是否定它（之前在对话里把它说成“否定”是错的）。`ISO_COND_ARITH=1`（算术平均，接近总 κ 的 Kirchhoff 值）的结果见下。
- N1000（每个格子减半）：脉冲 rms 32 → 26 m/s（−19 %），频率加倍；平均通量变 0.4 %。按 70 kK 时约 50 m 的稳态宽度、再按质量通量缩放，100 kK 时物理宽度约 25–30 m，所以要分辨它需要约 25 m 的格子，比 N500 细约 10 倍。
- 数据：`outputs/model_column/hot_wall_resolution/tr_pulses_metrics.json`（mac-mini）；图和动画见 `docs/viz_catalog/model_column.md` 的 “Outer pressure closure rerun and hot-wall TR pulses” 一节。
- `ISO_COND_ARITH=1`（算术平均，底部面 κ ≈ 0.36）：脉冲 rms 31 → **5.9 m/s**（小 5 倍），顶格速度摆动 44 → 7.6 m/s，平均顶面通量 1.149 vs 1.141 × 10⁻⁸（0.7 %）。三种规则下脉冲随底部面 κ 单调变化（0.064 → 114，0.17 → 31，0.36 → 6 m/s）：原因是调和平均卡住了流进 TR 底格子的热流。效果远大于网格减半（−20 %）。尚未验证：22 kK 参考算例、高压角落、网格收敛。
