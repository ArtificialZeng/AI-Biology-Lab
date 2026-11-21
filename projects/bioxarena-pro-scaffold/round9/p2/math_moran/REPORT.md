# M1 任务报告：星形图中性 Moran 一致时间精确公式证明

## 结论摘要

Brewster, Svoboda, Roscow, Chatterjee, Tkadlec & Nowak（*PNAS Nexus* 2025, `doi:10.1093/pnasnexus/pgaf252`; `arXiv:2503.09841v2`）在正文公式 (9) 中猜想：在 $N = n+1$ 个顶点（$1$ 个中心、$n \ge 1$ 片叶子、无自环）的星形图 $S_n$ 上，从全员互异类型出发，在中性 birth-death（`bd`）更新下的期望一致时间（吸收时间 / 总溯祖时间）为
$$T_N^{\mathrm{bd}} = n^3 - n^2 + n \cdot H_n.$$

**结论：该猜想为真。** 我们通过中性对偶（后向合并随机游走）将指数规模状态空间精确强合并压缩至 $2n - 1$ 个瞬态状态 $(k, c)$，并给出对任意叶子对称中性更新规则 $(p_{\mathrm{in}}, p_{\mathrm{out}})$（满足 $p_{\mathrm{in}} > 0,\, p_{\mathrm{out}} \ge 0,\, n(p_{\mathrm{in}}+p_{\mathrm{out}})=1$）的统一精确主定理：
$$T_N(p_{\mathrm{in}}, p_{\mathrm{out}}) = n H_n + n(n-1)\frac{p_{\mathrm{out}}}{p_{\mathrm{in}}}.$$
由此直接推得：
1. **Birth-death（`bd`）**：$T_N^{\mathrm{bd}} = n^3 - n^2 + n H_n$（证明原猜想公式 (9)）；
2. **Death-birth（`db`）**：$T_N^{\mathrm{db}} = n H_n + n - 1 = (N-1)H_{N-1} + N - 2$；
3. **$\lambda$-混合更新（以概率 $\lambda$ 用 `bd`、$1-\lambda$ 用 `db`）**：$T_N^{(\lambda)} = n H_n + n(n-1)\frac{n\lambda + 1 - \lambda}{\lambda + n(1-\lambda)}$，特别地 $\lambda=1/2$ 时 $T_N^{(1/2)} = n^2 - n + n H_n$。

完整数学证明见 [proof.md](proof.md)，项目概览与复现指令见 [README.md](README.md)，文献来源记录见 [sources.md](sources.md)。

---

## 任务 1：原文定义核对

通过 `arXiv:2503.09841v2` 源码（`main.tex`、`si/proofs/star.tex`、`si/proofs/well_mixed.tex`）及 *PNAS Nexus* 正式全文（Europe PMC `PMC12363668`）与 Crossref 元数据核对，确认原文精确设定如下：

| 设定项 | 原文出处 | 核对结果 |
|---|---|---|
| 顶点与 $n, N$ 关系 | `main.tex` 第 536–538 行 | 中心顶点 $c$ 度为 $n$，外围叶子 $v_1,\dots,v_n$ 度为 $1$，总顶点数 $N = n + 1$（$n \ge 1$）。 |
| 自环与重边 | `SI Section 1` | 简单连通无向图，**无自环**、无重边（仅完全图研究特例时另加自环）。 |
| 初始条件 | `SI Section 1`、Methods Eq. (15) | $|\tau| = |V| = N$ 且 $X_0(V) = \tau$，即初始时全部 $N$ 个顶点类型互异。 |
| `bd` 更新法则 | `main.tex` Eq. (1)、Methods | 均匀选母体 $B \in V$（概率 $1/N$），再均匀选其邻点 $D \in \Gamma(B)$（概率 $1/\deg(B)$）被替代；中性进化全员适应度相等。 |
| 时间单位 | `SI Section 1`、Methods Eq. (15) | 离散更新步数（每步执行一次出生-死亡替代），直至全图仅余 $1$ 种类型为止的期望步数 $T_N$。 |
| 公式 (9) 原文 | `main.tex` 第 556–560 行 | `T_N = n^3 - n^2 + n \cdot H_n`，其中 $H_n = \sum_{k=1}^n 1/k$；原文 SI Section 5 仅证明了 $\Theta(N^3)$（`bd`）与 $\Theta(N \log N)$（`db`）阶数。 |

---

## 任务 2：文献查新核验

在 2026-11-21 通过 OpenAlex、Crossref、Europe PMC 与 arXiv API 完成检索与全文核查（检索范围局限性见 [sources.md](sources.md)）：

1. **原论文核验**：OpenAlex `W4413186912`、Crossref `10.1093/pnasnexus/pgaf252`（2025-07-31 发表，*PNAS Nexus* 4(8):pgaf252）、`arXiv:2503.09841v2`（2025-09-16 更新）、Europe PMC `PMC12363668` / `PMID:40838023` 一致。
2. **OpenAlex 全部 6 篇施引文献逐一核查**：
   - `W7154604876`（Svoboda 等，*Nat. Commun.* 2026，`10.1038/s41467-026-71777-2`）：研究空间适应度梯度对双类型固定概率的影响；
   - `W7163063073`（*J. Theor. Biol.* 2026，`10.1016/j.jtbi.2026.112517`）：研究有亲代无关突变下的稳态等位基因分布；
   - `W7115587857`（Brewster, Huang, Mitzenmacher & Nowak，*PLOS Comput. Biol.* 2026，`10.1371/journal.pcbi.1014829` / `PMC13645234` / `arXiv:2511.18252v3`）：全文核查确认其给出的是混合更新下星图/圈图双类型**固定概率**精确公式及一般图固定时间界，未解 $N$-类型一致时间公式 (9)；
   - `W4414010756`（Sharma 等，bioRxiv 2025，`10.1101/2025.09.04.674237`）：研究空间结构放大突变负荷；
   - `W7162518942`（Brewster, Lippner, Tkadlec & Nowak，2026，`arXiv:2605.26411v3`）：全文核查确认研究单突变体固定位置分布；
   - `W7202117459`（Pecho 等，2026，`arXiv:2608.06552v1`）：全文核查确认研究 replacer 表型的演化稳定性。
3. **Europe PMC 与 arXiv 补充检索**：Europe PMC `PMID:40838023` 显示 $2$ 条施引（均已含于 OpenAlex 上述 6 篇中）；arXiv 定向标题/摘要与作者检索（`Brewster`、`Tkadlec`、`Svoboda`）截至 2026-11 未检索到证明公式 (9) 的预印本。

---

## 任务 3 & 任务 5：证明思路与推广结果（含 `db`）

1. **中性对偶（[proof.md](proof.md) Lemma 1）**：单步替代 $(B \to D)$ 定义随机亲代映射 $\phi_{(B,D)} : V \to V$。因初态 $X_0$ 为单射且 $(\phi_t)_{t \ge 0}$ 独立同分布，时间反转 $(\phi_0, \dots, \phi_{t-1}) \stackrel{d}{=} (\phi_{t-1}, \dots, \phi_0)$ 给出前向一致时间 $\tau_{\mathrm{cons}}$ 与后向合并随机游走 $S_{m+1} = \phi_m(S_m)$（$S_0 = V$）的吸收时间 $\tau_{\mathrm{coal}} = \inf\{m \ge 0 : |S_m| = 1\}$ **同分布**，故 $T_N = \mathbb{E}[\tau_{\mathrm{coal}}]$。
2. **精确状态压缩（[proof.md](proof.md) Lemma 2）**：投影 $\pi(S) = (k, c)$（$k = |S \cap [n]|$ 为叶上谱系数，$c = \mathbf{1}[0 \in S]$ 为中心是否含谱系）构成强可合并马尔可夫链，仅含 $2n - 1$ 个瞬态与 $2$ 个吸收态 $(0, 1), (1, 0)$。因 $p_{\mathrm{in}} > 0$，吸收期望有限且线性方程组唯一可解。
3. **递推消元与裂项求和（[proof.md](proof.md) Theorem 1）**：
   - 由 $(k, 0)$（$k \ge 2$）得 $t(k, 0) = \frac{1}{k p_{\mathrm{in}}} + t(k-1, 1)$；
   - 边界 $(1, 1)$ 给出 $t(1, 1) = n + \frac{n(n-1)}{2}\frac{p_{\mathrm{out}}}{p_{\mathrm{in}}}$；
   - 代入 $(k, 1)$（$2 \le k \le n$）方程后，利用 $p_{\mathrm{in}} + p_{\mathrm{out}} = \frac{1}{n}$，$(k,1)$ 与 $(k-1,1)$ 的系数精确化简为相同的 $\frac{k}{n}$，得到精确差分
     $$t(k, 1) - t(k-1, 1) = \frac{n}{k} + \frac{n(n+1)r}{k(k+1)} \qquad \left(r := \frac{p_{\mathrm{out}}}{p_{\mathrm{in}}}\right).$$
   - 对所有 $1 \le k \le n$ 与 $2 \le k \le n$ 分别求得封闭通解：
     $$t(k, 1) = n H_k + \frac{n r (n k - 1)}{k + 1}, \qquad t(k, 0) = n H_{k-1} + \frac{n + n^2 r (k - 1)}{k},$$
     令 $k = n$ 即得 $T_N = n H_n + n(n - 1) r$。
4. **推论与结构附注**：
   - **猜想无反例**，对所有 $n \ge 1$ 成立：$T_N^{\mathrm{bd}} = n^3 - n^2 + n H_n$；
   - **`db` 精确公式**：$T_N^{\mathrm{db}} = n H_n + n - 1$（仅比一般图优惠券收集下界 $n H_n$ 多出 $n - 1$ 步）；
   - **$\lambda$-混合精确公式**：$T_N^{(\lambda)} = n H_n + n(n-1)\frac{n\lambda + 1 - \lambda}{\lambda + n(1-\lambda)}$；
   - **附加恒等式**：在 `bd` 更新且 $n \ge 2$ 时，$t_{\mathrm{bd}}(n, 0) = t_{\mathrm{bd}}(n, 1) = n^3 - n^2 + n H_n$（所有叶子溯祖时间等于全体顶点溯祖时间）。

---

## 任务 4：双独立检查器与故意改坏对照

全部验证脚本位于 `checks/`，不依赖缓存或中间状态文件，可在 `python3 -B checks/verify.py` 与 `python3 -B -O checks/verify.py` 下零差异复现；以下所有数值均直接出自 [checks/verification_summary.json](checks/verification_summary.json)。

### 检查器 1：精确有理数马尔可夫链求解器（[checks/exact_markov_checker.py](checks/exact_markov_checker.py)）

采用纯标准库 `fractions.Fraction` 高斯消元精确求解 $(I - P_{\mathrm{trans}})\tau = \mathbf{1}$，转移矩阵严格从更新规则定义独立构造，并检验每行和为 `1` 与吸收可达性：

1. **完全无约化 $N^N$ 顶点标号前向链**（直接实现原文 Methods Eq. (15)，$n=1,2,3$ 共 $9$ 组：$\lambda \in \{1, 0, 1/2\}$）：
   - `bd`：$n=1$（1 瞬态）得 `1`；$n=2$（7 瞬态）得 `7`；$n=3$（85 瞬态）得 `47/2`；
   - `db`：$n=1,2,3$ 分别得 `1`、`4`、`15/2`；$\lambda=1/2$ 分别得 `1`、`5`、`23/2`。全部零误差一致。
2. **仅做类型同构标号约化的前向链**（仍逐一保留全部 $n$ 个叶子顶点位置，不做叶子对称压缩，$n=1\dots 5$，$\lambda \in \{1, 0, 1/3, 1/2\}$ 共 $20$ 组）：
   - `bd`：$n=4$（37 瞬态）得 `169/3`，$n=5$（151 瞬态）得 `1337/12`；
   - `db`：$n=4, 5$ 分别得 `34/3`、`185/12`；$\lambda=1/3$ 得 `49/3`、`3187/132`。全部零误差一致。
3. **完整 $2^{n+1}-1$ 子集后向合并链 + 强可合并性检验**（$n=1\dots 6$，$\lambda \in \{1, 0, 1/2\}$ 共 $18$ 组）：
   - 对每个子集 $S \subseteq V$ 校验投影 $\pi(S)=(k,c)$ 的转移概率精确匹配强可合并条件，并验证全部瞬态子集（$n=6$ 时 $120$ 个瞬态）的精确吸收时间与 $t(k, c)$ 通解一致（`bd` $n=6$ 得 `1947/10`，`db` $n=6$ 得 `197/10`）。
4. **压缩 $(k, c)$ 链精确求解**（$n=1\dots 20$，$\lambda \in \{1, 0, 1/4, 1/2, 3/4\}$ 共 $100$ 组）：
   - 验证全部瞬态 $t(k,1), t(k,0)$ 与通解完全相等（如 $n=20$ 时 `bd` 为 `29766233275/3879876`，`db` 为 `352893319/3879876`），并验证 $n \ge 2$ 时 `t_bd(n,0) == t_bd(n,1)`。

### 检查器 2：独立 SymPy 符号推演校验 + 字面前向蒙特卡洛（[checks/independent_mc_and_symbolic.py](checks/independent_mc_and_symbolic.py)）

1. **SymPy 1.13.3 符号检验**：把 $n, k, r, \lambda$ 设为符号、$H(\cdot)$ 设为抽象调和函数（仅代入 $H(1)=1$ 与 $H(k)=H(k-1)+1/k$），检验 $(k,0)$、内点 $(k,1)$、边界 $(n,1)$、边界 $(1,1)$、基础情形 $n=1$、裂项差分、以及 `bd`/`db`/$\lambda$-混合特化与 $t_{\mathrm{bd}}(n,0)-t_{\mathrm{bd}}(n,1)$ 共 $8$ 项符号残差全部严格化简为 `0`。
2. **字面前向蒙特卡洛模拟**（完全独立于对偶与 $(k,c)$ 递推，逐脉冲执行单步出生-死亡替代）：共 $17$ 组设定（`bd` $n \in \{2,3,4,5,8,12\}$，`db` $n \in \{2,3,4,5,8,12,20,30\}$，$\lambda=1/2$ $n \in \{3,5,10\}$，每组 $3{,}000\text{–}15{,}000$ 次独立试验）。全部 $17$ 组均通过 $|z| < 4.0$ 检验，最大 $|z| = 2.316$（出现在 `db`, $n=3$：模拟均值 `7.4175 ± 0.0356` vs 精确值 `7.5`）。代表性结果：
   - `bd`, $n=8$（5,000 次）：`472.304 ± 4.064` vs 精确 `469.742857`（$z = 0.630$）；
   - `bd`, $n=12$（3,000 次）：`1602.894 ± 16.7187` vs 精确 `1621.238528`（$z = -1.097$）；
   - `db`, $n=30$（10,000 次）：`148.6548 ± 0.5646` vs 精确 `148.849614`（$z = -0.345$）。

### 故意改坏对照（[checks/negative_controls.py](checks/negative_controls.py)）

共构造 $8$ 项故意改坏的对抗/错误对照，全部被对应检查器成功捕获并拒绝：

| 对照名称 | 改坏内容 | 捕获检查器 | 观测误差 / 拒绝指标 |
|---|---|---|---|
| `off_by_one_N_instead_of_n` | 误将 $n$ 写成 $N=n+1$ | 精确链 + SymPy + MC | $n=3$ 精确差 `197/6`；$n=5$ MC $z = -83.08$ |
| `truncated_harmonic_H_{n-1}` | 调和级数漏最后一项 $H_{n-1}$ | 精确链 + SymPy + MC | $n=3$ 精确差 `-1`；$n=3$ MC $z = 7.01$ |
| `plus_one_additive_corruption` | 正确公式加常数 `+1` | 精确链 + SymPy + MC | $n=3$ 精确差 `1`；$n=2$ MC $z = -27.42$ |
| `log_instead_of_harmonic` | 误以 $\ln n$ 代替 $H_n$ | 精确链 + MC | $n=3$ 差 `-2.2042`；$n=3$ MC $z = 18.23$ |
| `self_loop_center_n3` | $n=3$ 中心顶点加自环 | 无约化精确前向链 | 期望时间变为 `472/15`（相对 `47/2` 差 `+239/30`） |
| `self_loop_all_n3` | $n=3$ 全顶点加自环 | 无约化精确前向链 | 期望时间变为 `280/9`（相对 `47/2` 差 `+137/18`） |
| `corrupted_kc_kernel_n4_eps_1_100` | 保行和前提下扰动 $(4,0)$ 转移概率 $1/100$ | 精确 $(k,c)$ 求解器 | `T_5` 变为 `3533/63`（偏离 `169/3` 差 `-16/63`） |
| `swapped_bd_vs_db_n4` | 用 `db` 模拟校验 `bd` 公式 | 精确链 + MC | $n=4$ MC 均值 `11.2146` vs `169/3`，$z = -727.34$ |

*(注：在初始设计 `log_instead_of_harmonic` 对照时曾试以 `n=8` 跑 `bd` 蒙特卡洛，发现仅得 $z=2.18$——这是因为 $n(H_n-\ln n)=O(n)$ 相对于 `bd` 吸收时间 $\Theta(n^3)$ 的标准差过小；改用 `n=3` 增大信噪比后得到 $z=18.23$ 稳定拒斥。)*

---

## 文件清单与复现

- [README.md](README.md)：项目摘要、定理与复现入口
- [proof.md](proof.md)：完整数学证明（中性对偶、强可合并压缩、递推精确解与 `bd`/`db`/$\lambda$-混合推论）
- [sources.md](sources.md)：原文定义核对与有界文献查新记录
- [REPORT.md](REPORT.md)：本任务报告
- [checks/common.py](checks/common.py)：精确分数算术、线性方程组求解与公式定义
- [checks/exact_markov_checker.py](checks/exact_markov_checker.py)：检查器 1（4 类精确有理数马尔可夫链核验）
- [checks/independent_mc_and_symbolic.py](checks/independent_mc_and_symbolic.py)：检查器 2（SymPy 符号检验 + 前向蒙特卡洛模拟）
- [checks/negative_controls.py](checks/negative_controls.py)：$8$ 项故意改坏对照
- [checks/verify.py](checks/verify.py)：主复现入口（普通模式与 `-O` 优化模式均保留全部检查）
- [checks/verification_summary.json](checks/verification_summary.json)：确定性校验结果档案
