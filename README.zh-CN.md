# AI Biology Lab｜AI 数学生物学与计算生物学实验室

**用 AI 研究生化振荡、分子开关和反应网络，公开别人能够检查的证明、代码和证据。**

[English](README.md) · [首个课题](projects/mixed-phosphorylation/README.md) · [研究路线](docs/RESEARCH_ROADMAP.md) · [媒体资料](docs/MEDIA_BRIEF.zh-CN.md)

由 [Zijian Zeng / ArtificialZeng](https://github.com/ArtificialZeng) 发起。我们希望做出值得科技媒体报道的科学成果，以明确的问题、可复现的结果和外部验证来支撑传播。

## 项目目录

| 项目 | 当前交付 | 证据范围 |
| --- | --- | --- |
| [混合机制生化振荡](projects/mixed-phosphorylation/README.md) | 精确失稳界、吸引周期轨道及非零周期修正 | Linux区间证书与独立精确复算；优先权和外审待确认 |
| [Paper-Code Auditor](projects/paper-code-auditor/README.md) | 复现入口缺陷、历史源码恢复与补丁 | 离线CPU检查通过；未完成训练结果复现或作者确认 |
| [ORD Data Auditor](projects/ord-data-auditor/README.md) | 标识符校验器、固定样本审计与校验位补丁 | 已复跑42数据集／20,720反应；不自动改写化合物身份 |
| [Protein Metric Auditor](projects/protein-metric-auditor/README.md) | NDCG行顺序敏感性全量审计 | 217实验／95配置；副指标排名变动，未推翻主榜单 |
| [MEMOTE方向审计](projects/memote-direction-auditor/README.md) | 等价反应写法反例与两行候选修复 | 真实完整ATP检测入口验证通过；待维护者确认 |
| [边界中继与孤立平衡分支](projects/relay-isola-auditor/README.md) | 给定单参数族中的精确质量作用反例 | 新上下文AI范围复核与两套精确计算；人类外审和优先权未确认 |
| [区室分配 beta 区间](projects/compartment-beta-gap/README.md) | 7/8 < beta < 1 正概率爆炸定理 | 独立 AI 数学复核及有限精确检查通过；待专家与优先性确认；未声称几乎必然爆炸 |
| [免疫衰减下的双毒株共存](projects/waning-coexistence-auditor/README.md) | 完整模型的共存存在性与条件逆命题，允许任意非负免疫衰减 | 新上下文AI数学审查与两套精确复算；广义猜想和根计数仍未解决 |

上游核验申请（按顺序提交，尚待维护者审核）：[Paper-Code PR #1](https://github.com/ziul-bio/DMS-SiteEffect-PLM/pull/1) · [ORD PR #1165](https://github.com/open-reaction-database/ord-schema/pull/1165)。

## 新结果：稳定振荡的周期仍有可量化偏差

固定原文速率、P=8、S=40，已证明：每个充分大的有限激酶总量K对应的守恒类内都有吸引周期轨道。其局部分支在模型时间下满足

$$T(K)=T_\infty+\frac{T_1}{K}+O(K^{-2}),\qquad T_1\in[7.65771717,7.68396573].$$

T1严格为正，说明周期从上方趋近极限，并非完全不变。结论经过Linux区间积分、新上下文数学审查及独立有理数复算；没有给出明确K0或余项常数。[完整证明与复现](projects/mixed-phosphorylation/validated-cycle/README.md) · [中文研究简报](projects/mixed-phosphorylation/validated-cycle/docs/BRIEF.zh-CN.md)。质性鲁棒性已有直接相关会议研究，优先权尚未确认；完整原猜想仍未解决。

## 首个结果

在 Conradi、Mincheva、Shiu 的混合机制双磷酸化模型中，固定原文 Table 1 速率、总磷酸酶 8、总底物 40，我们给出精确证明：

> 对每个有限的总激酶浓度 K ≥ 4660166/209569（约 22.2369），唯一正平衡点在六维守恒类中恰有两个不稳定、四个稳定特征值。因此，对所有 K ≥ 23 都成立。

证明包含 **125 个严格为正的有理数 Bernstein 系数**，覆盖整个无界参数区间。23 是充分界，尚未声称是最优阈值。

- [五页结果稿](projects/mixed-phosphorylation/manuscript/main.pdf)
- [完整证明](projects/mixed-phosphorylation/evidence/finite_kinase_proof.md)
- [机器可检查证书](projects/mixed-phosphorylation/evidence/finite_kinase_sign_certificate.json)

这份早期精确证书处理原文猜想6.2的第一项；单凭平衡点不稳定不能推出周期轨道。上面的新工作另行给出了周期轨道及周期修正证明，振幅的严格修正与完整猜想仍未解决。已有直接相关的 2023/2025 会议工作，新颖性与优先权待原作者确认。当前没有外部人类同行审稿或细胞实验。

![浓度变化的数值结果](projects/mixed-phosphorylation/evidence/concentration_effects.png)

图为浮点数值观察，不是周期轨道存在性的严格证明。

## 一键复核

Python 3.12 / 3.13；精确验证只依赖 SymPy，无需 GPU 或模型 API。

```sh
git clone https://github.com/ArtificialZeng/AI-Biology-Lab.git
cd AI-Biology-Lab
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-exact.txt
python scripts/reproduce.py verify
```

验证在临时副本中执行，不改动归档证据。详细命令见[复现说明](docs/REPRODUCIBILITY.md)。

## 接下来真正要做出的成果

1. 算出明确的有限K阈值及可用的余项常数。
2. 严格认证振幅修正，并澄清标量投影的不变性。
3. 扩展经过认证的参数区域。
4. 与相关研究核对优先权，取得外部人类复核。

这些是待完成目标。[研究路线](docs/RESEARCH_ROADMAP.md)列出验收条件；[外部验证说明](docs/VALIDATION.md)列出专家可以复核的问题。AI 参与推导和审查的范围见 [AI 使用说明](docs/AI_USE.md)。

## 2026-10-26：开放问题研究

新增两个精确研究包：[分舱模型猜想端点反例](projects/directed-cycle-identifiability/README.md)与[双磷酸化循环八参数对称子类的稳定性分类](projects/dual-futile-symmetry/README.md)。前者针对猜想的边界条件，原文主定理不受影响；后者留下完整十六参数猜想未解。优先权和原作者确认仍待核验。[筛选与范围记录](docs/OPEN_RESEARCH_2026-10-26.md)。
