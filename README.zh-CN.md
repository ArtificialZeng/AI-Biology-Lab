# AI Biology Lab｜AI 数学生物学与计算生物学实验室

**用 AI 研究生化振荡、分子开关和反应网络，公开别人能够检查的证明、代码和证据。**

[English](README.md) · [首个课题](projects/mixed-phosphorylation/README.md) · [研究路线](docs/RESEARCH_ROADMAP.md) · [媒体资料](docs/MEDIA_BRIEF.zh-CN.md)

由 [Zijian Zeng / ArtificialZeng](https://github.com/ArtificialZeng) 发起。我们希望做出值得科技媒体报道的科学成果，以明确的问题、可复现的结果和外部验证来支撑传播。

## 项目目录

| 项目 | 当前交付 | 证据范围 |
| --- | --- | --- |
| [混合机制生化振荡](projects/mixed-phosphorylation/README.md) | 有限激酶失稳界与精确证书 | 数学部分结果；新颖性与外审待确认 |
| [Paper-Code Auditor](projects/paper-code-auditor/README.md) | 复现入口缺陷、历史源码恢复与补丁 | 离线CPU检查通过；未完成训练结果复现或作者确认 |
| [ORD Data Auditor](projects/ord-data-auditor/README.md) | 标识符校验器、固定样本审计与校验位补丁 | 已复跑42数据集／20,720反应；不自动改写化合物身份 |

上游核验申请（按顺序提交，尚待维护者审核）：[Paper-Code PR #1](https://github.com/ziul-bio/DMS-SiteEffect-PLM/pull/1) · [ORD PR #1165](https://github.com/open-reaction-database/ord-schema/pull/1165)。

## 首个结果

在 Conradi、Mincheva、Shiu 的混合机制双磷酸化模型中，固定原文 Table 1 速率、总磷酸酶 8、总底物 40，我们给出精确证明：

> 对每个有限的总激酶浓度 K ≥ 4660166/209569（约 22.2369），唯一正平衡点在六维守恒类中恰有两个不稳定、四个稳定特征值。因此，对所有 K ≥ 23 都成立。

证明包含 **125 个严格为正的有理数 Bernstein 系数**，覆盖整个无界参数区间。23 是充分界，尚未声称是最优阈值。

- [五页结果稿](projects/mixed-phosphorylation/manuscript/main.pdf)
- [完整证明](projects/mixed-phosphorylation/evidence/finite_kinase_proof.md)
- [机器可检查证书](projects/mixed-phosphorylation/evidence/finite_kinase_sign_certificate.json)

**当前只处理原文猜想 6.2 的第一项。** 平衡点不稳定不能直接推出存在稳定周期轨道。周期、振幅与浓度变化目前有数值证据；完整猜想仍未解决。已有直接相关的 2023/2025 会议工作，新颖性与优先权待原作者确认。当前没有外部人类同行审稿或细胞实验。

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

1. 严格认证吸引周期轨道及其稳定性。
2. 证明所有充分大有限 K 下周期振荡持续存在。
3. 给出周期、振幅对酶浓度变化的严格误差界。
4. 解释哪些参数区域支持这种鲁棒性，在哪里失效。

这些是待完成目标。[研究路线](docs/RESEARCH_ROADMAP.md)列出验收条件；[外部验证说明](docs/VALIDATION.md)列出专家可以复核的问题。AI 参与推导和审查的范围见 [AI 使用说明](docs/AI_USE.md)。
