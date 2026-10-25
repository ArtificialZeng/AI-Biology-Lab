# Paper-Code Auditor：可复跑的复现入口审计与修复

本项目检查一项已发表蛋白质研究的公开复现路线，交付固定版本证据、离线CPU检查器，以及可应用的修复与历史源码恢复补丁。这是软件复现工程成果，尚无证据表明论文的科学结论因此改变。

审计对象为 [DMS-SiteEffect-PLM](https://github.com/ziul-bio/DMS-SiteEffect-PLM)，固定提交 `e168a2c825be0c3e28c300b8f6277dda1d999c4d`，检查日期2026年10月26日。

## 问题与交付

| 问题 | 验证结果 | 修复范围 |
| --- | --- | --- |
| cellular入口循环了未定义的数组 | 干净Bash控制流正常退出、调度0次；修正后计数198次 | 修正数组名；198是33数据集×2策略×3重复的调度计数，没有启动训练 |
| 两个wrapper调用已删除的partial trainer | 固定版本完整文件清单确认缺失；公开Git历史可恢复它及配套模型模块 | 从删除前版本成对恢复原源码，保留partial模型实现，使用模块启动方式 |
| 当前LoRA的site分支传入错误参数 | 实际setup函数稳定触发TypeError | 改为函数原有的test_size参数；三个种子的验证行数为300、300、299，位点不重叠 |
| 历史运行配置的GPU编号不一致 | wrapper只暴露一张GPU，trainer却选择逻辑编号1 | 改为逻辑编号0，保留用户对可见GPU的设置；未进行GPU运行验证 |

原partial trainer与配套模型来自公开提交 `e0c2d6ae06903cb356abb7e1b579058115e5623a`。恢复方案没有拿LoRA替代partial训练方法。LoRA参数修复是一项独立局部修复；原partial文件自身的位点切分调用本来正确。

## 如何核验

仓库根目录执行：

```sh
python -m pip install -r projects/paper-code-auditor/requirements.txt
python projects/paper-code-auditor/audit.py
```

检查器验证源文件哈希、控制流、实际数据处理AST、补丁可应用性、恢复来源、模块调用与参数兼容性。数据夹具只保留真实实验表的行标识和突变标识，不含序列或实验表型数值。没有加载模型、下载权重或启动训练。

最初怀疑Lasso与LoRA使用不公平的不同划分，但三组种子测试否定了这一怀疑。本项目保留这个负对照，不宣称发现了不公平比较。

## 结论边界

这些证据支持固定公开版本存在复现入口问题，以及补丁在检查范围内有效。公开历史恢复不等于证明原作者实际使用了这些精确脚本。兼容训练依赖、真实GPU训练和论文指标复现仍未完成；替代的作者自定义checkpoint不在此次恢复范围内。

已在实验室负责人授权后于2026年10月26日提交[上游PR #1](https://github.com/ziul-bio/DMS-SiteEffect-PLM/pull/1)，等待维护者核验。目前没有原作者背书、上游接受或独立外部人类复核。公开记录中未找到相关报告，不构成穷尽优先权证明。不据此宣称新的生物学发现、模型性能反转或完整复现了论文。

[完整说明](../README.md) · [恢复补丁](../patches/reproduction-fixes.patch) · [来源与哈希](../evidence/restoration_provenance.json) · [记录的验证结果](../evidence/verification.json) · [供上游审核的材料](UPSTREAM_REVIEW.md)
