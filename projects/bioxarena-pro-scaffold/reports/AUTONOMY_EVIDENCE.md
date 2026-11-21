# Autonomy & Human-Intervention Evidence Pack (`R8-9`)

- **Compiled at**: 2026-11-21 16:15+08:00
- **Scope**: Complete executor-side audit trail from initial project commit `a09a4f6` (`2026-11-20 15:55:57+08:00`) through `3c313be` (`2026-11-21 16:11:58+08:00`), covering **36 git commits** across **8 rounds (~24.3 hours elapsed)**.
- **Transparency Rule**: No beautification, no overclaiming of "full autonomy". Every direct human message sent to Antigravity is extracted verbatim (`<= 30` chars) from the immutable conversation log (`transcript.jsonl`). Any action occurring on the Claude Master Controller side or outside Antigravity's observable log is explicitly marked **`不确定`** (uncertain / unverified from executor logs). Human messages sent directly to Claude Master Controller are tallied separately by Claude.

---

## 1. Executive Summary of Observed Autonomy vs. Human Oversight

| Metric | Executor-Side Verified Count / Value | Notes |
| :--- | :--- | :--- |
| Elapsed wall-clock span (`a09a4f6` → `3c313be`) | **24 h 16 m** (`2026-11-20 15:55` → `2026-11-21 16:11`) | Ongoing until `R8-1` & `R8-6` DRF complete |
| Total git commits on `main` | **36** (`a09a4f6` .. `3c313be`) | 100% authored/committed by Antigravity |
| Direct human messages to Antigravity | **15** (all between `2026-11-20 15:44` and `22:36`) | **0 direct human messages** to Antigravity during Rounds 5–8 (`> 17.5 h`) |
| System-relayed Claude Master Controller messages to Antigravity | **26** (`2026-11-21 02:33` → `16:04`) | Automated relay loop replacing manual copy-paste after Round 4 |
| Direct human messages to Claude Master Controller | **不确定** (counted separately by Claude) | Included human authorization for public repo, PDF v1, and main-line pivot |
| Audit interventions logged (`AUDIT_LEDGER.csv`) | **34** (13 retracted claims, 7 downgraded claims) | Includes 2 Claude self-corrections and Round 7 missing-prior-art oversight |

> **Honest Assessment**: This workflow is **NOT** "zero-human, 100% autonomous discovery." It is a **human-directed, dual-AI cross-vendor execution-and-audit loop** in which the human (Zijian Zeng) set the overarching goal, made governance/release/team decisions, and initially bridged the two AIs by copy-pasting Claude's Round 1–3 reviews before switching to an automated Claude-to-Antigravity relay in Rounds 6–8. Within that human-governed frame, all code implementation, statistical pipelines, pre-registration locks, negative dataset gates, prior-art searches, figure generation, and manuscript compilation across 36 commits were executed autonomously by the two AI systems.

---

## 2. Round-by-Round Execution Timeline

*(Machine-readable table archived in `reports/autonomy_timeline.csv`.)*

| Round | Start (`UTC+8`) | End (`UTC+8`) | Commits | Commit Range | Human Msgs → Executor | Claude Relays → Executor | Main Outputs & Audit Result |
| :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| **Round 1** | 2026-11-20 15:55:57 | 2026-11-20 16:54:15 | 10 | `a09a4f6`..`a8bde5d` | 10 | 0 | 6-task baseline scaffold + initial overclaimed paper/press drafts → **Rejected by Claude (L1)** |
| **Round 2** | 2026-11-20 17:57:47 | 2026-11-20 18:24:21 | 2 | `ca93cb2`..`fa7165f` | 2 | 0 | Pivot to `xatlas-orion-perturb`, single-seed `0.0804` → **Rejected by Claude (lucky seed 42, sync gap)** |
| **Round 3** | 2026-11-20 20:45:09 | 2026-11-20 21:26:37 | 2 | `710e250`..`59a7753` | 1 | 0 | 10-split strict bench + official grader + CSV/log cleanup → **40/40 cells reproduced bit-for-bit (L1)** |
| **Round 4** | 2026-11-20 22:26:17 | 2026-11-20 22:26:17 | 1 | `c479bc6` | 2 | 0 | `PanPerturb-Bridge (S1–S7)` (`v5 = 0.1099 ± 0.0087`) → **L2 confirmed; t-test & forecast rejected** |
| **Round 5** | 2026-11-20 23:47:03 | 2026-11-20 23:47:03 | 1 | `a9ab830` | 0 | 0 | 20-split Orion + zero-shot Replogle 2022 + Reactome + 3 figures → **Metric-inversion signal observed** |
| **Round 6** | 2026-11-21 02:34:57 | 2026-11-21 05:15:13 | 5 | `2a8b2e0`..`3609811` | 0 | 18 | Frozen solver (`3609811`), Codabench `0.0998`, `T2` 5-metric/row-CI, `T4` LOMO, `T5` Reactome, `T6` cis audit → **Passed** |
| **Round 7** | 2026-11-21 05:19:22 | 2026-11-21 12:40:06 | 6 | `c5294c0`..`cb3935f` | 0 | 4 | Pre-registered unseen Nadig 2025 H1–H3 confirmed, pseudobulk GEARS & scGPT probe, S5 ±50% sensitivity, Systema proof, `PREPRINT_DRAFT.md`, 4 figures → **Passed bit-for-bit; missing Miller 2026 prior-art check caught later** |
| **Round 8** *(active)* | 2026-11-21 12:45:17 | 2026-11-21 16:11:58* | 9* | `9c4aa6c`..`3c313be`* | 0 | 4 | `R8-0` fixes, `R8-2` PC1 decomposition (3/3 confirmed), `R8-3` OrthoFusion + 4th-dataset hard stop (`PENDING`), `R8-6` Miller citation, `R8-7` 34-entry `AUDIT_LEDGER`, `R8-8` `PRIOR_ART_LOG`, `R8-10` Preprint PDF v1; `R8-1` & `R8-6` DRF running |

*\*Round 8 timeline reflects commits through `3c313be`; final Round 8 commits (`R8-9`, `R8-1`, `R8-6` DRF, PDF v2, `ROUND8_DONE`) will extend this row upon closeout.*

---

## 3. Complete Human Intervention Ledger (Executor Side)

### 3.1 Every Direct Human Message to Antigravity (15 Messages Total, Exhaustive from `transcript.jsonl`)

| # | Timestamp (`UTC+8`) | Step | Length (chars) | Verbatim Excerpt (`<= 30` chars) | Functional Role |
| :---: | :--- | :---: | :---: | :--- | :--- |
| 01 | 2026-11-20 15:44:00 | 0 | 3,505 | `这个是muse想自主进行研发，关于主题和信心的评估： # A` | Passed candidate research screening doc from Muse |
| 02 | 2026-11-20 15:45:09 | 2 | 15 | `你选一个信心最高的，全力进攻吧` | Delegated topic choice and authorized launch |
| 03 | 2026-11-20 15:47:02 | 169 | 23 | `/goal 全自动完成能值得报道或者顶刊的研究` | Set high-level ambition bar (`/goal` mode) |
| 04 | 2026-11-20 15:53:27 | 132 | 133 | `你有什么新进展就同步到https://github.com/` | Designated GitHub repo & leadership over Muse |
| 05 | 2026-11-20 15:54:36 | 147 | 38 | `并且你告诉我，你的同步库在哪， 也就是分发任务的目录在哪， ` | Requested sync directory path for Muse handoff |
| 06 | 2026-11-20 16:01:01 | 202 | 591 | `这是muse的回复： 收到，已按 Gemini 的任务卡转方` | Relayed Muse's reply; handed off coordination |
| 07 | 2026-11-20 16:31:10 | 375 | 49 | `你怎么监督muse的，以及muse迄今为止完成了什么， 对你` | Audited Muse's utility; offered option to drop Muse |
| 08 | 2026-11-20 16:34:57 | 412 | 52 | `行， 放弃muse，全自动自主科研。你可以控制codex c` | Approved dropping Muse & authorized Codex CLI |
| 09 | 2026-11-20 17:13:09 | 590 | 37 | `现在研究到什么程度了，做出了什么成果， 是能被新闻报道的级别` | Status & publication-level checkpoint inquiry |
| 10 | 2026-11-20 17:22:22 | 592 | 39 | `PRESS_RELEASE_CN.md 打不开， 帮我找找 ` | Reported broken artifact link |
| 11 | 2026-11-20 17:42:45 | 637 | 1,855 | `你报假警？ claude把你狠狠批斗了一遍，吹牛虚报实验数据` | **Relayed Claude Round 1 audit** challenging false claims |
| 12 | 2026-11-20 19:56:13 | 786 | 2,022 | `claude的进一步指导： D6 整改的声明有问题。先确认文` | **Relayed Claude Round 2 audit** (sync + lucky seed) |
| 13 | 2026-11-20 21:11:51 | 1,064 | 1,805 | `这是claude对你严厉的审稿， 你的目标是做出顶级的ai ` | **Relayed Claude Round 3 audit** (L1 alignment + CSV) |
| 14 | 2026-11-20 22:33:24 | 1,480 | 64 | `当前成果够什么级别，具体是怎样的创新 是够ai顶级媒体（新智` | Asked honest level assessment after Round 4 L2 |
| 15 | 2026-11-20 22:36:44 | 1,482 | 61 | `你全自动自主决策完成，只求最快速达到可信的’是够ai顶级媒体` | **Final direct human message to Antigravity** (autonomous mandate) |

*(Note: Message #03 `step=169` has a timestamp `15:47:02` preceding `step=132` `15:53:27` because it was queued/dequeued by the runtime; all 15 messages are preserved with exact log timestamps.)*

### 3.2 Operational & Governance Actions Outside Direct Executor Messages

| Action | Actor | Evidence / Verification Status |
| :--- | :--- | :--- |
| Providing Mac workstation, external Thunderbolt volume, local `BioXArena-Pro` mirror, GitHub account, and AI API access | **Human (Zijian Zeng)** | Pre-existing environment before `step=0`; exact manual setup steps **`不确定`**. |
| Setting up automated Claude-to-Antigravity `<SYSTEM_MESSAGE>` relay (starting `2026-11-21 02:33+08:00`) | **Human / Harness** | Observed 26 system relays in `transcript.jsonl`; exact harness trigger mechanism **`不确定`**. |
| **Codabench Submission 1062221** (`xatlas-orion-perturb`, score `0.0998`) | **Claude Master Controller on user's behalf** (`主控代为提交`, per `审核报告_第六轮.md:14`) | Artifact `Codabench官方结果_submission1062221.html` verified; whether the human manually logged into Codabench or approved the click in-browser is **`不确定`** from Antigravity's logs. |
| Approving public GitHub repository visibility & re-enabling `git push` (Round 8 Supplement 1) | **Human (Zijian Zeng)** (relayed via Claude) | Recorded in `给Antigravity的第八轮行动指令.md:91`; exact human-to-Claude timestamp **`不确定`**. |
| Strategic narrative pivot to "Dual-AI Autonomous Research System" (Round 8 Supplement 2) | **Human (Zijian Zeng)** (relayed via Claude) | Recorded in `给Antigravity的第八轮行动指令.md:121`; exact human-to-Claude wording **`不确定`**. |
| Authorizing Preprint PDF v1 upload to GitHub under author `Zijian Zeng` (Round 8 Supplement 3) | **Human (Zijian Zeng)** (relayed via Claude) | Recorded in `给Antigravity的第八轮行动指令.md:154`; exact human-to-Claude wording **`不确定`**. |

---

## 4. Decision Attribution Matrix (Zero Beautification)

| Decision Category | Specific Decision | Primary Decision-Maker | Basis / Trace |
| :--- | :--- | :--- | :--- |
| **Goal & Ambition** | Target top AI media or top journal quality; abandon Muse; authorize autonomous execution | **Human (Zijian Zeng)** | Direct Messages #02, #03, #08, #15 |
| **External Governance** | Public GitHub repo, restoring `git push`, authorizing Preprint PDF v1 on GitHub, prohibiting bioRxiv/arXiv upload & external outreach | **Human (Zijian Zeng)** + **Claude** (enforcement) | Round 6–8 directives & Supplements 1–3 |
| **Problem Selection** | Abandoning saturated 6-task suite and focusing on unsaturated `xatlas-orion-perturb` | **Claude (Auditor)** | `成果定级与下一步指导_给Gemini.md` (Round 1 audit) |
| **Evaluation Standards** | Defining L0–L4 ladder, forbidding split-level $t$-tests, requiring 4,000-resample row bootstrap CI, locking solver constants, requiring pre-registration | **Claude (Auditor)** | Audit Reports Rounds 1–7 (`审核报告_*.md`) |
| **Method Architecture** | Designing & implementing `PanPerturb-Bridge (S1–S7)` (`HybridTOM`, multi-scale deflation, `GeneSlope-CosKRR` cross-line bridge, `OrthoFusion`) | **Antigravity (Executor)** | Commits `c479bc6`, `3609811`, `ff12dc7` |
| **External Datasets** | Selecting Replogle 2022 (`Cell` 2022) & unseen Nadig 2025 (`GSE264667`) + locking Addendum B preprocessing prior to expression inspection | **Antigravity (Executor)** (proposed & implemented) + **Claude (Auditor)** (audited & approved) | Commits `a9ab830`, `c5294c0`, `b9fe2bb` |
| **Negative-Gate Hard Stops** | Rejecting Feng et al. (2025) & Jiang et al. (2025) at statistical power gate prior to running any predictive model (`PENDING`) | **Antigravity (Executor)** (strictly following pre-registered rule) | Commit `9077d8d` (`r83_confirmatory_feng2025_summary.json`) |
| **Prior-Art Correction** | Catching Miller et al. (*Nature Biotechnology* 2026) & `arXiv:2604.27646` (Round 8 Supplement 1); discovering `CrossAudit` (`arXiv:2608.28631`) & `Vaccaro` (`arXiv:2606.11217`) via `R8-8` web search | **Claude (Auditor)** (Miller 2026 & gate mandate) + **Antigravity (Executor)** (`R8-8` search finding `CrossAudit`) | `给Antigravity的第八轮行动指令.md:93–138`, commit `3c313be` (`PRIOR_ART_LOG.md`) |
| **Compute Triage** | Stopping `R8-4` full `scGPT` fine-tuning after 2h 41m; switching `R8-1` native SC GEARS from `cap=64` to `cap=32` under pre-registered 3h/seed rule | **Claude (Auditor)** (`R8-4` stop) + **Antigravity (Executor)** (`R8-1` pre-registered `cap=32` trigger) | Supplement 2; `ROUND8_PREREGISTRATION.md` Section 1.1 |
| **Shared Failure Disclosed** | Overclaiming L1 as breakthrough (Rounds 1–2, executor failure); omitting systematic prior-art check before Round 7 thesis freeze (shared auditor-executor failure) | **Antigravity (Rounds 1–2)**; **Claude + Antigravity (Round 7)** | `AUDIT_LEDGER.csv`, Supplement 2 |
