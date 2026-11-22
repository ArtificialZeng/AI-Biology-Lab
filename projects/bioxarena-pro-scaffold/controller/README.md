# Controller files (Claude, master controller)
Lives in the repo on the internal disk; no dependency on the external sync drive (moved 2026-11-22).
- `agy_env.sh` — source it to get agentapi env for the current Antigravity daemon (auto-detects the newest hub_*.json port).
- `watchers/` — background watchers that nudge idle agents and wake the controller at gates.
- `checks/` — the controller's independent re-computation scripts (dev and held-out).
- `directives/` — controller directives to Antigravity.

## Permanent Governance Rule (User Decision — 2026-11-22)
- **媒体稿、新闻稿、投稿推介一律不进公开仓库**（根目录 `.gitignore` 已屏蔽 `*MEDIA_BRIEF*` 与 `*媒体*`），**执行代理（主线及全部 A-T 代理）一律不写、不发**。
- 只有主控按用户明确要求起草，经用户审核同意后才能发布。
- 执行代理不得自行开启新方向；任何新方向提议须先写入 `status/` 等待主控批准。
