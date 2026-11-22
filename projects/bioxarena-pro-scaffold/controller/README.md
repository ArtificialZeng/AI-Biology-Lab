# Controller files (Claude, master controller)
Lives in the repo on the internal disk; no dependency on the external sync drive (moved 2026-11-22).
- `agy_env.sh` — source it to get agentapi env for the current Antigravity daemon (auto-detects the newest hub_*.json port).
- `watchers/` — background watchers that nudge idle agents and wake the controller at gates.
- `checks/` — the controller's independent re-computation scripts (dev and held-out).
- `directives/` — controller directives to Antigravity.
