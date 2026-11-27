#!/usr/bin/env bash
# Controller driver: run both parties' Round 17 synthetic modes as black boxes (CPU, zero real cells), then the
# controller's independent recompute (r17_check.py) on each export. It never reads round17/src/*.py or
# round17/audit/*.py and never writes under round17/results/.
#
# Usage (from anywhere):
#   controller/checks/r17_check_synthetic.sh <OUT_DIR> [N_BOOT=400] [N_BOOT_BONF56=2000] [N_PERM=1000]
# Outputs: <OUT_DIR>/{main_smoke,at3_smoke}/ (party exports), <OUT_DIR>/{main,at3}_{run.log,check.txt,check.json}.
# Exit status: 0 only if every party run that completed gives R17_CHECK_MATCH.
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="$ROOT/round16/.venv/bin/python"
OUT="${1:?usage: r17_check_synthetic.sh <OUT_DIR> [N_BOOT] [N_BOOT_BONF56] [N_PERM]}"
NB="${2:-400}"; NB56="${3:-2000}"; NP="${4:-1000}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
case "$OUT" in "$ROOT"/round17/results*) echo "refusing to write under round17/results" >&2; exit 2;; esac
cd "$ROOT" || exit 2

"$PY" -I controller/checks/r17_check.py --self-test || exit 1
had_results=0; [ -e "$ROOT/round17/results" ] && had_results=1

status=0
run_party() {  # name, command...
  local name="$1"; shift
  echo "== $name synthetic run" >&2
  if "$@" > "$OUT/${name}_run.log" 2>&1; then
    echo "== $name controller recompute" >&2
    "$PY" -I controller/checks/r17_check.py "$OUT/${name}_smoke" --quiet --out "$OUT/${name}_check.json" \
      | tee "$OUT/${name}_check.txt"
    [ "${PIPESTATUS[0]}" -eq 0 ] || status=1
  else
    echo "== $name synthetic run FAILED (see $OUT/${name}_run.log); not checked" | tee "$OUT/${name}_check.txt"
    status=1
  fi
}

run_party main "$PY" round17/src/run_eval_r17.py --synthetic-smoke --smoke-out-dir "$OUT/main_smoke" --device cpu \
  --n-boot "$NB" --n-boot-bonf56 "$NB56" --n-perm "$NP" --n-splits-diag 20
run_party at3 "$PY" round17/audit/run_audit_r17.py --phase synthetic_e2e --cohort all --device cpu \
  --output-dir "$OUT/at3_smoke" --n-boot "$NB" --n-boot-bonf56 "$NB56" --n-perm "$NP" --n-splits-diag 20

if [ "$had_results" -eq 0 ] && [ -e "$ROOT/round17/results" ]; then
  echo "WARNING: a synthetic run created round17/results" >&2; status=1
fi
exit $status
