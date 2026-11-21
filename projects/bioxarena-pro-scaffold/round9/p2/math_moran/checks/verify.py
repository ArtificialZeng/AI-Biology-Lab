"""Master verification entry point for star graph neutral Moran consensus time."""
import sys
sys.dont_write_bytecode = True

import argparse
import json
from pathlib import Path
import time

from common import require
from exact_markov_checker import run_all_exact_checks
from independent_mc_and_symbolic import run_all_checker2
from negative_controls import run_negative_controls

HERE = Path(__file__).resolve().parent
SUMMARY_PATH = HERE / "verification_summary.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify star Moran consensus-time proofs and controls.")
    parser.add_argument(
        "--write-summary",
        action="store_true",
        help="Write checks/verification_summary.json (default is read-only verification).",
    )
    args = parser.parse_args()

    t0 = time.time()
    c1 = run_all_exact_checks()
    c2 = run_all_checker2()
    nc = run_negative_controls()
    elapsed = round(time.time() - t0, 3)

    summary = {
        "checker1_exact_markov": c1,
        "checker2_symbolic_and_mc": c2,
        "negative_controls": nc,
    }

    if args.write_summary:
        SUMMARY_PATH.write_text(json.dumps(summary, indent=2) + "\n")
    elif SUMMARY_PATH.is_file():
        recorded = json.loads(SUMMARY_PATH.read_text())
        require(recorded == summary, "Live verification results do not match verification_summary.json")

    max_z = max(abs(r["z_score"]) for r in c2["monte_carlo"])
    print(
        f"ALL CHECKS PASS in {elapsed}s: "
        f"Checker1(unreduced={len(c1['unreduced_forward_checks'])}, "
        f"canonical={len(c1['canonical_forward_checks'])}, "
        f"subset={len(c1['subset_coalescing_checks'])}, "
        f"kc={c1['lumped_kc_checks_count']}); "
        f"Checker2(sympy=PASS, mc={len(c2['monte_carlo'])}, max|z|={max_z:.3f}); "
        f"NegativeControls={nc['negative_controls_passed']} caught."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
