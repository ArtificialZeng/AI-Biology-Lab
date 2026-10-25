#!/usr/bin/env bash
# Reuse the already fixed CAPD build from the neighboring cycle checker.
set -euo pipefail
sensitivity_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
capd_build_dir="${CAPD_BUILD:-$sensitivity_root/../build-capd}"
cmake -S "$sensitivity_root" -B "$sensitivity_root/build-checker" \
 -DCAPD_BUILD="$capd_build_dir" -DCMAKE_BUILD_TYPE=Release
cmake --build "$sensitivity_root/build-checker" --parallel "${BUILD_JOBS:-4}"
"$sensitivity_root/build-checker/verify_period_sensitivity" 1e-9 30
