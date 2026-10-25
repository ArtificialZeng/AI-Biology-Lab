#!/usr/bin/env bash
# Fixed-source local build; no global install and no publication.
set -euo pipefail
research_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
capd_revision=37057e5899bccaa6cd5701057b7e6ea5e232ef8b
capd_source="$research_root/CAPD"
if [[ ! -d "$capd_source/.git" ]]; then
 git clone --no-checkout https://github.com/CAPDGroup/CAPD.git "$capd_source"
 git -C "$capd_source" checkout --detach "$capd_revision"
fi
[[ "$(git -C "$capd_source" rev-parse HEAD)" == "$capd_revision" ]] || {
 echo "CAPD source revision mismatch" >&2; exit 1;
}
if [[ "$(uname -s)" == Darwin ]]; then
 interval_backend=NATIVE
 echo "Mac run is development evidence; Linux x86 reproduction is required."
else
 interval_backend=FILIB
fi
cmake -S "$capd_source" -B "$research_root/build-capd" \
 -DCAPD_INTERVAL_TYPE="$interval_backend" -DCAPD_BUILD_TESTS=OFF \
 -DCAPD_BUILD_EXAMPLES=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build "$research_root/build-capd" --parallel "${BUILD_JOBS:-4}"
cmake -S "$research_root" -B "$research_root/build-checker" \
 -DCAPD_BUILD="$research_root/build-capd" -DCMAKE_BUILD_TYPE=Release
cmake --build "$research_root/build-checker" --parallel "${BUILD_JOBS:-4}"
"$research_root/build-checker/verify_cycle" 1e-8 30
