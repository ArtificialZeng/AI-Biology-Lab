#!/bin/bash
# Controller: independent recompute of the Addendum B (Occam) primary results from a party's exported arrays.
# Scope: primary axis rhj_6983_axis, baseline "reconstructed" (the H-I1_Occam / H-I2_Occam primary arms).
# For each target line, the party's observed.npz and preds_Occam_rhj_6983_axis_reconstructed_{matched,mis_<s>}.npz
# are linked into the layout the unchanged controller Round 16 scorer reads (r16_check.py: Occam arrays under the
# fixed_index_zero names), scored with seed 20261601 / n_boot 4000, and compared with the party's summary_addb.json:
# Occam matched r_orth_cal (estimate + 97.5% CI), Occam - FCT (H-I2), per-source paired contrasts and categories,
# and the H-I1_Occam line verdict. |diff| <= 1e-5 for numbers, identical calls.
# Usage: controller/checks/r16_addb_check.sh <party_results_dir> <out_dir>
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cd "$P" || exit 9
RES=$1; OUT=$2; C=controller/checks; mkdir -p "$OUT"
for ln in a549_ko dld1 cd4t_rest gse295214_h1 feng_ipsc; do
  A="$RES/arrays/$ln"; [ -d "$A" ] || { echo "MISSING $A"; continue; }
  read seen unseen < <(python3 -I -c "
import json; s=json.load(open('$RES/summary_addb.json'))['lines']['$ln']['strata_sources']
print(','.join(s['seen_by_checkpoint']), ','.join(s['unseen_by_checkpoint']))")
  L="$OUT/layout/$ln"; rm -rf "$L"; mkdir -p "$L"
  ln -s "$P/$A/observed.npz" "$L/observed.npz"
  ln -s "$P/$A/preds_Occam_rhj_6983_axis_reconstructed_matched.npz" "$L/preds_STATE_SE_fixed_index_zero_matched.npz"
  for s in ${seen//,/ } ${unseen//,/ }; do
    ln -s "$P/$A/preds_Occam_rhj_6983_axis_reconstructed_mis_$s.npz" "$L/preds_STATE_SE_fixed_index_zero_mis_$s.npz"
  done
  python3 -I $C/r16_check.py "$L" "$seen" "$unseen" "$OUT/$ln.json" 4000 20261601 > /dev/null || echo "r16_check failed on $ln"
done
python3 -I - "$RES" "$OUT" <<'EOF'
import json, sys
res, out = sys.argv[1:3]
S = json.load(open(f"{res}/summary_addb.json"))["lines"]
worst, bad, n = 0.0, [], 0
def num(tag, a, b):
    global worst, n
    n += 1; d = abs(float(a) - float(b)); worst = max(worst, d)
    if not d <= 1e-5: bad.append(f"{tag}: party={a} controller={b} diff={d:.3g}")
for ln, L in S.items():
    try: c = json.load(open(f"{out}/{ln}.json"))
    except FileNotFoundError: bad.append(f"{ln}: no controller recompute"); continue
    P = L["primary"]; m = P["Occam_matched_orth"]
    num(f"{ln}.matched.est", m["estimate"], c["arms"]["matched"]["point"])
    num(f"{ln}.matched.lo", m["ci975_lo"], c["arms"]["matched"]["ci975"][0]); num(f"{ln}.matched.hi", m["ci975_hi"], c["arms"]["matched"]["ci975"][1])
    f = P["Occam_vs_FCT"]
    num(f"{ln}.vsFCT.est", f["estimate"], c["H_I2"]["point"])
    num(f"{ln}.vsFCT.lo", f["ci975_lo"], c["H_I2"]["ci975"][0]); num(f"{ln}.vsFCT.hi", f["ci975_hi"], c["H_I2"]["ci975"][1])
    n += 1
    if P["H_I2_Occam"]["classification"] != c["H_I2"]["call"]:
        bad.append(f"{ln}.H_I2: party={P['H_I2_Occam']['classification']} controller={c['H_I2']['call']}")
    for s, ps in P["Occam_paired_vs_mismatched"]["per_source"].items():
        num(f"{ln}.paired.{s}.est", ps["estimate"], c["pairs"][s]["point"])
        num(f"{ln}.paired.{s}.lo", ps["ci975_lo"], c["pairs"][s]["ci975"][0])
        cat_p = P["H_I1_Occam"]["per_source"][s]["category"]; n += 1
        if cat_p != c["pairs"][s]["call"]: bad.append(f"{ln}.H_I1.{s}: party={cat_p} controller={c['pairs'][s]['call']}")
    n += 1
    if P["H_I1_Occam"]["line_summary"] != c["H_I1_line"]:
        bad.append(f"{ln}.H_I1 line: party={P['H_I1_Occam']['line_summary']} controller={c['H_I1_line']}")
    print(f"{ln:13s} Occam {c['arms']['matched']['point']:+.4f}  Occam-FCT {c['H_I2']['point']:+.4f} {c['H_I2']['call']:22s} H-I1 {c['H_I1_line']}")
print(json.dumps({"n_checks": n, "n_fail": len(bad), "worst": worst}))
for b in bad: print("FAIL", b)
print("R16_ADDB_CONTROLLER_MATCH" if not bad else "R16_ADDB_CONTROLLER_DIFF")
EOF
