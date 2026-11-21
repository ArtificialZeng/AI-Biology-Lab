"""Minimal reproduction for scPerturBench issue SCPB-03:
1. Perturbation_generalization/Chemical/myPRnet.py overwrites split_ood_multi_task3
   with split_ood_multi_task2 (lines 30-31, 37-38, 64-65), causing seed 3 to evaluate
   only the intersection of seed 2 and seed 3 OOD perturbations in all_dataset_chemical.csv.
2. Perturbation_generalization/Chemical/myCPA_chemical.py (lines 103-106) uses wrong
   obs column 'split' (instead of split_key='split_ood_multi_task{seed}') and performs
   an AnnData slice view assignment.
"""
import os
import re
import pandas as pd

SCPB_ROOT = os.environ.get(
    "SCPERTURB_ROOT",
    "/tmp/round9_S4/repos/scPerturBench/scPerturBench-main",
)

print("=== 1. Inspect myPRnet.py bug & prove impact in all_dataset_chemical.csv ===")
with open(f"{SCPB_ROOT}/Perturbation_generalization/Chemical/myPRnet.py") as f:
    pr_lines = f.readlines()
for ln in [30, 31, 37, 38, 64, 65]:
    print(f"  myPRnet.py:{ln}: {pr_lines[ln - 1].rstrip()}")

chem = pd.read_csv(f"{SCPB_ROOT}/Results/Chemical_perturbation/all_dataset_chemical.csv")
for ds in ["sciplex3_A549", "sciplex3_K562", "sciplex3_MCF7", "sciplex3_comb"]:
    s2 = set(chem[(chem["method"] == "trainMean") & (chem["DataSet"] == ds) & (chem["seed"] == 2)]["perturb"])
    s3 = set(chem[(chem["method"] == "trainMean") & (chem["DataSet"] == ds) & (chem["seed"] == 3)]["perturb"])
    pr3 = set(chem[(chem["method"] == "PRnet") & (chem["DataSet"] == ds) & (chem["seed"] == 3)]["perturb"])
    print(f"  {ds:14s}: |seed2|={len(s2):3d}, |seed3|={len(s3):3d}, |s2 & s3|={len(s2 & s3):2d}, |PRnet_seed3|={len(pr3):2d}, exact_match={pr3 == (s2 & s3)}")

print("\n=== 2. Inspect myCPA_chemical.py vs fixed myCPA.py ===")
with open(f"{SCPB_ROOT}/Perturbation_generalization/Chemical/myCPA_chemical.py") as f:
    cpa_lines = f.readlines()
for ln in [74, 98, 103, 104, 105, 106, 109]:
    print(f"  myCPA_chemical.py:{ln}: {cpa_lines[ln - 1].rstrip()}")
