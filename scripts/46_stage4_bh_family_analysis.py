"""
Script 46: Complete 96-Test Benjamini-Hochberg Family Audit & Step-by-Step Derivation
Rules:
- Read directly from results/tables/stage4_primary_edger_pseudobulk.csv
- Print CSV path, mtime, SHA256
- Sort 96 tests by P-value
- Show raw BH ratio (p * 96 / rank) and monotonic step-up enforcement
- Explicitly trace COL8A1 (FDR 0.1287) and THBS3 (FDR 0.2103)
- Assert that computed FDR matches CSV BH_FDR column identically
- Seed 42, no causal language.
"""
import os
import hashlib
import datetime
import numpy as np
import pandas as pd

np.random.seed(42)

print("=" * 80)
print("STAGE 4: BENJAMINI-HOCHBERG MULTIPLE TESTING FAMILY AUDIT (m = 96 TESTS)")
print("=" * 80)

primary_csv = "results/tables/stage4_primary_edger_pseudobulk.csv"
st = os.stat(primary_csv)
mtime_utc = datetime.datetime.fromtimestamp(st.st_mtime, tz=datetime.timezone.utc).isoformat()
with open(primary_csv, "rb") as f:
    sha256 = hashlib.sha256(f.read()).hexdigest()

print(f"Target CSV:   {primary_csv}")
print(f"Mtime (UTC):  {mtime_utc}")
print(f"SHA256:       {sha256}")
print(f"Size (Bytes): {st.st_size}")

df = pd.read_csv(primary_csv)

# Filter to 96 tests passing evidence rule
pass_ev = df[(df["Evidence_Status"] == "PASS") & df["PValue"].notna()].copy()
m = len(pass_ev)
print(f"\nTotal tests in BH Family (passing descriptive evidence rule): m = {m}")
assert m == 96, f"Assertion failed: expected m=96 tests, found {m}"

# Sort by P-value ascending
pass_ev = pass_ev.sort_values("PValue").reset_index(drop=True)
pass_ev["rank"] = np.arange(1, m + 1)

# Step 1: Raw Benjamini-Hochberg ratio: p * m / rank
pass_ev["raw_BH_ratio"] = pass_ev["PValue"] * m / pass_ev["rank"]

# Step 2: Monotonic step-up enforcement: min_{k >= i} (raw_BH_ratio_k)
raw_ratios = pass_ev["raw_BH_ratio"].values
monotonic_fdr = np.minimum.accumulate(raw_ratios[::-1])[::-1]
# Cap at 1.0
monotonic_fdr = np.minimum(monotonic_fdr, 1.0)
pass_ev["computed_BH_FDR"] = monotonic_fdr

# Assert match with CSV's BH_FDR column
diff = np.abs(pass_ev["computed_BH_FDR"].values - pass_ev["BH_FDR"].values)
max_diff = np.max(diff)
print(f"Maximum difference between manual step-by-step formula and CSV BH_FDR: {max_diff:.2e}")
assert max_diff < 1e-6, "Assertion failed: computed FDR does not match CSV BH_FDR column!"
print("Assertion PASSED: The FDR column in stage4_primary_edger_pseudobulk.csv reproduces the BH formula exactly.")

# Save sorted table
out_csv = "results/tables/stage4_primary_96_bh_family_sorted.csv"
pass_ev.to_csv(out_csv, index=False)
st_out = os.stat(out_csv)
mtime_out = datetime.datetime.fromtimestamp(st_out.st_mtime, tz=datetime.timezone.utc).isoformat()
with open(out_csv, "rb") as f:
    sha_out = hashlib.sha256(f.read()).hexdigest()

print(f"\nSaved full sorted 96-test table to: {out_csv}")
print(f"Mtime (UTC): {mtime_out}")
print(f"SHA256:      {sha_out}")

# Print trace for COL8A1 and THBS3
print("\n" + "=" * 80)
print(">>> DETAILED MATHEMATICAL TRACE FOR COL8A1 AND THBS3:")
print("=" * 80)

# COL8A1
col8a1_row = pass_ev[
    (pass_ev["Gene"] == "COL8A1") &
    (pass_ev["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") &
    (pass_ev["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].iloc[0]

rk_col = int(col8a1_row["rank"])
pval_col = col8a1_row["PValue"]
raw_ratio_col = col8a1_row["raw_BH_ratio"]
fdr_col = col8a1_row["BH_FDR"]

print(f"1. COL8A1 [Stromal, LV_UF vs LV_NOT_UF]:")
print(f"   Rank i:              {rk_col}")
print(f"   Raw P-value:         {pval_col:.8f}")
print(f"   Raw BH ratio:        P * 96 / {rk_col} = {pval_col:.8f} * 96 / {rk_col} = {raw_ratio_col:.6f}")
print(f"   Subsequent ratios:   min_{{k>={rk_col}}} (P_k * 96 / k) = {fdr_col:.6f}")
print(f"   Reported BH-FDR:     {fdr_col:.4f} (or 0.1287)")

# THBS3
thbs3_row = pass_ev[
    (pass_ev["Gene"] == "THBS3") &
    (pass_ev["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") &
    (pass_ev["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].iloc[0]

rk_thbs = int(thbs3_row["rank"])
pval_thbs = thbs3_row["PValue"]
raw_ratio_thbs = thbs3_row["raw_BH_ratio"]
fdr_thbs = thbs3_row["BH_FDR"]

print(f"\n2. THBS3 [Stromal, LV_UF vs LV_NOT_UF]:")
print(f"   Rank i:              {rk_thbs}")
print(f"   Raw P-value:         {pval_thbs:.8f}")
print(f"   Raw BH ratio:        P * 96 / {rk_thbs} = {pval_thbs:.8f} * 96 / {rk_thbs} = {raw_ratio_thbs:.6f}")
print(f"   Subsequent ratios:   min_{{k>={rk_thbs}}} (P_k * 96 / k) = {fdr_thbs:.6f}")
print(f"   Reported BH-FDR:     {fdr_thbs:.4f} (or 0.2103)")

# Print full 96-test table
print("\n" + "=" * 80)
print(">>> FULL 96-TEST TABLE SORTED BY P-VALUE (VERBATIM PROGRAMMATIC PRINT):")
print("=" * 80)

display_cols = ["rank", "Cell_Type", "Contrast", "Gene", "PValue", "raw_BH_ratio", "BH_FDR"]
print(pass_ev[display_cols].to_string(index=False))

print("\nBH family audit completed successfully.")
print("=" * 80)
