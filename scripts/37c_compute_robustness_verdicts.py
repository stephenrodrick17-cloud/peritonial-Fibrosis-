"""
Script 37c: Compute Stage 4 Robustness Verdicts from edgeR Pseudobulk Outputs
Rules:
(1) Primary FDR < 0.05
(2) Sign of log2FC unchanged in sensitivity analyses (b), (c), (d), and (e)
(3) Sign of log2FC unchanged in every leave-one-donor-out run (16/16)
Otherwise "not robust".
Pairs with INSUFFICIENT DATA reported as INSUFFICIENT DATA.
"""
import pandas as pd
import numpy as np

# Load primary results
df_prim = pd.read_csv("results/tables/stage4_primary_edger_pseudobulk.csv")
# Load sensitivity results
df_sens = pd.read_csv("results/tables/stage4_sensitivity_edger_pseudobulk.csv")

# Filter stromal primary contrast: LV_UF vs LV_NOT_UF
st_prim = df_prim[
    (df_prim["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") &
    (df_prim["Contrast"] == "LV_UF_vs_LV_NOT_UF")
].copy()

hub_genes = ["COL11A1", "THBS3", "COL3A1", "FN1", "COL8A1", "VCAN", 
             "EDIL3", "LOX", "INHBA", "ISM1", "COMP"]

def get_sign(val):
    if pd.isna(val):
        return "NA"
    return "+" if val > 0 else ("-" if val < 0 else "0")

rows = []
for g in hub_genes:
    p_row = st_prim[st_prim["Gene"] == g]
    if len(p_row) == 0:
        continue
    
    p_row = p_row.iloc[0]
    ev = p_row["Evidence_Status"]
    lfc = p_row["log2FC"]
    pval = p_row["PValue"]
    fdr = p_row["BH_FDR"]
    prim_sign = get_sign(lfc)
    
    if ev == "INSUFFICIENT DATA":
        rows.append({
            "Gene": g,
            "Evidence": "INSUFFICIENT DATA",
            "Primary_log2FC": round(lfc, 4) if pd.notna(lfc) else np.nan,
            "Primary_P": round(pval, 4) if pd.notna(pval) else np.nan,
            "Primary_FDR": round(fdr, 4) if pd.notna(fdr) else np.nan,
            "Sign_Primary": prim_sign,
            "Sign_b_noUF3": "—",
            "Sign_c_noDbl": "—",
            "Sign_d_ge700g": "—",
            "Sign_e_ge50cells": "—",
            "Sign_f_HVG": "—",
            "LODO_Consistent": "—",
            "Verdict": "INSUFFICIENT DATA"
        })
        continue
    
    # Get sensitivities
    def get_sens_sign(aname):
        sub = df_sens[(df_sens["Analysis"] == aname) & (df_sens["Gene"] == g) & (df_sens["Contrast"] == "LV_UF_vs_LV_NOT_UF")]
        if len(sub) > 0 and pd.notna(sub["log2FC"].iloc[0]):
            return get_sign(sub["log2FC"].iloc[0])
        return "NA"
    
    s_b = get_sens_sign("sens_b_no_LV_UF3")
    s_c = get_sens_sign("sens_c_no_scDblFinder")
    s_d = get_sens_sign("sens_d_ge700genes")
    s_e = get_sens_sign("sens_e_ge50cells")
    s_f = get_sens_sign("sens_f_top2000_hvg_tmm")
    
    # LODO consistency
    lodo_subs = df_sens[
        df_sens["Analysis"].str.startswith("lodo_out_") & 
        (df_sens["Gene"] == g) & 
        (df_sens["Contrast"] == "LV_UF_vs_LV_NOT_UF")
    ]
    lodo_signs = [get_sign(val) for val in lodo_subs["log2FC"] if pd.notna(val)]
    n_match = sum(1 for s in lodo_signs if s == prim_sign)
    lodo_consistent = f"{n_match}/{len(lodo_signs)}"
    
    # Robustness rule:
    # 1. Primary FDR < 0.05
    # 2. Sign unchanged in (b), (c), (d), (e)
    # 3. Sign unchanged in all LODO runs (16/16)
    c1 = (pd.notna(fdr) and fdr < 0.05)
    c2 = (s_b == prim_sign and s_c == prim_sign and s_d == prim_sign and s_e == prim_sign)
    c3 = (n_match == len(lodo_signs) and len(lodo_signs) == 16)
    
    verdict = "supported" if (c1 and c2 and c3) else "not robust"
    
    rows.append({
        "Gene": g,
        "Evidence": "PASS",
        "Primary_log2FC": round(lfc, 4),
        "Primary_P": f"{pval:.4e}" if pd.notna(pval) else "NA",
        "Primary_FDR": f"{fdr:.4e}" if pd.notna(fdr) else "NA",
        "Sign_Primary": prim_sign,
        "Sign_b_noUF3": s_b,
        "Sign_c_noDbl": s_c,
        "Sign_d_ge700g": s_d,
        "Sign_e_ge50cells": s_e,
        "Sign_f_HVG": s_f,
        "LODO_Consistent": lodo_consistent,
        "Verdict": verdict
    })

df_rob = pd.DataFrame(rows)
df_rob.to_csv("results/tables/stage4_stromal_robustness_verdicts.csv", index=False)

print("=" * 80)
print("STAGE 4: STROMAL ROBUSTNESS VERDICTS (Primary contrast: LV_UF vs LV_NOT_UF)")
print("=" * 80)
print(df_rob.to_string(index=False))
