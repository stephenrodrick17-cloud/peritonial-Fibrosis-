# scripts/audit_checks/check_directive2.py
import pandas as pd

prim = pd.read_csv("results/tables/stage4_primary_edger_pseudobulk.csv")
sens = pd.read_csv("results/tables/stage4_sensitivity_edger_pseudobulk.csv")

p_st = prim[(prim["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") & 
            (prim["Contrast"] == "LV_UF_vs_LV_NOT_UF")].set_index("Gene")

s_b = sens[(sens["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") & 
           (sens["Contrast"] == "LV_UF_vs_LV_NOT_UF") & 
           (sens["Analysis"] == "sens_b_no_LV_UF3")].set_index("Gene")

s_e = sens[(sens["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") & 
           (sens["Contrast"] == "LV_UF_vs_LV_NOT_UF") & 
           (sens["Analysis"] == "sens_e_ge50cells")].set_index("Gene")

s_c = sens[(sens["Cell_Type"] == "stromal / mesothelial-lineage (unresolved)") & 
           (sens["Contrast"] == "LV_UF_vs_LV_NOT_UF") & 
           (sens["Analysis"] == "sens_c_no_scDblFinder")].set_index("Gene")

print("=== DONOR COUNTS PER GROUP ===")
print(f"Primary:                     LV_UF = {p_st['Donors_LV_UF'].iloc[0]}, LV_NOT_UF = {p_st['Donors_LV_NOT_UF'].iloc[0]}")
print(f"No LV_UF-3 (sens_b):         LV_UF = {s_b['Donors_LV_UF'].iloc[0]}, LV_NOT_UF = {s_b['Donors_LV_NOT_UF'].iloc[0]}")
print(f">= 50 cells/donor (sens_e):   LV_UF = {s_e['Donors_LV_UF'].iloc[0]}, LV_NOT_UF = {s_e['Donors_LV_NOT_UF'].iloc[0]}")
print(f"scDblFinder-removed (sens_c): LV_UF = {s_c['Donors_LV_UF'].iloc[0]}, LV_NOT_UF = {s_c['Donors_LV_NOT_UF'].iloc[0]}")

hub_order = ["VCAN", "EDIL3", "COL8A1", "THBS3", "FN1", "COMP", "COL3A1", "LOX", "INHBA", "COL11A1", "ISM1"]

# Let's inspect the columns
runs = [
    ("Primary (m=98 family FDR)", p_st),
    ("No LV_UF-3 (sens_b)", s_b),
    (">= 50 cells/donor (sens_e)", s_e),
    ("scDblFinder-removed (sens_c)", s_c)
]

for title, df_run in runs:
    print(f"\n--- {title} ---")
    print(f"{'Gene':8s} | {'log2FC':7s} | {'PValue':7s} | {'FDR (m=98)':10s} | {'Evidence_Status':18s}")
    for g in hub_order:
        r = df_run.loc[g]
        fdr_val = f"{r['BH_FDR']:.4f}" if ('BH_FDR' in r and pd.notna(r['BH_FDR'])) else "NA"
        print(f"{g:8s} | {r['log2FC']:7.4f} | {r['PValue']:7.4f} | {fdr_val:10s} | {r['Evidence_Status']:18s}")
