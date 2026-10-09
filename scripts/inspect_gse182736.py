"""
Inspect GSE182736_All_Expressed_miRNA.xlsx to find actual header row and data structure
"""
import pandas as pd, numpy as np

XLSX = r"d:\Peritoneal Project\data\raw\GSE182736_All_Expressed_miRNA.xlsx"

# Read without assuming header
df_raw = pd.read_excel(XLSX, sheet_name="Sheet1", header=None)
print(f"Raw shape: {df_raw.shape}")
print("\nFirst 20 rows, first 10 cols:")
print(df_raw.iloc[:20, :10].to_string())
print("\nRows 20-40, first 10 cols:")
print(df_raw.iloc[20:40, :10].to_string())

# Look for rows that look like sample IDs / column headers
for i in range(min(30, len(df_raw))):
    row = df_raw.iloc[i].astype(str)
    if any(re.search(r'GSM|sample|control|UF|Case|Saline|PD|fold|logFC|p.val|adj', x, re.I)
           for x in row):
        print(f"\nRow {i} (looks like header or data):")
        print(row.tolist())

import re
