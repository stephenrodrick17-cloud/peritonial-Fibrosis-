"""
audit/generate_verbatim_tables.py
Generates audit/verbatim_tables.txt:
For every CSV file in results/tables/, prints:
- File name
- File modification timestamp (UTC)
- DataFrame shape (rows x columns)
- Exact verbatim df.to_string()
"""
import os
import glob
import time
import pandas as pd

table_files = sorted(glob.glob("results/tables/*.csv"))
output_path = "audit/verbatim_tables.txt"

with open(output_path, "w", encoding="utf-8") as out:
    out.write("=" * 80 + "\n")
    out.write("AUDIT: VERBATIM CONTENTS OF ALL RESULTS TABLES\n")
    out.write(f"Generated at: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}\n")
    out.write(f"Total CSV tables: {len(table_files)}\n")
    out.write("=" * 80 + "\n\n")
    
    for tf in table_files:
        fname = os.path.basename(tf)
        mtime = time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime(os.path.getmtime(tf)))
        df = pd.read_csv(tf)
        
        out.write("-" * 80 + "\n")
        out.write(f"FILE: {fname}\n")
        out.write(f"PATH: {tf}\n")
        out.write(f"MODIFIED (UTC): {mtime}\n")
        out.write(f"SHAPE: {df.shape[0]} rows x {df.shape[1]} columns\n")
        out.write("-" * 80 + "\n")
        out.write(df.to_string(index=False))
        out.write("\n\n")

print(f"Successfully generated {output_path} ({len(table_files)} tables).")
