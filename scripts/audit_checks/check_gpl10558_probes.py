# scripts/audit_checks/check_gpl10558_probes.py
import gzip
import pandas as pd

annot_file = "data/raw/GPL10558.annot.gz"

probes_to_check = [
    "ILMN_1739794", # CD3E
    "ILMN_2396444", # CD14 (Set A)
    "ILMN_1740015", # CD14 (Set B)
    "ILMN_2134453", # FCGR3B (Set A)
    "ILMN_1728639", # FCGR3B (Set B)
    "ILMN_1687301", # VCAN
    "ILMN_1778237", # FN1 (top-table)
    "ILMN_2366463", # FN1 (MaxMean)
]

found_rows = []
with gzip.open(annot_file, "rt", encoding="utf-8", errors="replace") as f:
    # Skip until !platform_table_begin
    for line in f:
        if line.startswith("!platform_table_begin"):
            break
    
    header = f.readline().strip().split("\t")
    id_idx = header.index("ID")
    symbol_idx = header.index("Gene symbol")
    title_idx = header.index("Gene title")
    entrez_idx = header.index("Gene ID")
    accession_idx = header.index("GenBank Accession")
    
    for line in f:
        if line.startswith("!platform_table_end"):
            break
        parts = line.strip().split("\t")
        if len(parts) > id_idx:
            probe_id = parts[id_idx]
            if probe_id in probes_to_check:
                symbol = parts[symbol_idx] if len(parts) > symbol_idx else ""
                title = parts[title_idx] if len(parts) > title_idx else ""
                entrez = parts[entrez_idx] if len(parts) > entrez_idx else ""
                acc = parts[accession_idx] if len(parts) > accession_idx else ""
                found_rows.append({
                    "Probe_ID": probe_id,
                    "Gene_Symbol": symbol,
                    "Gene_Title": title,
                    "Entrez_ID": entrez,
                    "GenBank_Accession": acc
                })

df = pd.DataFrame(found_rows)
print(df.to_string(index=False))
