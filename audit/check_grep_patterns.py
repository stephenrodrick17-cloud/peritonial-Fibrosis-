import os
import re

patterns = [
    '1365', '1,365', '1598', '1,598', '1526', '1,526',
    'fabricat', 'HPA0', 'CAB0', 'CHEMBL', 'DrugBank', 'DB1', 'PMID',
    'AUC = 1.000', 'AUC=1.0', 'AUC = 1.0', 'AUC=1.000'
]

repo_root = os.getcwd()
results = {p: [] for p in patterns}

for root, dirs, files in os.walk(repo_root):
    if any(skip in root for skip in ['.git', '.system_generated', 'venv', '.venv', 'scratch']):
        continue
    for f in files:
        if f.endswith(('.png', '.pdf', '.jpg', '.jpeg', '.gz', '.rds', '.RData')):
            continue
        filepath = os.path.join(root, f)
        relpath = os.path.relpath(filepath, repo_root)
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as fp:
                for line_no, line in enumerate(fp, 1):
                    line_lower = line.lower()
                    for pat in patterns:
                        if pat.lower() in line_lower:
                            results[pat].append((relpath, line_no, line.strip()))
        except Exception as e:
            pass

print('=== COMPREHENSIVE GREP REPORT ===')
for pat, matches in results.items():
    print(f'\n>>> PATTERN: "{pat}" (Total matches: {len(matches)})')
    non_tsvs = [m for m in matches if not m[0].endswith('.tsv')]
    print(f'  Non-TSV matches: {len(non_tsvs)}')
    for m in non_tsvs[:25]:
        print(f'    [{m[0]}:{m[1]}] {m[2][:120]}')
