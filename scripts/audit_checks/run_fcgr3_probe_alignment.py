# scripts/audit_checks/run_fcgr3_probe_alignment.py
import os
import gzip
from Bio import Align

# 1. Output directory
os.makedirs("results/audit", exist_ok=True)
out_file = "results/audit/fcgr3_alignment_raw.txt"

# 2. Read probe sequences from GPL10558.annot.gz
annot_path = "data/raw/GPL10558.annot.gz"
probes_to_find = ["ILMN_1703679", "ILMN_2112580", "ILMN_1728639", "ILMN_2134453"]
probe_data = {}

with gzip.open(annot_path, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!platform_table_begin"):
            header = f.readline().strip().split("\t")
            break
    id_col = header.index("ID")
    sym_col = header.index("Gene symbol")
    acc_col = header.index("GenBank Accession")
    seq_col = header.index("Platform_SEQUENCE")
    
    for line in f:
        if line.startswith("!platform_table_end"):
            break
        parts = line.strip().split("\t")
        pid = parts[id_col]
        if pid in probes_to_find:
            probe_data[pid] = {
                "probe_id": pid,
                "symbol": parts[sym_col],
                "accession_annot": parts[acc_col],
                "sequence": parts[seq_col],
                "length": len(parts[seq_col])
            }

# 3. Read transcripts
transcripts = {}
for acc, fname in [("NM_000569", "data/raw/NM_000569.fasta"), ("NM_000570", "data/raw/NM_000570.fasta")]:
    with open(fname, "r") as f:
        head = f.readline().strip()
        seq = "".join(l.strip().upper() for l in f if l.strip())
        transcripts[acc] = {
            "header": head,
            "accession_version": head.split()[0].lstrip(">"),
            "sequence": seq,
            "length": len(seq)
        }

# Helper for reverse complement
def rev_comp(s):
    tab = str.maketrans("ACGTN", "TGCAN")
    return s.translate(tab)[::-1]

# 4. Pairwise Aligner Configuration
aligner = Align.PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.open_gap_score = -2.0
aligner.extend_gap_score = -1.0

scoring_params = {
    "mode": aligner.mode,
    "match_score": aligner.match_score,
    "mismatch_score": aligner.mismatch_score,
    "open_gap_score": aligner.open_gap_score,
    "extend_gap_score": aligner.extend_gap_score
}

lines = []
lines.append("================================================================================")
lines.append("FCGR3 PROBE COMPREHENSIVE LOCAL ALIGNMENT AUDIT (STAGE 6)")
lines.append("================================================================================")
lines.append(f"Aligner: Bio.Align.PairwiseAligner (Biopython)")
lines.append(f"Scoring Parameters: {scoring_params}\n")

lines.append("--- REFERENCE TRANSCRIPTS LOADED ---")
for acc, tinfo in transcripts.items():
    lines.append(f"Accession & Version: {tinfo['accession_version']}")
    lines.append(f"Header: {tinfo['header']}")
    lines.append(f"Length: {tinfo['length']} bp\n")

lines.append("--- PROBES READ FROM GPL10558.annot.gz ---")
for pid in probes_to_find:
    pinfo = probe_data[pid]
    lines.append(f"Probe ID: {pinfo['probe_id']} | GPL10558 Symbol: {pinfo['symbol']} | Accession: {pinfo['accession_annot']}")
    lines.append(f"Sequence: {pinfo['sequence']} (Length: {pinfo['length']} bp)\n")

lines.append("================================================================================")
lines.append("PAIRWISE ALIGNMENT RESULTS ACROSS BOTH STRANDS (+ AND -) AND BOTH TRANSCRIPTS")
lines.append("================================================================================\n")

for pid in probes_to_find:
    pinfo = probe_data[pid]
    pseq = pinfo["sequence"]
    lines.append(f"################################################################################")
    lines.append(f"PROBE: {pid} (Target: {pinfo['symbol']}, Annot Acc: {pinfo['accession_annot']})")
    lines.append(f"Sequence (5'->3'): {pseq}")
    lines.append(f"################################################################################\n")
    
    for t_acc in ["NM_000569", "NM_000570"]:
        tinfo = transcripts[t_acc]
        tseq = tinfo["sequence"]
        t_ver = tinfo["accession_version"]
        
        for strand_name, qseq in [("+ (Sense)", pseq), ("- (Reverse Complement)", rev_comp(pseq))]:
            alignments = aligner.align(tseq, qseq)
            best_aln = alignments[0]
            score = best_aln.score
            
            # Extract coordinates
            t_start = best_aln.coordinates[0][0]
            t_end = best_aln.coordinates[0][-1]
            q_start = best_aln.coordinates[1][0]
            q_end = best_aln.coordinates[1][-1]
            
            # Count matches, mismatches, gaps
            # Formatted alignment string
            aln_str = format(best_aln)
            
            # Count exact matches, mismatches, indels in the aligned block
            t_aligned_sub = best_aln[0]
            q_aligned_sub = best_aln[1]
            matches = sum(1 for a, b in zip(t_aligned_sub, q_aligned_sub) if a == b and a != "-")
            mismatches = sum(1 for a, b in zip(t_aligned_sub, q_aligned_sub) if a != b and a != "-" and b != "-")
            gaps_target = sum(1 for a in t_aligned_sub if a == "-")
            gaps_query = sum(1 for b in q_aligned_sub if b == "-")
            total_gaps = gaps_target + gaps_query
            
            lines.append(f"Target: {t_ver} ({'FCGR3A' if '569' in t_acc else 'FCGR3B'}, Length: {tinfo['length']} bp) | Strand: {strand_name}")
            lines.append(f"  Alignment Score: {score}")
            lines.append(f"  Transcript Coordinates (0-based): {t_start}..{t_end} (1-based: {t_start + 1}..{t_end})")
            lines.append(f"  Probe Subsequence Aligned: {q_start}..{q_end} of {len(qseq)} bp")
            lines.append(f"  Matches: {matches} | Mismatches: {mismatches} | Gaps (Indels): {total_gaps} (Target gaps: {gaps_target}, Query gaps: {gaps_query})")
            lines.append("  Alignment Diagram:")
            for line_aln in aln_str.strip().split("\n"):
                lines.append("    " + line_aln)
            lines.append("")

output_text = "\n".join(lines)
with open(out_file, "w", encoding="utf-8") as f:
    f.write(output_text)

print(f"Alignment report successfully written to {out_file} ({len(lines)} lines).")
