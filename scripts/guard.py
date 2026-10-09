"""
Hub Gene Blind Guard System
Strict code-level enforcement to prevent accessing any of the 11 hub genes during Stage 3 scRNA-seq processing.
"""
import os
import datetime

BLOCKED_HUB_GENES = {
    "ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", 
    "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"
}

LOG_FILE = "audit/guard_audit.log"

def guard_check(query_genes, stage="unspecified"):
    """
    Checks whether any query gene matches one of the 11 blocked hub genes.
    Raises RuntimeError if a violation is detected.
    Logs all clean checks to audit/guard_audit.log.
    """
    os.makedirs("audit", exist_ok=True)
    if isinstance(query_genes, str):
        query_genes = [query_genes]
    
    violations = [str(g).upper() for g in query_genes if str(g).upper() in BLOCKED_HUB_GENES]
    
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    if violations:
        msg = f"[GUARD VIOLATION] {timestamp} | Stage: {stage} | Hub genes attempted: {violations}"
        with open(LOG_FILE, "a") as f:
            f.write(msg + "\n")
        raise RuntimeError(
            f"[CODE-LEVEL GUARD TRIGGERED]: Access to hub genes {violations} is strictly "
            f"prohibited during Stage 3 (hub-gene-blind scRNA processing)!"
        )
    else:
        msg = f"[GUARD PASS] {timestamp} | Stage: {stage} | Queried {len(query_genes)} genes. 0 hub genes."
        with open(LOG_FILE, "a") as f:
            f.write(msg + "\n")
        return True
