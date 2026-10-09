"""
Finish downloading GSE248762_RAW_clean.tar
"""
import os
import time
import subprocess

tar_path = "data/raw/GSE248762_RAW_clean.tar"
url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248762/suppl/GSE248762_RAW.tar"
expected_size = 1014528000

while True:
    current_size = os.path.getsize(tar_path) if os.path.exists(tar_path) else 0
    if current_size >= expected_size:
        print(f"DOWNLOAD COMPLETE! Total size: {current_size:,} bytes")
        break
    pct = (current_size / expected_size) * 100
    print(f"Resuming curl from byte {current_size:,} / {expected_size:,} ({pct:.1f}%)...")
    
    cmd = [
        "curl.exe", "-C", "-",
        "--connect-timeout", "30",
        "--max-time", "120",
        "-o", tar_path,
        url
    ]
    res = subprocess.run(cmd)
    sz = os.path.getsize(tar_path) if os.path.exists(tar_path) else 0
    if sz >= expected_size:
        print(f"Download finished! Total bytes: {sz:,}")
        break
    print("Transfer paused, resuming in 1 second...")
    time.sleep(1)
