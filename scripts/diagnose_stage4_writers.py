import glob
import os
import re
import hashlib
import datetime

def get_file_info(filepath):
    abs_path = os.path.abspath(filepath)
    if not os.path.exists(abs_path):
        return None
    stat = os.stat(abs_path)
    mtime_utc = datetime.datetime.fromtimestamp(stat.st_mtime, tz=datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
    h = hashlib.sha256()
    with open(abs_path, 'rb') as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return {
        "path": abs_path,
        "mtime": mtime_utc,
        "size": stat.st_size,
        "sha256": h.hexdigest()
    }

def main():
    print("=" * 100)
    print("STAGE 4 CSV FILES ON DISK NOW")
    print("=" * 100)
    stage4_files = sorted(glob.glob('results/tables/stage4_*.csv'))
    for f in stage4_files:
        info = get_file_info(f)
        print(f"File:   {os.path.basename(f)}")
        print(f"  Path:   {info['path']}")
        print(f"  mtime:  {info['mtime']}")
        print(f"  Size:   {info['size']} bytes")
        print(f"  SHA256: {info['sha256']}")
        print()

    print("=" * 100)
    print("SEARCHING FOR WRITERS IN ALL SCRIPTS")
    print("=" * 100)
    script_files = sorted(glob.glob('scripts/*.*'))
    
    # We want to check every stage4 csv
    stage4_basenames = [os.path.basename(f) for f in stage4_files]
    
    # Map from stage4 basename to list of writers
    writers = {b: [] for b in stage4_basenames}
    
    for sf in script_files:
        if not (sf.endswith('.py') or sf.endswith('.R') or sf.endswith('.sh')):
            continue
        s_info = get_file_info(sf)
        with open(sf, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
            for line_no, line in enumerate(lines, 1):
                for b in stage4_basenames:
                    if b in line:
                        # Check if it writes
                        if any(w in line for w in ['to_csv', 'write.csv', 'fwrite', 'write_csv']):
                            writers[b].append({
                                "script": sf,
                                "line_no": line_no,
                                "line": line.strip(),
                                "script_mtime": s_info['mtime']
                            })

    for b, w_list in writers.items():
        print(f"Target CSV: {b}")
        if not w_list:
            print("  (No direct writer line matched with filename string literal, check scripts that construct path)")
        for w in w_list:
            print(f"  Writer Script: {w['script']} (mtime: {w['script_mtime']})")
            print(f"    Line {w['line_no']}: {w['line']}")
        print()

if __name__ == '__main__':
    main()
