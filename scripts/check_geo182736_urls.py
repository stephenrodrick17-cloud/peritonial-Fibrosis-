"""Check GSE182736 URL availability and fetch via NCBI eutils/GEO API"""
import requests, os

urls = [
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182736&targ=self&form=text&view=quick",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182736/soft/GSE182736_family.soft.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182736/matrix/GSE182736_series_matrix.txt.gz",
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&id=GSE182736",
]

for url in urls:
    try:
        r = requests.head(url, timeout=30, allow_redirects=True)
        cl = r.headers.get("Content-Length", "?")
        print(f"{r.status_code} | {cl:>12} bytes | {url[:100]}")
    except Exception as e:
        print(f"ERROR: {e} | {url[:100]}")
