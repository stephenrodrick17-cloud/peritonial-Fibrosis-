# scripts/audit_checks/verify_dois.py
import urllib.request
import json
import sys

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

dois = {
    "Charoentong 2017 (Cell Reports)": "10.1016/j.celrep.2016.12.019",
    "Bindea 2013 (Immunity)": "10.1016/j.immuni.2013.10.003",
    "Naba 2012 (Mol Cell Proteomics)": "10.1074/mcp.M111.014647",
    "Naba 2016 (Matrix Biology)": "10.1016/j.matbio.2015.06.003"
}

print("=== VERIFYING DOIS VIA DOI.ORG CONTENT NEGOTIATION ===")
for ref_key, doi in dois.items():
    url = f"https://doi.org/{doi}"
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/vnd.citationstyles.csl+json", "User-Agent": "Antigravity/1.0 (mailto:audit@example.com)"}
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            title = data.get("title", "No title")
            container = data.get("container-title", "No container")
            author = data.get("author", [{}])[0].get("family", "Unknown")
            issued = data.get("issued", {}).get("date-parts", [[None]])[0][0]
            print(f"\n[{ref_key}]")
            print(f"  First Author: {author}")
            print(f"  DOI:          {doi}")
            print(f"  URL:          {url}")
            print(f"  Resolved:     SUCCESS (HTTP {resp.status})")
            print(f"  Title:        {title}")
            print(f"  Journal:      {container}")
            print(f"  Year:         {issued}")
    except Exception as e:
        print(f"\n[{ref_key}]")
        print(f"  DOI:          {doi}")
        print(f"  Failed:       {e}")
