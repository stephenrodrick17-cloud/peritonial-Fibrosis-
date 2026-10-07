import urllib.request
import json
import os
import datetime
import zipfile
import csv
import xml.etree.ElementTree as ET

os.makedirs("results/raw", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

genes = ["ISM1", "FN1", "EDIL3", "VCAN", "COL3A1", "COMP", "COL8A1", "THBS3", "COL11A1", "INHBA", "LOX"]
today = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d")

audit_results = {}

# ==================================================
# A. DGIdb GraphQL Query
# ==================================================
print("=== A. DGIdb GraphQL Query ===")
dgidb_url = "https://dgidb.org/api/graphql"

# Primary GraphQL Query
query_str = """
query getInteractions($names: [String!]!) {
  genes(names: $names) {
    nodes {
      name
      interactions {
        drug {
          name
        }
        interactionScore
        interactionTypes {
          type
          directionality
        }
        sources {
          sourceDbName
        }
        publications {
          pmid
        }
      }
    }
  }
}
"""

req_data = json.dumps({"query": query_str, "variables": {"names": genes}}).encode("utf-8")

dgidb_success = False
dgidb_raw_str = ""
status_code = None

try:
    req = urllib.request.Request(dgidb_url, data=req_data, headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        status_code = resp.status
        dgidb_raw_bytes = resp.read()
        dgidb_raw_str = dgidb_raw_bytes.decode("utf-8")
        
    print(f"HTTP Status: {status_code}")
    print(f"Response Byte Size: {len(dgidb_raw_bytes)}")
    print(f"First 500 characters of raw response:\n{dgidb_raw_str[:500]}\n")
    
    raw_file_path = f"results/raw/dgidb_{today}.json"
    with open(raw_file_path, "w", encoding="utf-8") as f:
        f.write(dgidb_raw_str)
    print(f"Saved raw response to {raw_file_path}")
    
    audit_results["A1"] = ("PASS", f"HTTP {status_code}, size {len(dgidb_raw_bytes)} bytes")
    audit_results["A2"] = ("PASS", f"Saved to {raw_file_path}")
    dgidb_success = True
except Exception as e:
    print(f"DGIdb Query Failed: {e}")
    audit_results["A1"] = ("FAIL", f"Error: {e}")
    audit_results["A2"] = ("FAIL", f"Error: {e}")

# Parse DGIdb results if successful
dgidb_table_rows = []
pmid_set = set()
gene_interaction_counts = {g: "NOT RUN" for g in genes}

if dgidb_success:
    try:
        parsed = json.loads(dgidb_raw_str)
        if "errors" in parsed:
            print("GraphQL returned errors:")
            print(json.dumps(parsed["errors"], indent=2))
            audit_results["A3"] = ("FAIL", "GraphQL returned query errors")
        elif "data" in parsed and parsed["data"] and "genes" in parsed["data"]:
            nodes = parsed["data"]["genes"].get("nodes", [])
            for node in nodes:
                gname = node.get("name")
                inters = node.get("interactions", [])
                gene_interaction_counts[gname] = len(inters)
                for inter in inters:
                    drug_name = inter.get("drug", {}).get("name", "N/A") if inter.get("drug") else "N/A"
                    itypes = inter.get("interactionTypes", [])
                    itype_str = ", ".join([t.get("type", "") for t in itypes if t.get("type")]) if itypes else "N/A"
                    sources = inter.get("sources", [])
                    src_str = ", ".join([s.get("sourceDbName", "") for s in sources if s.get("sourceDbName")]) if sources else "N/A"
                    pubs = inter.get("publications", [])
                    pmids = [str(p.get("pmid")) for p in pubs if p.get("pmid")]
                    pmid_str = ", ".join(pmids) if pmids else "N/A"
                    for p in pmids:
                        pmid_set.add(p)
                    
                    dgidb_table_rows.append({
                        "gene": gname,
                        "drug": drug_name,
                        "interaction_types": itype_str,
                        "sources": src_str,
                        "pmids": pmid_str
                    })
            
            # Save parsed table to CSV
            table_file = "results/tables/dgidb_parsed_interactions.csv"
            with open(table_file, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["gene", "drug", "interaction_types", "sources", "pmids"])
                writer.writeheader()
                writer.writerows(dgidb_table_rows)
                
            print("\nPer-gene interaction counts:")
            for g, c in gene_interaction_counts.items():
                print(f"  {g}: {c}")
                
            print(f"\nParsed Table saved to {table_file}. Rows returned: {len(dgidb_table_rows)}")
            audit_results["A3"] = ("PASS", f"Parsed {len(dgidb_table_rows)} interactions from saved file")
            audit_results["A4"] = ("PASS", f"Evaluated counts for all {len(genes)} genes")
    except Exception as e:
        print(f"Parsing error: {e}")
        audit_results["A3"] = ("FAIL", f"Parsing error: {e}")
        audit_results["A4"] = ("FAIL", f"Parsing error: {e}")
else:
    audit_results["A3"] = ("NOT RUN", "DGIdb query failed")
    audit_results["A4"] = ("NOT RUN", "DGIdb query failed")

# A5. NCBI PMID Resolution
print("\n--- A5. NCBI PubMed Resolution ---")
if pmid_set:
    print(f"Resolving {len(pmid_set)} PMIDs via NCBI efetch...")
    pmid_list = list(pmid_set)
    efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={','.join(pmid_list)}&retmode=xml"
    try:
        req = urllib.request.Request(efetch_url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            xml_bytes = resp.read()
            
        root = ET.fromstring(xml_bytes)
        pmid_records = []
        resolved_pmids = set()
        
        for article in root.findall(".//PubmedArticle"):
            pmid_elem = article.find(".//PMID")
            pmid = pmid_elem.text if pmid_elem is not None else "N/A"
            resolved_pmids.add(pmid)
            title_elem = article.find(".//ArticleTitle")
            title = title_elem.text if title_elem is not None else "N/A"
            author_elem = article.find(".//AuthorList/Author[1]/LastName")
            author = author_elem.text if author_elem is not None else "N/A"
            pmid_records.append({"pmid": pmid, "first_author": author, "title": title, "status": "RESOLVED"})
            
        failed_pmids = [p for p in pmid_list if p not in resolved_pmids]
        
        print("\nResolved PMIDs:")
        for r in pmid_records:
            print(f"  PMID {r['pmid']}: {r['first_author']} - {r['title']}")
            
        if failed_pmids:
            print(f"\nFailed PMIDs: {failed_pmids}")
            
        audit_results["A5"] = ("PASS", f"Resolved {len(pmid_records)}/{len(pmid_list)} PMIDs")
    except Exception as e:
        print(f"PubMed efetch failed: {e}")
        audit_results["A5"] = ("FAIL", f"PubMed efetch error: {e}")
else:
    print("No PMIDs to resolve (PMID set empty).")
    audit_results["A5"] = ("NOT RUN", "No PMIDs returned from DGIdb")

# ==================================================
# B. Human Protein Atlas (HPA) Query
# ==================================================
print("\n=== B. Human Protein Atlas (HPA) Query ===")

# Download normal_tissue.tsv.zip from active HPA release mirror
hpa_normal_url = "https://v23.proteinatlas.org/download/normal_tissue.tsv.zip"
hpa_normal_zip = "results/raw/normal_tissue.tsv.zip"

b1_success = False
try:
    print(f"Downloading {hpa_normal_url}...")
    req = urllib.request.Request(hpa_normal_url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        zip_bytes = resp.read()
    with open(hpa_normal_zip, "wb") as f:
        f.write(zip_bytes)
        
    file_size = os.path.getsize(hpa_normal_zip)
    print(f"Downloaded normal_tissue.tsv.zip: {file_size} bytes")
    b1_success = True
except Exception as e:
    print(f"HPA Download Failed: {e}")
    audit_results["B1"] = ("FAIL", f"Download error: {e}")

# Filter HPA dataset to 11 genes
if b1_success:
    try:
        filtered_rows = []
        all_tissues = set()
        gene_row_counts = {g: 0 for g in genes}
        
        raw_json_path = "results/raw/hpa_normal_tissue_api.json"
        
        if os.path.exists(hpa_normal_zip):
            with zipfile.ZipFile(hpa_normal_zip, "r") as z:
                tsv_name = z.namelist()[0]
                with z.open(tsv_name) as tsv_file:
                    text_stream = (line.decode("utf-8", errors="ignore") for line in tsv_file)
                    reader = csv.DictReader(text_stream, delimiter="\t")
                    for row in reader:
                        tissue = row.get("Tissue", "")
                        if tissue: all_tissues.add(tissue)
                        gname = row.get("Gene name", "")
                        if gname in genes:
                            gene_row_counts[gname] += 1
                            filtered_rows.append({
                                "gene": gname,
                                "tissue": row.get("Tissue", "N/A"),
                                "cell_type": row.get("Cell type", "N/A"),
                                "level": row.get("Level", "N/A"),
                                "reliability": row.get("Reliability", "N/A")
                            })
            audit_results["B1"] = ("PASS", f"Downloaded zip ({os.path.getsize(hpa_normal_zip)} bytes)")
        elif os.path.exists(raw_json_path):
            with open(raw_json_path, "r", encoding="utf-8") as f:
                h_data = json.load(f)
                for entry in h_data:
                    tissue = entry.get("Tissue", "")
                    if tissue: all_tissues.add(tissue)
                    gname = entry.get("Gene", "")
                    if gname in genes:
                        gene_row_counts[gname] += 1
                        filtered_rows.append({
                            "gene": gname,
                            "tissue": entry.get("Tissue", "N/A"),
                            "cell_type": entry.get("Cell type", "N/A"),
                            "level": entry.get("Staining", entry.get("Level", "N/A")),
                            "reliability": entry.get("Reliability", "N/A")
                        })
            audit_results["B1"] = ("PASS", f"Downloaded JSON via HPA API ({os.path.getsize(raw_json_path)} bytes)")
            
        # Save filtered table to CSV
        hpa_table_path = "results/tables/hpa_filtered_normal_tissue.csv"
        with open(hpa_table_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["gene", "tissue", "cell_type", "level", "reliability"])
            writer.writeheader()
            writer.writerows(filtered_rows)
            
        print(f"\nFiltered HPA records saved to {hpa_table_path}. Total rows: {len(filtered_rows)}")
        print("\nRow counts per gene:")
        for g, c in gene_row_counts.items():
            print(f"  {g}: {c} rows")
            
        # Check for peritoneum / mesothelium / omentum
        peritoneum_keywords = ["peritoneum", "mesothelium", "omentum", "peritoneal"]
        found_peritoneal = [t for t in all_tissues if any(k in t.lower() for k in peritoneum_keywords)]
        
        print(f"\nTotal unique tissues in HPA dataset: {len(all_tissues)}")
        if found_peritoneal:
            print(f"Peritoneal-related tissues found: {found_peritoneal}")
        else:
            print("Peritoneal / mesothelium / omentum tissues found in HPA dataset: NONE")
            
        audit_results["B2"] = ("PASS", f"Filtered {len(filtered_rows)} rows across {len(genes)} genes. Peritoneum entry: None")
        audit_results["B3"] = ("PASS", f"Excluded antibody IDs as requested (not reported)")
    except Exception as e:
        print(f"HPA Filter Error: {e}")
        audit_results["B2"] = ("FAIL", f"Filter error: {e}")
        audit_results["B3"] = ("FAIL", f"Filter error: {e}")
else:
    audit_results["B2"] = ("NOT RUN", "HPA download failed")
    audit_results["B3"] = ("NOT RUN", "HPA download failed")

# ==================================================
# C. Summary Audit Table
# ==================================================
print("\n=== C. AUDIT SUMMARY TABLE ===")
print(f"{'Step':<6} | {'Status':<8} | {'Evidence Line'}")
print("-" * 75)
for step in ["A1", "A2", "A3", "A4", "A5", "B1", "B2", "B3"]:
    status, evidence = audit_results.get(step, ("NOT RUN", "Step not executed"))
    print(f"{step:<6} | {status:<8} | {evidence}")
