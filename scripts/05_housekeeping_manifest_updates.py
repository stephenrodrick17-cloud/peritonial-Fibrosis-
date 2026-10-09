"""
Script 05: Housekeeping Manifest Updates and README Verification (H2, H3, H4)
Updates provenance/data_manifest.json with:
- Detailed per-sample confirmation for GSE130387 (Mus musculus biological tissue with Affymetrix synthetic spike-in probes on GPL17107)
- Explicit sample-size and study design notation for GSE121372 (n=1 per condition/timepoint, DESCRIPTIVE ONLY)
- Platform annotation updates for GSE62928 (GPL13158)
"""

import json
import os

def update_manifest():
    manifest_file = 'provenance/data_manifest.json'
    if not os.path.exists(manifest_file):
        raise FileNotFoundError(f"{manifest_file} not found")
        
    with open(manifest_file, 'r') as f:
        manifest = json.load(f)
        
    # H2: GSE130387 Species & Platform details
    if 'GSE130387' in manifest:
        manifest['GSE130387']['biological_species'] = 'Mus musculus'
        manifest['GSE130387']['species_explanation'] = (
            "Per-sample GSM metadata confirms biological source is 100% Mus musculus "
            "(C57BL/6 mouse peritoneal tissue: 3 PDF-induced fibrosis vs 3 saline controls). "
            "The NCBI taxon annotation lists 'synthetic construct; Mus musculus' because the "
            "Affymetrix GeneChip miRNA 4.0 Array (GPL17107) includes synthetic oligonucleotide "
            "spike-in and normalization probes alongside biological miRNA probes."
        )
        manifest['GSE130387']['curation_status'] = "VERIFIED MOUSE PERITONEAL TISSUE (Cross-Species Mapping Required)"

    # H3: GSE121372 Un-replicated Timecourse details
    if 'GSE121372' in manifest:
        manifest['GSE121372']['replicate_structure'] = 'n=1 per condition per timepoint (No biological replicates)'
        manifest['GSE121372']['statistical_designation'] = 'DESCRIPTIVE ONLY: Fold changes, no P-values reported as inference'
        manifest['GSE121372']['sample_details'] = [
            {"gsm": "GSM3433436", "title": "HPMC_control_6hr", "condition": "Control", "timepoint": "6h", "n": 1},
            {"gsm": "GSM3433437", "title": "HPMC_TGF-b1_6hr", "condition": "TGF-beta1", "timepoint": "6h", "n": 1},
            {"gsm": "GSM3433438", "title": "HPMC_control_24hr", "condition": "Control", "timepoint": "24h", "n": 1},
            {"gsm": "GSM3433439", "title": "HPMC_TGF-b1_24hr", "condition": "TGF-beta1", "timepoint": "24h", "n": 1}
        ]

    # H4: GSE62928 Platform
    if 'GSE62928' in manifest:
        manifest['GSE62928']['platform_gpl'] = '13158'
        manifest['GSE62928']['platform_title'] = 'Affymetrix HT HG-U133 Plus PM Array Plate (GPL13158)'

    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)
    print("Updated provenance/data_manifest.json with H2, H3, and H4 annotations.")

if __name__ == '__main__':
    update_manifest()
