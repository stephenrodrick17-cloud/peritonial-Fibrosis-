# Cohort provenance and patient-overlap assessment

This assessment records the recruiting-site, country, and recruitment-period
information available in the linked paper methods and GEO records. “Not
reported/verified” means the reviewed sources did not establish the fact; a
publication year, ethics approval number, author affiliation, or GEO submission
date was not treated as a recruitment year. No patient-level identifiers are
available across all four cohorts, so no pair is declared formally independent.
This document does not change the preregistered independence rule.

## Cohort-level evidence

| GEO cohort | Linked paper and source | Recruiting center and country | Recruitment years | Overlap evidence |
|---|---|---|---|---|
| **GSE125498** | Parikova et al., PMID [31424258](https://pubmed.ncbi.nlm.nih.gov/31424258/), DOI [10.33549/physiolres.934158](https://doi.org/10.33549/physiolres.934158); GEO [GSE125498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125498). | The GEO contact and PubMed affiliation identify the Institute for Clinical and Experimental Medicine (IKEM), Prague, Czech Republic. The paper compares peritoneal-effluent cells from short- and long-term PD patients. GEO does not give a more specific sample-collection site. | **Not reported in the reviewed GEO record or PubMed abstract.** The study publication is from 2019; this is not taken as a recruitment date. | No statement of overlap with the other three cohorts found in the reviewed GEO record or linked abstract. |
| **GSE130888** | Si et al., PMID [31167927](https://pubmed.ncbi.nlm.nih.gov/31167927/), DOI [10.1126/scitranslmed.aav5341](https://doi.org/10.1126/scitranslmed.aav5341); GEO [GSE130888](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130888). | **China.** The linked paper’s affiliations include the Third and Sixth Affiliated Hospitals of Sun Yat-sen University in Guangzhou and Shunde Hospital of Southern Medical University in Foshan. GEO describes effluent and normal-peritoneum specimens but does not map samples to a specific hospital. Thus a single recruiting center cannot be assigned from the reviewed sources. | **Not reported in the reviewed GEO record or PubMed abstract.** The 2019 article publication date is not a recruitment date. | No patient overlap with the other cohorts is stated in the reviewed GEO record or abstract. Sample-level recruitment hospital and years remain unresolved. |
| **GSE248762** | Diao et al., PMID [39687014](https://pubmed.ncbi.nlm.nih.gov/39687014/), PMCID [PMC11647153](https://pmc.ncbi.nlm.nih.gov/articles/PMC11647153/), DOI [10.1016/j.isci.2024.111383](https://doi.org/10.1016/j.isci.2024.111383); GEO [GSE248762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248762). | The full-text methods identify patient-derived PD effluent and regular PD-center follow-up at the **First Affiliated Hospital of Sun Yat-sen University**, Guangzhou, Guangdong, China. | **2021–2022** (June 1, 2021 to June 1, 2022, as stated in the paper methods). | The paper describes 16 patients in its three PD-vintage groups. No overlap with the other GEO cohorts is stated; GEO sample records do not expose reusable patient identifiers. |
| **GSE92455** | Ruiz-Carpio et al., PMID [28327551](https://pubmed.ncbi.nlm.nih.gov/28327551/), PMCID [PMC5361179](https://pmc.ncbi.nlm.nih.gov/articles/PMC5361179/), DOI [10.1038/srep44941](https://doi.org/10.1038/srep44941); GEO subseries [GSE92453](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92453) (ex-vivo phenotype array) and [GSE92454](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92454) (in-vitro treatment array), under [GSE92455](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92455). | **Spain (Madrid region).** The methods identify ethics approvals at Centro de Biología Molecular “Severo Ochoa” and Hospital Universitario Puerta de Hierro Majadahonda. The paper explicitly names Hospital La Paz for a separate longitudinal effluent validation series; it does not unambiguously attribute the 16-patient ex-vivo microarray samples to that hospital. Therefore, the specific recruiting center for the array cohort is **not resolved** by the reviewed methods. | **Not reported for the ex-vivo microarray cohort.** The paper reports a 2010 ethics protocol identifier, which does not establish recruitment years. | No overlap with the other cohorts is stated. The ex-vivo array has 17 sample records from 16 patients, so records are not all independent participants. |

## Pairwise assessment

The labels below use the preregistered categories. “Different centers” is based
on the documented institutional sites/countries, not on patient identifiers;
it does not mean overlap was formally ruled out.

| Pair | Label | Basis |
|---|---|---|
| GSE125498 – GSE130888 | **different centers (patient overlap unlikely, not formally confirmed)** | IKEM, Prague, Czech Republic versus Chinese institutions; different countries. |
| GSE125498 – GSE248762 | **different centers (patient overlap unlikely, not formally confirmed)** | IKEM, Prague, Czech Republic versus First Affiliated Hospital of Sun Yat-sen University, Guangzhou, China. |
| GSE125498 – GSE92455 | **different centers (patient overlap unlikely, not formally confirmed)** | Prague, Czech Republic versus Madrid region, Spain. |
| GSE130888 – GSE248762 | **unknown** | Both are associated with Sun Yat-sen University in Guangzhou, but GSE130888 does not assign its samples to a specific hospital; shared institution/city means overlap cannot be dismissed from metadata alone. |
| GSE130888 – GSE92455 | **different centers (patient overlap unlikely, not formally confirmed)** | Chinese institutions versus Madrid region, Spain. |
| GSE248762 – GSE92455 | **different centers (patient overlap unlikely, not formally confirmed)** | Guangzhou, China versus Madrid region, Spain. |

**Conclusion:** Independence is not confirmed across all cohorts. In particular,
GSE130888 versus GSE248762 remains **unknown** because the former’s patient
recruiting hospital is not identified in the reviewed metadata and both studies
are associated with Sun Yat-sen University in Guangzhou. Patient-level matching
is not possible from the GEO sample records reviewed. This uncertainty remains
subject to the preregistered rule.

## Source-access limitations

The full texts of the GSE125498, GSE248762, and GSE92455 papers were reviewed
for methods. GSE125498 methods describe long-dwell effluent collection from 33
stable patients under an IKEM ethics protocol, but give no calendar recruitment
years. The reviewed GSE130888 GEO record and PubMed abstract identify sample
types and author institutions but do not state recruitment years or the
specimen-collection hospital; the full article methods were not accessible from
the reviewed source endpoints. These details are marked unresolved rather than
inferred from publication dates or affiliations.
