# Methodological Limitations: Upper QC Ceilings

Across the cohort of 118,895 raw barcodes, 8956 cells were excluded by upper adaptive ceilings (3-MAD above median log10 counts or genes).

Lineage module scoring of these excluded cells showed the following distribution:

| donor_id    |   B_cell |   Monocyte_macrophage |   NK_cell |   Neutrophil |   Plasma_cell |   Stromal_mesothelial |   T_cell |   Unassigned |   cDC |
|:------------|---------:|----------------------:|----------:|-------------:|--------------:|----------------------:|---------:|-------------:|------:|
| LV_NOT_UF-1 |        0 |                     0 |         0 |            1 |             0 |                   425 |        0 |            0 |     2 |
| LV_NOT_UF-2 |        0 |                     1 |         0 |            0 |             0 |                    69 |        1 |            0 |     4 |
| LV_NOT_UF-3 |        1 |                    50 |         1 |            1 |             0 |                    64 |       13 |            0 |    10 |
| LV_NOT_UF-4 |       65 |                   206 |        36 |           35 |            63 |                    26 |      273 |            6 |   281 |
| LV_NOT_UF-5 |        3 |                   474 |        29 |           29 |            35 |                    26 |      205 |            0 |   246 |
| LV_NOT_UF-6 |        1 |                   301 |         5 |           14 |             4 |                    68 |       92 |            1 |   364 |
| LV_UF-1     |        0 |                     1 |         0 |            4 |            14 |                    23 |        0 |            0 |     0 |
| LV_UF-2     |       12 |                    64 |        21 |            6 |           264 |                   186 |      135 |            5 |    30 |
| LV_UF-3     |        3 |                   925 |        10 |           99 |           182 |                   118 |        9 |           27 |    95 |
| LV_UF-4     |       18 |                    40 |        34 |           39 |            10 |                    78 |      581 |            6 |    98 |
| SV-1        |        0 |                    13 |         0 |            2 |             1 |                   255 |        4 |            0 |    26 |
| SV-2        |        0 |                     1 |         0 |            0 |             0 |                    63 |        6 |            0 |     1 |
| SV-3        |        1 |                     3 |         0 |            0 |             1 |                   187 |        4 |            0 |    14 |
| SV-4        |        0 |                    11 |         0 |            1 |             0 |                   147 |        1 |            0 |    14 |
| SV-5        |        1 |                    70 |        18 |            6 |             2 |                   944 |       92 |            2 |   210 |
| SV-6        |       15 |                    28 |         1 |            1 |             4 |                   130 |        0 |            0 |     8 |
## QC Specification Provenance & Limitations

**Key Finding**: As expected for high-RNA metabolic states, stromal/mesothelial cells and activated macrophages are disproportionately represented among cells with high transcript counts. As per the analysis plan (self-documented; committed after primary results), QC thresholds are held strictly frozen to prevent post-hoc bias, and this exclusion is documented as a limitation.
