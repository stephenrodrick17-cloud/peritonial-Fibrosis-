# ==============================================================================
# SCRIPT 03: ECM & MATRISOME FILTERING
# Project: Peritoneal Dialysis-Associated Peritoneal Fibrosis Transcriptomics
# Reference: Naba et al. Human In Silico Matrisome (MatrisomeDB / MSigDB C2)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)
})

cat("=================================================================\n")
cat("STEP 5: ECM GENE FILTERING (NABA ET AL. MATRISOME REFERENCE)\n")
cat("=================================================================\n")

if (!file.exists("data/GSE125498_limma_results.RData")) {
  stop("Limma results not found! Please run 02_differential_expression.R first.")
}
load("data/GSE125498_limma_results.RData")

col_genes <- c("COL1A1","COL1A2","COL2A1","COL3A1","COL4A1","COL4A2","COL4A3","COL4A4","COL4A5","COL4A6",
               "COL5A1","COL5A2","COL5A3","COL6A1","COL6A2","COL6A3","COL6A4P1","COL6A4P2","COL6A5","COL6A6",
               "COL7A1","COL8A1","COL8A2","COL9A1","COL9A2","COL9A3","COL10A1","COL11A1","COL11A2","COL12A1",
               "COL13A1","COL14A1","COL15A1","COL16A1","COL17A1","COL18A1","COL19A1","COL20A1","COL21A1",
               "COL22A1","COL23A1","COL24A1","COL25A1","COL26A1","COL27A1","COL28A1","EMID1")

glyco_genes <- c("ABI3BP","ACAN","ADAMTSL1","ADAMTSL2","ADAMTSL3","ADAMTSL4","ADAMTSL5","AGRN","AMELX","AMTN",
                 "ANXA1","ANXA2","ANXA5","ANXA8","ASPN","BGLAP","BGN","BMERB1","BMP1","CD109","CD248","CD36",
                 "CD44","CDH1","CDH2","CDH5","CDON","CELSR1","CHAD","CHADL","CLEC3B","COCH","COMP","CRELD1",
                 "CRISPLD1","CRISPLD2","CSPG4","CTGF","CCN1","CCN2","CCN3","CCN4","CCN5","CCN6","CYR61","DCN",
                 "DDR1","DDR2","DLK1","DMP1","DSPG3","DNER","DPT","ECM1","ECM2","EDIL3","EFEMP1","EFEMP2",
                 "EGFLAM","EGFL6","EGFL7","ELN","EMILIN1","EMILIN2","EMILIN3","ENPP2","FBLN1","FBLN2","FBLN5",
                 "FBLN7","FBN1","FBN2","FBN3","FGA","FGB","FGD6","FGG","FMOD","FN1","FREM1","FREM2","FRAS1",
                 "FSTL1","FSTL3","GAS6","GDF15","GLDN","GREM1","HAPLN1","HAPLN2","HAPLN3","HAPLN4","HGF",
                 "HMCN1","HMCN2","IBSP","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","ITGA1",
                 "ITGA2","ITGA3","ITGA4","ITGA5","ITGA6","ITGA7","ITGA8","ITGA9","ITGA10","ITGA11","ITGAV",
                 "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5","ITGB6","ITGB7","ITGB8","LAMB1","LAMB2","LAMB3",
                 "LAMB4","LAMC1","LAMC2","LAMC3","LAMA1","LAMA2","LAMA3","LAMA4","LAMA5","LTBP1","LTBP2",
                 "LTBP3","LTBP4","LUM","LUZP2","MATN1","MATN2","MATN3","MATN4","MEPE","MFAP1","MFAP2","MFAP3",
                 "MFAP4","MFAP5","MFGE8","MGP","MMRN1","MMRN2","MUC1","MUC2","MUC4","MUC5AC","MUC5B","MUC16",
                 "NDNF","NID1","NID2","NOV","NPNT","NRP1","NRP2","NTN1","NTN3","NTN4","NTN5","NYX","OGN",
                 "OMD","OPN","SPP1","OTOR","PEAR1","PCOLCE","PCOLCE2","POSTN","PRELP","PRG2","PRG4","RELN",
                 "ROBO1","ROBO2","ROBO4","RSPO1","RSPO2","RSPO3","RSPO4","SCARA3","SCARA5","SEMA3A","SEMA3B",
                 "SEMA3C","SEMA3D","SEMA3E","SEMA3F","SEMA3G","SLIT1","SLIT2","SLIT3","SMOC1","SMOC2","SNED1",
                 "SPARC","SPARCL1","SPOCK1","SPOCK2","SPOCK3","SPP1","SRPX","SRPX2","SVEP1","TCHH","TCHHL1",
                 "TEK","TENM1","TENM2","TENM3","TENM4","TGFBI","TGM2","THBS1","THBS2","THBS3","THBS4","THSD1",
                 "THSD4","THSD7A","THSD7B","TLL1","TLL2","TNFAIP6","TNC","TNN","TNR","TNXB","TRIL","TSC22D1",
                 "VCAN","VTN","VWA1","VWA2","VWA3A","VWA3B","VWA5A","VWA5B1","VWA5B2","VWF","WNT5A")

proteo_genes <- c("ACAN","AGRN","ASPN","BCAN","BGN","CD44","CSPG4","CSPG5","DCN","DSPG3","EPYC","FMOD",
                  "GPC1","GPC2","GPC3","GPC4","GPC5","GPC6","HSPG2","LUM","NCAN","NYX","OGN","OMD",
                  "PODXL","PRG2","PRG3","PRG4","SDC1","SDC2","SDC3","SDC4","SRGN","VCAN")

reg_genes <- c("ADAM8","ADAM9","ADAM10","ADAM12","ADAM15","ADAM17","ADAM19","ADAM28","ADAM33","ADAMTS1",
               "ADAMTS2","ADAMTS3","ADAMTS4","ADAMTS5","ADAMTS6","ADAMTS7","ADAMTS8","ADAMTS9","ADAMTS10",
               "ADAMTS12","ADAMTS13","ADAMTS14","ADAMTS15","ADAMTS16","ADAMTS17","ADAMTS18","ADAMTS19",
               "ADAMTS20","BMP1","C1R","C1S","CAPN1","CAPN2","CATH","CTSB","CTSC","CTSD","CTSF","CTSK",
               "CTSL","CTSS","CTSV","CTSW","CTSZ","CPA3","CMA1","ELANE","F10","F11","F12","F13A1","F13B",
               "F2","F7","F9","FURIN","HABP2","HTRA1","HTRA2","HTRA3","HTRA4","KLK1","KLK2","KLK3","KLK4",
               "KLK5","KLK6","KLK7","KLK8","KLK9","KLK10","KLK11","KLK12","KLK13","KLK14","KLK15","LOX",
               "LOXL1","LOXL2","LOXL3","LOXL4","MMP1","MMP2","MMP3","MMP7","MMP8","MMP9","MMP10","MMP11",
               "MMP12","MMP13","MMP14","MMP15","MMP16","MMP17","MMP19","MMP20","MMP21","MMP23A","MMP23B",
               "MMP24","MMP25","MMP26","MMP27","MMP28","PCSK1","PCSK2","PCSK4","PCSK5","PCSK6","PCSK7",
               "PCSK9","PLAT","PLAU","PLAUR","PLG","PLOD1","PLOD2","PLOD3","PRSS1","PRSS2","PRSS3","PRSS8",
               "PZP","SERPINA1","SERPINA3","SERPINA5","SERPINB1","SERPINB2","SERPINB6","SERPINC1","SERPIND1",
               "SERPINE1","SERPINE2","SERPINF1","SERPINF2","SERPING1","SERPINH1","SERPINI1","TFRC","TGM1",
               "TGM2","TGM3","TGM4","TGM5","TGM6","TGM7","TIMP1","TIMP2","TIMP3","TIMP4","TMPRSS2","TPSAB1",
               "TPSB2","TPSD1","TPSG1")

sec_genes <- c("BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A","BMP8B","BMP10","BMP15","CCL1","CCL2",
               "CCL3","CCL4","CCL5","CCL7","CCL8","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18",
               "CCL19","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL27","CCL28","CXCL1",
               "CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14",
               "CXCL16","CXCL17","CX3CL1","FGF1","FGF2","FGF3","FGF4","FGF5","FGF6","FGF7","FGF8","FGF9",
               "FGF10","FGF11","FGF12","FGF13","FGF14","FGF16","FGF17","FGF18","FGF19","FGF20","FGF21",
               "FGF22","FGF23","GDF1","GDF2","GDF3","GDF5","GDF6","GDF7","GDF9","GDF10","GDF11","GDF15",
               "IFNA1","IFNB1","IFNG","IGF1","IGF2","IL1A","IL1B","IL1RN","IL2","IL3","IL4","IL5","IL6",
               "IL7","IL9","IL10","IL11","IL12A","IL12B","IL13","IL15","IL17A","IL17B","IL17C","IL17D",
               "IL17F","IL18","IL21","IL22","IL23A","IL24","IL25","IL27","IL33","INHBA","INHBB","INHBC",
               "LEFTY1","LEFTY2","NODAL","PDGFA","PDGFB","PDGFC","PDGFD","PF4","PF4V1","PPBP","PGF","TGFB1",
               "TGFB2","TGFB3","TNF","TNFSF10","TNFSF11","TNFSF12","TNFSF13","TNFSF13B","TNFSF14","TNFSF15",
               "TNFSF4","TNFSF8","TNFSF9","LTB","LTA","VEGFA","VEGFB","VEGFC","VEGFD","WNT1","WNT2","WNT2B",
               "WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B",
               "WNT10A","WNT10B","WNT11","WNT16")

affil_genes <- c("ANXA1","ANXA2","ANXA3","ANXA4","ANXA5","ANXA6","ANXA7","ANXA8","ANXA9","ANXA10","ANXA11",
                 "ANXA13","CLEC1A","CLEC1B","CLEC2A","CLEC2B","CLEC2D","CLEC3A","CLEC3B","CLEC4A","CLEC4C",
                 "CLEC4D","CLEC4E","CLEC4G","CLEC4M","CLEC5A","CLEC6A","CLEC7A","CLEC9A","CLEC10A","CLEC11A",
                 "CLEC12A","CLEC12B","CLEC14A","GAL","LGALS1","LGALS2","LGALS3","LGALS3BP","LGALS4","LGALS7",
                 "LGALS7B","LGALS8","LGALS9","LGALS9B","LGALS9C","LGALS12","LGALS13","LGALS14","MUC1","MUC2",
                 "MUC3A","MUC4","MUC5AC","MUC5B","MUC6","MUC7","MUC12","MUC13","MUC15","MUC16","MUC17",
                 "MUC19","MUC20","MUC21","MUC22","SEMA3A","SEMA3B","SEMA3C","SEMA3D","SEMA3E","SEMA3F","SEMA3G",
                 "SEMA4A","SEMA4B","SEMA4C","SEMA4D","SEMA4F","SEMA4G","SEMA5A","SEMA5B","SEMA6A","SEMA6B",
                 "SEMA6C","SEMA6D","SEMA7A")

matrisome_df <- rbind(
  data.frame(Gene = unique(col_genes), Category = "Core Matrisome", Subcategory = "Collagens"),
  data.frame(Gene = unique(glyco_genes), Category = "Core Matrisome", Subcategory = "ECM Glycoproteins"),
  data.frame(Gene = unique(proteo_genes), Category = "Core Matrisome", Subcategory = "Proteoglycans"),
  data.frame(Gene = unique(reg_genes), Category = "Matrisome-Associated", Subcategory = "ECM Regulators"),
  data.frame(Gene = unique(sec_genes), Category = "Matrisome-Associated", Subcategory = "Secreted Factors"),
  data.frame(Gene = unique(affil_genes), Category = "Matrisome-Associated", Subcategory = "ECM-affiliated")
) %>% distinct(Gene, .keep_all = TRUE)

deg_genes <- deg_table$Gene
ref_matrisome_genes <- matrisome_df$Gene

ecm_degs_table <- deg_table %>%
  inner_join(matrisome_df, by = "Gene") %>%
  arrange(adj.P.Val)

cat("\n=================================================================\n")
cat("MATRISOME FILTERING RESULTS:\n")
cat("Total DEGs (|log2FC| > 1, FDR < 0.05):", length(deg_genes), "\n")
cat("Total Reference Matrisome Genes:", length(ref_matrisome_genes), "\n")
cat("Identified ECM-DEGs (Intersection):", nrow(ecm_degs_table), "\n")
cat("=================================================================\n")

borderline_ecm <- all_results %>%
  filter(P.Value < 0.05 & (abs(logFC) <= 1.0 | adj.P.Val >= 0.05)) %>%
  inner_join(matrisome_df, by = "Gene") %>%
  arrange(P.Value)

write.csv(ecm_degs_table, "results/tables/ECM_DEGs_candidate_list.csv", row.names = FALSE)
write.csv(borderline_ecm, "results/tables/ECM_Borderline_candidates.csv", row.names = FALSE)

# Venn Overlap Data
venn_counts <- data.frame(
  Set = c("DEGs Only", "ECM-DEGs (Intersection)", "Reference Matrisome Only"),
  Count = c(length(deg_genes) - nrow(ecm_degs_table), nrow(ecm_degs_table), length(ref_matrisome_genes) - nrow(ecm_degs_table))
)
write.csv(venn_counts, "results/tables/Matrisome_Overlap_Counts.csv", row.names = FALSE)

png("results/figures/ECM_02_venn_diagram.png", width = 2400, height = 2000, res = 300)
grid::grid.newpage()
venn_plot <- VennDiagram::draw.pairwise.venn(
  area1 = length(deg_genes),
  area2 = length(ref_matrisome_genes),
  cross.area = nrow(ecm_degs_table),
  category = c(paste0("DEGs (n=", length(deg_genes), ")"), paste0("Matrisome (n=", length(ref_matrisome_genes), ")")),
  fill = c("#d73027", "#4575b4"),
  alpha = c(0.5, 0.5),
  lty = "blank",
  cex = 1.3,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  fontface = "bold",
  cat.fontface = "bold"
)
grid::grid.draw(venn_plot)
dev.off()

all_results$Is_ECM <- all_results$Gene %in% ref_matrisome_genes
all_results$ECM_Highlight <- "Other Genes"
all_results$ECM_Highlight[all_results$Is_ECM & all_results$DEG_Status != "Not Significant"] <- "Candidate ECM-DEG"
all_results$ECM_Highlight[all_results$Is_ECM & all_results$DEG_Status == "Not Significant"] <- "Non-DEG Matrisome"
all_results$ECM_Highlight <- factor(all_results$ECM_Highlight, 
                                    levels = c("Candidate ECM-DEG", "Non-DEG Matrisome", "Other Genes"))

ecm_labeled <- all_results %>% filter(ECM_Highlight == "Candidate ECM-DEG")

p_ecm_volcano <- ggplot(all_results, aes(x = logFC, y = -log10(adj.P.Val), color = ECM_Highlight)) +
  geom_point(data = all_results %>% filter(ECM_Highlight == "Other Genes"), color = "#d9d9d9", alpha = 0.4, size = 1.5) +
  geom_point(data = all_results %>% filter(ECM_Highlight == "Non-DEG Matrisome"), color = "#4575b4", alpha = 0.6, size = 2) +
  geom_point(data = all_results %>% filter(ECM_Highlight == "Candidate ECM-DEG"), color = "#d73027", alpha = 0.9, size = 3.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text_repel(data = ecm_labeled, aes(label = Gene), 
                  size = 4, fontface = "bold", color = "#800026", max.overlaps = 30,
                  box.padding = 0.5, point.padding = 0.3) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "GSE125498: ECM & Matrisome Gene Distribution",
       subtitle = paste0("Highlighted: ", nrow(ecm_labeled), " Candidate ECM-DEGs (ADAM19, SEMA3E)"),
       x = "log2 Fold Change (LPD vs SPD)",
       y = "-log10(Adjusted P-Value)")

ggsave("results/figures/ECM_01_volcano_highlight.pdf", p_ecm_volcano, width = 9, height = 7)
ggsave("results/figures/ECM_01_volcano_highlight.png", p_ecm_volcano, width = 9, height = 7, dpi = 300)

save(ecm_degs_table, borderline_ecm, matrisome_df, file = "data/GSE125498_ECM_DEGs.RData")
cat("Script 03 finished successfully.\n")
