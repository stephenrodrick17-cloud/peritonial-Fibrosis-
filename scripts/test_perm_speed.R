suppressPackageStartupMessages(library(limma))
expr <- read.csv('results/tables/GSE125498_full_expression_matrix.csv', row.names=1)
meta <- read.csv('results/tables/GSE125498_sample_metadata.csv')
matri <- read.csv('results/revision/matrisome_volcano/GSE125498_matrisome_volcano_source_data.csv')
matri_genes <- intersect(matri$Gene, rownames(expr))
X <- as.matrix(expr[, meta$sample_id])

cat("Starting 10-perm test on all 11,741 genes...\n")
t0 <- Sys.time()
set.seed(42)
for (i in 1:10) {
  grp <- sample(meta$group)
  design <- model.matrix(~ factor(grp, levels=c('SPD', 'LPD')))
  fit <- lmFit(X, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
}
t1 <- Sys.time()
cat("10 perms on all 11,741 genes took:", as.numeric(t1 - t0, units="secs"), "seconds\n")
