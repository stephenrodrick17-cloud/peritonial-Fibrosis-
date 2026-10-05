suppressPackageStartupMessages(library(edgeR))
dge <- DGEList(counts=matrix(rpois(100, 10), nrow=10), group=rep(1:2, each=5))
rownames(dge) <- paste0("gene_", 1:10)
dge <- calcNormFactors(dge)
fit <- glmQLFit(dge, model.matrix(~group, data=dge$samples))
qlf <- glmQLFTest(fit)
g <- "gene_2"
idx <- which(rownames(fit$counts) == g)
df_denom <- if (length(qlf$df.total) == 1) qlf$df.total else qlf$df.total[idx]
print(df_denom)
print(is.numeric(df_denom))
print(qt(0.975, df = df_denom))
