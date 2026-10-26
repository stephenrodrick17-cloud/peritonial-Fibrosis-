cat('R version: ', R.version.string, '\n')
pkgs <- c('limma', 'fgsea', 'clusterProfiler', 'msigdbr', 'GSEABase', 'ggplot2', 'data.table')
for (p in pkgs) {
  cat(p, ': ', requireNamespace(p, quietly=TRUE), '\n')
}
