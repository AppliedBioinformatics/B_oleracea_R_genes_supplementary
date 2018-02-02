# Compares the variable/core genes with in cluster/outside of cluster genes
cat('From Python:\nvariable, in cluster 87\ncore, in cluster 219\nvariable, not in cluster 18\ncore, not in cluster 111\n')

M <- as.table(rbind(c(18, 111), c(87, 219)))
dimnames(M) <- list(Status= c('not in cluster', 'in cluster'), PAV=c('variable','core'))
M
chisq.test(M)
chisq.test(M, correct=F)
