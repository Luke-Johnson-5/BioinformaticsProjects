library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds <- DESeq(dds)

res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
 


