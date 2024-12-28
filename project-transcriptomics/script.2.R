library("DESeq2")

#dds <- DESeqDataSet(se, design = ~ dex)
dds <- DESeqDataSet(se, design = ~ cell + dex)

dds <- dds[ rowSums(counts(dds)) > 1, ]
rld <- rlog(dds, blind=FALSE)



