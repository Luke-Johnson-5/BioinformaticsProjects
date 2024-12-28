#!/usr/local/bin/Rscript

#sampleTable <- read.csv("/courses/example-data/csv/sample_table.all.csv",row.names=1) 
sampleTable <- read.csv("/groups/bioc6242/data/csv/sample_table.reduced.csv",row.names=1) 

fileNames <- file.path("/groups/bioc6242/data/bam/", paste0(sampleTable$Run, "_all.bam"))

library("Rsamtools")

bamFiles <- BamFileList(fileNames, yieldSize=2000000)

library("GenomicFeatures")

txdb <- makeTxDbFromGFF("/groups/bioc6242/data/gtf/Homo_sapiens.GRCh37.75.gtf", format="gtf",circ_seqs=character())

ebg <- exonsBy(txdb, by="gene")

library("GenomicAlignments")

library("BiocParallel")

register(SerialParam())

se <- summarizeOverlaps(features=ebg, reads=bamFiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )

colData(se) <- DataFrame(sampleTable)

se$dex <- relevel(se$dex, "untrt")


library("DESeq2")

#dds <- DESeqDataSet(se, design = ~ cell + dex)
dds <- DESeqDataSet(se, design = ~ dex)

nrow(dds)

dds <- dds[ rowSums(counts(dds)) > 1, ]
  
nrow(dds)
 
rld <- rlog(dds, blind=FALSE)
  
head(assay(rld), 3)
 
sampleDists <- dist( t( assay(rld) ) )
  
library("pheatmap")
 
library("RColorBrewer")
 
sampleDistMatrix <- as.matrix( sampleDists )
 
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
  
colnames(sampleDistMatrix) <- NULL
 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
  
plotPCA(rld, intgroup = c("dex", "cell"))




dds <- DESeqDataSet(se, design = ~ cell + dex)

dds <- DESeq(dds)

res <- results(dds)

mcols(res, use.names=TRUE)

summary(res)

res.05 <- results(dds, alpha=.05)

table(res.05$padj < .05)

resLFC1 <- results(dds, lfcThreshold=1)

table(resLFC1$padj < 0.1)
 



#results(dds, contrast=c("cell", "N061011", "N61311"))
results(dds, contrast=c("cell", "N052611", "N080611"))




sum(res$pvalue < 0.05, na.rm=TRUE)

sum(!is.na(res$pvalue))

sum(res$padj < 0.1, na.rm=TRUE)

resSig <- subset(res, padj < 0.1)

head(resSig[ order(resSig$log2FoldChange), ])

head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])



topGene <- rownames(res)[which.min(res$padj)]

plotCounts(dds, gene=topGene, intgroup=c("dex"))



library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res),column="SYMBOL", keytype="ENSEMBL",multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
resOrdered <- res[order(res$padj),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrderedDF, file="results.csv")


