sampleTable <- read.csv("/groups/bioc6243/data/csv/sample_table.reduced.csv",row.names=1)
fileNames <- file.path("/groups/bioc6243/data/bam/", paste0(sampleTable$Run, ".Aligned.out.bam"))

library("Rsamtools")
bamFiles <- BamFileList(fileNames, yieldSize=2000000)

library("GenomicFeatures")
txdb <- makeTxDbFromGFF("/groups/bioc6243/data/gtf/Homo_sapiens.GRCh37.75.gtf", format="gtf",circ_seqs=character())
ebg <- exonsBy(txdb, by="gene")

library("GenomicAlignments")
library("BiocParallel")
register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamFiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )
colData(se) <- DataFrame(sampleTable)
se$dex <- relevel(factor(se$dex), "untrt")



