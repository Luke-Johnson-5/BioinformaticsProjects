#!/bin/sh

genomeDir="/groups/bioc6243/data/indexed-genome/"

fastqFileOne="/groups/bioc6243/data/fastq/SRR1039508_1.fastq"
fastqFileTwo="/groups/bioc6243/data/fastq/SRR1039508_2.fastq"
outPrefix="/groups/bioc6243/data/tmp/$USER-aln-SRR1039508."
logFile="/groups/bioc6243/data/tmp/$USER-log-SRR1039508.log"

module load star/2.7.10a 

STAR --genomeDir $genomeDir --readFilesIn $fastqFileOne $fastqFileTwo --runThreadN 1 --outFileNamePrefix $outPrefix  > $logFile




