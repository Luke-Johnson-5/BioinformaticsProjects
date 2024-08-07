#!/bin/sh

module load star/2.7.10a 

genomeDir="/groups/bioc6242/data/indexed-genome/"

STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeDir/Homo_sapiens.GRCh37.75.dna.chromosome.1.fa $genomeDir/Homo_sapiens.GRCh37.75.dna.chromosome.2.fa --runThreadN 8


