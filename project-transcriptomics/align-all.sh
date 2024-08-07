binDir="/softwares/STAR-2.5.2b/bin/Linux_x86_64"
genomeDir="/groups/bioc6243/data/indexed-genome/"
samDir="/groups/bioc6243/data/sam"
bamDir="/groups/bioc6243/data/bam"
fastqDir="/groups/bioc6243/data/fastq/"

module load star/2.7.10a
module load samtools/1.18

for runId in SRR1039512 SRR1039513 SRR1039516 SRR1039517 ; do
#for runId in SRR1039517 ; do
    STAR --genomeDir $genomeDir --readFilesIn $fastqDir/$runId""_1.fastq $fastqDir/$runId""_2.fastq  --runThreadN 8 --outFileNamePrefix $samDir/$runId. > align.$runId.log
    samtools  view -S -b $samDir/$runId".Aligned.out.sam" > $bamDir/$runId"_all.bam"
done


