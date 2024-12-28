#!/bin/sh

module load msgfplus

in_dir="/groups/bioc6243/data/msgfplus/"
out_dir="/groups/bioc6243/data/tmp/"

spectraFile="$in_dir/OR20070924_S_mix7_11.mzXML"
targetFile="$in_dir/18mix.fasta"
modFile="$in_dir/Mods.txt"

mzIdFile="$out_dir/$USER-output.mzid"
tsvFile="$out_dir/$USER-output.tsv"
paramSet="-t 10ppm -ti -1,2 -ntt 2 -tda 1 -inst 1"

module load msgfplus

echo $MSGFPlus

java -Xmx3500M -jar $MSGFPlus -s $spectraFile -d $targetFile $paramSet -mod $modFile -o $mzIdFile 

java -Xmx3500M -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -i $mzIdFile -o $tsvFile



