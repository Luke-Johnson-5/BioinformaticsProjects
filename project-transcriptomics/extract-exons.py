#!/usr/bin/python
import os,sys
from optparse import OptionParser
import csv
import operator
import collections




def main():
    
    
    # we are interested in MPV17 gene 
    gene_name = "MPV17"
    gene_id = "ENSG00000115204"
    gene_start = 27532360   
    gene_end = 27548547

    in_file = "/groups/bioc6243/data/gtf/Homo_sapiens.GRCh37.75.gtf"

    exon_dict = {}
    with open(in_file, 'r') as FR:
        reader = csv.reader(FR, delimiter='\t', quotechar='|')
        row_count = 0
        for row in reader:
            #chr_id,chr_pos,cigar,seq = row[2], row[3],row[5],row[9]
            if len(row) < 8:
                continue

            if row[2] == "exon" and row[8].find(gene_id) != -1:
                ann_dict = {}
                for ann in row[8].split(";"):
                    if ann == "":
                        continue
                    k, v = ann.strip().split(" ")[0], ann.strip().split(" ")[1]
                    ann_dict[k] = v
                start_pos = row[3]
                trs_id = ann_dict["transcript_id"]
                if trs_id not in exon_dict:
                    exon_dict[trs_id] = []
                exon_dict[trs_id].append(start_pos)

    for trs_id in exon_dict:
        print (trs_id, sorted(exon_dict[trs_id]))







if __name__ == '__main__':
	main()



