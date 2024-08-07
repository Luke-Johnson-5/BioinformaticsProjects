#!/usr/bin/python
import os,sys
from optparse import OptionParser
import csv
import operator
import collections


# given a chromosome position and CIGAR value,
# this function returns a list of splice junction
# positions in the chromosome

def get_splice_positions(chr_pos, cigar):
    pos = int(chr_pos)
    pos_list = []
    for s in cigar.split("N")[:-1]:
        l = 0
        for  v in s.split("M"):
            l += int(v) if v != "" else 0
        pos += l
        pos_list.append(pos)
    return pos_list




def main():
    
    
    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog 1.0")
    parser.add_option("-i","--in_file",action="store",dest="in_file",help="Input TSV file")

    (options,args) = parser.parse_args()
    for item in ([options.in_file]):
        if not (item):
            parser.print_help()
            sys.exit(0)

    in_file = options.in_file

    # we are interested in MPV17 gene 
    gene_name = "MPV17"
    gene_id = "ENSG00000115204"
    gene_start = 27532360   
    gene_end = 27548547

    # this dictionary will hold frequence (occurance count)
    # of splice junction positions
    splice_pos_freq = {}
    
    with open(in_file, 'r') as FR:
        reader = csv.reader(FR, delimiter='\t', quotechar='|')
        row_count = 0
        for row in reader:
            if row[0][0:1] in ["#","@"]:
                continue
            chr_id,chr_pos,cigar,seq = row[2], row[3],row[5],row[9]
            # we populate a list of conditions to select
            # alignment rows that indicate splicing junction
            cond_list = [chr_id == "2", int(chr_pos) >= gene_start, int(chr_pos) <= gene_end]
            cond_list.append(cigar.find("N") != -1)
            cond_list.append(cigar.find("S") == -1)
            if False not in cond_list:
                #print chr_pos, cigar
                # we call the get_splice_positions function 
                # passing chr_pos and CIGAR values to get the 
                # list of splice junction positions the CIGAR indicates
                splice_pos_list = get_splice_positions(chr_pos, cigar)
                for splice_pos in splice_pos_list:
                    if splice_pos not in splice_pos_freq:
                        splice_pos_freq[splice_pos] = 0
                    splice_pos_freq[splice_pos] += 1
            row_count += 1


    # here, we sort the splice_pos_freq dictionary by value
    sorted_splice_pos_freq = collections.OrderedDict(sorted(splice_pos_freq.items(), key=operator.itemgetter(1)))

    print ("slice_position", "frequency")
    for splice_pos in sorted_splice_pos_freq:
        print (splice_pos, sorted_splice_pos_freq[splice_pos])








if __name__ == '__main__':
	main()



