import os,sys
import string
from optparse import OptionParser
import glob
import json
import csv
from Bio import SeqIO
from Bio.Seq import Seq


__status__ = "Dev"
__version__ = "4.0"



def load_map_dict (map_file):
    
    map_dict = {}
    with open(map_file, 'r') as FR:
        for line in FR:
            k = line.strip().split(",")[0]
            v = line.strip().split(",")[1]
            if k not in map_dict:
                map_dict[k] = []
            map_dict[k].append(v)


    return map_dict



def load_counts():

    input_dir = data_dir + "/vcf/tcga/"
    pattern = input_dir + "*.vcf"
    file_list = glob.glob(pattern)
    count_two = {}
    count_one = {}
    file_count = 1
    for vcf_file in file_list:
        cancer_id = vcf_file.split("/")[-1].split("-")[0]
        if cancer_id not in count_one:
            count_one[cancer_id] = 0
        count_one[cancer_id] += 1
        file_count += 1
        with open(vcf_file, 'r') as FR:
            row_count = 0
            for line in FR:
                if line[0] == "#":
                    continue
                row_count += 1
                row = line.strip().split("\t")
                chr_id = row[0]
                chr_pos = int(row[1])
                ref_nt = row[3]
                alt_nt = row[4]
                pf = row[6]
                if pf != "PASS":
                    continue
                combo_id = "%s:%s:%s>%s" % (chr_id,chr_pos,ref_nt, alt_nt)
                if combo_id not in count_two:
                    count_two[combo_id] = {}
                if cancer_id not in count_two[combo_id]:
                    count_two[combo_id][cancer_id] = 0
                count_two[combo_id][cancer_id] += 1

    return count_one, count_two



##############################################
def main():
    
    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-c", "--cutoff", action = "store", dest = "cutoff", help = "Frequency cut off")
    parser.add_option("-r", "--rank", action = "store", dest = "rank", help = "Frequency rank")


    (options,args) = parser.parse_args()
    for file in ([options.cutoff]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    cut_off = float(options.cutoff)


    global data_dir
    data_dir = "/groups/bioc6243/data/"


    map_dict = {}
    map_dict["tcganames"] = load_map_dict(data_dir + "csv/tcganames.csv")


    count_one,count_two = load_counts()
    freq_dict = {}
    nd_dict = {}
    for combo_id in count_two:
        for cancer_id in count_two[combo_id]:
            n, d = count_two[combo_id][cancer_id], count_one[cancer_id]
            freq = "%5.3f" % (float(n)/float(d))
            if cancer_id not in freq_dict:
                freq_dict[cancer_id] = {}
                nd_dict[cancer_id] = {}
            freq_dict[cancer_id][combo_id] = float(freq)
            nd_dict[cancer_id][combo_id] = "%s,%s" % (n, d)
   
    print ("%s,%s,%s,%s,%s,%s" % ("frequency","subjects_positive","subjects_tested","mutation","cancer_id","cancer_name")) 
    for cancer_id in freq_dict:
        cancer_name = map_dict["tcganames"][cancer_id][0]
        n = 0
        for combo_id, freq in sorted(freq_dict[cancer_id].items(), key=lambda x: x[1], reverse=True):
            n += 1
            cond_list = [freq >= cut_off]
            if options.rank != None:
                rank = int(options.rank)
                cond_list.append(n == rank)
            if False not in cond_list:
                nd = nd_dict[cancer_id][combo_id]
                print ("%s,%s,%s,%s,%s" % (freq,nd,combo_id,cancer_id,cancer_name))





if __name__ == '__main__':
	main()


