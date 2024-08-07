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


def dump_formatted_cds (trs_id, pep_id, cds_seq, pep_seq):
        
    cds_pos = 0
    pep_pos = 0
    buffer = ">%s (%s)" % (trs_id, pep_id)
    while cds_pos < len(cds_seq) - 2:
        codon = cds_seq[cds_pos:cds_pos+3]
        aa = pep_seq[pep_pos] if pep_pos < len(pep_seq) else "*"
        nl = "\n" if pep_pos  % 10 == 0 else ""
        buffer += "%s%s(%s) " % (nl,codon,aa)
        cds_pos += 3
        pep_pos += 1
    return buffer.strip()




def main():
        
    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-t", "--trsid", action = "store", dest = "trsid", help = "transcsript ID")

    (options,args) = parser.parse_args()
    for file in ([options.trsid]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    trs_id = options.trsid
    
    global data_dir
    
    data_dir = "/groups/bioc6243/data/"
    cds_seq_file = data_dir + "fasta/Homo_sapiens.GRCh37.75.cds.all.fa"
    pep_seq_file = data_dir + "fasta/Homo_sapiens.GRCh37.75.pep.all.fa"
    #Load cds sequences
    cds_seq_hash = {}
    
    for record in SeqIO.parse(cds_seq_file, "fasta"):
        t_id = record.id.split(".")[0]
        cds_seq_hash[t_id] = record.seq.upper()
    
    pep_seq_hash = {}
    trsid2pepid = {}
    for record in SeqIO.parse(pep_seq_file, "fasta"):
        desc = record.description
        t_id = desc.split("transcript:")[1].split()[0].split(".")[0]
        pep_seq_hash[t_id] = str(record.seq.upper())
        trsid2pepid[t_id] =  record.id.split(".")[0]
    
    x = dump_formatted_cds (trs_id, trsid2pepid[trs_id], cds_seq_hash[trs_id], pep_seq_hash[trs_id])
    print (x)







if __name__ == '__main__':
	main()


