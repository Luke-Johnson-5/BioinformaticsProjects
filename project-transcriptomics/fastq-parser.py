import os,sys
import string
from optparse import OptionParser
from Bio import SeqIO
import re

__version__="1.0"
__status__ = "Dev"



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--fastq_file",action="store",dest="fastq_file",help="Input fastq_file")

    
    (options,args) = parser.parse_args()
    for file in ([options.fastq_file]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    fastq_file = options.fastq_file
        
    n = 0
    for record in SeqIO.parse(fastq_file, "fastq"):
        read_id = record.id
        read_seq = str(record.seq)
        read_qual = record.letter_annotations["phred_quality"]
        print read_id
        print read_seq
        print read_qual
        print "\n"
        n += 1



if __name__ == '__main__':
    main()








