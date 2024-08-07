import os,sys
import string
from optparse import OptionParser
from Bio import SeqIO

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

    #input fastq file  
    fastq_file = options.fastq_file
	
    #Our read sequences are 63 bases long
    read_len = 63

    #This list will contain sum of quality fo each position
    qual_bin = []
	
    #Here we Intialize qual_bin to contain 0
    for i in xrange (0,read_len):
        qual_bin.append(0)
    
    #Now we parse our fastq file to populate qual_bin
    for record in SeqIO.parse(fastq_file, "fastq"):	#For each read record
        for i in xrange (0,read_len):
            qual_bin[i] += int(record.letter_annotations["phred_quality"][i]) 

    #Now lets print out sum of quality for each position
    for i in xrange (0,read_len):
        print i, qual_bin[i]


if __name__ == '__main__':
    main()








