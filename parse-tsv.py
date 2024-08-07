#!/usr/bin/python
import os,sys
from optparse import OptionParser
import csv




###########################################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog 1.0")
    parser.add_option("-i","--in_file",action="store",dest="in_file",help="Input TSV file")


    #If user does not give all arguments, print message and exit
    (options,args) = parser.parse_args()
    for item in ([options.in_file]):
        if not (item):
	    parser.print_help()
            sys.exit(0)

    #Store user arguments in your own variables
    tsv_file = options.in_file
	
    #Parse TSV file
    with open(tsv_file, 'r') as FR:
        reader = csv.reader(FR, delimiter='\t', quotechar='|')
        row_count = 0
        for row in reader:
            if row[0][0:1] != "#" and float(row[13]) <= 0.001:
                print row[1], row[2],row[8], row[13], row[10]
            row_count += 1









if __name__ == '__main__':
    main()



