#!/usr/bin/python

import sys
import os
import math
import getopt

def usage():
    print """
    HLAminerToTsv -l <HLAminer_log> -c <HLAminer_csv> -o <output-dir>

    Requested Parameters:
        -l/--hlaminer-log      [string:    path to HLAminer_HPRA.log]
        -c/--hlaminer-csv      [string:    path to HLAminer_HPRA.csv]

    Optional Parameters:
        -o/--output-dir        [string:    path to output dir.  Default: ./  ]

    Version:                    1.0.0
          """

#parameters
input_log = ''
input_csv = ''
output_dir= ''

def get_default_output():
    global output_dir
    output_dir = "./" 

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hl:c:o:",["help",
                                                    "hlaminer-log=",
                                                    "hlaminer-csv=",
                                                    "output_dir="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        elif opt in ("-l","--hlaminer-log"):
            global input_log
            input_log = arg
        elif opt in ("-c","--hlaminer-csv"):
            global input_csv
            input_csv = arg
        elif opt in ("-o","--output_dir"):
            global output_dir 
            output_dir = arg

hla_allele = []

def get_hla_alleles(infile):
    f=open(infile,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            line=line.strip()
            tmp = line.split(",")
            if len(tmp)==7 and tmp[0]!="Gene":
                if tmp[3][len(tmp[3])-1] in ["P", "S", "L", "C", "A", "Q"]:
                   tmp[3]=tmp[3][0:len(tmp[3])-1]
                hla_allele.append(tmp)
    f.close();

top2 = []
top1 = []

def get_tops(infile):

    f=open(infile,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            line=line.strip()
            tmp = line.split(",")
            if len(tmp)==4 and tmp[0]!="Allele":
                hlahla="HLA-"+tmp[0][0]+tmp[0][2:len(tmp[0])];
                if hlahla[len(hlahla)-1] in ["P", "S", "L", "C", "A", "Q"]:
                    hlahla=hlahla[0:len(hlahla)-1]
                    top2.append(hlahla)
                    if openopen==1:
                        top1.append(hlahla)
            if line.count("Prediction")==1:
                openopen=1
            else:
                openopen=0
    f.close();

allele_records = {}

def get_allele_records():
    for t in range(len(hla_allele)):
        if t!=0:
            allele_records["HLA-" + hla_allele[t][0] + hla_allele[t][3][2:]]= hla_allele[t][4]+"\t"+ hla_allele[t][5]+"\t"+ hla_allele[t][6]

def getType(name):
    type="IN_LOG"
    if name in top2:
        type="IN_CSV"
    if name in top1:
        type="Top_1"
    return type

def print_result(output_dir):
    fo = open(output_dir+"/HLAminer_alleles.tsv", "wb")
    for t in range(len(hla_allele)):
        if t!=0:
            name="HLA-" + hla_allele[t][0] + hla_allele[t][3][2:]
            fo.write("%s\t%s\t%s\n" % (name,getType(name),allele_records[name]))
    fo.close()
    
def main(argv):
    get_default_output()
    getParameters(argv[1:])
    if input_log=='' or input_csv=='': 
        usage()
        exit(1);
    get_hla_alleles(input_log)
    get_tops(input_csv)
    get_allele_records()
    print_result(output_dir) 

if __name__ == '__main__':
    sys.exit(main(sys.argv))
