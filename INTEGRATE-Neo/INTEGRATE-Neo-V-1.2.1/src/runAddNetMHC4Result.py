#!/usr/bin/python

import sys
import os
import math
import getopt

def usage():
    print """
    runAddNetMHC4Result -b <original-bedpe> -h <hla-alleles-tsv> -f <netMHC4.0.out.append.txt> 
                        -m <min-score> -o <output-bedpe>

    Requested Parameters:
        -b/--bedpe-in          [string:    path to input bedpe             ]
        -a/--hla-allele        [string:    path to hla allele tsv file     ]
        -f/--netMHC4-result    [string:    path to netMHC4.0.out.append.txt]

    Optional Parameters:
        -m/--affinity-score    [value:     Default: 500.             min affinity score   ]
        -o/--output-bedpe      [string:    Default: result.bedpe.    path to outpout bedpe]        

    Version:                    1.0.0
          """

#parameters
inFile = ''               #input bedpe
inHLA = ''                #input hla alleles
inNetResult = ''          #netMHC4.0.out.append.txt
min_score=500
outFile = 'result.bedpe'  #result

net_num=3
hla_col_num=0


def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hb:a:f:m:o:",["help",
                                                   "bedpe-in=",
                                                   "hla-allele=",
                                                   "netMHC4-result=",
                                                   "affinity-score=",
                                                   "output-bedpe="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        elif opt in ("-b", "--bedpe-in"):
            global inFile
            inFile= arg
        elif opt in ("-a", "--hla-allele"):
            global inHLA
            inHLA = arg
        elif opt in ("-f","--netMHC4-result"):
            global inNetResult
            inNetResult = arg
        elif opt in ("-m","--affinity-score"):
            global min_score
            min_score = int(arg)
        elif opt in ("-o","--output-bedpe"):
            global outFile
            outFile = arg


records_hla = {}

def getHLARecord():
    f = open(inHLA,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            tmp=line.split("\t")
            tmp[len(tmp)-1]=tmp[len(tmp)-1][0:len(tmp[len(tmp)-1])-1]
            records_hla[tmp[0]]=tmp
            global hla_col_num
            hla_col_num=len(tmp)
    f.close()

records_bedpe = {}

def getBedRecord():
    f = open(inFile, "r")
    index=0
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            index=index+1
            records_bedpe[index]=line
    f.close()

atg = []

def assign():
    back_num=0
    fin = open(inNetResult,"r")
    while True:
        line=fin.readline()
        if line=="":
            break
        else:
            line=line.strip()
            if line.count("Protein")==1:
                protein=line.split(".")[0].split(" ")[1]
                for x in range(back_num):
                    full_rec=records_bedpe[int(protein)]
                    rc=full_rec.split("\t")
                    for t in range(len(rc)):
                        atg[len(atg)-1-x].append(rc[t])
                    atg[len(atg)-1-x].append(protein)
                back_num=0
            if line.count("HLA")>0 and line.count("NetMHC")==0 and line.count("peptide")==0:
                tmp = line.split()
                if len(tmp)>=14 and float(tmp[12])<=float(min_score):
                    hlahla=tmp[1][0:7]+":"+tmp[1][7:9]
                    ttt = []
                    ttt.append(tmp[1])
                    ttt.append(tmp[2])
                    ttt.append(tmp[12])
                    rc=records_hla[hlahla]
                    for t in range(len(rc)):
                        ttt.append(rc[t])
                    atg.append(ttt)
                    back_num=back_num+1
    fin.close()

    #for x in range(len(atg)):
    #    print '\t'.join(atg[x])

def getOutput():
    global atg
    atg_len=0
    if len(atg)>0:
        atg_len=len(atg[0])
    atg = sorted(atg, key = lambda x: (int(x[atg_len-1]),x[0],float(x[2])))
    pre_id='0'
    pre_hla='hla'
    s=net_num+hla_col_num
    atg2 = []
    for x in range(len(atg)):
        rc=atg[x]
        if rc[atg_len-1]!=pre_id or rc[0]!=pre_hla:
            aa = []
            for kk in range(10):
                aa.append(rc[s+kk])
            if rc[s+10][len(rc[s+10])-1]=='\n':
                aa.append(rc[s+10][0:len(rc[s+10])-1])
            else:
                aa.append(rc[s+10])
            aa.append(rc[1])
            aa.append(rc[2])
            for i in range(hla_col_num):
                aa.append(rc[i+3])
            #print '\t'.join(aa)
            atg2.append(aa)
        pre_id=rc[atg_len-1]
        pre_hla=rc[0]
    f=open(outFile,"w")
    atg2 = sorted(atg2, key = lambda x: (float(x[10])),reverse=True)
    for xx in range(len(atg2)):
        rc=atg2[xx]
        for yy in range(len(rc)-1):
            f.write("%s\t" % (rc[yy]))
        f.write("%s\n" % rc[len(rc)-1])
    f.close()        


def main(argv):
    getParameters(argv[1:])
    if inFile=='' or inHLA=='' or inNetResult=='':
        usage()
        exit(1)
    getHLARecord()
    getBedRecord()
    assign()
    getOutput()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
