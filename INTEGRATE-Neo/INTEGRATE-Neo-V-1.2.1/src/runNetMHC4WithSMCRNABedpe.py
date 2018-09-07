#!/usr/bin/python

import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT

def usage():
    print """
    python runNetMHC4WithSMCRNABedpe.py -a <hla-allele-file> -f <fusion-bedpe-file> -p <peptide-lengths>
                         -o <output-dir> -n <path-to-netMHC4.0> -v <available-alleles-file>
                         -x "<netMHC4.0-options>" -k

    Requested Parameters:
        -a/--hla-allele        [string:    path to HLA.allele.tsv         ]
        -f/--fusion-bedpe      [string:    path to SMC-RNA fusions.bedpe  ]
        -p/--peptide-lengths   [string:    comma sperated peptide lengths ]
        -v/--available-alleles [string:    path to allelelist             ]

    Optional Parameters:
        -o/--output-dir        [string:    path to output dir.  Default: ./        ]
        -k/--keep-tmp          [string:    keep the tmp folder.                    ]
        -n/--path-to-netMHC4   [string:    path to netMHC4.0.   Default: netMHC4.0 ]
        -x/--netMHC4.0-options ["string":  netMHC options.      Default: ""        ]

    netMHC4.0 options for -x:
        Please type netMHC4.0 to see after installation

    Version:                    1.0.0
          """

#parameters
hla_allele_file = ''
fusion_bedpe_file = ''
peptide_lengths = ''
output_dir = ''
netMHC_options = ''
path_to_netMHC4 = ''
avail_file = ''
is_rm_tmp=True

def setDefault():
    global output_dir
    output_dir = './'
    global path_to_netMHC4
    path_to_netMHC4 = 'netMHC4.0'

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"ha:f:p:v:o:kn:x:", ["help",
                                                             "hla-allele=",
                                                             "fusion-bedpe=",
                                                             "peptide-lengths=",
                                                             "available-alleles=",
                                                             "output-dir=",
                                                             "keep-tmp",
                                                             "hla-options=",
                                                             "path-to-netMHC4=",
                                                             "netMHC4-options="])
    except getopt.GetoptError:
        print "Except"
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        elif opt in ("-k","--keep-tmp"):
            global is_rm_tmp
            is_rm_tmp = False
        elif opt in ("-a", "--hla-allele"):
            global hla_allele_file
            hla_allele_file = arg
        elif opt in ("-f", "--fusion-bedpe"):
            global fusion_bedpe_file
            fusion_bedpe_file = arg
        elif opt in ("-p", "--peptide_lengths"):
            global peptide_lengths
            peptide_lengths = arg
        elif opt in ("-v", "--available-alleles"):
            global avail_file
            avail_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-n", "--path-to-netMHC4"):
            global path_to_netMHC4
            path_to_netMHC4 = arg
        elif opt in ("-x", "--netMHC4-options"):
            global netMHC4_options
            netMHC4_options = arg

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir( path, 0755 )

avail_allele = []

def get_avail_allele():
    f=open(avail_file, "r")
    f.readline()
    while True:
        line=f.readline()
        if line=="":
           break
        else:
           tmp=line.split("\t")
           avail_allele.append(tmp[0])
    f.close()

hla_allele_dic = {}

def get_allele_dic():
    f=open(hla_allele_file,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            tmp=line.split("\t")
            dicString=''
            tmp[len(tmp)-1]=tmp[len(tmp)-1][0:len(tmp[len(tmp)-1])-1]
            for x in range(len(tmp)):
                if x==1:
                    dicString=tmp[x]
                if x>1:
                    dicString=dicString+'\t'+tmp[x]
            hla_allele_dic[tmp[0]]=dicString
    f.close()

def get_allele_string():
    allele_str=''
    f=open(hla_allele_file,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            tmp=line.split("\t")
            tmp2=tmp[0].replace(":","")
            if tmp2 in avail_allele:
                if allele_str=='':
                    allele_str=tmp2
                else:
                    allele_str=allele_str+','+tmp2
    f.close()         
    return allele_str

fusion_records = {}

def get_fusion_records():
    numRecord=0
    f=open(fusion_bedpe_file,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            fusion_records[numRecord]=line
            numRecord=numRecord+1
    f.close()



def get_fasta(klen,outfilename):
    seqs = []
    f=open(fusion_bedpe_file,"r")
    numRecord=0
        
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            numRecord=numRecord+1
            tmp = line.split("\t")
            tmpLen = len(tmp)
            peps = tmp[tmpLen-5].split(",")
            pepLs = tmp[tmpLen-4].split(",")
            for x in range(len(peps)):
                pepLen=len(peps[x])
                pepL=int(pepLs[x])
                if klen-1 <= pepL:
                    start=pepL-klen+1
                    end=min(pepLen,pepL+klen-1)
                    seq = []
                    seq.append(str(numRecord))
                    seq.append(peps[x][start:end])
                    seqs.append(seq)
    f.close()

    f=open(outfilename,"w")

    for t in range(len(seqs)):
        f.write(">%s\n%s\n" % (seqs[t][0],seqs[t][1]))
    f.close()

ks = []

def get_ks():
    global peptide_lengths
    tmp=peptide_lengths.split(",")
    for x in range(len(tmp)):
        ks.append(tmp[x])    

def run_netMHC4():
    allele_strings=get_allele_string()    

    netMHC4_file=output_dir+'/netMHC4.0.out.append.txt'

    f=open(netMHC4_file,"w")
    f.close();

    for x in range(len(ks)):
        klen=int(ks[x])
        
        fasta_file=output_dir+'/tmp/netMHC4.0.k'+str(klen)+'.fasta'
        get_fasta(klen,fasta_file)
        
        cmd = path_to_netMHC4 +' -a '+ allele_strings + ' -l ' + str(klen) + ' ' + fasta_file;
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        f=open(netMHC4_file,"a")
        f.write(output)
        f.close()

def remove_tmp():
    if is_rm_tmp:
        cmd = 'rm -rf ' + output_dir +'/tmp'

def main(argv):
    setDefault()
    print hla_allele_file,fusion_bedpe_file,peptide_lengths,avail_file
    getParameters(argv[1:])
    print hla_allele_file,fusion_bedpe_file,peptide_lengths,avail_file
    if hla_allele_file=='' or fusion_bedpe_file=='' or peptide_lengths=='' or avail_file=='':
        usage()
        exit(1);
    global output_dir
    output_dir=os.path.realpath(output_dir);
    make_dir(output_dir)
    make_dir(output_dir+'/tmp')
    get_avail_allele()
    get_allele_dic()
    get_fusion_records()
    get_ks()
    run_netMHC4()
    remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
