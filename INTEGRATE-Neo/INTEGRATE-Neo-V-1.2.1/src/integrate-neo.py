#!/usr/bin/python

import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT

def usage():
    print """
    integrate-neo -s <setup.ini> -1 <fastq1> -2 <fastq2> -x <hlaminer-options> -t <hla-allele-tsv> 
                  -f <fusion-bedpe> -r <reference> -g <gene-model> -u <rule-file> -d <dinucleotides-file>
                  -l <peptide-lengths> -a <affinity-score> -y <netMH4-options> -o <output-dir> -k  

    Parameters:

        -s/--setup-file   [Optional]                    
                          [string:   path to setup file. Default: setup.ini in the same dir.                 ]

        -1/--fastq1       [Requested if -t not set] 
                          [string:   path to fastq1 to call hla alleles.                                     ]

        -2/--fastq2       [Requested if -t not set] 
                          [string:   path to fastq2 to call hla alleles.                                     ]

        -x/--hlaminer-opt [Optional]    
                          [string:   additional HLAminer options. encompass with "".                         ]

        -t/--hla-allele   [Optional, ignore -1,-2,-x if set ] 
                          [string:   path to tsv file for hla alleles.                                       ]

        -f/--fusion-bedpe [Requested]
                          [string:   path to file of fuions in SMC-RNA 11-column bedpe format.               ]   
  
        -r/--reference    [Requested]
                          [string:   Reference genome.                                                       ]

        -g/--gene-model   [Requested]
                          [string:   GenePred format annotation.                                             ]
        
        -u/--rule-file    [Optional]
                          [string:   path to rule.txt.  Default: rule.txt in the same dir.                   ]
        
        -d/--di-nucleo    [Optional]
                          [string:   path to difile.txt.  Default: difile.txt in the same dir.               ]

        -l/--pep-lengths  [Optional]
                          [string:   Comma separated string for lengths of peptides. Default: 8,9,10,11      ]

        -a/--affinity     [Optional]
                          [value:    affinity score to consider antigen.  Default: 500.                      ]

        -y/--netMHC4-opt  [Optional]
                          [string:   additional netMHC4 options. encompass with "".                          ]

        -o/--ouput-dir    [Optional]
                          [string:   output dir.   Default:  fusion_antigen_out                              ]

        -k/--keep-tmp     [Optional]
                          [if set, keeps the tmp folder under the output dir                                 ] 

    Version:                    1.2.1
          """

#parameters

#requested:
fastq1=''
fastq2=''
fusionBedpe=''
reference=''
geneModel=''

#optional:

setupFile='' #
hlaminerOpt=''
hlaAllele=''
ruleFile=''  #
di_file=''
pepLengths='8,9,10,11'
affinity=500
netMHC4Opt=''
outputDir='./fusion_antigen_out'
keepTmp=False

#from setup:
hla_abc_cds=''
hla_nom_p=''
hlaminer=''
perl='perl'
bwa='bwa'
netMHC4='netMHC4.0'
available_alleles_file=''

def initialSetupFile():
    cur=os.path.dirname(os.path.abspath(__file__))
    global setupFile
    setupFile=cur+"/"+"setup.ini"

def initialRuleFile():
    cur=os.path.dirname(os.path.abspath(__file__))
    global ruleFile
    ruleFile=cur+"/"+"rule.txt"
    global di_file
    di_file=cur+"/"+"difile.txt"


def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"h1:2:f:r:g:s:x:t:u:d:l:a:y:o:k",["help",
                                                          "fastq1=",
                                                          "fastq2=",
                                                          "fusion-bedpe=",
                                                          "reference=",
                                                          "gene-model=",
                                                          "setup-file=",
                                                          "hlaminer-opt=",
                                                          "hla-allele=",
                                                          "rule-file=",
                                                          "di-nucleo=",
                                                          "pep-lengths=",
                                                          "affinity=",
                                                          "netMHC4-opt=",
                                                          "ouput-dir=",
                                                          "keep-tmp="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        elif opt in ("-1", "--fastq1"):
            global fastq1
            fastq1 = arg
        elif opt in ("-2", "--fastq2"):
            global fastq2
            fastq2 = arg
        elif opt in ("-t", "hla-allele"):
            global hlaAllele
            hlaAllele = arg
        elif opt in ("-f", "--fusion-bedpe"):
            global fusionBedpe
            fusionBedpe = arg
        elif opt in ("-r", "--reference"):
            global reference
            reference = arg
        elif opt in ("-g", "--gene-model"):
            global geneModel
            geneModel = arg
        elif opt in ("-u", "--rule-file"):
            global ruleFile
            ruleFile = arg
        elif opt in ("-d", "--di-nucleo"):
            global di_file
            di_file = arg
        elif opt in ("-l", "--pep-lengths"):
            global pepLengths
            pepLengths = arg
        elif opt in ("-a", "--affinity"):
            global affinity
            affinity = arg
        elif opt in ("-x", "--hlaminer-opt"):
            global hlaminerOpt
            hlaminerOpt = arg
        elif opt in ("-y", "--netMHC4-opt"):
            global netMHC4Opt
            netMHC4Opt = arg
        elif opt in ("-o", "--output-dir"):
            global outputDir
            outputDir = arg
        elif opt in ("-k", "--keep-tmp"):
            global keepTmp
            keepTmp = True

def getFromSetup():
    global hla_abc_cds,hla_nom_p,hlaminer,perl,bwa,netMHC4,available_alleles_file
    f = open(setupFile,"r")
    while True:
        line=f.readline()
        if line=='':
            break
        else:
            tmp=line.split("\t")
            tmp1=tmp[1][0:len(tmp[1])-1]
            if tmp[0]=='hla_abc_cds' and tmp1!='':
                hla_abc_cds=tmp1
            if tmp[0]=='hla_nom_p' and tmp1!='':
                hla_nom_p=tmp1
            if tmp[0]=='HLAminer' and tmp1!='':
                hlaminer=tmp1
            if tmp[0]=='perl' and tmp1!='':
                perl=tmp1
            if tmp[0]=='bwa' and tmp1!='':
                bwa=tmp1
            if tmp[0]=='netMHC4' and tmp1!='':
                netMHC4=tmp1
            if tmp[0]=='available_alleles_file' and tmp1!='':
                available_alleles_file=tmp1
    f.close()
    if hla_abc_cds=='' or hla_nom_p=='' or hlaminer=='' or available_alleles_file=='':
        print "Please setup hla_abc_cds, hla_nom_p, HLAminer, and available_alleles_file in setup.ini"
        exit(1)

def main(argv):
    cur=os.path.dirname(os.path.abspath(__file__))
    python=sys.executable

    initialSetupFile()
    initialRuleFile()
    getParameters(argv[1:])
    getFromSetup()
    global hlaAllele
    if hlaAllele=='':
        if fastq1=='' or fastq2=='' or fusionBedpe=='' or reference=='' or geneModel=='':
            usage()
            exit(1)
    else:
        if fusionBedpe=='' or reference=='' or geneModel=='':
            usage()
            exit(1)

    print "[FA]creating and working and result directories..."
    if not os.path.exists(outputDir):
        os.mkdir(outputDir, 0755 )
    if not os.path.exists(outputDir+"/"+"tmp"):
        os.mkdir(outputDir+"/"+"tmp", 0755 )

    path, filename = os.path.split(fusionBedpe)
    bedpeAnnot=outputDir+"/"+filename+"."+"annot"

    print "[FA]annotating fusion bedpe file..."
    cmd = cur+"/"+"fusionBedpeAnnotator"+" --reference-file "+reference+" --gene-annotation-file "+geneModel+" --di-file "+di_file \
          +" --input-file "+fusionBedpe+" --output-file "+bedpeAnnot
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

    print "[FA]subsetting annotated fusion bedpe file..."
    bedpeSubset=outputDir+"/"+filename+"."+"subset"
    cmd = cur+"/"+"fusionBedpeSubsetter"+" --input-file "+bedpeAnnot+" --rule-file "+ruleFile+" --output-file "+bedpeSubset
    print cmd
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

    print "[FA]hla..."
    if hlaAllele=='':
        if True:    
            cmd = python+" "+cur+"/"+"runHLAminer.py"+" -a "+hla_abc_cds+" -n "+hla_nom_p+" -m "+hlaminer \
                +" -1 "+fastq1+" -2 "+fastq2 \
                +" -o "+outputDir \
                +" -p "+perl+" -b "+bwa
            if hlaminerOpt!='':
                cmd=cmd+" +x "+"\""+hlaminerOpt+"\""
            if keepTmp==True:
                cmd=cmd +" -k"
            print cmd
            p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            output = p.stdout.read()
            print output 
            cmd = python+" "+cur+"/"+"HLAminerToTsv.py"+" -l "+outputDir+"/"+"HLAminer_HPRA.log" \
                +" -c "+outputDir+"/"+"HLAminer_HPRA.csv" \
                +" -o "+outputDir
            print cmd
            p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            output = p.stdout.read()
            print output
            hlaAllele=outputDir+"/"+"HLAminer_alleles.tsv"

    print "[FA]NetHMC4"
    if True:
        cmd = python+" "+cur+"/"+"runNetMHC4WithSMCRNABedpe.py"+" -a "+hlaAllele+" -f "+bedpeSubset \
            +" -p "+pepLengths+" -o "+outputDir+" -n "+netMHC4+" -v "+available_alleles_file
        if netMHC4Opt!='':
            cmd=cmd+" +x "+"\""+netMHC4Opt+"\""
        if keepTmp==True:
            cmd=cmd +" -k"       
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        print output 
        cmd = python+" "+cur+"/"+"runAddNetMHC4Result.py"+" -b "+bedpeSubset+" -a "+hlaAllele \
            +" -f "+outputDir+"/"+"netMHC4.0.out.append.txt"+" -m "+str(affinity)+" -o "+outputDir+"/"+"result.bedpe"
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        print cmd
        output = p.stdout.read()
        print output

if __name__ == '__main__':
    sys.exit(main(sys.argv))
