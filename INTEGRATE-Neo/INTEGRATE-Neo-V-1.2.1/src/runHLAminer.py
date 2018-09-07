#!/usr/bin/python

import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT

def usage():
    print """
    runHLAminer -a <hla-abc-cds> -n <hla-nom-p> -m <path-to-hlaminer> -1 <input-file-1> -2 <input-file-2>
                -o <output-dir> -k -p <path-to-perl> -b <path-to-bwa> -x "<hlaminer-options>"

    Requested Parameters:
        -a/--hla-abc-cds       [string:    path to HLA-ABC-CDS.fasta]
        -n/--hla-nom-p         [string:    path to hla_nom_p.txt    ]
        -m/--path-to-hlaminer  [string:    path to HLAminer         ]
        -1/--input-file-1      [string:    input FASTA/FASTQ file 1 ]
        -2/--input-file-2      [string:    input FASTA/FASTQ file 2 ]
        
    Optional Parameters:
        -o/--output-dir        [string:    path to output dir.  Default: ./  ]
        -k/--keep-tmp          [string:    keep the tmp folder.              ]
        -p/--path-to-perl      [string:    path to perl.        Default: perl]
        -b/--path-to-bwa       [string:    path to bwa.         Default: bwa ]
        -x/--hla-options       ["string:"  HLAminer options.    Default: ""  ] 

    HLAminer options for -x:
        -i minimum % sequence identity...............<99>
        -q minimum log10 (phred-like) expect value...<30>
        -s minimum score.............................<1000>
        -n consider null alleles (1=yes/0=no)........<0>
        -l label (run name) -optional-
   
    Version:                    1.0.0
          """

#parameters
input_file_1 = ''
input_file_2 = ''
hla_abc_cds = ''
hla_nom_p = ''
output_dir = ''
path_to_hlaminer = ''
path_to_perl = ''
path_to_bwa = ''
is_rm_tmp=True
hla_options=''

def setDefault():
    global output_dir
    output_dir = './'
    global path_to_perl
    path_to_perl = 'perl'
    global path_to_bwa
    path_to_bwa = 'bwa'

def use_real_path():
    global hla_abc_cds
    path, filename = os.path.split(hla_abc_cds)
    path=os.path.dirname(os.path.realpath(hla_abc_cds));
    hla_abc_cds=path+'/'+filename
    global hla_nom_p
    path, filename = os.path.split(hla_nom_p)
    path=os.path.dirname(os.path.realpath(hla_nom_p));
    hla_nom_p=path+'/'+filename
    global output_dir
    output_dir=os.path.realpath(output_dir);


def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hk1:2:a:n:o:m:x:p:b:",["help",
                                                             "keep-tmp",
                                                             "input-file-1=",
                                                             "input-file-2=",
                                                             "hla-abc-cds=",
                                                             "hla-nom-p=",
                                                             "output-dir=",
                                                             "path-to-hlaminer=",
                                                             "hla-options=",
                                                             "path-to-perl=",
                                                             "path-to-bwa="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        elif opt in ("-k","--keep-tmp"):
            global is_rm_tmp
            is_rm_tmp = False
        elif opt in ("-1", "--input-file-1"):
            global input_file_1
            input_file_1 = arg
        elif opt in ("-2", "--input-file-2"):
            global input_file_2
            input_file_2 = arg
        elif opt in ("-a", "--hla-abc-cds"):
            global hla_abc_cds
            hla_abc_cds = arg
        elif opt in ("-n", "--hla-nom-p"):
            global hla_nom_p
            hla_nom_p = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-m", "--path-to-hlaminer"):
            global path_to_hlaminer
            path_to_hlaminer = arg
        elif opt in ("-p", "--path-to-perl"):
            global path_to_perl
            path_to_perl = arg
        elif opt in ("-b", "--path-to-bwa"):
            global path_to_bwa
            path_to_bwa = arg
        elif opt in ("-x", "--hla-options"):
            global hla_options
            hla_options = arg

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir( path, 0755 ) 
def run_bwa():
    cmd = path_to_bwa + ' aln -e 0 -o 0 -t 4 -f ' + output_dir + '/tmp/1.sai ' + hla_abc_cds + ' ' + input_file_1
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

    cmd = path_to_bwa + ' aln -e 0 -o 0 -t 4 -f ' + output_dir + '/tmp/2.sai ' + hla_abc_cds + ' ' + input_file_2
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

    cmd= path_to_bwa+ ' sampe -f ' + output_dir +'/tmp/all.sam ' + hla_abc_cds +' '+ output_dir +'/tmp/1.sai ' + output_dir +'/tmp/2.sai ' + input_file_1 + ' ' + input_file_2
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

def run_hla_miner():
    #cmd = 'cd '+output_dir
    #p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    #output = p.stdout.read()
    #print output
    cmd = path_to_perl +' '+ path_to_hlaminer +' -h '+ hla_abc_cds +' -p '+ hla_nom_p +' -a '+ output_dir +'/tmp/all.sam'
    cmd = cmd +' '+ hla_options
    p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

def remove_tmp():
    if is_rm_tmp:
        cmd = 'rm -rf ' + output_dir +'/tmp'

def main(argv):
    setDefault()
    getParameters(argv[1:])
    print hla_abc_cds,hla_nom_p,path_to_hlaminer,input_file_1,input_file_2
    if hla_abc_cds=='' or hla_nom_p=='' or path_to_hlaminer=='' or input_file_1=='' or input_file_2=='':
        usage()
        exit(1);
    print "before ...",output_dir
    use_real_path()
    print output_dir
    make_dir(output_dir)
    make_dir(output_dir+'/tmp')
    run_bwa()
    run_hla_miner()
    remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
