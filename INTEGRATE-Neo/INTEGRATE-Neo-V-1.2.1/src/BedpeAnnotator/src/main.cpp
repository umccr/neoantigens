/*
 * main.cpp
 *
 *  Created on: Jun 16, 2016
 *      Author: Jin Zhang
 */

#include "Bedpe.h"
#include "Reference.h"
#include "Gene2.h"
#include "Annotate.h"
#include "Dinucleo.h"
#include "AnnotateDi.h"
#include <iostream>
#include <getopt.h>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std;

map<int,char> intChar;
map<char,char> charChar;
map<string,char> tableAmino;

string version("0.2.2");

int usage()
{
    cerr<<endl;
    cerr<<"fusionBedpeAnnotator version "+version<<endl;
    cerr<<endl;
    cerr<<"Usage:"<<endl;
    cerr<<"fusionBedpeAnnotator --reference-file reference-file --gene-annotation-file annot-file --di-file di-file ";
    cerr<<"--input-file input-bedpe-file --output-file output-bedpe-file"<<endl;
    cerr<<endl;
    cerr<<"Required parameters:"<<endl;
    cerr<<"    -r/--reference-file"<<endl;
    cerr<<"    -g/--gene-annotation-file"<<endl;
    cerr<<"    -d/--di-file"<<endl;
    cerr<<"    -i/--input-file"<<endl;
    cerr<<"    -o/--output-file"<<endl;
    cerr<<endl;
    cerr<<"Options:"<<endl;
    cerr<<"    -r/--reference-file          <string>    [ fasta                                            ]"<<endl;
    cerr<<"    -g/--gene-annotation-file    <string>    [ 11 columns             refer to INTEGRATE 0.2.5  ]"<<endl;
    cerr<<"    -i/--input-file              <string     [ bedpe                  refer to SMC-RNA          ]"<<endl;
    cerr<<"    -d/--di-file                 <string>    [ 2 columns              dinucleo canonical splice ]"<<endl;
    cerr<<"    -o/--output-file             <string>    [ bedpe                                            ]"<<endl;
    cerr<<endl;
    return 0;
}


int main(int argc, char * argv[])
{
    InitialIntChar();
    int c;
    int option_index = 0;

    string input_file="";
    string output_file="";
    string gene_file="";
    string reference_file="";
    string di_file="";

    static struct option long_options[] = {
        {"input-file",            required_argument, 0,  'i' },
        {"output-file",           required_argument, 0,  'o' },
        {"gene-annotation-file",  required_argument, 0,  'g' },
        {"reference-file",        required_argument, 0,  'r' },
        {"di-file",               required_argument, 0,  'd' },
        {"help",                  no_argument,       0,  'h' },
        {0, 0, 0, 0}
    };

    while(1)
    {
        c = getopt_long(argc, argv, "i:o:g:d:r:h",
                 long_options, &option_index);
        if (c==-1)
        {
            break;
        }
        switch(c)
        {
            case 'h':
                usage();
                exit(0);
            case 'r':
                reference_file=optarg;
                break;
            case 'g':
                gene_file=optarg;
                break;
            case 'd':
                di_file=optarg;
                break;
            case 'o':
                output_file=optarg;
                break;
            case 'i':
                input_file=optarg;
                break;
            default:
                break;
        }

    }

    if (reference_file.compare("")==0 || gene_file.compare("")==0 || output_file.compare("")==0 || input_file.compare("")==0 || di_file.compare("")==0 )
    {
        usage();
        exit(1);    
    }

    cerr<<"loading gene annotation file..."<<endl;    

    Gene g;
    g.loadGenesFromFile((char*)gene_file.c_str());
    g.setGene();    

    cout<<"loading reference..."<<endl;
    Reference ref;
    ref.setIsInt(0);
    ref.test((char*)reference_file.c_str());

    cerr<<"loading bedpe file..."<<endl;

    Bedpe bpe;
    bpe.loadFromFile((char*)input_file.c_str());

    cerr<<"annotate..."<<endl;

    Annotate annot;
    annot.assignJunctions(g,ref,bpe);
    annot.getJunctionType(g,ref,bpe);
    annot.getJunctionPeptide(g,ref);
    annot.annotate(bpe,g,ref);

    Dinucleo di;   
    di.loadFromFile((char*)di_file.c_str());

    AnnotateDi annotDi;
    annotDi.annotateDi(bpe,di,ref);
 
    bpe.printBedpe((char*)output_file.c_str());

    return 0;
}
