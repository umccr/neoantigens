/****
  Jin Zhang, Jun 23, 2016
*****/

#include <getopt.h>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std;  

#include "CutterByRule.h"

map<int,char> intChar;
map<char,char> charChar;
map<string,char> tableAmino;

string version("0.1.0");

int usage()
{
    cerr<<endl;
    cerr<<"subsetter version "+version<<endl;
    cerr<<endl;
    cerr<<"Usage:"<<endl;
    cerr<<"subsetter --input-file input-file --rule-file rule-file --output-file output-file ";
    cerr<<endl;
    cerr<<"Required parameters:"<<endl;
    cerr<<"    -i/--input-file"<<endl;
    cerr<<"    -r/--rule-file"<<endl;
    cerr<<"    -o/--output-file"<<endl;
    cerr<<endl;
    return 0;
}


int main(int argc, char * argv[])
{
    int c;
    int option_index = 0;

    string input_file="";
    string output_file="";
    string rule_file="";

    static struct option long_options[] = {
        {"input-file",            required_argument, 0,  'i' },
        {"output-file",           required_argument, 0,  'o' },
        {"rule-file",             required_argument, 0,  'r' },
        {"help",                  no_argument,       0,  'h' },
        {0, 0, 0, 0}
    };

    while(1)
    {
        c = getopt_long(argc, argv, "i:o:r:h",
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
            case 'i':
                input_file=optarg;
                break;
            case 'o':
                output_file=optarg;
                break;
            case 'r':
                rule_file=optarg;
                break;
            default:
                break;
        }

    }

    if(rule_file.compare("")==0 || input_file.compare("")==0 || output_file.compare("")==0)
    {
        usage();
        exit(1);
    }


    string line;
    ifstream myfile ((char *)rule_file.c_str());
    getline (myfile,line);
    std::vector<std::string> tmp = my_split(line, '\t');
    string ruleString=tmp[1];
    myfile.close();


    CutterByRule cbr;
    Bedpe in;
    in.loadFromFile((char*)input_file.c_str());


    Bedpe out;
    string tmpString=output_file+".tmp";
    cbr.cut(in, out, ruleString, tmpString);
    out.printBedpe((char *)output_file.c_str());
        

    return 0;
}
