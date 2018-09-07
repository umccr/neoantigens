/*
 * Util.cpp
 *
 *  Created on: Apr 28, 2013
 *      Author: jinzhang
 */

#include "Util.h"


/*
 * count the length of chars in a file
 *
 *
 */


uint32_t getFilelength(char *file)
{
	struct stat filestatus;
    stat( file, &filestatus );
    uint32_t length = filestatus.st_size;
    return length;
}


/*
 * Read a block of chars from file till a \n of file
 * return the actual length;
 *
 */

int readBlock(char * block, int length, FILE *infile)
{
	fread(block,1,length,infile);
	if(block[length-1]!='\n' && block[length-1]!=EOF)
	{
		fgets(block+length,1024,infile);
		length+=strlen(block+length);
	}


	return length;
}


int InitialIntChar()
{
	intChar.insert(pair<int,char>(1,'A'));
	intChar.insert(pair<int,char>(2,'C'));
	intChar.insert(pair<int,char>(4,'G'));
	intChar.insert(pair<int,char>(8,'T'));
	intChar.insert(pair<int,char>(15,'N'));

	charChar.insert(pair<char,char>('A','T'));
	charChar.insert(pair<char,char>('C','G'));
	charChar.insert(pair<char,char>('G','C'));
	charChar.insert(pair<char,char>('T','A'));
	charChar.insert(pair<char,char>('N','N'));
    
    
    tableAmino.insert(pair<string,char>("GCT",'A'));
    tableAmino.insert(pair<string,char>("GCC",'A'));
    tableAmino.insert(pair<string,char>("GCA",'A'));
    tableAmino.insert(pair<string,char>("GCG",'A'));
    tableAmino.insert(pair<string,char>("CGT",'R'));
    tableAmino.insert(pair<string,char>("CGC",'R'));
    tableAmino.insert(pair<string,char>("CGA",'R'));
    tableAmino.insert(pair<string,char>("CGG",'R'));
    tableAmino.insert(pair<string,char>("AGA",'R'));
    tableAmino.insert(pair<string,char>("AGG",'R'));
    tableAmino.insert(pair<string,char>("AAT",'N'));
    tableAmino.insert(pair<string,char>("AAC",'N'));
    tableAmino.insert(pair<string,char>("GAT",'D'));
    tableAmino.insert(pair<string,char>("GAC",'D'));
    tableAmino.insert(pair<string,char>("TGT",'C'));
    tableAmino.insert(pair<string,char>("TGC",'C'));
    tableAmino.insert(pair<string,char>("CAA",'Q'));
    tableAmino.insert(pair<string,char>("CAG",'Q'));
    tableAmino.insert(pair<string,char>("GAA",'E'));
    tableAmino.insert(pair<string,char>("GAG",'E'));
    tableAmino.insert(pair<string,char>("GGT",'G'));
    tableAmino.insert(pair<string,char>("GGC",'G'));
    tableAmino.insert(pair<string,char>("GGA",'G'));
    tableAmino.insert(pair<string,char>("GGG",'G'));
    tableAmino.insert(pair<string,char>("CAT",'H'));
    tableAmino.insert(pair<string,char>("CAC",'H'));
    tableAmino.insert(pair<string,char>("ATT",'I'));
    tableAmino.insert(pair<string,char>("ATC",'I'));
    tableAmino.insert(pair<string,char>("ATA",'I'));
    tableAmino.insert(pair<string,char>("TTA",'L'));
    tableAmino.insert(pair<string,char>("TTG",'L'));
    tableAmino.insert(pair<string,char>("CTT",'L'));
    tableAmino.insert(pair<string,char>("CTC",'L'));
    tableAmino.insert(pair<string,char>("CTA",'L'));
    tableAmino.insert(pair<string,char>("CTG",'L'));
    tableAmino.insert(pair<string,char>("AAA",'K'));
    tableAmino.insert(pair<string,char>("AAG",'K'));
    tableAmino.insert(pair<string,char>("ATG",'M'));
    tableAmino.insert(pair<string,char>("TTT",'F'));
    tableAmino.insert(pair<string,char>("TTC",'F'));
    tableAmino.insert(pair<string,char>("CCT",'P'));
    tableAmino.insert(pair<string,char>("CCC",'P'));
    tableAmino.insert(pair<string,char>("CCA",'P'));
    tableAmino.insert(pair<string,char>("CCG",'P'));
    tableAmino.insert(pair<string,char>("TCT",'S'));
    tableAmino.insert(pair<string,char>("TCC",'S'));
    tableAmino.insert(pair<string,char>("TCA",'S'));
    tableAmino.insert(pair<string,char>("TCG",'S'));
    tableAmino.insert(pair<string,char>("AGT",'S'));
    tableAmino.insert(pair<string,char>("AGC",'S'));
    tableAmino.insert(pair<string,char>("ACT",'T'));
    tableAmino.insert(pair<string,char>("ACC",'T'));
    tableAmino.insert(pair<string,char>("ACA",'T'));
    tableAmino.insert(pair<string,char>("ACG",'T'));
    tableAmino.insert(pair<string,char>("TGG",'W'));
    tableAmino.insert(pair<string,char>("TAT",'Y'));
    tableAmino.insert(pair<string,char>("TAC",'Y'));
    tableAmino.insert(pair<string,char>("GTT",'V'));
    tableAmino.insert(pair<string,char>("GTC",'V'));
    tableAmino.insert(pair<string,char>("GTA",'V'));
    tableAmino.insert(pair<string,char>("GTG",'V'));
    tableAmino.insert(pair<string,char>("TAA",'X'));//X for stop
    tableAmino.insert(pair<string,char>("TGA",'X'));
    tableAmino.insert(pair<string,char>("TAG",'X'));
    
	return 0;
}

char getCharComp(char reada)
{
	return charChar[reada];

}



char getCharA(int reada)
{
	return intChar[reada];
}












char getAmino(string a)
{
    return tableAmino[a];
}






int getPeptide(vector<char> & seq5p, vector<char> & seq, int start_pos, vector<char> & peptide, int & full, int & left)
{
    full=0;
    for(int i=start_pos-1;i<seq5p.size()-2;i+=3)
    {
        full++;
    }
    
    left=seq5p.size()-3*full-start_pos+1;

 
    for(int i=start_pos-1;i<seq.size()-2;i+=3)
    {
        string a;
        for(int j=0;j<3;j++)
            a=a+seq[i+j];
        char am=getAmino(a);
        if(am!='X')
            peptide.push_back(am);
        else
            break;
    }
    return 0;
}



std::vector<std::string> &my_split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> my_split(const std::string &s, char delim) {
       std::vector<std::string> elems;
       my_split(s, delim, elems);
       return elems;
}

double f_score(double sensitivity, double precision)
{
    return 2.0*(sensitivity*precision)/(sensitivity+precision);
}

string my_db2string(double db)
{
    std::ostringstream strs;
    strs << db;
    return strs.str();
}

int get_pseudo_counts(string & input, pseudo_counts_t & pct)
{
   pct.t_t=0;pct.f_t=0;pct.truth_t=0;pct.t_g=0;pct.f_g=0;pct.truth_g=0;
   vector<string> numStrs=my_split(input,',');
   if(numStrs.size()!=6)
   {
       cerr<<"We need 6 numbers for pseudo counts."<<endl;
       exit(0);
   }
   int t_t,f_t,truth_t,t_g,f_g,truth_g;
   t_t=atoi(numStrs[0].c_str());
   f_t=atoi(numStrs[1].c_str());
   truth_t=atoi(numStrs[2].c_str());
   t_g=atoi(numStrs[3].c_str());
   f_g=atoi(numStrs[4].c_str());
   truth_g=atoi(numStrs[5].c_str());
   
   if(t_t>truth_t)
   {
       cerr<<"For pseudo counts: True positive fusions transcripts should be less than or equal to transcripts in truth."<<endl;
       exit(0);
   }
   if(t_g>truth_g)
   {
       cerr<<"For pseudo counts: True positive fusions should be less than or equal to fusions in truth."<<endl;
       exit(0);
   }
   if(t_t<t_g)
   {
       cerr<<"For pseudo counts: True positive fusions transcripts should be greater than or equal to true positive fusions ."<<endl;
       exit(0);
   }
   if(f_t<f_g)
   {
       cerr<<"For pseudo counts: False positive fusions transcripts should be greater than or equal to false positive fusions ."<<endl;
       exit(0);
   }
   if(truth_t<truth_g)
   {
       cerr<<"For pseudo counts: Fusions transcripts in truth should be greater than or equal to fusions in truth."<<endl;
       exit(0);
   }

   pct.t_t=t_t;
   pct.f_t=f_t;
   pct.truth_t=truth_t;

   pct.t_g=t_g;
   pct.f_g=f_g;
   pct.truth_g=truth_g;


   return 0;
}



