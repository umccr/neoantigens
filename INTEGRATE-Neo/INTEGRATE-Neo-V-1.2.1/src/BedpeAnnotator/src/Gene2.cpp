/* 
 *  Created on: Mar 4, 2013
 *      Author: Jin Zhang
 */

#include "Gene2.h"



Gene::Gene() {
        // TODO Auto-generated constructor stub
}

Gene::~Gene() {
        // TODO Auto-generated destructor stub
}

bool myGeneSortFunc(gene_t i, gene_t j)
{
        if(i.chr.compare(j.chr)<0)
        {
                return true;
        }
        else if(i.chr.compare(j.chr)==0)
        {
                if(i.leftLimit<j.leftLimit)
                        return true;
                else
                        return false;
        }
        else
                return false;
}

bool myTransSortFunc(transcript_t i, transcript_t j)
{
        if(i.name2.compare(j.name2)<0)
                return true;
        else if(i.name2.compare(j.name2)==0)
        {
                if(i.chr.compare(j.chr)<0)
                {
                        return true;
                }
                else if(i.chr.compare(j.chr)==0)
                {
                        if(i.txStart<j.txStart)
                                return true;
                        else
                                return false;
                }
                else
                        return false;

        }
        else
                return false;
}

int Gene::loadGenesFromFile(char* file) {

	uint32_t length=getFilelength(file);
	
        FILE *infile;

	infile = fopen(file, "r");
	if (!infile) {
		cerr<<"Couldn't open file for reading: "<<file<<"."<<endl;
		exit(1);
	}

	char * fileContent;

	try
	{
		fileContent= new char [length];
	}
	catch(exception& e)
	{
		cerr << "Trying to allocate Memory to load Genes. exception caught: " << e.what() << endl;
		return 1;
	}


	size_t result = fread (fileContent,1,length,infile);
	if (result != length)
	{
		cerr << "Fail to read the genes file"<<endl;
		exit (1);
	}

	char * line_buffer=fileContent;
	char* nextN=NULL;

	char* p=NULL;
	char* NC=NULL;
	char intChar [1024];

 	//int bin;
 	char nameC [1024];
 	char chromC [1024];
 	char strandC [1024];
 	uint32_t txStart;
 	uint32_t txEnd;
 	uint32_t cdsStart;//Dec 7, 2015
 	uint32_t cdsEnd;
 	int exonCount;
 	char exonStartsC [1000000]; //should be ok
 	char exonEndsC [1000000];	  //should be ok
 	//int score;
 	char name2C[1024];
 	char * chr;
 	string chrStr;

	int strand;

	nextN=strchr(line_buffer,'\n');
	line_buffer=nextN+1;



	int num=0;
	while (1) {
		nextN=strchr(line_buffer,'\n');

		//sscanf(line_buffer,"%d %s %s %s %u %u %u %u %d %s %s %d %s", &bin, nameC, chromC, strandC, &txStart, &txEnd,
		//	    			&cdsStart, &cdsEnd, &exonCount, exonStartsC, exonEndsC, &score, name2C);

		nextN[0]='\0';
		int numnum=sscanf(line_buffer,"%s\t%s\t%s\t%u\t%u\t%u\t%u\t%d\t%s\t%s\t%s", nameC, chromC, strandC, &txStart, &txEnd, &cdsStart, &cdsEnd, &exonCount, exonStartsC, exonEndsC, name2C);
		if(numnum!=11)
		{
			cerr<<"error loading genes at: "<<line_buffer<<endl;
			cerr<<"From 0.3.0, INTEGRATE also use cdsStart and cdsEnd, check your annotation file should have 11 columns."<<endl;
                	exit(1);
		}


	   	uint32_t * ps=new uint32_t [exonCount];
	   	uint32_t * pe=new uint32_t [exonCount];

	   	p=exonStartsC;
	   	for(int i=0;i<exonCount;i++)
	    {
	    	NC=strchr(p, ',');
	    	strncpy(intChar, p, NC-p);
	    	intChar[NC-p]='\0';
	    	ps[i]=atol(intChar);
	    	p=NC+1;
	    }

	   	p=exonEndsC;
	   	for(int i=0;i<exonCount;i++)
	   	{
	   		NC=strchr(p, ',');

	   		strncpy(intChar, p, NC-p);
	   		intChar[NC-p]='\0';
	   		pe[i]=atol(intChar);
	   		p=NC+1;
	    }

	   	if(strcmp(strandC,"+")==0)
	   		strand=0;
	    else
	    	strand=1;

	   	transcript_t tt;
	   	//tt.bin=bin;
	    tt.name=string(nameC);
	   	if(strstr(chromC,"chr"))
	   	{
	   		chr=chromC+3;
	    }
	   	else
	   		chr=chromC;
	    chrStr=string(chr);
	    tt.chr=chrStr;
            tt.strand=strand;
	    tt.txStart=txStart;
	    tt.txEnd=txEnd;
	    tt.cdsStart=cdsStart;//Dec 7, 2015
	    tt.cdsEnd=cdsEnd;
	    tt.exonCount=exonCount;
	    tt.exonStarts=ps;
	    tt.exonEnds=pe;
	    //tt.score=score;
	    tt.name2=string(name2C);

            //Dec 2015, let us not use the transcripts truncated in coding region
            int isRealAdd=1;
            if(tt.cdsStart!=tt.cdsEnd && ((tt.txStart==tt.cdsStart)||(tt.txEnd==tt.cdsEnd)))
                  isRealAdd=0;
            if(isRealAdd==1)
	    {
                  transcripts.push_back(tt);
	          num++;
            }

	    if(nextN-fileContent >= length-1)
	    {
	    	cerr<<num<<" complete transcripts loaded."<<endl;
	    	break;
	    }
	    else
	    {
	    	line_buffer=nextN+1;
	    }

	}
	sort(transcripts.begin(),transcripts.end(),myTransSortFunc);
	return 0;

}

int Gene::setGene() {


	int fid=0;
	if(transcripts.size()<1)
	{
		cerr<<"No gene annotation"<<endl;
		exit(1);
	}
	else
	{
		gene_t gt;
		gt.name2=transcripts[0].name2;
		gt.strand=transcripts[0].strand;
		gt.transIds.push_back(0);
		gt.leftLimit=transcripts[0].txStart;
		gt.rightLimit=transcripts[0].txEnd;
		gt.chr=transcripts[0].chr;
		gt.fakeId=-1;
		genes.push_back(gt);
	}

	for(int i=1;i<transcripts.size();i++)
	{
		if(transcripts[i].name2.compare(transcripts[i-1].name2)!=0)
		{
			gene_t gt;
			gt.name2=transcripts[i].name2;
			gt.strand=transcripts[i].strand;
			gt.transIds.push_back(i);
			gt.leftLimit=transcripts[i].txStart;
			gt.rightLimit=transcripts[i].txEnd;
			gt.chr=transcripts[i].chr;
			gt.fakeId=-1;
			genes.push_back(gt);
		}
		else
		{

			if(transcripts[i].chr!=transcripts[i-1].chr)
			{
				gene_t gt;
				gt.name2=transcripts[i].name2;
				gt.strand=transcripts[i].strand;
				gt.transIds.push_back(i);
				gt.leftLimit=transcripts[i].txStart;
				gt.rightLimit=transcripts[i].txEnd;
				gt.chr=transcripts[i].chr;
				if(genes[genes.size()-1].fakeId!=-1)
					gt.fakeId=fid++;
				else
				{
					genes[genes.size()-1].fakeId=fid++;
					gt.fakeId=fid++;
				}
				genes.push_back(gt);
				continue;
			}
			else
			{
				if(transcripts[i].txStart > genes[genes.size()-1].rightLimit)
				{
					gene_t gt;
					gt.name2=transcripts[i].name2;
					gt.strand=transcripts[i].strand;
					gt.transIds.push_back(i);
					gt.leftLimit=transcripts[i].txStart;
					gt.rightLimit=transcripts[i].txEnd;
					gt.chr=transcripts[i].chr;
					if(genes[genes.size()-1].fakeId!=-1)
						gt.fakeId=fid++;
					else
					{
						genes[genes.size()-1].fakeId=fid++;
						gt.fakeId=fid++;
					}
					genes.push_back(gt);
					continue;
				}
			}

			//if(transcripts[i].txStart<genes[genes.size()-1].leftLimit)//this is not possible now
			//{
			//	genes[genes.size()-1].leftLimit=transcripts[i].txStart;
			//}
			if(transcripts[i].txEnd>genes[genes.size()-1].rightLimit)
			{
				genes[genes.size()-1].rightLimit=transcripts[i].txEnd;
			}
			genes[genes.size()-1].transIds.push_back(i);
		}
	}


	sort(genes.begin(),genes.end(),myGeneSortFunc);
	return 0;
}

int Gene::isInGene(string chr, uint32_t pos, vector<int>& geneIds) {
	gene_t dumbGene;
	dumbGene.chr=chr;
	dumbGene.leftLimit=pos;

	vector<gene_t>::iterator up=upper_bound(genes.begin(),genes.end(),dumbGene,myGeneSortFunc);
	if (up-genes.begin()==0)
		return 0;
	up--;
	while(up-genes.begin()>=0 && (*up).chr.compare(chr)==0 && pos > (*up).leftLimit)
	{
		if((*up).leftLimit < pos && (*up).rightLimit > pos)
		{
			geneIds.push_back(up-genes.begin());
		}
		up--;
	}
	if(geneIds.size()>0)
	{
		return 1;
	}
	else
		return 0;

}


gene_t * Gene::getGene(int index) {
	return &(genes[index]);
}


int Gene::getStrand(int geneId) {
	return genes[geneId].strand;
}

string Gene::getName2(int geneId) {
        return genes[geneId].name2;
}


uint32_t Gene::getLimitLeft(int geneId) {
	return genes[geneId].leftLimit;
}

uint32_t Gene::getLimitRight(int geneId) {
	return genes[geneId].rightLimit;
}

int Gene::getIndex(string name, vector<int> & ids) {
	for(int i=0;i<genes.size();i++)
	{
		if(genes[i].name2.compare(name)==0)
		{
			ids.push_back(i);
		}
	}
	return 0;
}


int Gene::getBestDiff(int gid, int pos, int isbkLeft)
{
    int best=1000000000;
    for(int i=0;i<genes[gid].transIds.size();i++)
    {
        int tranId=genes[gid].transIds[i];
        int count=transcripts[tranId].exonCount;
        for(int j=0;j<count;j++)
        {
            int pp1=transcripts[tranId].exonStarts[j];
            int pp2=transcripts[tranId].exonEnds[j];
            if(isbkLeft==1 && abs(pos-pp1)<best)
                best=abs(pos-pp1);
            if(isbkLeft==0 && abs(pos-pp2)<best)
                best=abs(pos-pp2);
        }
    }
    return best;
}

int Gene::getBestDiff2(vector<int> & gids, int pos, int isbkLeft)
{
    int best=1000000000;
    for(int x=0;x<gids.size();x++)
    {
        int gid=gids[x];
        for(int i=0;i<genes[gid].transIds.size();i++)
        {
            int tranId=genes[gid].transIds[i];
            int count=transcripts[tranId].exonCount;
            for(int j=0;j<count;j++)
            {
                int pp1=transcripts[tranId].exonStarts[j];
                int pp2=transcripts[tranId].exonEnds[j];
                if(isbkLeft==1 && abs(pos-pp1)<best)
                    best=abs(pos-pp1);
                if(isbkLeft==0 && abs(pos-pp2)<best)
                    best=abs(pos-pp2);
            }
        }
    }
    return best;
}


int Gene::isAt5p(int gid, int isbkLeft)
{
    int is5p;
    if((isbkLeft==1 && genes[gid].strand==0) || (isbkLeft==0 && genes[gid].strand==1))
        is5p=0;
    else
        is5p=1;
    return is5p;
}

int Gene::isAt5pSL(int strand, int isbkLeft)
{
    int is5p;
    if((isbkLeft==1 && strand==0) || (isbkLeft==0 && strand==1))
        is5p=0;
    else
        is5p=1;
    return is5p;
}


int Gene::getCodingAndBaseLeft(int tranId, int exonNum, int isbkLeft, int & isCoding, int & baseLeft)
{
//cout<<"name "<<transcripts[tranId].name<<endl;
//cout<<"isbkLeft "<<isbkLeft<<endl;
    transcript_t tt=transcripts[tranId];
    if(tt.cdsStart==tt.cdsEnd)
        isCoding=0;
    else
        isCoding=1;
    baseLeft=-1;
    if(isCoding==1)
    {
        int len_seq=0;
	int diff=0;
	int number;

        if(transcripts[tranId].strand==0 && isbkLeft==0)
		number=exonNum;
	if(transcripts[tranId].strand==0 && isbkLeft==1)
		number=exonNum-1;
        if(transcripts[tranId].strand==1 && isbkLeft==1)
		number=transcripts[tranId].exonCount-exonNum;
        if(transcripts[tranId].strand==1 && isbkLeft==0)
                number=transcripts[tranId].exonCount-exonNum+1;
//cout<<"number = "<<number<<endl;
        for(int i=0;i<number;i++)
        {
            len_seq+=tt.exonEnds[i]-tt.exonStarts[i];
            if(tt.exonEnds[i] < tt.cdsStart)
            {
//cout<<tt.exonEnds[i]<<" XXX "<<tt.exonStarts[i]<<endl;
                diff+=tt.exonEnds[i]-tt.exonStarts[i];
//cout<<"diff = "<<diff<<endl;
            }
            else if (tt.exonEnds[i] >= tt.cdsStart && tt.exonStarts[i] <= tt.cdsStart)
            {
               diff+=tt.cdsStart-tt.exonStarts[i];
//cout<<tt.cdsStart<<" YYY "<<tt.exonStarts[i]<<endl;
//cout<<"diff = "<<diff<<endl;
            }

        }
        if(len_seq-diff>=0 && tt.cdsEnd > tt.exonEnds[number-1])
        {
//cout<<"len_seq = "<<len_seq<<" diff = "<<diff<<endl;
            baseLeft=(len_seq-diff)%3;
//cout<<"baseLeft "<<baseLeft<<endl;
        }
        if(baseLeft!=-1 && isbkLeft==1)//this !=-1 is OK
        {
            baseLeft=(3-baseLeft)%3;
//cout<<"baseLeft change to "<<baseLeft<<endl;
        }

        if(transcripts[tranId].strand==0 && isbkLeft==0 && tt.cdsStart > tt.exonEnds[number-1])
            baseLeft=-1;
	if(transcripts[tranId].strand==1 && isbkLeft==1 && tt.cdsEnd < tt.exonStarts[number-1])
            baseLeft=-1;
        if(transcripts[tranId].strand==0 && isbkLeft==0 && tt.cdsEnd < tt.exonEnds[number-1])
            baseLeft=-2;
        if(transcripts[tranId].strand==1 && isbkLeft==1 && tt.cdsStart > tt.exonStarts[number-1])
            baseLeft=-2;
//cout<<"baseLeft final to "<<baseLeft<<endl;
    }
    return 0;
}

int Gene::getCumu5pNT(int tranId, int exonNum, int isbkLeft, int & codingNT)
{
    transcript_t tt=transcripts[tranId];
    int isCoding;
    if(tt.cdsStart==tt.cdsEnd)
        isCoding=0;
    else
        isCoding=1;
    
    if(isCoding==1)//transcript coding, not necessarily junction
    {
	int diff=0;
	int number;
        int len_seq=0;        

        if(transcripts[tranId].strand==0 && isbkLeft==0)
        {
            number=exonNum;
            for(int i=0;i<number;i++)
            {
                len_seq+=tt.exonEnds[i]-tt.exonStarts[i];
                if(tt.exonEnds[i] < tt.cdsStart)
                {
                    diff+=tt.exonEnds[i]-tt.exonStarts[i];
                }
                else if (tt.exonEnds[i] >= tt.cdsStart && tt.exonStarts[i] <= tt.cdsStart)
                {
                    diff+=tt.cdsStart-tt.exonStarts[i];
                }
            }
            codingNT=len_seq-diff;
        }

        if(transcripts[tranId].strand==1 && isbkLeft==1)
        {
            number=transcripts[tranId].exonCount-exonNum;
            for(int i=number-1;i>=0;i--)
            {
                len_seq+=tt.exonEnds[i]-tt.exonStarts[i];
                if(tt.exonStarts[i] > tt.cdsEnd)
                {
                    diff+=tt.exonEnds[i]-tt.exonStarts[i];
                }
                else if (tt.exonStarts[i] <= tt.cdsEnd && tt.exonEnds[i] >= tt.cdsEnd)
                {
                    diff+=tt.exonEnds[i]-tt.cdsEnd;
                }
            }
            codingNT=len_seq-diff;
        }
    }
    return 0;
}


junction_t assign_junction(int gId, int is5p, string chr, int strand, int pos1, int pos2, int exonNum, int is_coding, int baseLeft, string name, int coding_start, int cumu5pNT)
{

    junction_t jt;
    jt.gId=gId;
    jt.is5p=is5p;
    jt.isCoding=is_coding;
    jt.coding_start=coding_start;
    jt.coding_left=baseLeft;
    jt.chr=chr;
    jt.strand=strand;
    jt.pos1=pos1;
    jt.pos2=pos2;
    jt.exonNum=exonNum;
    jt.name=name;
    jt.cumu5pNT=cumu5pNT;
    return jt;
}


int Gene::getBestExon2(int gid, int pos, int isbkLeft, vector<junction_t> & juncs)
{

    int is5p, strand, pos1, pos2, exonNum;
    int is_coding;
    int baseLeft;
    string chr;
    string name;

    int best=getBestDiff(gid,pos,isbkLeft);
    is5p=isAt5p(gid,isbkLeft);

    chr=genes[gid].chr;
    strand=genes[gid].strand;

    for(int i=0;i<genes[gid].transIds.size();i++)
    {
        int tranId=genes[gid].transIds[i];
        int count=transcripts[tranId].exonCount;

        int seq_len=0;

        for(int j=0;j<count;j++)
        {

            int pp1=transcripts[tranId].exonStarts[j];
            int pp2=transcripts[tranId].exonEnds[j];

            seq_len+=pp2-pp1;

            if(isbkLeft==1)
            {
                if(abs(pos-pp1)==best)
                {
                    name=transcripts[tranId].name;
                    pos1=pp1;
                    pos2=pp2;
                    if(strand==0)
                    {
                        exonNum=j+1;
                    }
                    else
                    {
                        exonNum=count-j;
                    }
                    getCodingAndBaseLeft(tranId, exonNum, isbkLeft, is_coding, baseLeft);
                    int coding_start=0;
                    if(is5p==1 && is_coding==1 && baseLeft>=0)
                    {
                        if(transcripts[tranId].cdsEnd>=pos1 && transcripts[tranId].cdsEnd<=pos2)
                            coding_start=pos2-transcripts[tranId].cdsEnd+1;
                        else
                            coding_start=((pos2-pos1)-baseLeft)%3+1;
                    }
//cout<<"push baseLeft "<<baseLeft<<endl;
                    juncs.push_back(assign_junction(gid, is5p, chr, strand, pos1, pos2, exonNum, is_coding, baseLeft, name, coding_start,0));
                }

            }
            else
            {

                if(abs(pos-pp2)==best)
                {
                    name=transcripts[tranId].name;
                    pos1=pp1;
                    pos2=pp2;
                    if(strand==0)
                    {
                        exonNum=j+1;
                    }
                    else
                    {
                        exonNum=count-j;
                    }

                    getCodingAndBaseLeft(tranId, exonNum, isbkLeft, is_coding, baseLeft);
                    int coding_start=0;
                    if(is5p==1 && is_coding==1 && baseLeft>=0)
                    {
                        if(transcripts[tranId].cdsStart>=pos1 && transcripts[tranId].cdsStart<=pos2)
                            coding_start=transcripts[tranId].cdsStart+1-pos1;
                        else
                            coding_start=((pos2-pos1)-baseLeft)%3+1;
                    }
//cout<<"push2 baseLeft "<<baseLeft<<endl;
                    juncs.push_back(assign_junction(gid, is5p, chr, strand, pos1, pos2, exonNum, is_coding, baseLeft, name, coding_start,0));
                }
            }

        }
    }


    return 0;

}

int Gene::getBestExon3(vector<int> & gids, int pos, int isbkLeft, vector<junction_t> & juncs)
{
//cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^ "<<endl;
int best=getBestDiff2(gids,pos,isbkLeft);
if(best>1)//////
  return 0;

for(int x=0;x<gids.size();x++)
{

    int gid=gids[x];
//  int is5p, strand; change to using transcript
    int pos1, pos2, exonNum;    
    int is_coding;
    int baseLeft;
    string chr;
    string name;
//cout<<"^^^^^"<<gid<<" "<<isbkLeft<<endl;
//    is5p=isAt5p(gid,isbkLeft);
//cout<<"is5p="<<is5p<<endl;
    chr=genes[gid].chr;
//    strand=genes[gid].strand;

    for(int i=0;i<genes[gid].transIds.size();i++)
    {
        int is5p, strand;
  
        int tranId=genes[gid].transIds[i];
        int count=transcripts[tranId].exonCount;
         
        strand=transcripts[tranId].strand;
        is5p=isAt5pSL(strand, isbkLeft);

        int seq_len=0;

        for(int j=0;j<count;j++)
        {

            int pp1=transcripts[tranId].exonStarts[j];
            int pp2=transcripts[tranId].exonEnds[j];

            seq_len+=pp2-pp1;

            if(isbkLeft==1)
            {
                if(abs(pos-pp1)==best)
                {
                    name=transcripts[tranId].name;
                    pos1=pp1;
                    pos2=pp2;
                    if(strand==0)
                    {
                        exonNum=j+1;
                    }
                    else
                    {
                        exonNum=count-j;
                    }
                    getCodingAndBaseLeft(tranId, exonNum, isbkLeft, is_coding, baseLeft);
                    int cumu5pNT=0;
                    getCumu5pNT(tranId, exonNum, isbkLeft, cumu5pNT);
                    int coding_start=0;
                    if(is5p==1 && is_coding==1 && baseLeft>=0)
                    {
                        if(transcripts[tranId].cdsEnd>=pos1 && transcripts[tranId].cdsEnd<=pos2)
                            coding_start=pos2-transcripts[tranId].cdsEnd+1;
                        else
                            coding_start=((pos2-pos1)-baseLeft)%3+1;
                    }
//cout<<"push baseLeft "<<baseLeft<<endl;

//cout<<gid<<"\t"<<is5p<<"\t"<<chr<<"\t"<<strand<<"\t"<<pos1<<"\t"<<pos2<<"\t"<<exonNum<<"\t"<<is_coding<<"\t"<<baseLeft<<"\t"<<name<<"\t"<<coding_start<<"\t"<<cumu5pNT<<"11111111111111111111111"<<endl;
                    juncs.push_back(assign_junction(gid, is5p, chr, strand, pos1, pos2, exonNum, is_coding, baseLeft, name, coding_start,cumu5pNT));
                }

            }
            else
            {

                if(abs(pos-pp2)==best)
                {
                    name=transcripts[tranId].name;
                    pos1=pp1;
                    pos2=pp2;
                    if(strand==0)
                    {
                        exonNum=j+1;
                    }
                    else
                    {
                        exonNum=count-j;
                    }

                    getCodingAndBaseLeft(tranId, exonNum, isbkLeft, is_coding, baseLeft);
                    int cumu5pNT=0;
                    getCumu5pNT(tranId, exonNum, isbkLeft, cumu5pNT);
                    int coding_start=0;
                    if(is5p==1 && is_coding==1 && baseLeft>=0)
                    {
                        if(transcripts[tranId].cdsStart>=pos1 && transcripts[tranId].cdsStart<=pos2)
                            coding_start=transcripts[tranId].cdsStart+1-pos1;
                        else
                            coding_start=((pos2-pos1)-baseLeft)%3+1;
                    }
//cout<<"push2 baseLeft "<<baseLeft<<endl;
//cout<<gid<<"\t"<<is5p<<"\t"<<chr<<"\t"<<strand<<"\t"<<pos1<<"\t"<<pos2<<"\t"<<exonNum<<"\t"<<is_coding<<"\t"<<baseLeft<<"\t"<<name<<"\t"<<coding_start<<"\t"<<cumu5pNT<<"22222222222222222222222"<<endl;
                    juncs.push_back(assign_junction(gid, is5p, chr, strand, pos1, pos2, exonNum, is_coding, baseLeft, name, coding_start, cumu5pNT));
                }
            }

        }
    }

}////////
    return 0;

}
