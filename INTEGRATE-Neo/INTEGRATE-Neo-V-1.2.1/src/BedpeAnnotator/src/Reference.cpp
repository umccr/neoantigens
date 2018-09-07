/*
 * Reference.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: jinzhang
 *      Edit and change: Jan 29, 2013
 */

#include "Reference.h"


map<char,int> mask1;

int initialMask()
{
/*
	mask1['a']=0;
	mask1['c']=1;
	mask1['g']=2;
	mask1['t']=3;
	mask1['n']=4;
	mask1['A']=0;
	mask1['C']=1;
	mask1['G']=2;
	mask1['T']=3;
	mask1['N']=4;
*/
	mask1['A']=0;
        mask1['B']=4;
        mask1['C']=1;
        mask1['D']=4;
        mask1['E']=4;
        mask1['F']=4;
        mask1['G']=2;
        mask1['H']=4;
        mask1['I']=4;
        mask1['J']=4;
        mask1['K']=4;
        mask1['L']=4;
        mask1['M']=4;
        mask1['N']=4;
        mask1['O']=4;
        mask1['P']=4;
        mask1['Q']=4;
        mask1['R']=4;
        mask1['S']=4;
        mask1['T']=3;
        mask1['U']=4;
        mask1['V']=4;
        mask1['W']=4;
        mask1['X']=4;
        mask1['Y']=4;
	mask1['Z']=4;

        mask1['a']=0;
        mask1['b']=4;
        mask1['c']=1;
        mask1['d']=4;
        mask1['e']=4;
        mask1['f']=4;
        mask1['g']=2;
        mask1['h']=4;
        mask1['i']=4;
        mask1['j']=4;
        mask1['k']=4;
        mask1['l']=4;
        mask1['m']=4;
        mask1['n']=4;
        mask1['o']=4;
        mask1['p']=4;
        mask1['q']=4;
        mask1['r']=4;
        mask1['s']=4;
        mask1['t']=3;
        mask1['u']=4;
        mask1['v']=4;
        mask1['w']=4;
        mask1['x']=4;
        mask1['y']=4;
	mask1['z']=4;



	return 0;
}

inline int getMask(char a)
{
	return mask1[a];
}


map<int, char> charMap;

int initialMap()
{
	charMap[0]='A';
	charMap[1]='C';
	charMap[2]='G';
	charMap[3]='T';
	charMap[4]='N';
	return 0;
}

char getCharFromBinary(int a)
{
	return charMap[a];
}




Reference::Reference() {
	// TODO Auto-generated constructor stub
	reference=NULL;
	referenceC=NULL;
	numOfRefSeq=0;
	refLength=0;
	initialMap();
	initialMask();
}

int last; // the number of chars in the last uint32_t;
uint32_t numInRef; //the number of chars in reference;


/*
 *parse a block from the fasta file and store in bits in Reference.reference.
 *update names and coordinates.
 *
 */


int Reference::parseBlock(char * block, uint32_t & c_number, uint32_t length, int final)
{

    char * blockRM;
    if(isInt==1)
    	blockRM=new char [length+1024];
    char* line_buffer=block;
    char* p=NULL;
    char* nextN=NULL;
    int localc=0;
    while (1)
    {
    	nextN=strchr(line_buffer,'\n');

        if (line_buffer[0]=='>')
        {
            p=NULL;
            char tmp[nextN-line_buffer+2];
            tmp[nextN-line_buffer+1]='\0';
            memcpy(tmp,line_buffer,nextN-line_buffer+1);
            for(int i=0;i<strlen(tmp);i++)
            {
            	if(tmp[i]==' ' || tmp[i]=='\t' || tmp[i]=='\n')
            	{
            		p=tmp+i;
            		break;
            	}
            }
            char aa[1000];
            aa[p-tmp-1]='\0';

            strncpy(aa,tmp+1,p-tmp-1);

            char *trunkStr;

	    	if(strstr(aa,"chr"))
	    	{
	    		trunkStr=aa+3;
	    	}
	    	else
	    		trunkStr=aa;

            nameList.push_back(string(trunkStr));
            posListLeft.push_back(c_number+1);
            numOfRefSeq++;

            if(numOfRefSeq>1)
            {
            	cout<<"\t"<<posListLeft.back()-1<<endl;
            	posListRight.push_back(posListLeft.back()-1);
            	cout<<aa<<"\t"<<posListLeft[numOfRefSeq-1];
            }
            else
            {
            	cout<<aa<<"\t"<<posListLeft[0];
            }
            nameIndex[aa]=numOfRefSeq-1;
        }
        else
        {
            int len=nextN-line_buffer;//length of seq
            if(isInt==1)
	    {
	    	memcpy(blockRM+localc,line_buffer,len);
            }
	    else
	    {
		memcpy(referenceC+c_number,line_buffer,len);
  	    }
	    //addRef(line_buffer,len);
            c_number+=len;
            localc+=len;
            /*
            cout<<"pushed"<<endl;
            for(int k=0;k<len;k++)
            {
            	cout<<line_buffer[k];
            }
            cout<<endl;
             */
        }

        if(nextN-block==length-1)
        {
        	break;
        }
        else
        {
        	line_buffer=nextN+1;
        }

    }
    if(final==1)
    {
    	//cout<<"in final!"<<endl;
        if(numOfRefSeq==1)
	{
	//	cout<<"\t"<<posListLeft.back()-1<<endl;
               posListRight.push_back(posListLeft.back()-1);
        //       cout<<nameList[0]<<"\t"<<posListLeft[numOfRefSeq-1];
	}



        posListRight[numOfRefSeq-1]=c_number;
        refLength=posListRight[numOfRefSeq-1];
        cout<<"\t"<<posListRight[numOfRefSeq-1]<<endl;
    }
    if(isInt==1)
    {
	addRef(blockRM, localc);
    	delete [] blockRM;
    }
    return 0;
}





/*
 * add the len of chars to bits into reference
 *
 */

int Reference::addRef(char* seq, uint32_t len) {

	//cout<<"in addRef"<<endl;
/*
	for(int k=0;k<len;k++)
	{
		cout<<seq[k];
	}
	cout<<endl;
*/

	reference[numInRef]=0;
	for(int i=0;i<len;i++)
	{

		if(last==8 && i+7<len)
		{


			numInRef++;

			reference[numInRef]=(getMask(seq[i])<<28)|(getMask(seq[i+1])<<24)|(getMask(seq[i+2])<<20)|(getMask(seq[i+3])<<16)|(getMask(seq[i+4])<<12)|(getMask(seq[i+5])<<8)|(getMask(seq[i+6])<<4)|(getMask(seq[i+7]));

			last=8;

			i=i+7;
			continue;
		}

		uint32_t number=getMask(seq[i]);

		if(last==8)
		{
			last=0;
			numInRef++;
		}
		number=number<<4*(8-last-1);
		reference[numInRef]=reference[numInRef]|number;
		last++;

		/*
		if(last==8)
		{
			cout<<reference[numInRef]<<endl;
		}
		*/
	}
	return 0;
}

uint32_t blockSize=10000000;

int Reference::loadReference(char * filename, uint32_t length)
{
    last=0;
    numInRef=0;

    FILE *infile;

    infile = fopen(filename, "r");
    if (!infile) {
        cout<<"Couldn't open file for reading: "<<filename<<"."<<endl;
        exit(1);
    }



    char * block=new char [blockSize+1024];

    uint32_t c_number=0;

    uint32_t f_number=0;

    while(f_number!=length)
    {
    	uint32_t readed;
    	if(length-f_number>blockSize)
    	{
    		readed=readBlock(block, blockSize, infile);
    	}
    	else
    	{
    		readed=readBlock(block, length-f_number, infile);
    	}
    	f_number+=readed;
    	//cout<<f_number<<" "<<length<<endl;
		
	
    	if(f_number==length)
    	{
    		//cout<<"pass 11"<<endl;
    		parseBlock(block, c_number, readed,1);
    	}
    	else
    	{
    		//cout<<"pass 00"<<endl;
    		parseBlock(block, c_number, readed,0);
    	}
    }
/*
    if(f_number==length)
    {
    	uint32_t readed=readBlock(block, length-f_number, infile);
    	f_number+=readed;
    	parseBlock(block, c_number, readed,1);
    	cout<<"pass1"<<endl;
        parseBlock(block, c_number, readed, 1);
    }
*/
    delete [] block;
	return 0;
}


int Reference::loadRef(char* fileName) {

	// if exist reference, remove it;
	if(reference!=NULL)
	{
		try
		{
			delete [] reference;
		}
		catch (exception& e)
		{
		    cerr << "Trying to delete reference in getRefFromFile. exception caught: " << e.what() << endl;
		    return 1;
		}
	}

        if(referenceC!=NULL)
        {
                try
                {
                        delete [] referenceC;
                }
                catch (exception& e)
                {
                    cerr << "Trying to delete reference in getRefFromFile. exception caught: " << e.what() << endl;
                    return 1;
                }
        }
	numOfRefSeq=0;

	nameList.clear();
	posListLeft.clear();
	posListRight.clear();

	

	

	//allocate memory
	uint32_t length=getFilelength(fileName);
	try
	{
		if(isInt==1)
		{
			reference= new uint32_t [length/8+1024];
		}
		else
		{
			referenceC= new char [length];
		}
	}
	catch(exception& e)
	{
		cerr << "Trying to new reference in getRefFromFile. exception caught: " << e.what() << endl;
		return 1;
	}


	loadReference(fileName,length);

    /*
	int aa=refLength/8;
	int bb=refLength%8;
	for(int i=0;i<aa;i++)
	{
		//cout<<"reference"<<i<<" "<<reference[i]<<endl;
		for(int k=8;k>0;k--)
		{
			int a=(15&(reference[i]>>4*(k-1)));
			cout<<getCharFromBinary(a);
		}
	}
	if(bb>0)
	{
		for(int k=8;k>8-bb;k--)
		{
			int a=(15&(reference[aa]>>4*(k-1)));
			cout<<getCharFromBinary(a);
		}
	}
	cout<<endl;
	*/
	return 0;

}







void Reference::test(char * filename) {
	loadRef(filename);
/*
	cout<<getRefChar(1)<<endl;
	cout<<getRefChar(2)<<endl;
	cout<<getRefChar(3)<<endl;

	cout<<getRefChar(getRefLength()-2)<<endl;
	cout<<getRefChar(getRefLength()-1)<<endl;
	cout<<getRefChar(getRefLength())<<endl;

*/
}





string Reference::getCharName(int i) {
	return this->nameList[i];
}

uint32_t Reference::getPosLeft(int i) {
	return this->posListLeft[i];
}

uint32_t Reference::getPosRight(int i) {
	return this->posListRight[i];
}



int Reference::getListSize() {
	return this->nameList.size();
}


int Reference::chrNameToIndex(string chrname)
{
	try
	{
		return nameIndex.at(chrname);
	}
	catch(exception& e)
	{
		cerr << "Inquire a not existed reference sequence in chrNameToIndex. exception caught: " << e.what() << endl;
		return 1;
	}

	return 0;
}

uint32_t Reference::to_ref_pos(string chrname,uint32_t poschr)
{

	int index=chrNameToIndex(chrname);
	if(index<0 || index>numOfRefSeq-1)
	{
		cerr << chrname<<" "<<poschr<<" in to_ref_pos, wrong index of chromosome provided"<<endl;
		exit(1);
	}
	return posListLeft[index]+poschr-1;

}

uint32_t Reference::to_ref_pos(int tid,uint32_t poschr)
{

	int index=tid;
	if(index<0 || index>numOfRefSeq-1)
	{
		cerr << tid<<" "<<poschr<<" in to_ref_pos, wrong index of chromosome provided"<<endl;
		exit(1);
	}
	return posListLeft[index]+poschr-1;

}



uint32_t Reference::to_chr_pos(uint32_t posref, string & chrname)
{
	vector<uint32_t>::iterator low,up;
	if(posref>refLength)
	{
		cerr<<"pos "<<posref<<" larger than ref Length in to_chr_p"<<endl;
		exit(1);
	}
	up= upper_bound (posListLeft.begin(), posListLeft.end(), posref);// up > posref
	low=lower_bound (posListLeft.begin(), posListLeft.end(), posref);//low >=posref


	int index;
	if(int(up-posListLeft.begin())-int(low-posListLeft.begin()+1)<=1)
	{
		index=low-posListLeft.begin()-1;
	}
	else
	{
		cerr<<"didn't finder where the poschr in to_chr_pos"<<endl;
		exit(1);
	}

	if(index>0)
	{
		chrname=nameList[index];
		return posref-posListRight[index-1];
	}
	else
	{
		chrname=nameList[0];
		return posref;

	}
}


uint32_t Reference::to_chr_pos(int& tid, uint32_t posref) {
	vector<uint32_t>::iterator low,up;
	if(posref>refLength)
	{
		cerr<<"pos "<<posref<<" larger than ref Length in to_chr_p"<<endl;
		exit(1);
	}
	up= upper_bound (posListLeft.begin(), posListLeft.end(), posref);// up > posref
	low=lower_bound (posListLeft.begin(), posListLeft.end(), posref);//low >=posref


	int index;
	if(int(up-posListLeft.begin())-int(low-posListLeft.begin()+1)<=1)
	{
		index=low-posListLeft.begin()-1;
	}
	else
	{
		cerr<<"didn't finder where the poschr in to_chr_pos"<<endl;
		exit(1);
	}

	if(index>0)
	{
		tid=index;
		return posref-posListRight[index-1];
	}
	else
	{
		tid=index;
		return posref;

	}
}





uint32_t Reference::getRefLength()
{
	return refLength;
}


char Reference::getRefChar(uint32_t pos) {
	//return reference[pos-1];

	if(isInt==1)
	{
		int aa=pos/8;
		int bb=pos%8;

		if(bb==0)
		{
			int a=15&reference[aa-1];
			return getCharFromBinary(a);
		}
		if(bb>0)
		{

			int a=(15&(reference[aa]>>4*(8-bb)));
			return getCharFromBinary(a);

		}

	}
	else
		return getCharFromBinary(getMask(referenceC[pos-1]));

}


Reference::~Reference() {
	// TODO Auto-generated destructor stub
	try
	{
		if(isInt==1)
			delete [] reference;
		else
			delete [] referenceC;
	}
	catch (exception& e)
	{
	    cerr << "Trying to delete reference in ~Reference. exception caught: " << e.what() << endl;
	}
}



int Reference::getBlock(uint32_t start, uint32_t end, char * block)
{
	if(isInt==1)
	{
		cerr<<"Reference is represented by integer, you are asking from char"<<endl;
		exit(0);
	}
	
	memcpy(block,referenceC+start-1,end-start+1);
	return 0;
}
