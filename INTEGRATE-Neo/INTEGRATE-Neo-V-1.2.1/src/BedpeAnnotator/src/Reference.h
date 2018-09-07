/*
 * Reference.h
 *
 *  Created on: Oct 13, 2011
 *      Author: jinzhang
 *  Edit and change: Jan 29, 2013
 *  Edit and change: Apr 27, 2013
 */

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <iostream>
#include <stdint.h>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <exception>
#include <typeinfo>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <map>

using namespace std;


#include "Util.h"



class Reference {
private:
	vector<string> nameList;//chr names with 'chr' removed.
	vector<uint32_t> posListLeft;//left pos of a ref seq
	vector<uint32_t> posListRight;//right pos of a ref seq
	map<string,int> nameIndex;//store the chrname and index
	int numOfRefSeq;//total num of ref seqs
	uint32_t * reference;//reference;
	uint32_t refLength;//length of reference
	char * referenceC;//reference in Char
	int isInt;


//	int loadReference(char * filename, uint32_t length);//function to load reference
//	//given chrname, return index in the vectors; this assume chrname is existed!!!


	int parseBlock(char * block, uint32_t & c_number,uint32_t length, int final);
	int	addRef(char * seq, uint32_t len);
	int chrNameToIndex(string chrname);
	int loadReference(char * filename, uint32_t length);//function to load reference

public:

	Reference();
	int loadRef(char * fileName);
	void test(char * file);//test

	string getCharName(int i);
	uint32_t getPosLeft(int i);
	uint32_t getPosRight(int i);

	int getListSize();


	uint32_t to_ref_pos(string chrname,uint32_t poschr);// given chr and pos on chr, return ref pos
	uint32_t to_ref_pos(int tid, uint32_t poschr);
	uint32_t to_chr_pos(uint32_t refpos, string & chrname);//given pos on ref, return chr pos, and given the chr name
	uint32_t to_chr_pos(int & tid, uint32_t refpos);
	uint32_t getRefLength();//return length of ref
	char getRefChar(uint32_t pos);

	int getBlock(uint32_t start, uint32_t end, char * block);


/*
	int getRefFromFile(char * file);//get ref from file
	int setChrUsed();//set chrUsed=1, if reference fasta use, for example chr1 not 1;
	int setChrByValue(int value);
	int getChrUsed();// return chrUsed
	uint32_t to_ref_pos(string chrname,uint32_t poschr);// given chr and pos on chr, return ref pos
	uint32_t to_chr_pos(uint32_t poschr, string & chrname);//given pos on ref, return chr pos, and given the chr name
	uint32_t getRefLength();//return length of ref
	char getRefChar(uint32_t pos);
	int getListSize();


	string getCharName(int i);
	uint32_t getPosLeft(int i);
	uint32_t getPosRight(int i);

	void test(char * file);//test

*/

	virtual ~Reference();

	int setIsInt(int isInt)
	{
		this->isInt=isInt;
		return 0;
	}

	int getIsInt()
	{
		return isInt;
	}


};



#endif /* REFERENCE_H_ */
