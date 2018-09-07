/*
 * Gene2.h
 *
 *  Created on: Mar 4, 2016
 *  Updated on: Oct 24, 2016
 *      Author: Jin Zhang
 */

#ifndef GENE_H_
#define GENE_H_

#include "MyTypes.h"
#include "Util.h"

#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>


using namespace std;


class Gene {
private:
	vector<gene_t> genes;
	vector<transcript_t> transcripts;

public:
	Gene();
	virtual ~Gene();


        /*load transcripts from input annotation file with 11 columns.*/
	int loadGenesFromFile(char * file);
	int setGene();



        /*given one coordinate, return all the genes the coordinate is in*/
	int isInGene(string chr, uint32_t pos, vector<int> & geneIds);
        
	gene_t * getGene(int index);

	int getStrand(int geneId);

	string getName2(int geneId);
		

	int getSize(){return genes.size();}

	uint32_t getLimitLeft(int geneId);
	uint32_t getLimitRight(int geneId);

        /*given hugo name of a gene, return the ids*/
	int getIndex(string name, vector<int> & ids);

        /*for choosing exons*/
        int getBestExon2(int gid, int pos, int isbkLeft, vector<junction_t> & juncs);
        int getBestDiff(int gid, int pos, int isbkLeft);
        int getBestExon3(vector<int> & gids, int pos, int isbkLeft, vector<junction_t> & juncs);
        int getBestDiff2(vector<int> & gids, int pos, int isbkLeft);

        int isAt5p(int gid, int isbkLeft);
        int isAt5pSL(int strand, int isbkLeft);
        int getCodingAndBaseLeft(int tranId, int exonNum, int isbkLeft, int & isCoding, int & baseLeft);
        int getCumu5pNT(int tranId, int exonNum, int isbkLeft, int & codingNT); 

};


#endif /* GENE_H_ */
