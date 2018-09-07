/*
 * Bedpe.h
 * Mar 4 2016
 * Jin Zhang  
 */

#ifndef BEDPE_H_
#define BEDPE_H_

#include "MyTypes.h"
#include "Util.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <algorithm>

using namespace std;

class Bedpe
{
    private:
        vector<bedpe_t> bedpevec;
    public:
        Bedpe(){};
        int loadFromFile(char * filename);
        int size(){return bedpevec.size();}
        bedpe_t getBedpe(int i){ return bedpevec[i];}
        int getPos(bedpe_t & tt, uint32_t & tt_pos5, uint32_t & tt_pos3);
        int printBedpe(char * file);
        int isBedpeSame(bedpe_t & a, bedpe_t & b);
        int uniq();
        int insertCol(int index, string str);
};


#endif 

