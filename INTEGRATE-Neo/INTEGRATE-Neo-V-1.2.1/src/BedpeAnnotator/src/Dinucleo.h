/*
 * Dinucleo.h
 * Aug 24 2016
 * Jin Zhang  
 */

#ifndef DINUCLEO_H_
#define DINUCLEO_H_

#include "MyTypes.h"
#include "Util.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <stdint.h>

using namespace std;

#include "Reference.h"


typedef struct
{
    char FiveOne;
    char FiveTwo;
    char ThreeOne;
    char ThreeTwo;
} splice_t;

class Dinucleo
{
    private:
        vector<splice_t> spliceVec;
    public:
        Dinucleo(){};
        int loadFromFile(char * filename);
        int size(){return spliceVec.size();}
        splice_t getSplice(int i){ return spliceVec[i];}
        int isCanoSplice(char fiveone, char fivetwo, char threeone, char threetwo);
        int isCanoSplice(string chr5, string strand5, uint32_t junc5, string chr3, string strand3, uint32_t junc3, Reference & ref);
};


#endif 

