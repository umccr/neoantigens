/*
 * Gene2.h
 *
 *  Created on: Jun 16, 2016
 *      Author: Jin Zhang
 */

#ifndef ANNOT_H_
#define ANNOT_H_

#include "MyTypes.h"
#include "Util.h"

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <string>
#include <stdint.h>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>


using namespace std;

#include "Bedpe.h"
#include "Reference.h"
#include "Gene2.h"


class Annotate {

public:

    int assignJunctions(Gene& g, Reference & ref, Bedpe & bdpe);
    int getJunctionType(Gene& g, Reference & ref, Bedpe & b);
    int getCanonical(fusion_junction_t & fjt, bedpe_t & b);
    int getJunctionPeptide(Gene& g, Reference & ref);
    int annotate(Bedpe & bedpe, Gene& g, Reference & ref);
};


#endif
