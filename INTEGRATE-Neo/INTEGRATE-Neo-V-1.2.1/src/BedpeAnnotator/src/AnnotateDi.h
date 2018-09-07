/*
 * AnnotateDi.h
 *
 *  Created on: Aug 24, 2016
 *      Author: Jin Zhang
 */

#ifndef ANNOT_DI_H_
#define ANNOT_DI_H_

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
#include "Dinucleo.h"

class AnnotateDi {

public:

    int getCanonicalDi(bedpe_t & bp, Dinucleo & di ,Reference & ref);
    int annotateDi(Bedpe & bedpe, Dinucleo & di, Reference & ref);
};


#endif
