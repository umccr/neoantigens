/*
 * CutterByRule.h
 * Mar 6 2016
 * Jin Zhang
 */

#ifndef CutterByRule_H_
#define CutterByRule_H_

#include "MyTypes.h"
#include "Bedpe.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>

using namespace std;

class CutterByRule
{
    public:
        CutterByRule(){};
        int cut(Bedpe & in, Bedpe & out, string rule, string file); 
};


#endif
