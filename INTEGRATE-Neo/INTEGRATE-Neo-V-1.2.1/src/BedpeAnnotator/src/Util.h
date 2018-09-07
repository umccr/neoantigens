/*
 * Util.h
 *
 *  Created on: Apr 28, 2013
 *      Updated on Mar 4, 2016
 *      Author: jinzhang
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <MyTypes.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include <string>
#include <cstring>
#include <typeinfo>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <sstream>


using namespace std;


uint32_t getFilelength(char *file);
int readBlock(char * block, int length, FILE *infile);


extern map<int,char> intChar;
extern map<char,char> charChar;

extern map<string,char> tableAmino;


int InitialIntChar();


char getCharComp(char reada);
char getCharA(int reada);

int getPeptide(vector<char> & seq5p, vector<char> & seq, int start_pos, vector<char> & peptide, int & full, int & left);

std::vector<std::string> &my_split(const std::string &s, char delim, std::vector<std::string> &elems) ;

std::vector<std::string> my_split(const std::string &s, char delim) ;

double f_score(double sensitivity, double precision);

string my_db2string(double db);//db for double

int get_pseudo_counts(string & input, pseudo_counts_t & pct);

string vecToString(string sep, vector<string> & strs);
vector<string> twoPairVecToString(vector<string> & strs1, vector<string> & strs2);

#endif /* UTIL_H_ */
