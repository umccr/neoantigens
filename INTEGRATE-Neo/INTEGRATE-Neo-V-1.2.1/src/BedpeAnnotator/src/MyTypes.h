/*
 * MyTypes.h
 *
 *  Created on: Mar 4, 2016
 *      Author: Jin Zhang
 */

#ifndef MYTYPES_H_
#define MYTYPES_H_

#include <string>
#include <stdint.h>
#include <vector>

using namespace std;



typedef struct
{
	//int bin;
	string name;
	string chr;
	int strand;
	uint32_t txStart;
	uint32_t txEnd;
	uint32_t cdsStart;//Dec 7, 2015 start to use cdsStart and cdsEnd for in/out-of-frame and peptide
	uint32_t cdsEnd;
	int exonCount;
	uint32_t * exonStarts;
	uint32_t * exonEnds;
	//int score;
	string name2;

} transcript_t;


typedef struct
{
	uint32_t start;
	uint32_t end;
} exon_t;

typedef struct
{
	string chr;
	uint32_t start;
	uint32_t end;
	int strand;

	vector<int> transIds;
	vector<int> exonIds;
	int geneId;

} exon_map_t;




typedef struct
{
	string chr;
	int strand;
	uint32_t leftLimit;
	uint32_t rightLimit;
	string name2;
	vector<int> transIds;

	int fakeId;//-1 for no complex genes at multiple locations. 0-upper means yes. so that you see diff locations of the same gene;

	vector<int> anchors;


} gene_t;

typedef struct
{
    string chr1;
    uint32_t start1;
    uint32_t end1;
    string chr2;
    uint32_t start2;
    uint32_t end2;
    
    string name;
    double score;
    string strand1;
    string strand2;
    
    vector<string> others; 

} bedpe_t;

typedef struct
{
    string name;
    
    int num_res_trans;
    int num_truth_trans;    

    string sensitivity_t;
    string precision_t;
    string f_t;
    
    int num_res_gene;
    int num_truth_gene;

    string sensitivity_g;
    string precision_g;
    string f_g;
    
} evaluate_t;

typedef struct
{
    vector<int> ids1;
    vector<int> ids2;
} set_pair_t;

typedef struct
{
    int t_t;
    int f_t;
    int truth_t;

    int t_g;
    int f_g;
    int truth_g;
} pseudo_counts_t;

typedef struct
{
            int gId;
            int is5p;
            int isCoding;
            int coding_start; // the first full codon starts at position from pos1 or pos2 // has meaning if is5p==true
            int coding_left;  // 0 1 or 2, base left for the last codon of 5p or 0 1 2 for 3p, the first not full codon
            string chr;
            int strand;
            int pos1;
            int pos2;
            string name;
            int exonNum;
            int cumu5pNT;

} junction_t;

typedef struct
{
    int fusion_id;
    junction_t p5;
    junction_t p3;
    vector<char> seq1;
    vector<char> seq2;
    int isCanonical;
    int isInframe;
    string p5Type;
    string p3Type;
    string peptide;
    int left;
    int full;
    int isCanonicalDi;

} fusion_junction_t;

#endif
