#include <iostream>
#include "global.h"
#include "label-map.hh"

#ifndef BMERGE_H_
#define BMERGE_H_

unsigned int validate_merge(string & infile, string & mergefile, LabelMap &lm, unsigned int option, unsigned int & total_trees);
void setup_bmerge(string infilename, string mergefile, int option);
void bmerge_union(string infilename, string mergefile, unsigned int offset, FILE * fout);
void bmerge_intersection(string infilename, string mergefile, unsigned int offset, FILE * fout);
void bmerge_setDifference(string infilename, string mergefile, unsigned int offset, FILE * fout);
string encode(unsigned int * new_vec, unsigned int newlength);
unsigned int decode(string encoded, unsigned int * found);
unsigned int count_ones(string bipart);
unsigned int get_bitstring_length(string bitstring);
void decode_bitstring(string bitstring, bool * bs, unsigned int maxLength);
string rle(bool * temp, unsigned int length);
void parse_and_get(string str, string check, unsigned int & var);
unsigned int get_ntrees(string str);
unsigned int get_unique(string str, unsigned int ntrees);
#endif
