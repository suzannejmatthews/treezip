#include <iostream>
#include "global.h"
#include "hashfunc.hh"
#include "hash.hh"
#include <vector>

#ifndef TMERGE_H_
#define TMERGE_H_

void test_equivalence(string infile, string mergefile);
void setup_tmerge(string infilename, string mergefile, string outfile, int option, bool quiet);
bool null_case(string infilename, string mergefile, string outfile, int option);
void tmerge_parse_and_load(string filename, //name of file
			   unsigned int size, //number of unique trees
			   unsigned int offset, //offset to be used
			   unsigned int ** inverted_index, //the inverted index
			   unsigned int * inverted_index_sizes, //corresponding sizes
			   vector<string> & all_bipart, //vector of all the bipartitions 
			   unsigned int nt, //total number of trees
			   unsigned int encode_old, //size of old encoding 
			   unsigned int encode_new, //size of new encoding
 			   vector < vector<string> > &all_branch, //stores all the branches
			   vector < unsigned int > &branch_sizes, //stores all the branches
			   int * true_ids, //stores the "true" id values 
			   vector< vector<unsigned int> > alldups, //stores duplicates info 
			   unsigned int totalTrees,
			   unsigned int trueOffset); //stores the total number of trees in both files 

void tmerge_union(HashRFMap tmerge_hash, vector<unsigned int> & want, int * true_ids, int * chaindup, bool quiet);
//void tmerge_union(HashRFMap tmerge_hash, vector<unsigned int> & want, int * true_ids, vector< vector<unsigned int> > &alldups);
void tmerge_intersection(HashRFMap tmerge_hash, vector<unsigned int> & want, unsigned int offset, bool quiet);
void tmerge_setDifference(HashRFMap tmerge_hash, vector<unsigned int> & want, unsigned int offset, bool quiet);
void getdups(string filename, unsigned int &nt, vector<string> &dups);
#endif
