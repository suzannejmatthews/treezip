//*****************************************************************/
/*
This file is part of TreeZip.

TreeZip is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

TreeZip is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TreeZip.  If not, see <http://www.gnu.org/licenses/>.

*/
/*****************************************************/
#include <iostream>
#include <limits.h>
#include "hash.hh"

#ifndef _GLOBAL_H_
#define _GLOBAL_H_

extern unsigned int C; 
extern unsigned int NEWSEED; 
extern unsigned NUM_TREES; // number of trees
extern unsigned NUM_TAXA; // number of taxa
extern bool WEIGHTED; // default: unweighted
extern bool RATE; // default: rf distance
extern unsigned int P;
extern unsigned int UNIQUE_B;
extern HashRFMap vvec_hashrf;

//for FastHashRF
extern bool COMPRESSION;
extern bool PRINT;
extern bool PRINT_L;
extern bool PRINT_HASH;
extern bool MULTIFURCATING;
extern bool HASH_INPUT;
extern bool DATA;
extern bool AMAZING;
extern bool UNIQUE;
extern bool MC;
extern bool HETERO;
extern bool HCHECK;

//timer functions
extern double update_time;
extern double build_time;
extern double clean_time;
extern double total_time;
extern double index_time;
extern double print_time;

//micc
extern unsigned int ALL_COUNT; //global count variable for thresholding
extern unsigned int ** helper;
//extern unsigned int ** full_index;
//extern unsigned int * index_sizes;
//extern vector < vector<unsigned int> > full_index; 
extern unsigned int * helper_sizes;
extern unsigned int ** hashtable;
extern unsigned int * hash_lengths;
extern unsigned int tree_counter;
extern unsigned int NUMBIPART;
extern unsigned int MAXVAL;

#endif
