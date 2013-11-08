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
#include "hashfunc.hh"
#include "hash.hh"

// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;

#define BITS                                                                           64
#define LEFT               								0
#define RIGHT              								1
#define ROOT               								2
#define LABEL_WIDTH        								3

// the c value of m2 > c*t*n in the paper
unsigned int C = 1000; 
unsigned int NEWSEED = 0; 
unsigned NUM_TREES = 0; // number of trees
unsigned NUM_TAXA  = 0; // number of taxa
bool WEIGHTED  = false; // unweighted
bool RATE  = false;
unsigned int P = 50; //percent compression or "thresholding" to be done
unsigned int UNIQUE_B = 0; //percent compression or "thresholding" to be done
HashRFMap vvec_hashrf;

//for FastHashRF
bool COMPRESSION  = false;
bool PRINT = false;
bool PRINT_L = false;
bool PRINT_HASH  = false;
bool MULTIFURCATING = false;
bool HASH_INPUT = false;
bool DATA = false;
bool AMAZING  = false;
bool UNIQUE = false;
bool MC = false;
bool HETERO = false;
bool HCHECK = false;
double update_time = 0;
double build_time = 0;
double clean_time = 0;
double total_time = 0;
double index_time = 0;

double print_time = 0;
unsigned int ALL_COUNT =0; //stores the count for all the trees
unsigned int ** helper;
unsigned int * helper_sizes;
unsigned int ** hashtable;
unsigned int * hash_lengths;
unsigned int tree_counter = 0;
unsigned int NUMBIPART = 0;
unsigned int MAXVAL = 0;
