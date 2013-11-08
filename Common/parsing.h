//*****************************************************************/
/*
This is Tree*Zip, a compression software for phylogenetic trees. 
It is based on HashBase, which is in turn based off of HashRF and HashCS,
developed by SeungJin Sul

(c) 2009 Tree*Zip: Suzanne Matthews
(c) 2009 HashBase: Suzanne Matthews
(c) 2009 HashRF : SeungJin Sul 
(c) 2009 HashCS : SeungJin Sul

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

#ifndef PARSING_H_
#define PARSING_H_

#include "label-map.hh"
#include "global.h"

extern "C" {
#include <newick.h>
}

void GetTaxaLabels(NEWICKNODE *node, LabelMap &lm);
void initialize_hashtable(unsigned long long &M1, unsigned long long &M2);

bool *
dfs_compute_hash(
  NEWICKNODE* startNode,
  LabelMap &lm,
  HashRFMap &vvec_hashrf,
  unsigned treeIdx,
  unsigned &numBitstr,
  unsigned long long m1,
  unsigned long long m2, 
  NEWICKNODE* parent);

bool *
dfs_compute_hash_mc(
  NEWICKNODE* startNode,
  LabelMap &lm,
  HashRFMap &vvec_hashrf,
  unsigned treeIdx,
  unsigned &numBitstr,
  unsigned long long m1,
  unsigned long long m2, 
  NEWICKNODE* parent);

bool * 
dfs_compute_hash_unrooted2( 
  NEWICKNODE* startNode, 
  LabelMap &lm,
  HashRFMap &vvec_hashrf,
  unsigned treeIdx,
  unsigned &numBitstr,
  unsigned long long m1,
  unsigned long long m2,
  bool &processed,
  NEWICKNODE* parent);

bool * 
dfs_compute_hash_unrooted_mc( 
  NEWICKNODE* startNode, 
  LabelMap &lm,
  HashRFMap &vvec_hashrf,
  unsigned treeIdx,
  unsigned &numBitstr,
  unsigned long long m1,
  unsigned long long m2,
  bool &processed,
  NEWICKNODE* parent);

string itos(int i);	

int toDecimal(bool * temp);

void toBinary(bool * output, int input);
void encode(unsigned long long x);
void encode2(unsigned int x);
void encode3(unsigned int & x, bool contin, FILE *fp);

LabelMap collect_labels(string file1);
LabelMap collect_labels_light(string file1);

#endif
