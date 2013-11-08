//*****************************************************************/
//
// hash.cc -- Based on code written by Seung-Jin Sul
// Copyright (C) 2010 Suzanne J. Matthews
//             Department of Computer Science & Engineering
//             Texas A&M University
//             Contact: sjm@cse.tamu.edu
// Copyright (C) 2006-2009 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS IMPLEMENTATION
//		HashRFMap: Class for hashing bit
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details (www.gnu.org).
//
//*****************************************************************/

#include "hash.hh"
#include "hashfunc.hh"
#include <stdlib.h>

void 
HashRFMap::uhashfunc_init(
	unsigned int t, 
	unsigned int n, 
	unsigned int c)
{
	_HF.UHashfunc_init(t, n, c);	
}

void 
HashRFMap::uhashfunc_init_unique(
	unsigned int t, 
	unsigned int n, 
	unsigned int k, 
	unsigned int c,
	unsigned int newseed)
{
  _HF.UHashfunc_init_unique(t, n, k, c, newseed);	
}

void 
HashRFMap::uhashfunc_init(
	unsigned int t, 
	unsigned int n, 
	unsigned int c, 
	int32 newseed)
{
  _HF.UHashfunc_init(t, n, c, newseed);	
}

//Las-Vegas Randomized Hashing algorithm -- implemented by SJM
void 
HashRFMap::hashing_bs_lv(
	unsigned int treeIdx, 
	unsigned int numTaxa,
	unsigned long long hv1,
	unsigned long long hv2,
	float dist,
	bool w_option, 
	bool flip,
	bool * bs) 
{	
  ///////////////////////////////
  // double linked list
  ///////////////////////////////
  unsigned sizeVec = _hashtab2[hv1].size();
  unsigned int maxLength = 0;
  if (flip){
    for (unsigned int i = 0; i < numTaxa; i++){
      if (!bs[i])
	maxLength = i;
    }
  }
  else{
    for (unsigned int i = 0; i < numTaxa; i++){
      if (bs[i])
	maxLength = i;
    }
  }
  maxLength++;
  bool *tempbs = new bool[maxLength];
  if (flip){
    for (unsigned int i = 0; i < maxLength; i++)
      tempbs[i] = !bs[i];
  }
  else{
    for (unsigned int i = 0; i < maxLength; i++)
      tempbs[i] = bs[i];
  }
  if (sizeVec > 0) { //if there is a bipartition
    bool found = false;
    for (unsigned int i=0; i<sizeVec; ++i) {			
      if (_hashtab2[hv1][i]._hv2 == hv2) {
	found = true;
        for (unsigned int j = 0; j < maxLength; j++){
          if (tempbs[j] != _hashtab2[hv1][i]._bs[j]){  //if bitstrings are different
            found = false; //say it's no longer found
	    printf("bitstrings don't match! treeidx is: %u. faulty position is: %u\n", treeIdx, j);
	    exit(0);
	  }
        }
        if (found){ //if bitstrings match true*/
	  _hashtab2[hv1][i]._vec_treeidx.push_back(treeIdx);
	  if (w_option) //weighted option
	    _hashtab2[hv1][i]._vec_dist.push_back(dist);
	  //if (add)
	    _hashtab2[hv1][i]._count++;
	    delete [] tempbs;
	  break;
	}
      }
    }
    if (!found) { //if the tree is not found
      TREEIDX_STRUCT_T bk2;
      bk2._hv2 = hv2;
      bk2._vec_treeidx.push_back(treeIdx);
      bk2._bs = tempbs;
      bk2._bs_size = maxLength;
      //if (add)
	bk2._count = 1;
      if (w_option)
	bk2._vec_dist.push_back(dist);
      _hashtab2[hv1].push_back(bk2);
    }
  }
  else if (sizeVec == 0) {//if there is a bipartition, but sizeVec is 0
      TREEIDX_STRUCT_T bk2; //create a new bucket, and add it to it
      bk2._hv2 = hv2;
      bk2._bs = tempbs;
      bk2._bs_size = maxLength;
      //if (add)
	bk2._count = 1;
      bk2._vec_treeidx.push_back(treeIdx);
      if (w_option)
	bk2._vec_dist.push_back(dist);
      _hashtab2[hv1].push_back(bk2);
  }
}

// Montecarlo version -- does not store or check bitstrings 
void 
HashRFMap::hashing_bs_mc(
	unsigned int treeIdx, 
	unsigned int numTaxa,
	unsigned long long hv1,
	unsigned long long hv2,
	float dist,
	bool w_option) 
{	
  ///////////////////////////////
  // double linked list
  ///////////////////////////////
  unsigned sizeVec = _hashtab2[hv1].size();
  if (sizeVec > 0) { //if there is a bipartition
    bool found = false;
    for (unsigned int i=0; i<sizeVec; ++i) {			
      if (_hashtab2[hv1][i]._hv2 == hv2) {
	found = true;
	_hashtab2[hv1][i]._vec_treeidx.push_back(treeIdx);
	if (w_option) //weighted option
	  _hashtab2[hv1][i]._vec_dist.push_back(dist);
	_hashtab2[hv1][i]._count++;
	  break;
      }
    }
    if (!found) { //if the tree is not found
      TREEIDX_STRUCT_T bk2;
      bk2._hv2 = hv2;
      bk2._vec_treeidx.push_back(treeIdx);
      bk2._count = 1;
      if (w_option)
	bk2._vec_dist.push_back(dist);
      _hashtab2[hv1].push_back(bk2);
    }
  }
  else if (sizeVec == 0) {//if there is a bipartition, but sizeVec is 0
      TREEIDX_STRUCT_T bk2; //create a new bucket, and add it to it
      bk2._hv2 = hv2;
      bk2._count = 1;
      bk2._vec_treeidx.push_back(treeIdx);
      if (w_option)
	bk2._vec_dist.push_back(dist);
      _hashtab2[hv1].push_back(bk2);
  }
}

//hashing function for unique function -- implemented by SJM
void 
HashRFMap::hashing_bs_unique(
	unsigned int treeIdx, 
	unsigned int unique,
	unsigned long long hv1,
	unsigned long long hv2,
	bool *bs
) {	
  ///////////////////////////////
  // double linked list
  ///////////////////////////////
  unsigned sizeVec = _hashtab2[hv1].size();
   bool *tempbs = new bool[unique];
    for (unsigned int i = 0; i < unique; i++)
      tempbs[i] = bs[i];

  if (sizeVec > 0) { //if there is a bipartition
    bool found = false;
    for (unsigned int i=0; i<sizeVec; ++i) {			
      if (_hashtab2[hv1][i]._hv2 == hv2) {
	//check bitstring
	found = true;
        for (unsigned int j = 0; j  < unique; j++){
	  if (tempbs[j] != _hashtab2[hv1][i]._bs[j]) //if bitstrings are different
	    found = false; //say it's no longer found
        }
	if (found){
	  _hashtab2[hv1][i]._vec_treeidx.push_back(treeIdx);
	  _hashtab2[hv1][i]._count++;
	    delete [] tempbs;
	}
	break;
      }
    }
    if (!found) { //if the tree is not found
      TREEIDX_STRUCT_T bk2;
      bk2._hv2 = hv2;
      bk2._vec_treeidx.push_back(treeIdx);
      bk2._count = 1;
      bk2._bs = tempbs;
      _hashtab2[hv1].push_back(bk2);
    }
  }
  else if (sizeVec == 0) {//if there is a bipartition, but sizeVec is 0
    TREEIDX_STRUCT_T bk2; //create a new bucket, and add it to it
    bk2._hv2 = hv2;
    bk2._count = 1;
    bk2._bs = tempbs;
    bk2._vec_treeidx.push_back(treeIdx);
    _hashtab2[hv1].push_back(bk2);
  }
}

/*
void 
HashRFMap::hashrfmap_clear() {
  for (unsigned int i = 0; i < _hashtab2.size(); i++){
    for (unsigned int j = 0; j < _hashtab2[i].size(); j++){
      if (_hashtab2[i][j]._bs){
	//_hashtab2[i][j]._vec_treeidx.clear();
	//_hashtab2[i][j]._vec_dist.clear();
	delete _hashtab2[i][j]._bs;
	_hashtab2[i][j]._bs = NULL;
      }
    }
    //_hashtab2[i].clear();
  }
  //  _hashtab2.clear();
}
*/

