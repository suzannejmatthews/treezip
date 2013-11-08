//*****************************************************************/
/*
This is TreeZip, a compression software for phylogenetic trees. 
It is based on HashRF and HashCS, developed by SeungJin Sul

(c) 2010 TreeZip: Suzanne Matthews
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
/*****************************************************/

#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <bitset>
#include <string>
#include <sstream>
#include <algorithm>
#include "label-map.hh"
#include "hashfunc.hh"
#include "hash.hh"
#include "SCTree.h"
#include "parsing.h"
#include "global.h"

using namespace std;


#define add_of(c,a,b) ({ \
  typeof(a) __a=a; \
  typeof(b) __b=b; \
  (__b)<1 ? \
    ((__MIN(typeof(c))-(__b)<=(__a)) ? assign(c,__a+__b):1) : \
    ((__MAX(typeof(c))-(__b)>=(__a)) ? assign(c,__a+__b):1); \
})

#define __HALF_MAX_SIGNED(type) ((type)1 << (sizeof(type)*8-2))
#define __MAX_SIGNED(type) (__HALF_MAX_SIGNED(type) - 1 + __HALF_MAX_SIGNED(type))
#define __MIN_SIGNED(type) (-1 - __MAX_SIGNED(type))

#define __MIN(type) ((type)-1 < 1?__MIN_SIGNED(type):(type)0)
#define __MAX(type) ((type)~__MIN(type))

#define assign(dest,src) ({ \
  typeof(src) __x=(src); \
  typeof(dest) __y=__x; \
  (__x==__y && ((__x<1) == (__y<1)) ? (void)((dest)=__y),0:1); \
})

void
GetTaxaLabels(
  NEWICKNODE *node,
  LabelMap &lm)
{
  if (node->Nchildren == 0) {
    string temp(node->label);
    lm.push(temp);
  }
  else
    for (int i=0;i<node->Nchildren;++i)
      GetTaxaLabels(node->child[i], lm);
}

void PopulateTaxaLabels( NEWICKNODE *node, LabelMap &lm) {
  if (node->Nchildren == 0) { //if leaf node
    string temp(node->label);
    lm.add(temp);
  }
  else
    for (int i=0;i<node->Nchildren;++i)
      PopulateTaxaLabels(node->child[i], lm);
}

void initialize_hashtable(unsigned long long &M1, unsigned long long &M2) { 
  assert(NUM_TREES != 0 );
  //generate random numbers
  if (NEWSEED != 1000) { 
    vvec_hashrf.uhashfunc_init(NUM_TREES, NUM_TAXA, C, NEWSEED);
  }
  else { 
    vvec_hashrf.uhashfunc_init(NUM_TREES, NUM_TAXA, C);
  }
  M1 = vvec_hashrf._HF.getM1();
  M2 = vvec_hashrf._HF.getM2();

  //cout << "M1 is: " << M1 << endl;
  //cout << "M2 is: " << M2 << endl;
  vvec_hashrf._hashtab2.resize(M1);
  /*cout << "A1 random numbers:" << endl;
  for (unsigned int i = 0; i < NUM_TAXA; ++i){
    cout << "A1[" << i << "]= " << vvec_hashrf._HF.getA1(i) << endl;
  }
  cout << endl << endl;
  cout << "A2 random numbers:" << endl;
  for (unsigned int i = 0; i < NUM_TAXA; ++i){
    cout << "A2[" << i << "]= " << vvec_hashrf._HF.getA2(i) << endl;
  }
  cout << endl;*/
  
  /*vector<unsigned long long> random_nums;                                                                              
  //cout << "A1 random numbers:" << endl;                                                                                   
  unsigned long long tempr;                                                                                               
  for (unsigned int i = 0; i < NUM_TAXA; i++){                                                             
    tempr = vvec_hashrf._HF.getA1(i);                                                                                  
    //cout << "A1[" << i << "]= " << tempr << endl;                                                                        
    random_nums.push_back(tempr);                                                                                        
  }                                                                                                                      
  sort(random_nums.begin(), random_nums.end());                                                                          
  unsigned int real_unique = 0;                                                                                          
  tempr = -1;                                                                                                            
  for (unsigned int i = 0; i < NUM_TAXA; i++){                                                                
    if (tempr != random_nums[i]){                                                                                        
      real_unique++;                                                                                                     
      tempr = random_nums[i];                                                                                            
    }                                                                                                                     
  }                                                                                                                       
  cout << "number of unique random numbers (A1): " << real_unique << endl;                                               
  random_nums.clear();                                                                                                
  
  //cout << "A2 random numbers:" << endl;                                                                                   
  for (unsigned int i = 0; i < NUM_TAXA; i++){                                                                 
    tempr = vvec_hashrf._HF.getA2(i);                                                                                    
    //cout << "A2[" << i << "]= " << tempr << endl;                                                                        
    random_nums.push_back(tempr);                                                                                        
  }                                                                                                                       
  sort(random_nums.begin(), random_nums.end());                                                                           
  tempr = -1;                                                                                                             
  real_unique = 0;                                                                                                        
  for (unsigned int i = 0; i < NUM_TAXA; i++){                                                                 
    if (tempr != random_nums[i]){                                                                                         
      real_unique++;                                                                                                      
      tempr = random_nums[i];                                                                                            
    }                                                                                                                    
  }                                         
  cout << "number of unique random numbers (A2): " << real_unique << endl;
  exit(0);*/
}
 
//--- dfs traversal (rooted)
bool *
dfs_compute_hash(
  NEWICKNODE* startNode,
  LabelMap &lm,
  HashRFMap &vvec_hashrf,
  unsigned treeIdx,
  unsigned &numBitstr,
  unsigned long long m1,
  unsigned long long m2,
  NEWICKNODE* parent)
{
  if (HETERO && !HCHECK)
    return NULL;
  // If the node is leaf node, just set the place of the taxon name in the bit string to '1' 
  //bitstrings here are not used to compute hash values -- however, we do collect them.
  startNode->myparent = parent;
  //fprintf(stderr, "AT NODE: "); //debugging code
  if (startNode->Nchildren == 0) { // leaf node
    string temp(startNode->label);
    unsigned int idx = lm[temp];
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i  =0; i < NUM_TAXA; i++)
      bs[i] = 0;
    bs[idx] = true;
    ///debugging code
    /*cout << "leaf node bitstring is: " << endl;
       for (unsigned int i = 0; i < NUM_TAXA; i++)
       cout << bs[i];
       cout << endl;*/
    // Implicit BPs /////////////////////////
    // Set the hash values for each leaf node.
     unsigned long long temp1 = vvec_hashrf._HF.getA1(idx);
     unsigned long long temp2 = vvec_hashrf._HF.getA2(idx);
     temp1 = temp1 % m1;
     temp2 = temp2 % m2;
     startNode->hv1 = temp1;
     startNode->hv2 = temp2;
     

    int numOnes = 0;
    for (unsigned int i= 0; i < NUM_TAXA; i++)
      numOnes+=bs[i];

    if (numOnes > 1 || numOnes == 0) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring creation error\n";
      exit(0);
    }
    float dist = 0.0;
    if (WEIGHTED || HETERO){ 
      //if we are collecting branch lengths
      if (WEIGHTED)
	dist = startNode->weight; //grab the value
      vvec_hashrf.hashing_bs_lv(treeIdx, NUM_TAXA, temp1, temp2, dist, WEIGHTED,false, bs);
    }
    return bs;
  }
  else {
    bool *ebs[startNode->Nchildren];
    ///fprintf(stderr, "I am an internal node!\n"); //debugging code
    for (int i=0; i<startNode->Nchildren; ++i) {
      ebs[i] = dfs_compute_hash(startNode->child[i], lm, vvec_hashrf, treeIdx, numBitstr, m1, m2, startNode);
    }
    
    bool *bs = new bool[NUM_TAXA];
    for (unsigned int i  =0; i < NUM_TAXA; i++)
      bs[i] = 0;

    ///debugging code
    /*fprintf(stderr, "the bitstrings I'm going to be comparing are:\n");
    for (int i=0; i<startNode->Nchildren; ++i) {
      fprintf(stderr, "%s\n", (startNode->child[i])->label);
      }
      fprintf(stderr, "\n");*/

    for (int i=0; i<startNode->Nchildren; ++i) {
      //cerr << "bitstring: " << *(ebs[i]) << endl; //debugging code
 	if (ebs[i]) {
	  for (unsigned int j = 0; j < NUM_TAXA; j++)
	    bs[j] |= ebs[i][j];
	  delete [] ebs[i];
	  ebs[i] = NULL;
	}
	else {
	  if (!HETERO){
	    cout << "ERROR: null bitstring\n";
	    exit(0);
	  }
	}
    }
    ///debugging code
    /*fprintf(stderr, "new bitstring is:\n");
    for (unsigned int j = 0; j < NUM_TAXA; j++){
      cerr << bs[j];
    }
    fprintf(stderr, "\n");*/


   int numOnes = 0;
   for (unsigned int i= 0; i < NUM_TAXA; i++)
     numOnes+=bs[i];

    if (numOnes < 1) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring OR error\n";
      exit(0);
    }

    // For weighted RF
    float dist = 0.0;
    if (WEIGHTED)
      dist = startNode->weight;
    else
      dist = 1;
    
    ++numBitstr;
    
    // Implicit BPs ////////////
    // After an internal node is found, compute the hv1 and hv2
    unsigned long long temphv1=0;
    unsigned long long temphv2=0;
    
    for (int i=0; i<startNode->Nchildren; ++i) {    	
      unsigned long long t1 = temphv1;
      unsigned long long t2 = temphv2;
      unsigned long long h1 = startNode->child[i]->hv1;
      unsigned long long h2 = startNode->child[i]->hv2;
      // Check overflow  
      if ( add_of(temphv1, t1, h1) ) {
	cout << "ERROR: ullong add overflow!!!\n"; 
	cout << "t1=" << t1 << " h1=" << h1 << " t1+h1=" << t1+h1 << endl;
	exit(0);
      }
      if ( add_of(temphv2, t2, h2) ) {
	cout << "ERROR: ullong add overflow!!!\n"; 
	cout << "t2=" << t2 << " h2=" << h2 << " t2+h2=" << t2+h2 << endl;
	exit(0);
      }
    }
    
    unsigned long long temp1 = temphv1 % m1;
    unsigned long long temp2 = temphv2 % m2;
    startNode->hv1 = temp1;
    startNode->hv2 = temp2;
    
    // Store bitstrings in hash table
    if (parent != NULL){
      vvec_hashrf.hashing_bs_lv(treeIdx, NUM_TAXA, startNode->hv1, startNode->hv2, dist, WEIGHTED, false, bs);   // without TYPE-III using n-bits (Hash-RF)
    }
    return bs;
  } //end else
}

bool *
dfs_compute_hash_mc(
  NEWICKNODE* startNode,
  LabelMap &lm,
  HashRFMap &vvec_hashrf,
  unsigned treeIdx,
  unsigned &numBitstr,
  unsigned long long m1,
  unsigned long long m2,
  NEWICKNODE* parent)
{
  // If the node is leaf node, just set the place of the taxon name in the bit string to '1' 
  //bitstrings here are not used to compute hash values -- however, we do collect them.
  if (HETERO && !HCHECK)
    return NULL;
  startNode->myparent = parent;
  ///fprintf(stderr, "AT NODE: ");
  if (startNode->Nchildren == 0) { // leaf node
    string temp(startNode->label);
    unsigned int idx = lm[temp];
    ///fprintf(stderr, "At leaf node: %s\n", startNode->label);
     bool * bs = new bool[NUM_TAXA];
     for (unsigned int i  =0; i < NUM_TAXA; i++)
       bs[i] = 0;
     bs[idx] = true;

     /*(cout << "leaf node bitstring is: " << endl;
       for (unsigned int i = 0; i < NUM_TAXA; i++)
       cout << bs[i];
       cout << endl;*/
    // Implicit BPs /////////////////////////
    // Set the hash values for each leaf node.
     unsigned long long temp1 = vvec_hashrf._HF.getA1(idx);
     unsigned long long temp2 = vvec_hashrf._HF.getA2(idx);
     temp1 = temp1 % m1;
     temp2 = temp2 % m2;
     startNode->hv1 = temp1;
     startNode->hv2 = temp2;
     

    int numOnes = 0;
    for (unsigned int i= 0; i < NUM_TAXA; i++)
      numOnes+=bs[i];

    if (numOnes > 1 || numOnes == 0) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring creation error\n";
      exit(0);
    }
    float dist = 0.0;
    if (WEIGHTED){ 
      //if we are collecting branch lengths
      dist = startNode->weight; //grab the value
      vvec_hashrf.hashing_bs_mc(treeIdx, NUM_TAXA, temp1, temp2, dist, WEIGHTED);
    }
    return bs;
  }
  else {
    bool *ebs[startNode->Nchildren];
    ///fprintf(stderr, "I am an internal node!\n");
    for (int i=0; i<startNode->Nchildren; ++i) {
      ebs[i] = dfs_compute_hash_mc(startNode->child[i], lm, vvec_hashrf, treeIdx, numBitstr, m1, m2, startNode);
    }
    
    bool *bs = new bool[NUM_TAXA];
    for (unsigned int i  =0; i < NUM_TAXA; i++)
      bs[i] = 0;

    ///fprintf(stderr, "the bitstrings I'm going to be comparing are:\n");
    /*for (int i=0; i<startNode->Nchildren; ++i) {
      fprintf(stderr, "%s\n", (startNode->child[i])->label);
      }
      fprintf(stderr, "\n");*/

    for (int i=0; i<startNode->Nchildren; ++i) {
      //cerr << "bitstring: " << *(ebs[i]) << endl;
 	if (ebs[i]) {
	  for (unsigned int j = 0; j < NUM_TAXA; j++)
	    bs[j] |= ebs[i][j];
	  delete [] ebs[i];
	  ebs[i] = NULL;
	}
	else {
	  if (!HETERO){
	    cout << "ERROR: null bitstring\n";
	    exit(0);
	  }
	}
    }

    ///fprintf(stderr, "new bitstring is:\n");
    ///for (unsigned int j = 0; j < NUM_TAXA; j++){
    ///  cerr << bs[j];
    ///}
    ///fprintf(stderr, "\n");


   int numOnes = 0;
   for (unsigned int i= 0; i < NUM_TAXA; i++)
     numOnes+=bs[i];

    if (numOnes < 1) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring OR error\n";
      exit(0);
    }

    // For weighted RF
    float dist = 0.0;
    if (WEIGHTED)
      dist = startNode->weight;
    else
      dist = 1;
    
    ++numBitstr;
    
    // Implicit BPs ////////////
    // After an internal node is found, compute the hv1 and hv2
    unsigned long long temphv1=0;
    unsigned long long temphv2=0;
    
    for (int i=0; i<startNode->Nchildren; ++i) {    	
      unsigned long long t1 = temphv1;
      unsigned long long t2 = temphv2;
      unsigned long long h1 = startNode->child[i]->hv1;
      unsigned long long h2 = startNode->child[i]->hv2;
      // Check overflow  
      if ( add_of(temphv1, t1, h1) ) {
	cout << "ERROR: ullong add overflow!!!\n"; 
	cout << "t1=" << t1 << " h1=" << h1 << " t1+h1=" << t1+h1 << endl;
	exit(0);
      }
      if ( add_of(temphv2, t2, h2) ) {
	cout << "ERROR: ullong add overflow!!!\n"; 
	cout << "t2=" << t2 << " h2=" << h2 << " t2+h2=" << t2+h2 << endl;
	exit(0);
      }
    }
    
    unsigned long long temp1 = temphv1 % m1;
    unsigned long long temp2 = temphv2 % m2;
    startNode->hv1 = temp1;
    startNode->hv2 = temp2;
    
    // Store bitstrings in hash table
    if (parent != NULL){
      vvec_hashrf.hashing_bs_mc(treeIdx, NUM_TAXA, startNode->hv1, startNode->hv2, dist, WEIGHTED);   // without TYPE-III using n-bits (Hash-RF)
    }
    return bs;
  } //end else
}

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
  NEWICKNODE* parent)
{

  if (HETERO && !HCHECK)
   return NULL;
  // If the node is leaf node, just set the place of the taxon name in the bit string to '1'
  // and push the bit string into stack 
  ///fprintf(stderr, "AT NODE: ");
  startNode->myparent = parent;
  if (startNode->Nchildren == 0) { // leaf node
    string temp(startNode->label);
    unsigned int idx = lm[temp];
    ///fprintf(stderr, "At leaf node: %s\n", startNode->label);
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i =0; i < NUM_TAXA; i++)
      bs[i] = 0;
    bs[idx] = true;

    // Implicit BPs /////////////////////////
    // Set the hash values for each leaf node.
    unsigned long long temp1 = vvec_hashrf._HF.getA1(idx);
    unsigned long long temp2 = vvec_hashrf._HF.getA2(idx);
    startNode->hv1 = temp1;
    startNode->hv2 = temp2;

    temp1 = temp1 % m1;
    temp2 = temp2 % m2;
    int numOnes = 0;
    for (unsigned int i  =0; i < NUM_TAXA; i++)
      numOnes += bs[i];
    if (numOnes > 1 || numOnes == 0) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring creation error\n";
      exit(0);
    }

    return bs;
  }
  else { //if we are not at a leaf node
     //perform DFSn all the children
    bool * ebs[startNode->Nchildren];
    ///fprintf(stderr, "I am an internal node!\n");
    for (int x = 0; x < startNode->Nchildren; ++x) { 
      ebs[x] = dfs_compute_hash_unrooted2(startNode->child[x],lm, vvec_hashrf, treeIdx, numBitstr, m1, m2, processed, startNode);
    }
            
    // At this point, we find a bipartition. 
    // Thus, OR the bitstrings and make a bit string for the bipartition
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    ///fprintf(stderr, "the bitstrings I'm going to be comparing are:\n");
    //for (int i=0; i<startNode->Nchildren; ++i) {
     // fprintf(stderr, "%s\n", (startNode->child[i])->label);
     // }
      //fprintf(stderr, "\n");

    for (int i=0; i<startNode->Nchildren; ++i) {
      //cerr << "bitstring: " << *(ebs[i]) << endl;
 	if (ebs[i]) {
	  for (unsigned int j =0; j < NUM_TAXA; j++)
	    bs[j] |= ebs[i][j];
	  delete ebs[i];
	  ebs[i] = NULL;
	}
	else {
	  if (!HETERO){
	    cout << "ERROR: null bitstring\n";
	    exit(0);
	  }
	}
    }
    
    //fprintf(stderr, "new bitstring is:\n");
    //for (unsigned int j = BITSETSZ-1; j > BITSETSZ-1-NUM_TAXA; j--){
     // cerr << (*bs)[j];
      //}
     // fprintf(stderr, "\n");

    int numOnes = 0;
    for (unsigned int i  =0; i < NUM_TAXA; i++)
      numOnes += bs[i];
    if (numOnes < 1) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring OR error\n";
      exit(0);
    }

    unsigned long long temp1 = 0;
    unsigned long long temp2 = 0;
    
    for (int i = 0; i < startNode->Nchildren; ++i) { 
      temp1 += startNode->child[i]->hv1;
      temp2 += startNode->child[i]->hv2;
    }
    startNode->hv1 = temp1 % m1;
    startNode->hv2 = temp2 % m2;

    float dist = 0.0;
    if (WEIGHTED)
      dist = startNode->weight;
    else
      dist = 1;
  
    bool toprocess = true;
    if ( (parent == NULL) || ( (parent->myparent == NULL) && (parent->Nchildren == 2) && (processed) ) ){
      ///fprintf(stderr, "we don't collect this bitstring!\n");
      toprocess = false;
    }
    else{
      if ( (parent->myparent == NULL) && (parent->Nchildren == 2) && (!processed)){ 
	  if (bs[0] != 0){ //don't process if it's oriented wrong
	    toprocess = false;
	    ///fprintf(stderr, "we don't collect this bitstring.. it's oriented wrong!\n");
	  }
	  else{
	    processed = true;
	  }
      }
      else{
	if ( bs[0] != 0 ){
	  ///fprintf(stderr, "need to flip bitstring!\n");
	  toprocess = false;

	  //calculate new h1 h2 values	    
	  unsigned long long temp1 = 0;
	  unsigned long long temp2 = 0;
	  unsigned int idx = 0;
	  for (unsigned int j = 0; j < NUM_TAXA; j++){
	    if ( !bs[j] ) {
	      idx = j;
	      temp1 += vvec_hashrf._HF.getA1(idx);
	      temp2 += vvec_hashrf._HF.getA2(idx);
	    }
	  }
	  temp1 = temp1 % m1;
	  temp2 = temp2 % m2;
	  if (numBitstr < NUM_TAXA-3){
	    ///fprintf(stderr, "adding flipped bitstring to hash table\n");
	    vvec_hashrf.hashing_bs_lv(treeIdx, NUM_TAXA, temp1, temp2, dist, WEIGHTED, true, bs); //flip bitstring
	    ++numBitstr;
	  }
	}
      }
    }    
    if (toprocess){
      if (numBitstr < NUM_TAXA-3){
	///fprintf(stderr, "adding bitstring to hash table!\n");
	vvec_hashrf.hashing_bs_lv(treeIdx, NUM_TAXA, startNode->hv1, startNode->hv2, dist, WEIGHTED, false, bs);   //no flipping needed
	++numBitstr;
      }
    }
    return bs;
  }

}

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
  NEWICKNODE* parent)
{
  // If the node is leaf node, just set the place of the taxon name in the bit string to '1'
  // and push the bit string into stack 
  ///fprintf(stderr, "AT NODE: ");
  startNode->myparent = parent;
  if (startNode->Nchildren == 0) { // leaf node
    string temp(startNode->label);
    unsigned int idx = lm[temp];
    ///fprintf(stderr, "At leaf node: %s\n", startNode->label);
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i =0; i < NUM_TAXA; i++)
      bs[i] = 0;
    bs[idx] = true;

    // Implicit BPs /////////////////////////
    // Set the hash values for each leaf node.
    unsigned long long temp1 = vvec_hashrf._HF.getA1(idx);
    unsigned long long temp2 = vvec_hashrf._HF.getA2(idx);
    startNode->hv1 = temp1;
    startNode->hv2 = temp2;

    temp1 = temp1 % m1;
    temp2 = temp2 % m2;
    int numOnes = 0;
    for (unsigned int i  =0; i < NUM_TAXA; i++)
      numOnes += bs[i];
    if (numOnes > 1 || numOnes == 0) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring creation error\n";
      exit(0);
    }

    return bs;
  }
  else { //if we are not at a leaf node
     //perform DFSn all the children
    bool * ebs[startNode->Nchildren];
    ///fprintf(stderr, "I am an internal node!\n");
    for (int x = 0; x < startNode->Nchildren; ++x) { 
      ebs[x] = dfs_compute_hash_unrooted_mc(startNode->child[x],lm, vvec_hashrf, treeIdx, numBitstr, m1, m2, processed, startNode);
    }
            
    // At this point, we find a bipartition. 
    // Thus, OR the bitstrings and make a bit string for the bipartition
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    ///fprintf(stderr, "the bitstrings I'm going to be comparing are:\n");
    //for (int i=0; i<startNode->Nchildren; ++i) {
     // fprintf(stderr, "%s\n", (startNode->child[i])->label);
     // }
      //fprintf(stderr, "\n");

    for (int i=0; i<startNode->Nchildren; ++i) {
      //cerr << "bitstring: " << *(ebs[i]) << endl;
 	if (ebs[i]) {
	  for (unsigned int j =0; j < NUM_TAXA; j++)
	    bs[j] |= ebs[i][j];
	  delete ebs[i];
	  ebs[i] = NULL;
	}
	else {
	  cout << "ERROR: null bitstring\n";
	  exit(0);
	}
    }
    
    //fprintf(stderr, "new bitstring is:\n");
    //for (unsigned int j = BITSETSZ-1; j > BITSETSZ-1-NUM_TAXA; j--){
     // cerr << (*bs)[j];
      //}
     // fprintf(stderr, "\n");

    int numOnes = 0;
    for (unsigned int i  =0; i < NUM_TAXA; i++)
      numOnes += bs[i];
    if (numOnes < 1) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring OR error\n";
      exit(0);
    }

    unsigned long long temp1 = 0;
    unsigned long long temp2 = 0;
    
    for (int i = 0; i < startNode->Nchildren; ++i) { 
      temp1 += startNode->child[i]->hv1;
      temp2 += startNode->child[i]->hv2;
    }
    startNode->hv1 = temp1 % m1;
    startNode->hv2 = temp2 % m2;

    float dist = 0.0;
    if (WEIGHTED)
      dist = startNode->weight;
    else
      dist = 1;
  
    bool toprocess = true;
    if ( (parent == NULL) || ( (parent->myparent == NULL) && (parent->Nchildren == 2) && (processed) ) ){
      ///fprintf(stderr, "we don't collect this bitstring!\n");
      toprocess = false;
    }
    else{
      if ( (parent->myparent == NULL) && (parent->Nchildren == 2) && (!processed)){ 
	  if (bs[0] != 0){ //don't process if it's oriented wrong
	    toprocess = false;
	    ///fprintf(stderr, "we don't collect this bitstring.. it's oriented wrong!\n");
	  }
	  else{
	    processed = true;
	  }
      }
      else{
	if ( bs[0] != 0 ){
	  ///fprintf(stderr, "need to flip bitstring!\n");
	  toprocess = false;

	  //calculate new h1 h2 values	    
	  unsigned long long temp1 = 0;
	  unsigned long long temp2 = 0;
	  unsigned int idx = 0;
	  for (unsigned int j = 0; j < NUM_TAXA; j++){
	    if ( !bs[j] ) {
	      idx = j;
	      temp1 += vvec_hashrf._HF.getA1(idx);
	      temp2 += vvec_hashrf._HF.getA2(idx);
	    }
	  }
	  temp1 = temp1 % m1;
	  temp2 = temp2 % m2;
	  if (numBitstr < NUM_TAXA-3){
	    ///fprintf(stderr, "adding flipped bitstring to hash table\n");
	    vvec_hashrf.hashing_bs_mc(treeIdx, NUM_TAXA, startNode->hv1, startNode->hv2, dist, WEIGHTED);   // without TYPE-III using n-bits (Hash-RF)
	    ++numBitstr;
	  }
	}
      }
    }    
    if (toprocess){
      if (numBitstr < NUM_TAXA-3){
	///fprintf(stderr, "adding bitstring to hash table!\n");
	vvec_hashrf.hashing_bs_mc(treeIdx, NUM_TAXA, startNode->hv1, startNode->hv2, dist, WEIGHTED);   //no flipping needed
	++numBitstr;
      }
    }
    return bs;
  }

}

//---convert int to string---
string itos(int i)	
{
  stringstream s;
  s << i;
  return s.str();
}


//--------- auxilliary function(1) for bitstring encoding----------
int toDecimal(bool * temp){ 
  int sum = 0;
  int j = 1;
  for (short i = 19; i >=0; i--) {
    if (temp[i]) 
      sum += j;
     j*=2;
  }

  return sum;
}

//----------- auxilliary function(2) for bitstring encoding---------
void toBinary(bool *output, int input) {
  int j = 524288;
  int i = 0;
  //short temp = 0;

  while (i!=20){
    if ( (input - j) >= 0){
      output[i] = 1;
      input -= j;
    }
    else { 
      output[i] = 0;
    }
    i++;
    j/=2;
  }
}


void encode(unsigned long long x) {
   while (x > 127) {
     putchar((x & 127)|128);
     x >>= 7;
   }
   putchar(x);
}

void encode2(unsigned int x) {
  bool contin;
  unsigned int y = 0;

  while (x > 0){
    cout << "x is: " << x << endl;
    y = x;
    contin = false;
    if (y > 127)
      contin = true;
    while (y > 127)
      y = y/128;
    cout << "y is: " << y << endl;
    y <<=1;
    cout << "y is: " << y << endl;
    if (contin)
      y = y | 1;
    cout << "y is:" << y << endl;
    //printf("%x", y);
    x = x / 128;
    x = x << y;
    cout << "new x is:" << x << endl;
    x = x >> y;
  }
  //printf("\n");
}

void encode3(unsigned int &x, bool contin, FILE * fp){
  unsigned int y,z;
  if (x > 127){
    y =  x >> 7;
    encode3(y, 1, fp);
    y = x - (y << 7);
    z = y << 1;
    z = z | contin;
    //cout << z << endl;
    fprintf(fp, "%c", z);
  }
  else{
    z =  x << 1;
    z = z | contin;
    //cout << z << endl;
    fprintf(fp, "%c", z);
  }
}
	  
//------------ collect labels from a set of trees ---------
LabelMap collect_labels(string file1) {
  unsigned long ulLineCount, initial_ntaxa;
  ulLineCount = 0;  
  initial_ntaxa = 0;
  NEWICKTREE *newickTree = NULL;
  int err;
  FILE *fp;
  LabelMap lm, lm2;
  if (NUM_TREES == 0){ //first pass: get number of trees (needed for next op)
    ifstream fin;
    string line;
    fin.open(file1.c_str());
    if (!fin){
      cerr << "Cannot open file for reading!" << endl;
      exit(1);
    }
    while (!fin.eof()){
      getline(fin, line);
      if (fin.eof())
	break;
      ulLineCount++;
    }
    NUM_TREES = ulLineCount;
    fin.close();
  } 
  fp = fopen(file1.c_str(), "r"); //second pass: collect all labels in file
  if (!fp) { cout << "Fatal error: file open error\n"; exit(1); }
  for (unsigned int x = 0; x < NUM_TREES; x++){ 
    //cout << "x is: " << x << endl;
    newickTree = loadnewicktree2(fp, &err);
    if (!newickTree) {
      switch (err) {
      case -1:
	printf("Out of memory\n");
	break;
      case -2:
	printf("parse error on line %u\n", x);
	break;
      case -3:
	printf("Can't load file\n");
	break;
      default:
	printf("Error %d\n", err);
      }
      exit(0);
    }
    /*****************************************************/
    //cout << "\n*** Collecting the taxon labels ***\n"; 
    /*****************************************************/
    PopulateTaxaLabels(newickTree->root, lm);
    killnewicktree(newickTree);
  }
  NUM_TAXA = lm.size();
  fprintf(stderr, "Found: %d taxa\n", (int)lm.size());
  fclose(fp);
  fp = fopen(file1.c_str(), "r"); //third pass: check for heterogeneous
  if (!fp) { cout << "Fatal error: file open error\n"; exit(1); }
  for (unsigned int x = 0; x < NUM_TREES; x++){ 
    newickTree = loadnewicktree2(fp, &err);
    if (!newickTree) {
      switch (err) {
      case -1:
	printf("Out of memory\n");
	break;
      case -2:
	printf("parse error on line %u\n", x);
	break;
      case -3:
	printf("Can't load file\n");
	break;
      default:
	printf("Error %d\n", err);
      }
      exit(0);
    }
    /*****************************************************/
    //cout << "\n*** Collecting the taxon labels ***\n"; 
    /*****************************************************/
    lm2.clear(); //get labels for each tree
    PopulateTaxaLabels(newickTree->root, lm2);
    //if the size of any tree does not match the total leaf set, we have
    //a heterogeneous collection of trees.
    if ( (lm2.size() != NUM_TAXA) && (HETERO!= true)){ 
      fprintf(stderr, "Detected heterogeneous collection of trees!\n");
      HETERO = true;
    }
    killnewicktree(newickTree);
  }
  fclose(fp);
  return lm;
}

LabelMap collect_labels_light(string file1) {
  unsigned long ulLineCount, initial_ntaxa;
  ulLineCount = 0;  
  initial_ntaxa = 0;
  NEWICKTREE *newickTree = NULL;
  int err;
  if (NUM_TREES == 0){
    ifstream fin;
    string line;
    fin.open(file1.c_str());
    if (!fin){
      cerr << "Cannot open file for reading!" << endl;
      exit(1);
    }
    while (!fin.eof()){
      getline(fin, line);
      if (fin.eof())
	break;
      ulLineCount++;
    }
    NUM_TREES = ulLineCount;
    fin.close();
  } //the reason we do this separately is because if we do it together, there is a seg-fault!
  FILE *fp;
  LabelMap lm;
  fp = fopen(file1.c_str(), "r");
  //for the quick version, we just get the labels from the first tree
  newickTree = loadnewicktree2(fp, &err);
  if (!newickTree) {
    switch (err) {
    case -1:
      printf("Out of memory\n");
      break;
    case -2:
      printf("parse error\n");
      break;
    case -3:
      printf("Can't load file\n");
      break;
    default:
      printf("Error %d\n", err);
    }
    exit(0);
  }
  /*****************************************************/
  //cout << "\n*** Collecting the taxon labels ***\n"; 
  /*****************************************************/
  GetTaxaLabels(newickTree->root, lm);
  NUM_TAXA = lm.size();
  killnewicktree(newickTree);
  fprintf(stderr, "Found: %d taxa\n", (int)lm.size());
  fclose(fp);
  return lm;
}

