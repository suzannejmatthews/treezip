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

*/
/*****************************************************/

#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <valarray>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits>
#include <bitset>

#include "label-map.hh"
#include "hashfunc.hh"
#include "hash.hh"
#include "SCTree.h"
#include "parsing.h"
#include "global.h"

// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;

/*bool mmap_sorter(const pair<unsigned unsigned> & a, const pair<unsigned, unsigned> &b) { 
  return a.first > b.first;
  }*/

bool vvec_sorter(const vector<SCNode *> & a, const vector<SCNode*> &b) { 
  return a.size() > b.size();
}

//in here for ONLY HashCS! Do not use this function for any other purpose!!!
//please use the other compute_tree function for general purpose things.
string compute_tree(
    LabelMap lm,
    vector< bool * > my_bs,
    vector< float > my_branches,
    unsigned id,
    bool branch) {

  vector<vector<SCNode*> > vvec_distinctClusters2;

  //specify timers
  struct timeval update_start;
  struct timeval update_end;
  struct timeval build_start;
  struct timeval build_end;
  struct timeval clean_start;
  struct timeval clean_end;
  struct timeval total_start;  
  struct timeval total_end;
  gettimeofday(&total_start, NULL);
  gettimeofday(&update_start, NULL);

  //update distinct clusters
  for (unsigned int i = 0; i < my_bs.size(); ++i) {
    vector<SCNode*> vec_nodes2;
    unsigned int lmIndex = 0;
    for (unsigned int j = 0; j < NUM_TAXA; j++) {
      if (my_bs[i][j]) {
	SCNode* aNode = new SCNode();
	aNode->name = lm.name(lmIndex);
	vec_nodes2.push_back(aNode);
      }
      lmIndex++;
    }
    vvec_distinctClusters2.push_back(vec_nodes2);
  }
    
  SCTree *scTree = new SCTree();
  bool addedRoot = false;
  unsigned int intNodeNum = 0;
  vector<SCNode*> vec_garbageCan; // for collecting pointers
 
  gettimeofday(&update_end, NULL);
  update_time += update_end.tv_sec - update_start.tv_sec + (update_end.tv_usec - update_start.tv_usec) / 1.e6;
  
  gettimeofday(&build_start, NULL);
  for (unsigned int pos = 0; pos < vvec_distinctClusters2.size(); ++pos) {
    if (!addedRoot) {
      // The first cluster has all the taxa.
      // This constructs a star tree with all the taxa.
      // 1. Dangle all the taxa as root's children by adjusting parent link.
      // 2. Push all the node* in the root's children
      // 3. Push all the nodes in the tree's nodelist.
      // 4. Push all the nodes' parent in the tree's parentlist.
      
      for (unsigned int i=0; i < vvec_distinctClusters2[pos].size(); ++i) {
	vvec_distinctClusters2[pos][i]->parent = scTree->root;
	scTree->root->children.push_back(vvec_distinctClusters2[pos][i]);
	scTree->nodelist.push_back(vvec_distinctClusters2[pos][i]);
	assert(scTree->nodelist[0]->name == "root");
	scTree->parentlist2.insert(map<string,int>::value_type(vvec_distinctClusters2[pos][i]->name, 0));
      }
      addedRoot = true;
    }
    else {
      // For the next biggest cluster,
      // 1. Find node list to move (= vvec_distinctClusters2[itr->second]) and
      //    Get the parent node of the to-go nodes.
      // 2. Make an internal node.
      // 3. Insert the node in the nodelist of the tree and update the parentlist accordingly
      // 4. Adjust the to-go nodes' parent link.
      // 5. Adjust the parent node's link to children (delete the moved nodes from children).
      
      // 1. --------------------------------------------------------------------------
      SCNode* theParent = NULL;
      theParent = scTree->nodelist[scTree->parentlist2[vvec_distinctClusters2[pos][0]->name]];
      assert(theParent != NULL);
      assert(theParent->name != "");
      
      // 2. --------------------------------------------------------------------------
      //			string newIntNodeName = "int" + itostr(intNodeNum, 10);
      string newIntNodeName = "int" + itos(intNodeNum);
      SCNode* newIntNode = new SCNode();
      vec_garbageCan.push_back(newIntNode);
      newIntNode->name = newIntNodeName;
      newIntNode->parent = theParent;
      if (WEIGHTED)
	newIntNode->bl = my_branches[pos];
      
      // 3. --------------------------------------------------------------------------
      assert(newIntNodeName.size() != 0);
      scTree->nodelist.push_back(newIntNode);
      assert(scTree->nodelist[scTree->nodelist.size()-1]->name == newIntNode->name);
      
      scTree->parentlist2.insert(map<string, unsigned>::value_type(newIntNodeName, scTree->nodelist.size()-1));
      
      for (unsigned int i=0; i<vvec_distinctClusters2[pos].size(); ++i) {
	// 4. --------------------------------------------------------------------------
	vvec_distinctClusters2[pos][i]->parent = newIntNode;
	
	// We have to update parentlist in the tree.
	assert(vvec_distinctClusters2[pos][i]->parent->name == scTree->nodelist[scTree->nodelist.size()-1]->name);
	
	scTree->parentlist2[vvec_distinctClusters2[pos][i]->name] = scTree->nodelist.size()-1;
	newIntNode->children.push_back(vvec_distinctClusters2[pos][i]);
	
	// 5. --------------------------------------------------------------------------
	// Delete the moved nodes from parent's children.
	vector<SCNode*>::iterator itr2;
	
	for (itr2 = theParent->children.begin(); itr2 != theParent->children.end(); ++itr2) {
	  if (vvec_distinctClusters2[pos][i]->name == (*itr2)->name) {
	    theParent->children.erase(itr2);
	    break;
	  }
	}
      }
      theParent->children.push_back(newIntNode);
      intNodeNum++;
    }
  }
  
  string mytree;

  if (branch && WEIGHTED)
    mytree = scTree->GetTreeString(true);
  else
    mytree = scTree->GetTreeString(false);
  //scTree->DrawOnTerminal(true);
 
  gettimeofday(&build_end, NULL);
  build_time += build_end.tv_sec - build_start.tv_sec + (build_end.tv_usec - build_start.tv_usec) / 1.e6;
  gettimeofday(&clean_start, NULL);
  
  for (unsigned int i = 0; i < vvec_distinctClusters2.size(); ++i) {
    for (unsigned int j=0; j<vvec_distinctClusters2[i].size(); ++j) {
      SCNode * tmp = vvec_distinctClusters2[i][j];
      delete tmp;
      vvec_distinctClusters2[i][j] = NULL;
    }
  }

  for (unsigned i=0; i<vec_garbageCan.size(); ++i) {
    if (vec_garbageCan[i]) {
      SCNode* temp = vec_garbageCan[i];
      delete temp;
      vec_garbageCan[i] = NULL;
    }
  }
  
  if (scTree->root) {
    SCNode * tmp = scTree->root;
    delete tmp;
    scTree->root = NULL;
  }
  scTree->nodelist.clear();
  scTree->parentlist2.clear();
  delete scTree;
  gettimeofday(&clean_end, NULL);
  clean_time += clean_end.tv_sec - clean_start.tv_sec + (clean_end.tv_usec - clean_start.tv_usec) / 1.e6;
 gettimeofday(&total_end, NULL);
 total_time += total_end.tv_sec - total_start.tv_sec + (total_end.tv_usec - total_start.tv_usec) / 1.e6;

  return mytree;
  }


string compute_tree(
    LabelMap lm,
    vector< bool * > my_bs,
    vector< float > my_branches,
    unsigned id,
    bool branch,
    vector<unsigned int> bs_sizes) {

  vector<vector<SCNode*> > vvec_distinctClusters2;

  //specify timers
  struct timeval update_start;
  struct timeval update_end;
  struct timeval build_start;
  struct timeval build_end;
  struct timeval clean_start;
  struct timeval clean_end;
  struct timeval total_start;  
  struct timeval total_end;
  gettimeofday(&total_start, NULL);
  gettimeofday(&update_start, NULL);

  //update distinct clusters
  for (unsigned int i = 0; i < my_bs.size(); ++i) {
    vector<SCNode*> vec_nodes2;
    unsigned int lmIndex = 0;
    unsigned int bs_size = bs_sizes[i];
    //cout << "bs_size is: " << bs_size << endl;
    for (unsigned int j = 0; j < bs_size; j++) {
      if (my_bs[i][j]) {
	SCNode* aNode = new SCNode();
	aNode->name = lm.name(lmIndex);
	vec_nodes2.push_back(aNode);
      }
      lmIndex++;
    }
    vvec_distinctClusters2.push_back(vec_nodes2);
  }
  

  /* for (unsigned int npos = 0; npos < vvec_distinctClusters2.size(); ++npos) { 
    for (unsigned int x = 0; x < vvec_distinctClusters2[npos].size(); ++x) { 
      cout << vvec_distinctClusters2[npos][x]->name << " ";
    }
    cout << endl;
  }
  exit(0);
  */
  //if (vvec_distinctClusters2.size() != (NUM_TAXA-2)){
  //  fprintf(stderr, "ERROR!: Tree %d has only %d bipartitions!!\n Exiting...\n", id, vvec_distinctClusters2.size());
  //  exit(0);
  //  }
  
  
  SCTree *scTree = new SCTree();
  bool addedRoot = false;
  unsigned int intNodeNum = 0;
  vector<SCNode*> vec_garbageCan; // for collecting pointers
  
  gettimeofday(&update_end, NULL);
  update_time += update_end.tv_sec - update_start.tv_sec + (update_end.tv_usec - update_start.tv_usec) / 1.e6;
  
  

  gettimeofday(&build_start, NULL);
  for (unsigned int pos = 0; pos < vvec_distinctClusters2.size(); ++pos) {
    if (!addedRoot) {
      // The first cluster has all the taxa.
      // This constructs a star tree with all the taxa.
      // 1. Dangle all the taxa as root's children by adjusting parent link.
      // 2. Push all the node* in the root's children
      // 3. Push all the nodes in the tree's nodelist.
      // 4. Push all the nodes' parent in the tree's parentlist.
      
      for (unsigned int i=0; i < vvec_distinctClusters2[pos].size(); ++i) {
	vvec_distinctClusters2[pos][i]->parent = scTree->root;
	scTree->root->children.push_back(vvec_distinctClusters2[pos][i]);
	scTree->nodelist.push_back(vvec_distinctClusters2[pos][i]);
	assert(scTree->nodelist[0]->name == "root");
	scTree->parentlist2.insert(map<string,int>::value_type(vvec_distinctClusters2[pos][i]->name, 0));
      }
      addedRoot = true;
    }
    else {
      // For the next biggest cluster,
      // 1. Find node list to move (= vvec_distinctClusters2[itr->second]) and
      //    Get the parent node of the to-go nodes.
      // 2. Make an internal node.
      // 3. Insert the node in the nodelist of the tree and update the parentlist accordingly
      // 4. Adjust the to-go nodes' parent link.
      // 5. Adjust the parent node's link to children (delete the moved nodes from children).
      
      // 1. --------------------------------------------------------------------------
      SCNode* theParent = NULL;
      theParent = scTree->nodelist[scTree->parentlist2[vvec_distinctClusters2[pos][0]->name]];
      assert(theParent != NULL);
      assert(theParent->name != "");
      
      // 2. --------------------------------------------------------------------------
      //			string newIntNodeName = "int" + itostr(intNodeNum, 10);
      string newIntNodeName = "int" + itos(intNodeNum);
      SCNode* newIntNode = new SCNode();
      vec_garbageCan.push_back(newIntNode);
      newIntNode->name = newIntNodeName;
      newIntNode->parent = theParent;
      if (WEIGHTED)
	newIntNode->bl = my_branches[pos];
      
      // 3. --------------------------------------------------------------------------
      assert(newIntNodeName.size() != 0);
      scTree->nodelist.push_back(newIntNode);
      assert(scTree->nodelist[scTree->nodelist.size()-1]->name == newIntNode->name);
      
      scTree->parentlist2.insert(map<string, unsigned>::value_type(newIntNodeName, scTree->nodelist.size()-1));
      
      for (unsigned int i=0; i<vvec_distinctClusters2[pos].size(); ++i) {
	// 4. --------------------------------------------------------------------------
	vvec_distinctClusters2[pos][i]->parent = newIntNode;
	
	// We have to update parentlist in the tree.
	assert(vvec_distinctClusters2[pos][i]->parent->name == scTree->nodelist[scTree->nodelist.size()-1]->name);
	
	scTree->parentlist2[vvec_distinctClusters2[pos][i]->name] = scTree->nodelist.size()-1;
	newIntNode->children.push_back(vvec_distinctClusters2[pos][i]);
	
	// 5. --------------------------------------------------------------------------
	// Delete the moved nodes from parent's children.
	vector<SCNode*>::iterator itr2;
	
	for (itr2 = theParent->children.begin(); itr2 != theParent->children.end(); ++itr2) {
	  if (vvec_distinctClusters2[pos][i]->name == (*itr2)->name) {
	    theParent->children.erase(itr2);
	    break;
	  }
	}
      }
      theParent->children.push_back(newIntNode);
      intNodeNum++;
    }
  }
  
  string mytree;

  if (branch && WEIGHTED)
    mytree = scTree->GetTreeString(true);
  else
    mytree = scTree->GetTreeString(false);
  //scTree->DrawOnTerminal(true);
 
  gettimeofday(&build_end, NULL);
  build_time += build_end.tv_sec - build_start.tv_sec + (build_end.tv_usec - build_start.tv_usec) / 1.e6;
  gettimeofday(&clean_start, NULL);
  
  for (unsigned int i = 0; i < vvec_distinctClusters2.size(); ++i) {
    for (unsigned int j=0; j<vvec_distinctClusters2[i].size(); ++j) {
      SCNode * tmp = vvec_distinctClusters2[i][j];
      delete tmp;
      vvec_distinctClusters2[i][j] = NULL;
    }
  }

  for (unsigned i=0; i<vec_garbageCan.size(); ++i) {
    if (vec_garbageCan[i]) {
      SCNode* temp = vec_garbageCan[i];
      delete temp;
      vec_garbageCan[i] = NULL;
    }
  }
  
  if (scTree->root) {
    SCNode * tmp = scTree->root;
    delete tmp;
    scTree->root = NULL;
  }
  scTree->nodelist.clear();
  scTree->parentlist2.clear();
  delete scTree;
  gettimeofday(&clean_end, NULL);
  clean_time += clean_end.tv_sec - clean_start.tv_sec + (clean_end.tv_usec - clean_start.tv_usec) / 1.e6;
 gettimeofday(&total_end, NULL);
 total_time += total_end.tv_sec - total_start.tv_sec + (total_end.tv_usec - total_start.tv_usec) / 1.e6;

  return mytree;
		
}
