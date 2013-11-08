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


#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits>
#include <bitset>
#include <map>
#include <utility>

#include "label-map.hh"

using namespace std;

#ifndef _BUILDTREE_H_
#define _BUILDTREE_H_

/*
string compute_tree( 
    //multimap<unsigned, unsigned, greater<unsigned> > mmap_cluster, 
    vector<pair<unsigned, unsigned> > mmap_cluster,
    vector<vector<SCNode*> > vvec_distinctClusters2,
    LabelMap lm,
    vector< bitset<BITSETSZ>* > my_bs,
    unsigned id);
 */

string compute_tree( 
    LabelMap lm,
    vector< bool* > my_bs,
    vector<float> my_branches,
    unsigned id,
    bool branch,
    vector<unsigned int> bs_sizes);

string compute_tree( 
    LabelMap lm,
    vector< bool* > my_bs,
    vector<float> my_branches,
    unsigned id,
    bool branch);
#endif 
