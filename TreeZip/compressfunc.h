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
#include <vector>

using namespace std;

#ifndef _COMPRESSFUNC_H_
#define _COMPRESSFUNC_H_
 
void compress(string file1, bool bUbid, bool quiet, vector<string> nextaxa, bool use);
void decompress(string file, int random_decompress, bool branch, bool calc_unique_trees, unsigned char consensus, bool quiet);
void sample_hashtable(string trzfile);
void convert_nexus(string in, string out);
bool is_nexus(string file);
void parse_nexus(string filename);
void decompress_nexus(string file, int random_decompress, bool branch, bool calc_unique_trees, unsigned char consensus, bool quiet);
void clear_hashtable();
#endif

