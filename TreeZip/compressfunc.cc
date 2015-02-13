//*****************************************************************/
/*
This is TreeZip, a compresson software for phylogenetic trees. 
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
#include <iomanip>
#include <sstream>
#include <unistd.h>
#include "RandomLib/Random.hpp"

#include "label-map.hh"
#include "hashfunc.hh"
#include "hash.hh"
#include "SCTree.h"
#include "parsing.h"
#include "global.h"
#include "buildtree.h"
#include "bmerge.h"
#include "tmerge.h"
#include "compressfunc.h"
#include "nexus_compress.hh"
#include "nexus_decompress.hh"
#include "normalizer.h"
#include <utility>

// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;

typedef struct uniq { 
  unsigned int tid;
  int pos;
} myuniq;

bool sort_bitstrings(const pair<unsigned, pair<string, string> > & a, const pair<unsigned, pair< string, string > > &b) { 
  if (a.first == b.first)
    return ( (a.second).first > (b.second).first);
  else
    return a.first > b.first;
}

void count_bipartitions(unsigned int & bipart_counter, unsigned long & tree_counter, vector< pair<unsigned, pair<string, string > > >  & ordered_bitstrings, LabelMap & lm, LabelMap & lm_sorted){
  stringstream ss;
  unsigned int my_tid;
  string bitstring;
  unsigned int loctax = 0;
  unsigned int maxLength;
  string mytax;
  bool * sortedbs = new bool[NUM_TAXA];
  for (unsigned int ls = 0; ls < NUM_TAXA; ls++)
    sortedbs[ls] = 0; //initialize everything to zero
  for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size();    
    if (sizeVec) {
      bipart_counter += sizeVec;
      for (unsigned int i=0; i<sizeVec; ++i) {	  
	unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
	tree_counter += sizeTreeIdx;
	for (unsigned int j = 0; j < sizeTreeIdx; j++){
	  my_tid = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	  //cout << my_tid << " ";
	  helper_sizes[my_tid]++;
	}
	//cout << " ";
	bool * temp = vvec_hashrf._hashtab2[hti][i]._bs;
	int num_ones = 0;
	ss.str("");
	unsigned int bs_size = vvec_hashrf._hashtab2[hti][i]._bs_size;
	maxLength = 0;
	for (unsigned int ls = 0; ls < bs_size; ls++){ //change order of bits
	  if (temp[ls]){
	    mytax = lm.name(ls);
	    loctax = lm_sorted[mytax];
	    sortedbs[loctax] = 1;
	    if (loctax > maxLength)
	      maxLength = loctax;
	    num_ones++;
	  }
	}
	maxLength++;
	string rle_bipart = rle(sortedbs, maxLength); //run-length encode sorted bipartition
	for (unsigned int ls = 0; ls < maxLength; ls++) //reset sortedbs  
	  sortedbs[ls] = 0;
	ss.str("");
	string loc;
	ss << hti << "." << i;
	loc = ss.str();
	//ss.str("");
	//ss << i;
	//loc = loc + "." + ss.str();
	ordered_bitstrings.push_back(make_pair(num_ones, make_pair(rle_bipart, loc)));
	//loc = "";
	ss.str("");   
      } //for everything in sizevec
    }
  }
}

string encode_branch_lengths(unsigned int x, unsigned int y, unsigned int sizeTreeIdx, unsigned int encode_size){
  string len, tmplen;
  unsigned int size = 0, strpad, tmpad;
  string prefixes, myprefix;
  stringstream ss, stm;
  short digit;
  ss << " ";
  //first see if there are branch lengths greater than one
  prefixes="";
  for (unsigned int a = 0; a < sizeTreeIdx; a++){
    //printf("%f\n", vvec_hashrf._hashtab2[x][y]._vec_dist[a]);
    if (vvec_hashrf._hashtab2[x][y]._vec_dist[a] >= 1){
      digit = (int)vvec_hashrf._hashtab2[x][y]._vec_dist[a];
      stm << vvec_hashrf._hashtab2[x][y]._vec_treeidx[a];
      stm >> myprefix;
      stm.clear();
      if (myprefix.size() < (encode_size*2)){
	strpad = (encode_size*2)-myprefix.size();
	while (strpad > 0){
	  myprefix = "0" + myprefix;
	  strpad--;
	}
      }     
      if (encode_size > 1){
	strpad = 0;
	while ((strpad < encode_size)){
	  tmplen = myprefix.substr(0,2);
	  tmpad = atoi(tmplen.c_str());
	  tmpad+=33;
	  //fprintf(fout, "%c", tmpad);
	  //printf("%c", tmpad);
	  ss << char(tmpad);
	  strpad++;
	  myprefix = myprefix.substr(2);
	}
      }
      else{
	tmpad = atoi(myprefix.c_str());
	tmpad+=33;
	//fprintf(fout, "%c", tmpad);
	//printf("%c", tmpad);
	ss << char(tmpad);
      }
      
      stm << digit;
      stm >> myprefix;
      stm.clear();
      if (digit > 9){	
	//fprintf(fout, "+%s+", myprefix.c_str());
	//printf("+%s+", myprefix.c_str());
	ss << "+" << myprefix << "+";
      }
      else
	//fprintf(fout, "%s", myprefix.c_str());
	//printf("%s", myprefix.c_str());
	ss << myprefix;
    }
  }
  //fprintf(fout, " ");
  //printf(" ");
  ss << " ";
  //now process rest of branch lengths, in order
  for (unsigned int a = 0; a < sizeTreeIdx; a++){
    stm << std::fixed << setprecision (6) << vvec_hashrf._hashtab2[x][y]._vec_dist[a];
    stm >> len;
    int pos = len.find_first_of("."); 
    len = len.substr(pos+1); //get rid of the prefix    
    //cout << len << " "; 
    size = len.size();
    unsigned int myi = 0;
    while (myi < 3){
      if (size>1){
	tmplen = len.substr(0,2);
	digit = atoi(tmplen.c_str());
	digit+=33;
	//fprintf(fout, "%c", digit);
	//printf("%c", digit);
	ss << char(digit);
	size-=2;
	len = len.substr(2);
      }
      else if(size==1){
	tmplen = len.substr(0,1);
	digit = atoi(tmplen.c_str());
	digit*=10;
	digit+=33;
	//fprintf(fout, "%c", digit);
	//printf("%c", digit);
	ss << char(digit);
	size-=1;
	len = len.substr(1);
      }
      else if (size < 1){
	digit = 33;
	//fprintf(fout, "%c", digit);
	//printf("%c", digit);
	ss << char(digit);
      }  
      myi++;
    }
    //printf("\n");    
    //cout << size << endl;
    stm.clear();   
  }
  //exit(0);
  return ss.str();
}

bool check_for_weights(string file1){
  ifstream fin;
  string tree;
  fin.open(file1.c_str());
  if (!fin){
    cerr << "cannot open file for reading!" << endl;
    exit(1);
  }
  getline(fin, tree);
  fin.close();
  int pos = tree.find_first_of(':');
  if (pos == -1)
    return false;
  return true;
}

bool is_nexus(string file1){
  ifstream fin;
  string nex;
  fin.open(file1.c_str());
  if (!fin){
    cerr << "cannot open file " << file1 << endl;
    exit(2);
  }
  getline(fin, nex);
  fin.close();
  size_t pos = nex.find("#NEXUS");
  if (pos != string::npos)
    return true;
  return false;
}

void convert_nexus(string in, string out){
  //convert the nexus file called 'in' to a trz file stored in 'out'
  ifstream fin;
  string line;
  int pos;
  ofstream fout;
  fin.open(in.c_str());
  if (!fin){
    cerr << "cannot open " << in << " for reading!" << endl;
    exit(1);
  }
  fout.open(out.c_str());
  if (!fout){
    cerr << "cannot open " << fout << " for writing!" << endl;
    exit(1);
  }
  bool print = false;
  bool first = false;
  while (!fin.eof()){
    getline(fin, line);
    if (fin.eof())
      break;
    if (line == "end;")
      break;
    if (!first){
      pos = line.find_first_of("=");
      if (pos != -1){
	print = true;
	first = true;
      }
      if (first)
	continue;
    }
    if (print)
      fout << line << endl;
  }
  fin.close();
  fout.close();
}

//this is the main function for parsing nexus files
//once this function is done, nexus_compress and nexus_decompress can be removed.
void print_out_temporary_file(string file, ofstream & fout){
  ifstream myfin;
  string myline;
  myfin.open(file.c_str());
  if (!myfin){
    cerr << "cannot find temporary file '" << file << "'!" << endl;
  }
  while (!myfin.eof()){
    getline(myfin, myline);
    if (myfin.eof())
      break;
    fout << myline << endl;
  }
  unlink(file.c_str());
  myfin.close();
}

int preprocessNEXUS(string file, string outfile){
  //preprocessing does some stuff to file to help with parsing
  //right now, all it does is remove comments
  ifstream fin;
  fin.open(file.c_str());
  if (!fin){
    cerr << "cannot open file for reading!" << endl;
    exit(1);
  }
  ofstream fout;
  fout.open(outfile.c_str());
  if (!fout){
    cerr << "cannot open outfile " << outfile << " for writing!" << endl;
    exit(1);
  }
  string line;
  bool incomment = false;
  while (!fin.eof()){
    getline(fin, line);
    if (fin.eof())
      break;

    int pos = line.find_first_of("[");
    if (pos != -1){
      if (pos > 0)
	line = line.substr(0, pos);
      incomment = true;
    }
    pos = line.find_first_of("]");
    if (pos != -1){
      assert(incomment==true);
      line = line.substr(pos+1);
      incomment=false;
    }
    if (!incomment){
      if (line != "")
	fout << line << endl;
    }
  }
  return 0;
}

void parse_nexus(string file){
  int success =-1;
  //first, normalize file for output -- uncomment this section for normalizing (NCL)
  //success = readFilepathAsNEXUS(file.c_str(), ".mytempnexus"); //NCL support
  success = preprocessNEXUS(file, ".mytempnexus");
  if (success != 0){
    cerr << "Error: Malformed NEXUS file!" << endl;
    unlink(".mytempnexus");
    exit(1);
    }
  ifstream fin;
  fin.open(".mytempnexus");
  if (!fin){
    cerr << "Cannot open temporary file .mytempnexus for reading!" << endl;
    exit(1);
  }
  string line, checkline, type;
  int pos = 0; 
  int pos2;
  bool etaxa, echar, eunalign, edata, etree, eset, eassume, edist, enote, epaup, emb, ecodon;
  short ntaxa, nchar, nunalign, ndata, ntree, nset, nassume, ndist, nnote, npaup, nmb, ncodon;
  etaxa = false;
  echar = false;
  eunalign = false;
  edata = false;
  etree = false;
  eset = false;
  eassume = false;
  edist = false;
  enote = false;
  epaup = false;
  emb = false;
  ecodon = false; 
  ntaxa = 0;
  nchar = 0;
  nunalign = 0;
  ndata = 0;
  ntree = 0;
  nset = 0;
  nassume = 0;
  ndist = 0;
  nnote = 0;
  npaup = 0;
  nmb = 0;
  ncodon = 0;
  while (!fin.eof()){
    getline(fin, line);
    if (fin.eof())
      break;
    checkline = line;
    transform(checkline.begin(), checkline.end(), checkline.begin(), ::toupper);
    pos = checkline.find("BEGIN");
    if (pos != -1){
      pos2 = line.find_first_of(" ");
      line = line.substr(pos2+1);
      pos2 = line.find_first_of(';');
      if (pos2 == -1){
	cerr << "malformed nexus file!" << endl;
	exit(1);
      }
      type = line.substr(0, pos2);
      transform(type.begin(), type.end(), type.begin(), ::tolower);
      if (type == "taxa"){
	taxa_compress(fin, etaxa, ntaxa);
      }
      else if (type == "trees"){
	trees_compress(fin, etree, ntree);
      }
      else if (type == "characters"){
	char_compress(fin, echar, nchar);
      }
      else if (type == "unaligned"){
	unalign_compress(fin, eunalign, nunalign);
      }
      else if (type == "data"){
	data_compress(fin, edata, ndata);
      }
      else if (type == "distances"){
	dist_compress(fin, edist, ndist);
      }
      else if (type == "codons"){
	codon_compress(fin, ecodon, ncodon);
      }
      else if (type == "sets"){
	set_compress(fin, eset, nset);
      }
      else if (type == "assumptions"){
	assume_compress(fin, eassume, nassume);
      }
      else if (type == "notes"){
	note_compress(fin, enote, nnote);
      }
      else if (type == "paup"){
	paup_compress(fin, epaup, npaup);
      }
      else if (type == "mrbayes"){
	mb_compress(fin, emb, nmb);
      }
      else {
	cerr << "unsupported NEXUS block!" << endl;
	exit(1);
      }
    }
    //ignore anything that is not one of these blocks
  }
  //now, check to see which ones are turned on, and output the blocks into the uniform TRZ file.
  string outfile = file+".trz";
  ofstream fout;
  ifstream myfin;
  string myline;
  fout.open(outfile.c_str());
  if (!fout){
    cerr << "cannot open TRZ file for final writing!" << endl;
    exit(1);
  }
  fout << "#NEXUS" << endl;
  if (etaxa){
    print_out_temporary_file(".nexus-taxa", fout);
  }
  if (echar){
    print_out_temporary_file(".nexus-chars", fout);
  }
  if (eunalign){
    print_out_temporary_file(".nexus-unalign", fout);
  }
  if (edata){
    print_out_temporary_file(".nexus-data", fout);
  }
  if (eset){
    print_out_temporary_file(".nexus-sets", fout);
  }
  if (ecodon){
    print_out_temporary_file(".nexus-codons", fout);
  }
  if (edist){
    print_out_temporary_file(".nexus-dist", fout);
  }
  if (eassume){
    print_out_temporary_file(".nexus-assume", fout);
  }
  if (enote){
    print_out_temporary_file(".nexus-notes", fout);
  }
  if (epaup){
    print_out_temporary_file(".nexus-paup", fout);
  }
  if (emb){ 
    print_out_temporary_file(".nexus-mb", fout);
  }
  if (etree){
    print_out_temporary_file(".nexus-trees", fout);
  }
  fout.close();
  unlink(".mytempnexus");
}

void decompress_nexus(string file, int random_decompress, bool branch, bool calc_unique_trees, unsigned char consensus, bool quiet){
  ifstream fin;
  fin.open(file.c_str());
  if (!fin){
    cerr << "cannot open file for reading!" << endl;
    exit(2);
  }
  string outfile = file.substr(0, file.size()-4);
  outfile+=".e";
  ofstream fout, temp;
  fout.open(outfile.c_str());
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(2);
  }
  string line;
  //read in #NEXUS block
  getline(fin, line); 
  if (line != "#NEXUS"){
    cerr << "not a nexus file!" << endl;
    exit(1);
  }
  fout << line << endl << endl;
  int pos;
  string myid, type; 
  while (!fin.eof()){
    getline(fin, line);
    if (fin.eof())
      break;
    pos = line.find_first_of(" ");
    type = line.substr(0, pos);
    line = line.substr(pos+1);
    cout << "type is: " << type << endl;
    if (type == "TAXAL"){
      fout << "BEGIN TAXA;" << endl;
      taxa_decompress(fin, fout, line);
      fout << "END;" << endl << endl;
    }
    else if (type == "TREES"){
      cout << "We found the trees block!" << endl;
      fout << "BEGIN TREES;" << endl;
      trees_decompress(fin, fout);
      fout << "END;" << endl << endl;
    }
    else if (type == "CHARACTERS"){
      stringstream ss(line);
      unsigned int type, lines;
      ss >> lines >> type;
      if (type == 0)
	fout << "BEGIN CHARACTERS;" << endl;
      else
	fout << "BEGIN DATA;" << endl;
      common_decompress(fin, fout, lines);
      fout << "END;" << endl << endl;
    }
    else if (type == "DIST"){
      unsigned int lines = atoi(line.c_str()); 
      fout << "BEGIN DIST;" << endl;
      dist_decompress(fin, fout, lines);
      fout << "END;" << endl;
    }
    else if ( (type == "UNALIGNED") || (type == "SETS") || (type == "CODONS") || 
	      (type == "NOTES") || (type == "ASSUME") || (type == "PAUP") || 
	      (type == "MRBAYES") ){
      unsigned int lines = atoi(line.c_str());
      if (type != "ASSUME")
	fout << "BEGIN " << type << ";" << endl;
      else
	fout << "BEGIN ASSUMPTIONS;" << endl;
      common_decompress(fin, fout, lines);
      fout << "END;" << endl << endl;
    }
    else{
      cerr << "unrecognized block! " << type << " Exiting..." << endl;
      exit(1);
    }
  }
  fin.close();
}

//compresses a tree file into compressed form
void compress(string file1, bool bUbid, bool quiet, vector<string> nexTaxa, bool use) {
  FILE *fp;
  NEWICKTREE *newickTree;
  int err;
  unsigned long long M1=0;
  unsigned long long M2=0;

  //int offset = 0; //if doing parallel, may need to pass this again as a parameter
  struct timeval label_start;
  struct timeval label_end;
  struct timeval compress_start;
  struct timeval compress_end;
  struct timeval processing_start;
  struct timeval processing_end;
  struct timeval dfs_start;
  struct timeval dfs_end;
  struct timeval kill_start;
  struct timeval kill_end;
  struct timeval init_start;
  struct timeval init_end;
  struct timeval load_start;
  struct timeval load_end;
  struct timeval hash_start;
  struct timeval hash_end;
  double dfs_time = 0;
  double kill_time = 0;
  double load_time = 0;
  //char encoded[11] = "ABCDEFGHIJ";

  gettimeofday(&compress_start, NULL);
  //check if it is a nexus file
  //if (is_nexus(file1)){
  // nexus_compress(file1, bUbid, quiet);
  // return;
  //}
  gettimeofday(&label_start, NULL);

  LabelMap lm = collect_labels(file1); //collect labels and number of trees
  ///cout << "Hetero is: " << HETERO << endl;
  gettimeofday(&label_end, NULL);
  double label_time = label_end.tv_sec - label_start.tv_sec + (label_end.tv_usec - label_start.tv_usec) / 1.e6;
  if (!quiet)
    fprintf(stderr,"\nLabel Collection time: %g seconds \n",label_time);
  gettimeofday(&processing_start, NULL);
  gettimeofday(&init_start, NULL);
  initialize_hashtable(M1, M2); //initialize contents of hashtable
  gettimeofday(&init_end, NULL);
  double init_time = init_end.tv_sec - init_start.tv_sec + (init_end.tv_usec - init_start.tv_usec) / 1.e6;
  if (!quiet)
    fprintf(stderr,"\ninitialize_hashtable(): %g seconds \n",init_time);
  //step 0: check for weights

  WEIGHTED = check_for_weights(file1);
  if (WEIGHTED){
    fprintf(stderr, "Automatic check: Trees are weighted!\n");
  }
  else{
    fprintf(stderr,"Automatic check: Trees are unweighted!\n");
  }
  //step 1: read in trees sequentially and hash each one
  fp = fopen(file1.c_str(), "r");


  if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}
  for (unsigned int treeIdx=0; treeIdx<NUM_TREES; ++treeIdx) {
    gettimeofday(&load_start, NULL);
    newickTree = loadnewicktree2(fp, &err);
    gettimeofday(&load_end, NULL);
    load_time += load_end.tv_sec - load_start.tv_sec + (load_end.tv_usec - load_start.tv_usec) / 1.e6;
    if (!newickTree) {
      switch (err) {
      case -1:
        printf("Out of memory\n");
        break;
      case -2:
        printf("parse error on line %u\n", treeIdx);
        break;
      case -3:
        printf("Can't load file\n");
        break;
      default:
        printf("Error %d\n", err);
      }
    }
    else {
      unsigned int numBitstr=0;
      gettimeofday(&dfs_start, NULL);
      //cout << "How about here? Iteration: " << treeIdx << endl; 
      dfs_compute_hash(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
      if (HETERO)
	break;
      gettimeofday(&dfs_end, NULL);
      gettimeofday(&kill_start, NULL);
      killnewicktree(newickTree);
      gettimeofday(&kill_end, NULL);
      dfs_time += dfs_end.tv_sec - dfs_start.tv_sec + (dfs_end.tv_usec - dfs_start.tv_usec) / 1.e6;
      kill_time += kill_end.tv_sec - kill_start.tv_sec + (kill_end.tv_usec - kill_start.tv_usec) / 1.e6;
    }
  }
  fclose(fp);
  gettimeofday(&processing_end, NULL);
  double processing_time = processing_end.tv_sec - processing_start.tv_sec + (processing_end.tv_usec - processing_start.tv_usec) / 1.e6;
  if (HETERO){
    HCHECK = true;
    cerr << "WARNING! Detected Heterogeneous trees! Restarting procedure..." << endl;
    clear_hashtable(); //clear hashtable
    lm.clear();
    lm = collect_labels(file1);
    initialize_hashtable(M1, M2); //initialize contents of hashtable
    fp = fopen(file1.c_str(), "r");
    if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}
    for (unsigned int treeIdx=0; treeIdx<NUM_TREES; ++treeIdx) {
      gettimeofday(&load_start, NULL);
      newickTree = loadnewicktree2(fp, &err);
      gettimeofday(&load_end, NULL);
      load_time += load_end.tv_sec - load_start.tv_sec + (load_end.tv_usec - load_start.tv_usec) / 1.e6;
      if (!newickTree) {
	switch (err) {
	case -1:
	  printf("Out of memory\n");
	  break;
	case -2:
	  printf("parse error on line %u\n", treeIdx);
	  break;
	case -3:
	  printf("Can't load file\n");
	  break;
	default:
	  printf("Error %d\n", err);
	}
      }
      else {
	unsigned int numBitstr=0;
	gettimeofday(&dfs_start, NULL);
	//cout << "How about here? Iteration: " << treeIdx << endl; 
	dfs_compute_hash(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
	gettimeofday(&dfs_end, NULL);
	gettimeofday(&kill_start, NULL);
	killnewicktree(newickTree);
	gettimeofday(&kill_end, NULL);
	dfs_time += dfs_end.tv_sec - dfs_start.tv_sec + (dfs_end.tv_usec - dfs_start.tv_usec) / 1.e6;
	kill_time += kill_end.tv_sec - kill_start.tv_sec + (kill_end.tv_usec - kill_start.tv_usec) / 1.e6;
      }
    }
  }
  if (!quiet){
    fprintf(stderr,"\nTree Processing time: %g seconds \n",processing_time);
    fprintf(stderr,"\nloadnewicktree2() time: %g seconds \n",load_time);
    fprintf(stderr,"\ndfs_compute_hash_w_bitstring() time: %g seconds \n",dfs_time);
    fprintf(stderr,"\nkillnewicktree() time: %g seconds \n",kill_time);
  }
  gettimeofday(&hash_start, NULL);
  unsigned long tree_counter = 0; //this determines if the trees are multifurcating or binary
  vector<short> n_biparts(NUM_TREES,0);

  float threshold = NUM_TREES - (float)(NUM_TREES*P)/100;
  //cout << "threshold is: " << threshold << endl;

  string outfile = file1 + ".trz";
  //if (!quiet)
    fprintf(stderr, "File to be written is: %s\n",outfile.c_str());

  FILE * fout;

  fout = fopen(outfile.c_str(), "wb");
  fprintf(fout, "TAXA ");

  LabelMap lm_sorted;
  string label;
  if (use && nexTaxa.size()){
    //this is for renaming the taxa for nexus compress
    int label_loc;
    string newlabel;
    for (unsigned int ls = 0; ls < lm.size(); ls++){
      label = lm.name(ls); //get the current label
      label_loc = atoi(label.c_str()); //get the location it represents
      newlabel = nexTaxa[label_loc-1]; //get the new label 
      lm.rename(label, newlabel, ls); //rename it
    }
  }
  for (unsigned int ls = 0; ls < lm.size(); ls++){
    label = lm.name(ls);
    lm_sorted.push(label);
  }
  lm_sorted.sortTaxa();
  if (use){
    fprintf(fout, "1");
    for (unsigned int i = 1; i < NUM_TAXA; i++)
      fprintf(fout, ":%d", i+1);
  }
  else{
    fprintf(fout, "%s", lm_sorted.name(0).c_str());
    for (unsigned int i = 1; i < NUM_TAXA; i++) { 
      fprintf(fout, ":%s", lm_sorted.name(i).c_str());
    }
  }
    fprintf(fout, "\n");
  fprintf(fout, "NTREES %d", NUM_TREES);
  if (HETERO)
    fprintf(fout, " H");
  fprintf(fout, "\n");
  vector< pair< unsigned, pair<string, string > > > ordered_bitstrings;
  //float * branch_lengths = NULL;
  unsigned int numtrees_strsize = 0;
  unsigned int encode_size = 0;
 
  if (WEIGHTED){ //if weighted
    stringstream tmpss;
    string numtrees_str;
    tmpss << NUM_TREES;
    tmpss >> numtrees_str;
    numtrees_strsize = numtrees_str.size();
    encode_size = (numtrees_strsize+1)/2;
  }

  unsigned int bipart_counter = 0;
  stringstream ss;

  //step2: allocate inverted index and helper_sizes array --> initialize the latter
  helper = (unsigned int **)malloc(NUM_TREES*sizeof(unsigned int*));
  helper_sizes = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  if ( (helper == NULL) || (helper_sizes == NULL) ){
    cerr << "cannot allocate inverted index arrays!" << endl;
    exit(2);
  }
  for (unsigned int i = 0; i < NUM_TREES; i++)
    helper_sizes[i] = 0;

  //step 3: count up the number of unique bipartitions, compute bitstring order and figure out inverted index dimensions
  count_bipartitions(bipart_counter, tree_counter, ordered_bitstrings, lm, lm_sorted);
  NUMBIPART = bipart_counter;

  unsigned int temp_size;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    temp_size = helper_sizes[i];
    helper[i] = (unsigned int *)malloc(temp_size * sizeof(unsigned int));
    helper_sizes[i] = 0;
  }
  
  //step 4: populate inverted index
  unsigned int bipart_id = 0;
  unsigned int item, loc;
  for (unsigned int hti = 0; hti < vvec_hashrf._hashtab2.size(); hti++) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size();    
    if (sizeVec) {
      for (unsigned int i = 0; i < sizeVec; i++) {	  
	unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
	for (unsigned int j  = 0; j < sizeTreeIdx; j++){
	  item = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	  //index[item].push_back(bipart_id);
	  loc = helper_sizes[item];
	  helper[item][loc] = bipart_id;
	  helper_sizes[item]++;
	}
	bipart_id++;
      } //for everything in sizevec
    }
  }

  if (!quiet)
    fprintf(stderr, "Rehashing trees to determine uniqueness\n");
  //step 5: rehash
  HashRFMap vvec_unique;
  vector<unsigned int> unique_ids;
  bool are_unique = false;
  unsigned int duplicates = 0;
  vvec_unique.uhashfunc_init_unique(NUM_TREES, NUM_TAXA, bipart_counter, C, NEWSEED);
  //new m1 and m2 values
  unsigned long long um1 = vvec_unique._HF.getM1();
  unsigned long long um2 = vvec_unique._HF.getM2();
  vvec_unique._hashtab2.resize(um1); //resize vvec_unique
  unsigned int locsize = 0;
  for (unsigned int i = 0; i < NUM_TREES; i++){ 
    unsigned int pos = 0;
    locsize = helper_sizes[i];
    unsigned long long temp1 = 0;
    unsigned long long temp2 = 0;
    bool * cbs = new bool[NUMBIPART];
    for (unsigned int j = 0; j < NUMBIPART; j++)
      cbs[j] = 0;
    for (unsigned int j = 0; j < locsize; j++){
      pos = helper[i][j];
      cbs[pos] = 1;
      temp1+= vvec_unique._HF.getA1(pos);
      temp2+= vvec_unique._HF.getA2(pos);
    }
    temp1 = temp1 % um1;
    temp2 = temp2 % um2;
    vvec_unique.hashing_bs_unique(i, NUMBIPART, temp1, temp2, cbs);
    delete [] cbs;
  }
  for (unsigned int hti=0; hti<vvec_unique._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_unique._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; ++i) {
	unique_ids.push_back(vvec_unique._hashtab2[hti][i]._vec_treeidx[0]);
	if (vvec_unique._hashtab2[hti][i]._vec_treeidx.size() > 1)
	  duplicates++;
      }
    } //end if sizevec
  } //end hashtable traversal
  sort(unique_ids.begin(), unique_ids.end());
  
  if (!quiet)
    fprintf(stderr, "Found %u unique trees\n", (unsigned int)unique_ids.size());
  if (unique_ids.size() == NUM_TREES){
    fprintf(stderr, "All trees are unique!\n");
    fprintf(fout, "UNIQUE_T all\n");
  }
  else{
    are_unique = true;
    fprintf(fout, "UNIQUE_T %d\n", (unsigned int)unique_ids.size());
    //free(new_vec);
  }
  
  //step 6: sort the bitstrings once:
  fprintf(fout, "NBIPARTITIONS %d\n", NUMBIPART);
  sort(ordered_bitstrings.begin(), ordered_bitstrings.end(), sort_bitstrings);
  
  vector< pair< unsigned int, pair<string, string> > >::iterator b_itr = ordered_bitstrings.begin();
  //DEBUG: Print out sorted order of the bitstrings:
  //vector<unsigned int> histogram(NUM_TAXA-1, 0);
  /*while (b_itr != ordered_bitstrings.end()) { 
   pair<string, string> tmp_pair = b_itr -> second;
   cout << b_itr->first << ": " << tmp_pair.first << endl;
   b_itr++;
   }*/
    //histogram[my_tmp_val]++;
  // cout << b_itr -> first << ":" << b_itr -> second << endl;
  //  b_itr++;
  // }
  //exit(0);
  /*cout << "print out histogram!" << endl;
  float temp_hist = 0.0;
  float denominator = 0.0;
  for (unsigned int i = 0; i < NUM_TAXA-1; ++i){
    if (histogram[i]){
      denominator = (float)NUM_TAXA/i;
      temp_hist = ceil(histogram[i]/denominator);
      cout << i << ": " << temp_hist << endl;
    }
  }
  exit(0);*/
  string tempval, digitstr;
  //short digit;
  stringstream stm;

  //step 7. in the order of the sorted bitstrings, compress bitstrings, associated tree ids, and branch lengths

  b_itr = ordered_bitstrings.begin();

  //create array to aid in tree compression
  int * myuniqhelper = NULL;
  if (are_unique){
    myuniqhelper = (int*)malloc(NUM_TREES*sizeof(int));
    if (myuniqhelper == NULL){
      cerr << "cannot allocate!" << endl;
      exit(2);
    }
    unsigned int myplace = 0;
    for (unsigned int i = 0; i < NUM_TREES; i++){
      if (i == unique_ids[myplace]){
	myuniqhelper[i] = myplace;
	myplace++;
      }
      else
	myuniqhelper[i] = -1;
    }

    //for (unsigned int i  = 0; i < NUM_TREES; i++)
    //  cout << i << " " << myuniqhelper[i] << endl;
    //exit(0);
  }


  vector<bool> check(NUM_TREES, false);
  vector<bool> checku(unique_ids.size(), false);
  unsigned int x = 0, y = 0, p = 0;
  string place, xval;
  unsigned long my_tree_count = 0;
  unsigned int deb = 0;
  string rle_bipart;
  float u_threshold = (float)unique_ids.size()/2;
  //now create the new labelMap
  while (b_itr != ordered_bitstrings.end()) {
    //cout << "We get here! " << deb << endl;
    deb++;
    pair<string, string> temp_pair = b_itr->second;
    rle_bipart = temp_pair.first;
    place = temp_pair.second;
    p = place.find_first_of(".");
    xval = place.substr(0,p);
    place = place.substr(p+1);
    x = atoi(xval.c_str());
    y = atoi(place.c_str());
    unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[x][y]._vec_treeidx.size();
    //cout << "sizeTreeIdx is: " << sizeTreeIdx << endl;
    my_tree_count += sizeTreeIdx;
    //bool * temp = vvec_hashrf._hashtab2[x][y]._bs;
    //unsigned int bs_size = vvec_hashrf._hashtab2[x][y]._bs_size;
    fprintf(fout, "%s ", rle_bipart.c_str());
    //tree id compression procedure (RLE)
    if( sizeTreeIdx < NUM_TREES) {	     
      //cout << "we get here!(-) deb is: " << deb << endl;
      unsigned int * new_vec;
      unsigned int newlength;
      if ( sizeTreeIdx > threshold) { //if it is greater than the threshold value, we're using thresholding
	newlength = NUM_TREES - sizeTreeIdx;
	new_vec = (unsigned int *)malloc(newlength * sizeof(unsigned int));
      }
      else{
	new_vec = (unsigned int *)malloc(sizeTreeIdx * sizeof(unsigned int));
	newlength = sizeTreeIdx;
      }
      unsigned int nx = 0;
      unsigned int true_count = 0;
      if (are_unique){	      
	for (unsigned int j =0; j < sizeTreeIdx; j++) { //first determine the ones that exist
	  int checker = vvec_hashrf._hashtab2[x][y]._vec_treeidx[j];
	  if (myuniqhelper[checker] >= 0){
	    nx = myuniqhelper[checker];
	    checku[nx] = true;
	    true_count++;
	  }
	  n_biparts[checker]++;
	} 
	nx = 0;
	if (true_count > u_threshold){ //compact
	  for (unsigned int j = 0; j < unique_ids.size(); j++) { //then store the ones that don't
	    if (checku[j] == false) { 
	      new_vec[nx] = j;
	      nx++;
	    }
	    else {
	      checku[j] = false;
	    }
	  }
	  string tid_encode = encode(new_vec, nx);
	  fprintf(fout, "-:%d:%s",nx, tid_encode.c_str());
	  free(new_vec);
	}
	else{ //don't compact
	  for (unsigned int j = 0; j < unique_ids.size(); j++) { //then store the ones that don't
	    if (checku[j] == true) { 
	      new_vec[nx] = j;
	      nx++;
	      checku[j] = false;
	    }
	  }
	  string tid_encode = encode(new_vec, nx);
	  fprintf(fout, "+:%d:%s",nx, tid_encode.c_str());
	  free(new_vec);
	}
      }
      else{ //if all are unique, follow the regular procedure
	nx = 0;
	if (sizeTreeIdx > threshold){
	  for (unsigned int j =0; j < sizeTreeIdx; j++) { //first determine the ones that exist
	    int checker = vvec_hashrf._hashtab2[x][y]._vec_treeidx[j];
	    check[checker] = true;
	    n_biparts[checker]++;
	  } 
	  for (unsigned int j = 0; j < NUM_TREES; ++j) { //then store the ones that don't
	    if (check[j] == false) { 
	      new_vec[nx] = j;
	      nx++;
	    }
	    else {
	      check[j] = false;
	    }
	  }
	  assert(newlength == nx);
	  string tid_encode = encode(new_vec, nx);
	  fprintf(fout, "-:%d:%s", nx, tid_encode.c_str());
	  free(new_vec);
	}
	else{
	  for (unsigned int i = 0; i < sizeTreeIdx; i++){
	    new_vec[i] = vvec_hashrf._hashtab2[x][y]._vec_treeidx[i];
	  }
	  nx = sizeTreeIdx;
	  string tid_encode = encode(new_vec, nx);
	  fprintf(fout, "+:%d:%s", nx, tid_encode.c_str());
	  free(new_vec);
	}
      }
    }
    else if (sizeTreeIdx == NUM_TREES) { //if sizeTreeIdx == NUM_TREES
      //cout << "we get here (==)! deb is: " << deb << endl;
      fprintf(fout, "-:0:"); 
      
      for (unsigned int j = 0; j < sizeTreeIdx; ++j) {	
	int temp  = vvec_hashrf._hashtab2[x][y]._vec_treeidx[j];
	n_biparts[temp]++;
      }
    }
    else { 
      cerr << "HUGE FREAKING ERROR!" << endl;
      cerr << "size is: " << sizeTreeIdx << " vs. NUM_TREES: " << NUM_TREES << endl;
      exit(3);
    }
    if (WEIGHTED){
      string branch_encode = encode_branch_lengths(x, y, sizeTreeIdx, encode_size);
      fprintf(fout, "%s", branch_encode.c_str());
    }
    fprintf(fout, "\n");
    //cout << endl;
    //printf("\n");

    b_itr++;
  }
  assert(bipart_counter == ordered_bitstrings.size());
  assert(tree_counter == my_tree_count);

  //step 8: Print out duplicate information:
  fprintf(fout, "DUPLICATES %d\n", duplicates);
  unsigned int * new_val = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  unsigned int newlength = 0;
  unsigned int tempsize = 0;
  unsigned int firstloc = 0;
  unsigned int myelem; 
  vector< vector<unsigned int> > all_dups;
  all_dups.resize(NUM_TREES);
  for (unsigned int hti=0; hti<vvec_unique._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_unique._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; ++i) {
	if (vvec_unique._hashtab2[hti][i]._vec_treeidx.size() > 1){
	  tempsize = vvec_unique._hashtab2[hti][i]._vec_treeidx.size();
	  firstloc = vvec_unique._hashtab2[hti][i]._vec_treeidx[0];
	  for (unsigned int j  =1; j < tempsize; j++){
	    myelem = vvec_unique._hashtab2[hti][i]._vec_treeidx[j];
	    //newlength++;
	    all_dups[firstloc].push_back(myelem);
	  }
	  //encode everything
	}
      }
    } //end if sizevec
  } //end hashtable traversal


  for (unsigned int i = 0; i < NUM_TREES; i++){
    tempsize = all_dups[i].size();
    if(tempsize){
      new_val[newlength] = i;
      newlength++;
      for (unsigned int j = 0; j < tempsize; j++){
	new_val[newlength] = all_dups[i][j];
	newlength++;
      }
      string dup_encode = encode(new_val, newlength);
      newlength = 0;
      fprintf(fout, "%s\n", dup_encode.c_str());
    }
  }
  
  fclose(fout);
  gettimeofday(&hash_end, NULL);
  double hash_time = hash_end.tv_sec - hash_start.tv_sec + (hash_end.tv_usec - hash_start.tv_usec) / 1.e6;
  printf("%s: %d\n", outfile.c_str(), (int)bipart_counter);
  if (!quiet)  {
    fprintf(stderr,"\nBipartition Table Processing time: %g seconds \n",hash_time);
    fprintf(stderr,"\nNumber of Unique Bipartitions: %d\n", (int)bipart_counter);
    fprintf(stderr,"Total Number of Bipartitions: %lu\n", tree_counter);
    if ( (!WEIGHTED && (tree_counter < ((NUM_TAXA-3)*NUM_TREES))) || (WEIGHTED && (tree_counter < 2*(NUM_TAXA-1)*NUM_TREES)) )  
      fprintf(stderr, "TYPE:multifurcating\n");
    else 
      fprintf(stderr, "TYPE:binary\n");
  }
  gettimeofday(&compress_end, NULL);
  double compress_time = compress_end.tv_sec - compress_start.tv_sec + (compress_end.tv_usec - compress_start.tv_usec) / 1.e6;
  if (!quiet)  
    fprintf(stderr,"\nCompress Function time: %g seconds \n",compress_time);
}

void populate_integrals(unsigned int * hold_integrals, string branches_int, unsigned int encode_size){
  unsigned int integral_size = branches_int.size();
  unsigned int intpos = 0;
  unsigned int treeloc, m;
  int treeval, n;
  unsigned char v;
  treeloc = 0;
  treeval = 0;
  m = 0;
  for (unsigned int i  =0; i < NUM_TREES; i++)
    hold_integrals[i] = 0;
  while (intpos < integral_size){
    while (m < encode_size){ //while we're less than the encoded size
      v = branches_int[intpos];
      treeloc+= v;
      treeloc-=33;
      m++;
      if (m < encode_size)
	treeloc *= 100;
      intpos++;
    }
    treeval = branches_int[intpos]; //next, get the tree value
    treeval-=48;
    intpos++;
    if (treeval < 0){
      n = 0;
      treeval = 0;
      while ( (n <= 9) && (n >= 0)){
	treeval*=10;
	n = branches_int[intpos];
	n-=48;
	treeval+= n;
	intpos++;
	n = branches_int[intpos];
	n-=48;
      }
      intpos++;
    }
    hold_integrals[treeloc] = treeval;
    m = 0;
    treeloc = 0;
    treeval = 0;
  }
}

void decompress_branch(unsigned int * hold_integrals, vector<unsigned int> my_set_of_ids, vector< vector<float> > & all_branches, string branches_frac){
  unsigned int myplace = 0;
  unsigned int numer = 0;
  unsigned int denom = 100;
  float weight = 0;
  unsigned char v;
  unsigned int val;
  unsigned int count = my_set_of_ids.size();
  if (count < NUM_TREES){
    for (unsigned int i = 0; i < count; ++i){
      unsigned int temp = my_set_of_ids[i];
      while (numer != 3){
	v = branches_frac[myplace];
	val = v;
	val -= 33;
	weight += (float)val/denom;
	denom*=100;
	myplace++;
	numer++;
      }
      weight+=hold_integrals[temp];
      //printf("%f\n", weight);
      hold_integrals[temp] = 0;
      all_branches[temp].push_back(weight);
      numer = 0;
      denom = 100;
      weight = 0;
    }
  }
  else{
    for (unsigned int i = 0; i < NUM_TREES; ++i){
      while (numer != 3){
	v = branches_frac[myplace];
	val = v;
	val -= 33;
	weight += (float)val/denom;
	denom*=100;
	myplace++;
	numer++;
      }
      weight+=hold_integrals[i];
      //printf("%f\n", weight);
      hold_integrals[i] = 0;
      all_branches[i].push_back(weight);
      numer = 0;
      denom = 100;
      weight = 0;
    }
  }
}

void parse_trz_file(vector< vector< bool *> > &all_bs, vector< vector<unsigned int> > & inverted_index, vector< bool* > & list_bs, vector< vector<float> > &all_branches, string file, LabelMap &lm, vector<bool> & is_dup, vector< vector<unsigned int> > &bs_sizes){  
  string mycount, ids, bipartition, line_type, str, taxa, treeline, bitstring, branches_int, branches_frac;
  unsigned int bipart_count, total_bipartitions, ntaxa, nbipart, num_unique, ndup_lines;
  vector<bool> check;
  vector<unsigned int> true_ids;
  vector< vector<unsigned int> > dups;
  bipart_count = 0;
  total_bipartitions = 0; 

  ifstream fin(file.c_str(), ios::binary);
  if (!fin) {
    cerr << "cannot open file!\n";
    exit(2);
  }
  
  //read in taxa
  getline(fin, str);  
  int pos = str.find_first_of(" ");
  line_type = str.substr(0, pos);
  if (line_type != "TAXA"){
    cerr << "Error! No taxa labels identified for file! Exiting..\n";
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  stringstream ss(str);
  ntaxa = 0;
  while(getline(ss, taxa, ':')){ 
    lm.push(taxa);
     ntaxa++;
  }
  NUM_TAXA = ntaxa;
  ss.clear();

  //read in number of trees
  getline(fin, str);
  NUM_TREES = get_ntrees(str);

  //read in number of unique trees
  getline(fin, str);
  num_unique = get_unique(str, NUM_TREES);

  all_bs.resize(NUM_TREES);
  bs_sizes.resize(NUM_TREES);
  inverted_index.resize(NUM_TREES);
  //list_bs.resize(num_unique);
  all_branches.resize(NUM_TREES);
  check.resize(NUM_TREES);
  true_ids.resize(num_unique);
  dups.resize(NUM_TREES);
  is_dup.resize(NUM_TREES);
  unsigned int * hold_integrals = NULL;
  //read in number of bipartitions to be read
  getline(fin, str);
  parse_and_get(str, "NBIPARTITIONS", nbipart);
  NUMBIPART = nbipart;
  //skip nbipart lines, get the duplicate information
  for (unsigned int i = 0; i < nbipart; i++)
    getline(fin, str);
  getline(fin, str); //now read the duplicates line
  parse_and_get(str, "DUPLICATES", ndup_lines);
  //cout << "ndup_lines is: " << ndup_lines << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++)
    is_dup[i] = 0;
  unsigned int decode_size = 0;
  unsigned int * found, dec_loc, dec_val;
  found = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  for (unsigned int i = 0; i < ndup_lines; i++){
    getline(fin, str); //read in the line
    pos = str.find_first_of("\n");
    str = str.substr(0, pos); //get rid of the newline
    decode_size = decode(str, found); //decode the line
    dec_loc = found[0];
    for (unsigned int j = 1; j < decode_size; j++){
      dec_val = found[j];
      dups[dec_loc].push_back(dec_val); //add the elements of the array to the associated dup structure
      is_dup[dec_val] = 1; //also set those locations in our bool structure to a 1
    }
  }

  //populate true_ids
  unsigned int tempj = 0;
  for (unsigned int i = 0; i < num_unique; i++){
    if (is_dup[tempj]){
      while (is_dup[tempj])
	tempj++;
    }
    true_ids[i] = tempj;
    tempj++;
  }

  //debug: print out the data structures
  /*cout << "true_ids: " << endl;
  for (unsigned int i  = 0; i < true_ids.size(); i++)
    cout << i << "-->" << true_ids[i] << endl;
  cout << endl;

  cout << "dups: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    cout << i << ": ";
    for (unsigned int j = 0; j < dups[i].size(); j++){
      cout << dups[i][j] << " ";
    }
    cout << endl;
    }*/
  fin.close(); //now close the file
  fin.open(file.c_str(), ios::binary); //reopen
  if (!fin) {
    cerr << "cannot open file!\n";
    exit(2);
  }
  for (unsigned int i  =0; i < 4; i++) //go back to where we need to be in order to start reading in the bipartitions
    getline(fin, str); 

  unsigned int counter = 0;
  unsigned int count = 0;
  unsigned int bipart_loc = 0;
  unsigned int numtrees_strsize, my_count;
  stringstream tmpss;
  string numtrees_str;
  unsigned int encode_size = 0;
  vector<unsigned int> my_set_of_ids;
  tmpss << NUM_TREES;
  tmpss >> numtrees_str;
  numtrees_strsize = numtrees_str.size();
  encode_size = (numtrees_strsize+1)/2;
  unsigned int maxLength = 0;
  while ( counter < nbipart) { 
    getline(fin, str);  
    pos = str.find_first_of(" ");
    bitstring = str.substr(0, pos); //contains bitstring
    str = str.substr(pos+1);
    pos = str.find_first_of("\n");
    treeline = str.substr(0, pos);

    //process bitstring first
    maxLength = get_bitstring_length(bitstring);//determine length of bitstring maxLength
    bool *bs = new bool[maxLength]; //allocate it to be maxLength	
    decode_bitstring(bitstring, bs, maxLength);

    //next, process tree line
    //first, determine the number of TIDs in the line:
    pos = treeline.find_first_of(":");
    line_type = treeline.substr(0, pos);
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(":");
    mycount = treeline.substr(0, pos);
    count = atoi(mycount.c_str());
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(" ");
    if (pos != -1){
      if (!WEIGHTED){
	fprintf(stderr, "\nFile has branch lengths! Weighted option: ON\n\n");	
	WEIGHTED = true;
	hold_integrals = (unsigned int*)malloc(NUM_TREES*sizeof(unsigned int));
	for (unsigned int x = 0; x < NUM_TREES; x++){
	  hold_integrals[x] = 0;
	}
      }
    }
    if (WEIGHTED){
      ids = treeline.substr(0, pos);
      treeline = treeline.substr(pos+1);
      pos = treeline.find_first_of(" ");
      branches_int = treeline.substr(0, pos); //integral portion of branches
      treeline = treeline.substr(pos+1);
      pos = treeline.find_first_of("\n");
      branches_frac = treeline.substr(0, pos); //fractional component of branches
    }
    else{
      pos = treeline.find_first_of("\n");
      ids = treeline.substr(0, pos);
    }

    //process tree ids
    if (count == 0) { 
      for (unsigned int b = 0; b < NUM_TREES; ++b) { 
	all_bs[b].push_back(bs); //push bipartition into everything
	bs_sizes[b].push_back(maxLength); //push bipartition length into everything
	inverted_index[b].push_back(bipart_loc);
      }
      //vec_bs.push_back(bs); //this is a strict consensus bipartition       
      if (WEIGHTED){
	if (branches_int != "")
	  populate_integrals(hold_integrals, branches_int, encode_size);
	my_set_of_ids.resize(NUM_TREES);

	decompress_branch(hold_integrals, my_set_of_ids, all_branches, branches_frac);
      }
      bipart_count++;
    }      
    else { //count != 0 (so we have tree ids to process)
      my_count = decode(ids, found);
      assert(my_count == count);

      if (line_type == "-") { //compressed line
	for (unsigned int i = 0; i < count; ++i) {
	  unsigned int temp = found[i];
	  check[temp] = true;
	}
	unsigned int true_id, sec_id;
	//cout << "we do get here though? num_unique is: " << num_unique << endl;
	for (unsigned int i = 0; i < num_unique; ++i) {
	  if (check[i] == false) {
	    true_id = true_ids[i];
	    my_set_of_ids.push_back(true_id);
	    all_bs[true_id].push_back(bs);
	    bs_sizes[true_id].push_back(maxLength);
	    inverted_index[true_id].push_back(bipart_loc);
	    if (dups[true_id].size() > 0){
	      for (unsigned int j = 0; j< dups[true_id].size(); j++){
		sec_id = dups[true_id][j];
		all_bs[sec_id].push_back(bs);
		bs_sizes[sec_id].push_back(maxLength);
		inverted_index[sec_id].push_back(bipart_loc);
		my_set_of_ids.push_back(sec_id);
	      }
	    }
	    //cout << "end of for stat" << endl;
	  }
	  else
	    check[i] = false;
	}
	if (WEIGHTED){	  
	  sort(my_set_of_ids.begin(), my_set_of_ids.end());
	  if (branches_int != "")
	    populate_integrals(hold_integrals, branches_int, encode_size);
	  decompress_branch(hold_integrals, my_set_of_ids, all_branches, branches_frac);
	}
      }
      else { //line is not compressed
	unsigned int true_id, sec_id;
	for (unsigned int i = 0; i < count; ++i) { 
	  unsigned int temp = found[i];
	  true_id = true_ids[temp];
	  all_bs[true_id].push_back(bs);
	  bs_sizes[true_id].push_back(maxLength);
	  inverted_index[true_id].push_back(bipart_loc);
	  my_set_of_ids.push_back(true_id);
	  if (dups[true_id].size() > 0){
	    for (unsigned int j = 0; j < dups[true_id].size(); j++){
	      sec_id = dups[true_id][j];
	      all_bs[sec_id].push_back(bs);
	      bs_sizes[sec_id].push_back(maxLength);
	      inverted_index[sec_id].push_back(bipart_loc);
	      my_set_of_ids.push_back(sec_id);
	    }
	  }
	} 
	if (WEIGHTED){
	  sort(my_set_of_ids.begin(), my_set_of_ids.end());	 
	  if (branches_int != "")
	    populate_integrals(hold_integrals, branches_int, encode_size);
	  decompress_branch(hold_integrals, my_set_of_ids, all_branches, branches_frac);
	}
      }
      bipart_count++;
    } //end if count != 0
    list_bs.push_back(bs); //push into the list 
    my_set_of_ids.clear();
    bipart_loc++;
    counter++;
  }

  free(found);
  if (bipart_count != counter) { 
    cerr << "ERROR! Bipartitions not processed correctly!" << endl;
    cout << "bipart_count: " << bipart_count << endl;
    cout << "counter: " << counter << endl;
    exit(1);
  }
  assert(lm.size()!=0); 

} 

void parse_trz_hcomp(string file, vector< vector<unsigned int> > & inverted_index, vector< vector<float> > &all_branches, LabelMap &lm, unsigned int hashtable_length, unsigned int * row){  
  string mycount, ids, bipartition, line_type, str, taxa, treeline, bitstring, branches_int, branches_frac;
  unsigned int bipart_count, total_bipartitions, ntaxa, nbipart, ndup_lines;
  vector<bool> check;
  vector<unsigned int> true_ids;
  vector< vector<unsigned int> > dups;
  unsigned int num_unique;
  vector <bool> is_dup;
  bipart_count = 0;
  total_bipartitions = 0; 

  ifstream fin(file.c_str(), ios::binary);
  if (!fin) {
    cerr << "cannot open file!\n";
    exit(2);
  }
  
  //read in taxa
  getline(fin, str);  
  int pos = str.find_first_of(" ");
  line_type = str.substr(0, pos);
  if (line_type != "TAXA"){
    cerr << "Error! No taxa labels identified for file! Exiting..\n";
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  stringstream ss(str);
  ntaxa = 0;
  while(getline(ss, taxa, ':')){ 
    lm.push(taxa);
     ntaxa++;
  }
  NUM_TAXA = ntaxa;
  ss.clear();

  //read in number of trees
  getline(fin, str);
  NUM_TREES = get_ntrees(str);

  //read in number of unique trees
  getline(fin, str);
  num_unique = get_unique(str, NUM_TREES);

  //resive data structures
  inverted_index.resize(NUM_TREES);
  all_branches.resize(NUM_TREES);
  check.resize(NUM_TREES);
  true_ids.resize(num_unique);
  dups.resize(NUM_TREES);
  is_dup.resize(NUM_TREES);
  unsigned int * hold_integrals = NULL;
  //read in number of bipartitions to be read
  getline(fin, str);
  parse_and_get(str, "NBIPARTITIONS", nbipart);
  NUMBIPART = nbipart;
  hashtable_length = NUMBIPART;

  //skip nbipart lines, get the duplicate information
  for (unsigned int i = 0; i < nbipart; i++)
    getline(fin, str);
  getline(fin, str); //now read the duplicates line
  parse_and_get(str, "DUPLICATES", ndup_lines);
  //cout << "ndup_lines is: " << ndup_lines << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++)
    is_dup[i] = 0;
  unsigned int decode_size = 0;
  unsigned int * found, dec_loc, dec_val;
  found = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  for (unsigned int i = 0; i < ndup_lines; i++){
    getline(fin, str); //read in the line
    pos = str.find_first_of("\n");
    str = str.substr(0, pos); //get rid of the newline
    decode_size = decode(str, found); //decode the line
    dec_loc = found[0];
    for (unsigned int j = 1; j < decode_size; j++){
      dec_val = found[j];
      dups[dec_loc].push_back(dec_val); //add the elements of the array to the associated dup structure
      is_dup[dec_val] = 1; //also set those locations in our bool structure to a 1
    }
  }

  //populate true_ids
  unsigned int tempj = 0;
  for (unsigned int i = 0; i < num_unique; i++){
    if (is_dup[tempj]){
      while (is_dup[tempj])
	tempj++;
    }
    true_ids[i] = tempj;
    tempj++;
  }

  //debug: print out the data structures
  /*cout << "true_ids: " << endl;
  for (unsigned int i  = 0; i < true_ids.size(); i++)
    cout << i << "-->" << true_ids[i] << endl;
  cout << endl;

  cout << "dups: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    cout << i << ": ";
    for (unsigned int j = 0; j < dups[i].size(); j++){
      cout << dups[i][j] << " ";
    }
    cout << endl;
    }*/
  fin.close(); //now close the file
  fin.open(file.c_str(), ios::binary); //reopen
  if (!fin) {
    cerr << "cannot open file!\n";
    exit(2);
  }
  for (unsigned int i  =0; i < 4; i++) //go back to where we need to be in order to start reading in the bipartitions
    getline(fin, str); 

  unsigned int counter = 0;
  unsigned int count = 0;
  unsigned int bipart_loc = 0;
  unsigned int numtrees_strsize, my_count;
  stringstream tmpss;
  string numtrees_str;
  unsigned int encode_size = 0;
  vector<unsigned int> my_set_of_ids;
  tmpss << NUM_TREES;
  tmpss >> numtrees_str;
  numtrees_strsize = numtrees_str.size();
  encode_size = (numtrees_strsize+1)/2;
  unsigned int mycounter = 0;
  hashtable = (unsigned int **)malloc(NUMBIPART*sizeof(unsigned int *)); // the hash table is NUMBIPART long 
  hash_lengths = (unsigned int *)malloc(NUMBIPART*sizeof(unsigned int)); //allocate corresponding array for this
  for (unsigned int i = 0; i < NUMBIPART; i++)
    hash_lengths[i] = 0;
  while ( counter < nbipart) { 
    getline(fin, str);  
    pos = str.find_first_of(" ");
    bitstring = str.substr(0, pos); //contains bitstring
    str = str.substr(pos+1);
    pos = str.find_first_of("\n");
    treeline = str.substr(0, pos);

    //process bitstring first
    bool *bs = new bool[NUM_TAXA];	
    decode_bitstring(bitstring, bs, NUM_TAXA);

    //next, process tree line
    //first, determine the number of TIDs in the line:
    pos = treeline.find_first_of(":");
    line_type = treeline.substr(0, pos);
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(":");
    mycount = treeline.substr(0, pos);
    count = atoi(mycount.c_str());
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(" ");
    if (pos != -1){
      if (!WEIGHTED){
	fprintf(stderr, "\nFile has branch lengths! Weighted option: ON\n\n");	
	WEIGHTED = true;
	hold_integrals = (unsigned int*)malloc(NUM_TREES*sizeof(unsigned int));
	for (unsigned int x = 0; x < NUM_TREES; x++){
	  hold_integrals[x] = 0;
	}
      }
    }
    if (WEIGHTED){
      ids = treeline.substr(0, pos);
      treeline = treeline.substr(pos+1);
      pos = treeline.find_first_of(" ");
      branches_int = treeline.substr(0, pos); //integral portion of branches
      treeline = treeline.substr(pos+1);
      pos = treeline.find_first_of("\n");
      branches_frac = treeline.substr(0, pos); //fractional component of branches
    }
    else{
      pos = treeline.find_first_of("\n");
      ids = treeline.substr(0, pos);
    }

    //process tree ids
    if (count == 0) { 
      ALL_COUNT++;
      hashtable_length--;
      for (unsigned int b = 0; b < NUM_TREES; ++b) { 
	inverted_index[b].push_back(bipart_loc);
      }
      //vec_bs.push_back(bs); //this is a strict consensus bipartition       
      if (WEIGHTED){
	if (branches_int != "")
	  populate_integrals(hold_integrals, branches_int, encode_size);
	my_set_of_ids.resize(NUM_TREES);
	decompress_branch(hold_integrals, my_set_of_ids, all_branches, branches_frac);
      }
      bipart_count++;
    }      
    else { //count != 0 (so we have tree ids to process)
      my_count = decode(ids, found);
      assert(my_count == count);

      if (line_type == "-") { //compressed line
	for (unsigned int i = 0; i < count; ++i) {
	  unsigned int temp = found[i];
	  check[temp] = true;
	}
	unsigned int true_id, sec_id;
	//cout << "we do get here though? num_unique is: " << num_unique << endl;
	for (unsigned int i = 0; i < num_unique; ++i) {
	  if (check[i] == false) {
	    true_id = true_ids[i];
	    my_set_of_ids.push_back(true_id);
	    inverted_index[true_id].push_back(bipart_loc);
	    if (dups[true_id].size() > 0){
	      for (unsigned int j = 0; j< dups[true_id].size(); j++){
		sec_id = dups[true_id][j];
		inverted_index[sec_id].push_back(bipart_loc);
		my_set_of_ids.push_back(sec_id);
	      }
	    }
	    //cout << "end of for stat" << endl;
	  }
	  else
	    check[i] = false;
	}
	//now, take care of of the bipartition associated with this
	sort(my_set_of_ids.begin(), my_set_of_ids.end());
	unsigned int mytotalsize = my_set_of_ids.size();
	if (mytotalsize > NUM_TREES/2){
	  vector<bool> tempcheck(NUM_TREES,0);
	  unsigned int temploc;
	  for (unsigned int j = 0; j < mytotalsize; j++){
	    temploc = my_set_of_ids[j];
	    tempcheck[temploc] = true;
	  }
	  unsigned int mynewsize = NUM_TREES - mytotalsize;
	  hashtable[mycounter] = (unsigned int *)malloc(mynewsize*sizeof(unsigned int));
	  hash_lengths[mycounter] = mynewsize;
	  unsigned int xpos = 0;
	  for (unsigned int j  = 0; j < NUM_TREES; j++){
	    if (tempcheck[temploc] == 0){
	      hashtable[mycounter][xpos] = j;
	      xpos++;
	      row[j]++;
	    }
	  }
	  ALL_COUNT++;
	}
	else {
	  hashtable[mycounter] = (unsigned int *)malloc(mytotalsize*sizeof(unsigned int));
	  hash_lengths[mycounter] = mytotalsize;
	  for (unsigned int j = 0; j < mytotalsize; j++){
	    hashtable[mycounter][j]  = my_set_of_ids[j];
	  }
	}

	if (WEIGHTED){	  
	  if (branches_int != "")
	    populate_integrals(hold_integrals, branches_int, encode_size);
	  decompress_branch(hold_integrals, my_set_of_ids, all_branches, branches_frac);
	}
	mycounter++;
      }
      else { //line is not compressed
	unsigned int true_id, sec_id;
	for (unsigned int i = 0; i < count; ++i) { 
	  unsigned int temp = found[i];
	  true_id = true_ids[temp];
	  inverted_index[true_id].push_back(bipart_loc);
	  my_set_of_ids.push_back(true_id);
	  if (dups[true_id].size() > 0){
	    for (unsigned int j = 0; j < dups[true_id].size(); j++){
	      sec_id = dups[true_id][j];
	      inverted_index[sec_id].push_back(bipart_loc);
	      my_set_of_ids.push_back(sec_id);
	    }
	  }
	} 
	sort(my_set_of_ids.begin(), my_set_of_ids.end());	 
	unsigned int mytotalsize = my_set_of_ids.size();
	if (mytotalsize > 1){
	  hashtable[counter] = (unsigned int *)malloc(mytotalsize*sizeof(unsigned int));
	  hash_lengths[counter] = mytotalsize;
	  for (unsigned int j = 0; j < mytotalsize; j++){
	    hashtable[counter][j]  = my_set_of_ids[j];
	  }
	}
	if (WEIGHTED){
	  if (branches_int != "")
	    populate_integrals(hold_integrals, branches_int, encode_size);
	  decompress_branch(hold_integrals, my_set_of_ids, all_branches, branches_frac);
	}
      }
      bipart_count++;
    } //end if count != 0
    my_set_of_ids.clear();
    bipart_loc++;
    counter++;
  }

  free(found);
  if (bipart_count != counter) { 
    cerr << "ERROR! Bipartitions not processed correctly!" << endl;
    cout << "bipart_count: " << bipart_count << endl;
    cout << "counter: " << counter << endl;
    exit(1);
  }
  assert(lm.size()!=0); 

} 

void parse_trz_hash(string file, vector< vector<unsigned int> > & inverted_index, vector< bool* > & list_bs, vector< vector<float> > &all_branches, LabelMap &lm, unsigned int & num_unique){  
  string mycount, ids, bipartition, line_type, str, taxa, treeline, bitstring, branches_int, branches_frac;
  unsigned int bipart_count, total_bipartitions, ntaxa, nbipart, ndup_lines;
  vector<bool> check;
  vector<unsigned int> true_ids;
  vector< vector<unsigned int> > dups;
  vector <bool> is_dup;
  bipart_count = 0;
  total_bipartitions = 0; 

  ifstream fin(file.c_str(), ios::binary);
  if (!fin) {
    cerr << "cannot open file!\n";
    exit(2);
  }
  
  //read in taxa
  getline(fin, str);  
  int pos = str.find_first_of(" ");
  line_type = str.substr(0, pos);
  if (line_type != "TAXA"){
    cerr << "Error! No taxa labels identified for file! Exiting..\n";
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  stringstream ss(str);
  ntaxa = 0;
  while(getline(ss, taxa, ':')){ 
    lm.push(taxa);
     ntaxa++;
  }
  NUM_TAXA = ntaxa;
  ss.clear();

  //read in number of trees
  getline(fin, str);
  NUM_TREES = get_ntrees(str);

  //read in number of unique trees
  getline(fin, str);
  num_unique = get_unique(str, NUM_TREES);

  //resize data structures
  inverted_index.resize(NUM_TREES);
  all_branches.resize(NUM_TREES);
  check.resize(NUM_TREES);
  true_ids.resize(num_unique);
  dups.resize(NUM_TREES);
  is_dup.resize(NUM_TREES);
  unsigned int * hold_integrals = NULL;
  //read in number of bipartitions to be read
  getline(fin, str);
  parse_and_get(str, "NBIPARTITIONS", nbipart);
  NUMBIPART = nbipart;
  hashtable = (unsigned int **)malloc(NUMBIPART*sizeof(unsigned int *)); // the hash table is NUMBIPART long 
  hash_lengths = (unsigned int *)malloc(NUMBIPART*sizeof(unsigned int)); //allocate corresponding array for this
  for (unsigned int i = 0; i < NUMBIPART; i++)
    hash_lengths[i] = 0;

  //skip nbipart lines, get the duplicate information
  for (unsigned int i = 0; i < nbipart; i++)
    getline(fin, str);
  getline(fin, str); //now read the duplicates line
  parse_and_get(str, "DUPLICATES", ndup_lines);
  //cout << "ndup_lines is: " << ndup_lines << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++)
    is_dup[i] = 0;
  unsigned int decode_size = 0;
  unsigned int * found, dec_loc, dec_val;
  found = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  for (unsigned int i = 0; i < ndup_lines; i++){
    getline(fin, str); //read in the line
    pos = str.find_first_of("\n");
    str = str.substr(0, pos); //get rid of the newline
    decode_size = decode(str, found); //decode the line
    dec_loc = found[0];
    for (unsigned int j = 1; j < decode_size; j++){
      dec_val = found[j];
      dups[dec_loc].push_back(dec_val); //add the elements of the array to the associated dup structure
      is_dup[dec_val] = 1; //also set those locations in our bool structure to a 1
    }
  }

  //populate true_ids
  unsigned int tempj = 0;
  for (unsigned int i = 0; i < num_unique; i++){
    if (is_dup[tempj]){
      while (is_dup[tempj])
	tempj++;
    }
    true_ids[i] = tempj;
    tempj++;
  }

  //debug: print out the data structures
  /*cout << "true_ids: " << endl;
  for (unsigned int i  = 0; i < true_ids.size(); i++)
    cout << i << "-->" << true_ids[i] << endl;
  cout << endl;

  cout << "dups: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    cout << i << ": ";
    for (unsigned int j = 0; j < dups[i].size(); j++){
      cout << dups[i][j] << " ";
    }
    cout << endl;
    }*/
  fin.close(); //now close the file
  fin.open(file.c_str(), ios::binary); //reopen
  if (!fin) {
    cerr << "cannot open file!\n";
    exit(2);
  }
  for (unsigned int i  =0; i < 4; i++) //go back to where we need to be in order to start reading in the bipartitions
    getline(fin, str); 

  unsigned int counter = 0;
  unsigned int count = 0;
  unsigned int bipart_loc = 0;
  unsigned int numtrees_strsize, my_count;
  stringstream tmpss;
  string numtrees_str;
  unsigned int encode_size = 0;
  vector<unsigned int> my_set_of_ids;
  tmpss << NUM_TREES;
  tmpss >> numtrees_str;
  numtrees_strsize = numtrees_str.size();
  encode_size = (numtrees_strsize+1)/2;

  while ( counter < nbipart) { 
    getline(fin, str);  
    pos = str.find_first_of(" ");
    bitstring = str.substr(0, pos); //contains bitstring
    str = str.substr(pos+1);
    pos = str.find_first_of("\n");
    treeline = str.substr(0, pos);

    //process bitstring first
    bool *bs = new bool[NUM_TAXA];	
    decode_bitstring(bitstring, bs, NUM_TAXA);

    //next, process tree line
    //first, determine the number of TIDs in the line:
    pos = treeline.find_first_of(":");
    line_type = treeline.substr(0, pos);
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(":");
    mycount = treeline.substr(0, pos);
    count = atoi(mycount.c_str());
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(" ");
    if (pos != -1){
      if (!WEIGHTED){
	fprintf(stderr, "\nFile has branch lengths! Weighted option: ON\n\n");	
	WEIGHTED = true;
	hold_integrals = (unsigned int*)malloc(NUM_TREES*sizeof(unsigned int));
	for (unsigned int x = 0; x < NUM_TREES; x++){
	  hold_integrals[x] = 0;
	}
      }
    }
    if (WEIGHTED){
      ids = treeline.substr(0, pos);
      treeline = treeline.substr(pos+1);
      pos = treeline.find_first_of(" ");
      branches_int = treeline.substr(0, pos); //integral portion of branches
      treeline = treeline.substr(pos+1);
      pos = treeline.find_first_of("\n");
      branches_frac = treeline.substr(0, pos); //fractional component of branches
    }
    else{
      pos = treeline.find_first_of("\n");
      ids = treeline.substr(0, pos);
    }

    //process tree ids
    if (count == 0) { 
      hashtable[counter] = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
      hash_lengths[counter] = NUM_TREES;
      for (unsigned int b = 0; b < NUM_TREES; ++b) { 
	hashtable[counter][b] = b;
	inverted_index[b].push_back(bipart_loc);
      }
      //vec_bs.push_back(bs); //this is a strict consensus bipartition       
      if (WEIGHTED){
	if (branches_int != "")
	  populate_integrals(hold_integrals, branches_int, encode_size);
	my_set_of_ids.resize(NUM_TREES);
	decompress_branch(hold_integrals, my_set_of_ids, all_branches, branches_frac);
      }
      bipart_count++;
    }      
    else { //count != 0 (so we have tree ids to process)
      my_count = decode(ids, found);
      assert(my_count == count);

      if (line_type == "-") { //compressed line
	for (unsigned int i = 0; i < count; ++i) {
	  unsigned int temp = found[i];
	  check[temp] = true;
	}
	unsigned int true_id, sec_id;
	//cout << "we do get here though? num_unique is: " << num_unique << endl;
	for (unsigned int i = 0; i < num_unique; ++i) {
	  if (check[i] == false) {
	    true_id = true_ids[i];
	    my_set_of_ids.push_back(true_id);
	    inverted_index[true_id].push_back(bipart_loc);
	    if (dups[true_id].size() > 0){
	      for (unsigned int j = 0; j< dups[true_id].size(); j++){
		sec_id = dups[true_id][j];
		inverted_index[sec_id].push_back(bipart_loc);
		my_set_of_ids.push_back(sec_id);
	      }
	    }
	    //cout << "end of for stat" << endl;
	  }
	  else
	    check[i] = false;
	}
	//now, take care of of the bipartition associated with this
	sort(my_set_of_ids.begin(), my_set_of_ids.end());
	unsigned int mytotalsize = my_set_of_ids.size();
	hashtable[counter] = (unsigned int *)malloc(mytotalsize*sizeof(unsigned int));
	hash_lengths[counter] = mytotalsize;
	for (unsigned int j = 0; j < mytotalsize; j++){
	  hashtable[counter][j]  = my_set_of_ids[j];
	}
	if (WEIGHTED){	  
	  if (branches_int != "")
	    populate_integrals(hold_integrals, branches_int, encode_size);
	  decompress_branch(hold_integrals, my_set_of_ids, all_branches, branches_frac);
	}
      }
      else { //line is not compressed
	unsigned int true_id, sec_id;
	for (unsigned int i = 0; i < count; ++i) { 
	  unsigned int temp = found[i];
	  true_id = true_ids[temp];
	  inverted_index[true_id].push_back(bipart_loc);
	  my_set_of_ids.push_back(true_id);
	  if (dups[true_id].size() > 0){
	    for (unsigned int j = 0; j < dups[true_id].size(); j++){
	      sec_id = dups[true_id][j];
	      inverted_index[sec_id].push_back(bipart_loc);
	      my_set_of_ids.push_back(sec_id);
	    }
	  }
	} 
	sort(my_set_of_ids.begin(), my_set_of_ids.end());	 
	unsigned int mytotalsize = my_set_of_ids.size();
	hashtable[counter] = (unsigned int *)malloc(mytotalsize*sizeof(unsigned int));
	hash_lengths[counter] = mytotalsize;
	for (unsigned int j = 0; j < mytotalsize; j++){
	  hashtable[counter][j]  = my_set_of_ids[j];
	}
	if (WEIGHTED){
	  if (branches_int != "")
	    populate_integrals(hold_integrals, branches_int, encode_size);
	  decompress_branch(hold_integrals, my_set_of_ids, all_branches, branches_frac);
	}
      }
      bipart_count++;
    } //end if count != 0
    list_bs.push_back(bs); //push into the list 
    my_set_of_ids.clear();
    bipart_loc++;
    counter++;
  }

  free(found);
  if (bipart_count != counter) { 
    cerr << "ERROR! Bipartitions not processed correctly!" << endl;
    cout << "bipart_count: " << bipart_count << endl;
    cout << "counter: " << counter << endl;
    exit(1);
  }
  assert(lm.size()!=0); 

} 

void sample_hashtable(string trzfile){
  //sample function that will print out a hash table -- indicates the use of the functions
  LabelMap lm; //stores the taxa labels
  vector< vector<unsigned int> > inverted_index; //stores the inverted index
  vector< vector<float> > all_branches; //stores the branches (in terms of inverted index)
  vector< bool * > list_bs; //list of the bipartitions, in bitstring form 
  unsigned int unique; //number of unique trees
  float percent_unique = 0.0;
  parse_trz_hash(trzfile, inverted_index, list_bs, all_branches, lm, unique);
  //print out some info about the dataset
  percent_unique = unique/NUM_TREES;
  cout << "File has " << NUM_TREES << " trees over " << NUM_TAXA << " taxa." << endl;
  cout << "Of those, " << unique << " trees are unique (" << percent_unique << "%)." << endl;
  cout << "File has " << NUMBIPART << " bipartitions." << endl;
  cout << "We show the hash table below, with the corresponding bitstring reps of the bipartitions:" << endl << endl;
  //keep in mind that hashtable and hash_lengths are global variables
  unsigned int mycurrsize = 0;
  for (unsigned int i = 0; i < NUMBIPART; i++){
    for (unsigned int j = 0; j < NUM_TAXA; j++){
      cout << list_bs[i][j];
    }
    cout << " --> [ ";
    mycurrsize = hash_lengths[i];
    for (unsigned int j = 0; j < mycurrsize; j++){
      cout << hashtable[i][j] << " ";
    }
    cout << "]" << endl;
  }
}

void parse_consensus(vector< bool* > & list_bs, vector<float> &list_branches, string file, LabelMap &lm, unsigned char consensus, vector<unsigned int> &bs_sizes){  
  string mycount, ids, bipartition, line_type, str, taxa, treeline, bitstring, branches_int, branches_frac;
  unsigned int bipart_count, total_bipartitions, ntaxa, nbipart;
  bipart_count = 0;
  total_bipartitions = 0; 

  ifstream fin(file.c_str(), ios::binary);
  if (!fin) {
    cerr << "cannot open file!\n";
    return;
  }
  
  getline(fin, str); //read in taxa
  int pos = str.find_first_of(" ");
  line_type = str.substr(0, pos);
  if (line_type != "TAXA"){
    cerr << "Error! No taxa labels identified for file! Exiting..\n";
    return;
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  stringstream ss(str);
  ntaxa = 0;
  while(getline(ss, taxa, ':')){ 
    lm.push(taxa);
     ntaxa++;
  }
  NUM_TAXA = ntaxa;
  ss.clear();

  getline(fin, str); //read in number of trees
  unsigned int unique_t = 0;
  NUM_TREES = get_ntrees(str);
  getline(fin, str); //get number of unique trees
  unique_t = get_unique(str, NUM_TREES);
  cout << "unique_t is: " << unique_t << endl;
  unsigned int * true_ids = NULL;
  unsigned int * dups = NULL;
  bool * is_dup;
  //unsigned int * hold_integrals;
  unsigned int * found = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  bool * check = (bool *)malloc(NUM_TREES*sizeof(bool));
  for (unsigned int i  = 0; i < NUM_TREES; i++)
    check[i] = 0;
  getline(fin, str); //read in number of bipartitions to be read
  parse_and_get(str, "NBIPARTITIONS", nbipart);
  NUMBIPART = nbipart;
  cout << "NUMBIPART is: " << NUMBIPART << endl;
  if (consensus == 'm'){
    for (unsigned int i  = 0; i < NUMBIPART; i++)
      getline(fin, str);

    getline(fin, str);
    unsigned int ndups = 0;
    parse_and_get(str, "DUPLICATES", ndups);
    true_ids = (unsigned int *)malloc(unique_t*sizeof(unsigned int));
    dups = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
    is_dup = (bool *)malloc(NUM_TREES*sizeof(unsigned int));

    for (unsigned int i  = 0; i < unique_t; i++)
      true_ids[i] = 0;

    for (unsigned int i  = 0; i < NUM_TREES; i++)
      dups[i] = 0;

    for (unsigned int i  =0; i < NUM_TREES; i++)
      is_dup[i] = 0;
    //process duplicates
 
    unsigned int dec_loc, decode_size, dec_val;
    for (unsigned int i = 0; i < ndups; i++){
      getline(fin, str); //read in the line
      pos = str.find_first_of("\n");
      str = str.substr(0, pos); //get rid of the newline
      decode_size = decode(str, found); //decode the line
      dec_loc = found[0];
      dups[dec_loc] = decode_size-1;
      for (unsigned int j = 1; j < decode_size; j++){
	dec_val = found[j];
	is_dup[dec_val] = 1; //also set those locations in our bool structure to a 1
      }                 
    }        

    
    //populate true_ids
    unsigned int tempj = 0;
    for (unsigned int i = 0; i < unique_t; i++){
      if (is_dup[tempj]){
	while (is_dup[tempj])
	  tempj++;
      }
      true_ids[i] = tempj;
      tempj++;
    }
   
    /*cout << "True ids: " << endl;
    for (unsigned int i = 0; i < unique_t; i++)
      cout << i << "--> " << true_ids[i] << endl;

    for (unsigned int i  =0; i < NUM_TREES; i++)
      cout << i <<  " = " << dups[i] << " duplicates" << endl; 
      exit(0);*/

    fin.close();
    fin.open(file.c_str());
    if (!fin) {
      cerr << "cannot open file!\n";
      return;
    }
    for (unsigned int i  =0; i < 4; i++) //now, get back to where we were
      getline(fin, str);
  }
  unsigned int counter = 0;
  unsigned int count = 0;
  //unsigned int num_trees, newcount, myplace, numer, denom, numtrees_strsize, val;
  //unsigned char v;
  float weight = 0;
  stringstream tmpss;
  string numtrees_str;
  unsigned int encode_size = 0;
  unsigned int numtrees_strsize = 0;
  unsigned int decode_size = 0;
  unsigned int threshold = (NUM_TREES+1)/2;
  unsigned int maxLength = 0;
  tmpss << NUM_TREES;
  tmpss >> numtrees_str;
  numtrees_strsize = numtrees_str.size();
  encode_size = (numtrees_strsize+1)/2;
  while ( counter < nbipart) { 
    getline(fin, str);  

    pos = str.find_first_of(" ");
    bitstring = str.substr(0, pos); //contains bitstring
    str = str.substr(pos+1);
    pos = str.find_first_of("\n");
    treeline = str.substr(0, pos);

    //process bitstring first
    maxLength = get_bitstring_length(bitstring);//determine length of bitstring maxLength
    bool *bs = new bool[maxLength]; //allocate it to be maxLength	
    decode_bitstring(bitstring, bs, maxLength);

    //next, process tree line
    //first, determine the number of TIDs in the line:
    string tids;
    pos = treeline.find_first_of(":");
    line_type = treeline.substr(0, pos);
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(":");
    mycount = treeline.substr(0, pos);
    count = atoi(mycount.c_str());
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(" ");
    tids = treeline.substr(0, pos);
    treeline = treeline.substr(pos+1);

    if (count == 0) { //if it is in all the trees
      list_bs.push_back(bs); //add  
      bs_sizes.push_back(maxLength);
      list_branches.push_back(1.0); //add support of 1.0
    }      
    else { //count != 0 (so we have tree ids to process)
      if (consensus == 'm') { //other lines are only processed when we we have majority consensus
	//decompress set of tree ids
	decode_size = decode(tids, found); //decode the line
	unsigned int trueid, myid;
	if (line_type == "+"){
	  for (unsigned int i = 0; i < decode_size; i++){
	    myid = found[i];
	    trueid = true_ids[myid];
	    count += dups[trueid];
	  } 
	}
	else{
	  //reverse the ids
	  for (unsigned int i = 0; i < decode_size; i++){
	    myid = found[i];
	    check[myid] = 1;
	  }
	  count  = unique_t - count;
	  for (unsigned int i = 0; i < unique_t; i++){
	    if (!check[i]){
	      trueid = true_ids[i];
	      count += dups[trueid];
	    }
	    else
	      check[i] = 0;
	  }
	  if (count >= threshold){
	    list_bs.push_back(bs);	
	    bs_sizes.push_back(maxLength);
	    weight = (float)(count)/NUM_TREES;
	    list_branches.push_back(weight);
	  }
	  //exit(0);
	} //end else line is compressed
      } //end consensus is m
    } //end else if count is greater than 0
    counter++;
    bipart_count++;
  }
  
  if (bipart_count != counter) { 
    cerr << "ERROR! Bipartitions not processed correctly!" << endl;
    cout << "bipart_count: " << bipart_count << endl;
    cout << "counter: " << counter << endl;
    exit(1);
  }
  assert(lm.size()!=0); 
} 

vector<unsigned int> unique(vector< vector<unsigned int> > inverted_index, 
		    vector< vector<float> >all_branches,
		    bool branch)
{
  fprintf(stderr, "Detecting unique trees...\n");
  
  HashRFMap vvec_unique;
  vector<unsigned int> unique_ids;
  unsigned int unique_count = 0;

  vvec_unique.uhashfunc_init_unique(NUM_TREES, NUM_TAXA, NUMBIPART, C, NEWSEED);
  cout << "NUM_TREES is: " << NUM_TREES << endl;
  cout << "NUMBIPART is: " << NUMBIPART << endl;

  //new m1 and m2 values
  unsigned long long um1 = vvec_unique._HF.getM1();
  unsigned long long um2 = vvec_unique._HF.getM2();
  
  vvec_unique._hashtab2.resize(um1); //resize vvec_unique

  unsigned int locsize = 0;

  for (unsigned int i = 0; i < NUM_TREES; i++){ 
    unsigned int pos = 0;
    locsize = inverted_index[i].size();
    unsigned long long temp1 = 0;
    unsigned long long temp2 = 0;
    bool * cbs = new bool[NUMBIPART];
    for (unsigned int j = 0; j < NUMBIPART; j++)
      cbs[j] = 0;
    for (unsigned int j = 0; j < locsize; j++){
      pos = inverted_index[i][j];
      cbs[pos] = 1;
      temp1+= vvec_unique._HF.getA1(pos);
      temp2+= vvec_unique._HF.getA2(pos);
    }
    temp1 = temp1 % um1;
    temp2 = temp2 % um2;
    vvec_unique.hashing_bs_unique(i, NUMBIPART, temp1, temp2, cbs);
  }

  if (!branch){
    for (unsigned int hti=0; hti<vvec_unique._hashtab2.size(); ++hti) {
      unsigned int sizeVec = vvec_unique._hashtab2[hti].size(); 
      if (sizeVec) {
	for (unsigned int i=0; i<sizeVec; ++i) {
	  unique_count++;	
	  unique_ids.push_back(vvec_unique._hashtab2[hti][i]._vec_treeidx[0]);
	}
      } //end if sizevec
    } //end hashtable traversal
    sort(unique_ids.begin(), unique_ids.end());
    fprintf(stderr, "Found %u unique trees\n", unique_count);
    return unique_ids;
  }
  else{
    bool found = false;
    vector<bool> candidates(NUM_TREES, false);
    for (unsigned int hti=0; hti<vvec_unique._hashtab2.size(); ++hti) {
      unsigned int sizeVec = vvec_unique._hashtab2[hti].size(); 
      if (sizeVec) {
	for (unsigned int i=0; i<sizeVec; ++i) {
	  unsigned int num_candidates = vvec_unique._hashtab2[hti][i]._vec_treeidx.size();
	  if (num_candidates > 1){
	    for (unsigned int j = 0; j < num_candidates; j++){
	      for (unsigned int k  = j+1; k < num_candidates; k++){
		for (unsigned int l = 0; l < NUMBIPART; l++){
		  if (all_branches[j][l] != all_branches[k][l])
		    break;
		  found = true; //it is a duplicate, it can be ignored
		}
		if (found){
		  candidates[k] = true; //update array location k to true, indicating its a duplicate
		}
	      }
	    }
	    for (unsigned int j = 0; j < num_candidates; j++){
	      if (candidates[j] ==false){
		unique_count++;
		unique_ids.push_back(vvec_unique._hashtab2[hti][i]._vec_treeidx[j]);
	      }
	      else{
		candidates[j] = false;
	      }
	    }
	  }
	  else{
	    unique_count++;	
	    unique_ids.push_back(vvec_unique._hashtab2[hti][i]._vec_treeidx[0]);
	  }
	}
      } //end if sizevec
    } //end hashtable traversal
    fprintf(stderr, "Found %u unique trees\n", unique_count);
    return unique_ids;
  }

}

void decompress(string file, int random_decompress, bool branch, bool calc_unique_trees, unsigned char consensus, bool quiet) { 
  LabelMap lm;
  vector< vector<bool* > > all_bs; //this stores all the bipartitions that are not in the "strict tree"
  vector< vector<unsigned int> > inverted_index, bs_sizes;
  vector<bool* > list_bs; 
  vector<float> list_branches;
  vector< vector<float> > all_branches; //this stores all the branches
  //vector<vector<SCNode*> > vvec_distinctClusters2;	
  vector<bool> check;
  string outfile;
  vector<unsigned int> my_unique;

  struct timeval parsing_start;
  struct timeval parsing_end;
  struct timeval building_start;
  struct timeval building_end;
  struct timeval compute_start;
  struct timeval compute_end;
  double compute_time = 0;
  bool calculate_consensus = false;
  gettimeofday(&parsing_start, NULL);
  string ext = file.substr(file.size()-4);
  if (ext != ".trz") {
    cerr << "Invalid extension. Please input a valid .trz file\n";
    exit(2);
  }
  if (is_nexus(file)){
    //nexus_decompress(file, random_decompress, branch, calc_unique_trees, consensus, quiet);
    decompress_nexus(file, random_decompress, branch, calc_unique_trees, consensus, quiet);
    return;
  }
  outfile = file.substr(0, file.size()-4);  
  outfile = outfile + ".e"; //temporary measure so original file does not get changed --omit when code is released

  if (consensus == 's' || consensus == 'm'){
    calculate_consensus = true;
  }
  if (calculate_consensus){
    vector<unsigned int> con_bs_sizes;
    parse_consensus(list_bs, list_branches, file, lm, consensus, con_bs_sizes);
    string consensus_tree;
    bool *bs2 = new bool[NUM_TAXA];
    for (unsigned int i =0; i < NUM_TAXA; i++)
      bs2[i] = 1;
    vector<bool* >::iterator bs_itr = list_bs.begin();
    list_bs.insert(bs_itr, bs2);
    vector<float>::iterator br_itr = list_branches.begin();
    list_branches.insert(br_itr, 0.0);
    vector<unsigned int>::iterator bz_itr = con_bs_sizes.begin();
    con_bs_sizes.insert(bz_itr, NUM_TAXA);
    if (consensus == 's')
      consensus_tree = compute_tree(lm, list_bs, list_branches, 0, false, con_bs_sizes);
    else{
      WEIGHTED = true;
      consensus_tree = compute_tree(lm, list_bs, list_branches, 0, true, con_bs_sizes);
    }
 
    ofstream fout;
    fout.open(outfile.c_str());
    if (!fout) { 
      cerr << "Error writing to file " << outfile << endl;
      return;
    }
    fout << consensus_tree << endl;
    fout.close();
    if (!quiet)
      fprintf(stderr, "Done calculating consensus. exiting...\n");
    return;
  }
  vector<bool> is_dup;
  parse_trz_file(all_bs, inverted_index, list_bs, all_branches, file, lm, is_dup, bs_sizes);
   gettimeofday(&parsing_end, NULL);
  double parsing_time = parsing_end.tv_sec - parsing_start.tv_sec + (parsing_end.tv_usec - parsing_start.tv_usec) / 1.e6;
  if (!quiet)
    fprintf(stderr,"\nFile Parsing time: %g seconds \n", parsing_time);
  
  /*cout << "printing out all branches data structure!" << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    mysize = all_branches[i].size();
    for (unsigned int j = 0; j < mysize; j++){
      cout << all_branches[i][j] << " ";
    }
    cout << endl;
    }*/
  //exit(0);
  /* unsigned int row_size;
 for (unsigned int i = 0; i < NUM_TREES; i++){
    row_size = inverted_index[i].size();
    for (unsigned int j =0; j < row_size; j++){
      cout << inverted_index[i][j] << " ";
    }
    cout << endl;
    }*/
  gettimeofday(&building_start, NULL);
  if (!quiet){
    fprintf(stderr, "Done processing compressed file. Inflating tree structures..\n");
    fprintf(stderr, "Number of taxa: %d\n",NUM_TAXA);
    fprintf(stderr, "Number of trees: %d\n",NUM_TREES);
  }
  
  //calculate the star bipartition
  vector<bool* >::iterator bs_itr;
  vector<unsigned int>::iterator bsz_itr;
  vector<float>::iterator br_itr;
  if (HETERO){
    cout << "NUM_TAXA is: **  " << NUM_TAXA << endl;
    bool *bs2 = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TREES; i++){
      for (unsigned int j = 0; j < NUM_TAXA; j++)
	bs2[j] = 0; //set to 0
      unsigned int bs_nbiparts = all_bs[i].size();
      for (unsigned int j =0; j < bs_nbiparts; j++){
	unsigned int bs_size = bs_sizes[i][j];
	for (unsigned int k = 0; k < bs_size; k++)
	  bs2[k] |= all_bs[i][j][k]; //OR operation
      }
      bool * mystar = new bool[NUM_TAXA]; //this will be inefficient for large, varied n.
      for (unsigned int j = 0; j < NUM_TAXA; j++)
	mystar[j] = bs2[j];
      bs_itr = all_bs[i].begin();
      bsz_itr = bs_sizes[i].begin();
      br_itr = all_branches[i].begin();
      all_bs[i].insert(bs_itr, mystar);
      bs_sizes[i].insert(bsz_itr, NUM_TAXA); //this will also need to be updated!
      all_branches[i].insert(br_itr, 0.0); // add a 0 as the value for the star bipart
    }
    delete [] bs2;
  }
  else{
    bool * mystar = new bool[NUM_TAXA]; 
    for (unsigned int j = 0; j < NUM_TAXA; j++)
      mystar[j] = 1; //for homogeneous collections, the star bipart is the same for all
    for (unsigned int i = 0; i < NUM_TREES; i++){
      bs_itr = all_bs[i].begin();
      bsz_itr = bs_sizes[i].begin();
      br_itr = all_branches[i].begin();
      all_bs[i].insert(bs_itr, mystar);
      bs_sizes[i].insert(bsz_itr, NUM_TAXA); //this will also need to be updated!
      all_branches[i].insert(br_itr, 0.0); // add a 0 as the value for the star bipart
    }
  }
  //exit(0); 

  //recreate the trees
  string tree;
  if (!quiet)
    fprintf(stderr, "Decompressed file located in: %s\n", outfile.c_str());
  ofstream fout;
  fout.open(outfile.c_str());
  if (!fout) { 
    cerr << "Error writing to file " << outfile << endl;
    return;
  }

  if (!random_decompress){
    if (calc_unique_trees){
      //my_unique = unique(inverted_index, all_branches, branch);
      gettimeofday(&compute_start, NULL); 
      for (unsigned int x = 0; x < NUM_TREES; x++){
	if (!is_dup[x]){
	  tree = compute_tree(lm, all_bs[x], all_branches[x], x, branch, bs_sizes[x]);
	  fout << tree << endl;
	}
      }
      gettimeofday(&compute_end, NULL);  
      compute_time += compute_end.tv_sec - compute_start.tv_sec + (compute_end.tv_usec - compute_start.tv_usec) / 1.e6;
    }
    else{ //if not calc unique trees
	for (unsigned int x = 0; x < NUM_TREES; ++x) {
	  gettimeofday(&compute_start, NULL); 
	  tree = compute_tree(lm, all_bs[x], all_branches[x], x, branch, bs_sizes[x]);
	  gettimeofday(&compute_end, NULL); 
	  compute_time += compute_end.tv_sec - compute_start.tv_sec + (compute_end.tv_usec - compute_start.tv_usec) / 1.e6;
	  fout << tree << endl;
	}
    } //end else not calc unique trees
  } //end if not random decompress
  else{ //if random decompress
    //cout << "random decompress:" << endl;
    check.resize(NUM_TREES, 0);
    unsigned int to_generate = NUM_TREES*random_decompress/100; //how many
    RandomLib::Random rnd;          // r created with random seed
    unsigned int tmp = NUM_TREES;
    unsigned int mynum, gencount;
    gencount = 0;
    cout << "generating " << to_generate << "random indices.." << endl;
    while(gencount < to_generate) {
      mynum = rnd.Integer<unsigned int>(tmp-1);
      if (!check[mynum]){
	gencount++;
	check[mynum] = true;
      }
    }                       
    cout << "done." << endl;
    cout << "rebuilding trees..\n" << endl;
    for (unsigned int x = 0; x < NUM_TREES; ++x) {
	if (check[x]){
	  gettimeofday(&compute_start, NULL); 
	  tree = compute_tree(lm, all_bs[x], all_branches[x], x, branch, bs_sizes[x]);
	  gettimeofday(&compute_end, NULL);  
	  compute_time += compute_end.tv_sec - compute_start.tv_sec + (compute_end.tv_usec - compute_start.tv_usec) / 1.e6;
	  fout << tree << endl;
	}
    }
  } //end else random decompress

  fout.close();
  gettimeofday(&building_end, NULL);
  double building_time = building_end.tv_sec - building_start.tv_sec + (building_end.tv_usec - building_start.tv_usec) / 1.e6;
  if (!quiet){
    fprintf(stderr,"\nTree Reconstruction time: %g seconds \n",building_time);
    fprintf(stderr,"\ncompute_tree() time: %g seconds \n",compute_time);
    fprintf(stderr,"     In compute_tree(): \n");
    fprintf(stderr,"\nUpdating Data Structures time: %g seconds \n",update_time);
    fprintf(stderr,"\nTree Building time: %g seconds \n", build_time);
    fprintf(stderr,"\nCleaning up time: %g seconds \n",clean_time);
    fprintf(stderr,"\nTotal time in function compute_tree(): %g seconds \n",total_time);
  }
}

void clear_hashtable(){
  NUM_TAXA = 0;
  NUM_TREES = 0;
  NUMBIPART = 0;
  for (unsigned int i = 0; i < vvec_hashrf._hashtab2.size(); i++){
    if (vvec_hashrf._hashtab2[i].size()){
      for (unsigned int j  = 0; j < vvec_hashrf._hashtab2[i].size(); j++){
	vvec_hashrf._hashtab2[i][j]._vec_treeidx.clear();
	vvec_hashrf._hashtab2[i][j]._vec_dist.clear();
	if (vvec_hashrf._hashtab2[i][j]._bs != NULL)
	  delete [] vvec_hashrf._hashtab2[i][j]._bs;
      }
    }
   vvec_hashrf. _hashtab2[i].clear();
  }
  vvec_hashrf._hashtab2.clear();
}
