#include <iostream>
#include <stdlib.h>
#include "global.h"
#include "bmerge.h"
#include "tmerge.h"
#include "label-map.hh"
#include <limits>
#include <sstream>
#include <algorithm>

char mycode[11] = "ABCDEFGHIJ";

bool sort_b(const pair<unsigned int, pair<string, unsigned int> > & a, const pair<unsigned int, pair<string, unsigned int> > &b) { 
    if (a.first == b.first)
      return ( (a.second).first > (b.second).first);
    else
      return a.first > b.first;
}

void test_equivalence(string infile, string mergefile){
  fprintf(stderr, "This function simply outputs if the two tree collections are equivalent.\n");
  ifstream fin, fin2;
  fin.open(infile.c_str());
  if (!fin){
    cerr << "cannot open file " << infile << " for reading!" << endl;
    exit(2);
  }
  fin2.open(mergefile.c_str());
  if (!fin2){
    cerr << "cannot open file " << mergefile << " for reading!" << endl;
    exit(2);
  }
  string line1, line2;
  unsigned int count = 1;
  while (!fin.eof() && !fin2.eof()){
    if (fin.eof() || fin2.eof())
      break;
    getline(fin, line1);
    getline(fin2, line2);
    if (line1 != line2){
      cout << "Files are not equivalent!" << endl;
      cout << "First difference found at line: " << count << endl;
      cout << line1 << " vs. " << line2 << endl;
      return;
    }
    count++;
  }
  cout << "Files are equivalent!" << endl;
}

void convert_bitstring(string bitstring, string & converted_bitstring){
  short bitstring_size = bitstring.size();
  unsigned int j = 0;
  int x = 0;
  int ypos = 0;                                                                                                            
  int y = 0;            
  converted_bitstring = ""; 
  stringstream ss;
  while ( x < bitstring_size ) {                                                                                           
    int val = bitstring[x];                                                                                                
    val-= 65;                                                                                                              
    if (val < 2){ 
      ss << val;
    }
    else{
      char type = bitstring[x];                                                                                            
      if (type == 'K')                                            
	j = 1;                                                                                                      
      if (type == 'L')                                                                                                     
	j = 0;
      ypos = x + 1;                                                                                                        
      val = -1;                                                                                                            
      y = 0;                                                                                                               
      while (val < 0) {                                                                                                    
	val = bitstring[ypos];                                                                                             
	val-=65;                                                                                                           
	if (val < 0) {                                                                                                     
	  ++ypos;                                                                                                          
	  ++y;                                                                                                             
	}                                                                                                                  
	else{                                                                                                              
	  ypos--;                                                                                                          
	}                                                                                                                  
	if (ypos == bitstring_size)  
	  break;                                                                                                           
      }                                                                                                                    
      assert(ypos > -1);                                                                                                   
      
      string amt = bitstring.substr((x+1),y);                                                                              
      x = ypos;                                                                                                            
      int iamt = atoi(amt.c_str());                                                                                        
      for (short z = 0; z < iamt; ++z) {                                                                                   
	ss << j;    
      }                                                                                                                    
    }                                                                                                                      
    ++x;                                                                                                                   
  }
  converted_bitstring = ss.str();
}

string compress_line(vector<unsigned int> tids, unsigned int nTrees){
  //takes a set of tree ids and compresses them
  unsigned int size, temp;
  string compressed;
  stringstream ss;
  if (tids.size() == nTrees){
    compressed = "-:0:";
    return compressed;
  }
  size = 0;
  unsigned int * new_vec = (unsigned int *)malloc(nTrees * sizeof(unsigned int));
  for (unsigned int i  =0; i < nTrees; i++)
    new_vec[i] = 0;

  if (tids.size() > nTrees/2){
    //need to reverse tree ids
    bool * check = (bool *)malloc(nTrees * sizeof(bool));
    for (unsigned int i =0; i < nTrees; i++)
      check[i] = 0;
    size = nTrees - tids.size();
    ss << "-:" << size << ":";
    for (unsigned int i  =0; i < tids.size(); i++){
      temp = tids[i];
      check[temp] = 1;
    }
    temp = 0;
    for (unsigned int i = 0; i < nTrees; i++){
      if (check[i] == 0){
	new_vec[temp] = i;
	temp++;
      }
    }
    free(check);
  }
  else{
    size = tids.size();
    ss << "+:" << size << ":";
    for (unsigned int i = 0; i < size; i++)
      new_vec[i] = tids[i];
  }
  compressed = encode(new_vec, size);
  ss << compressed;
  compressed = ss.str();
  ss.str("");
  free(new_vec);
  return compressed;
}

void  merge_bipartitions(string infilename, string mergefile, vector<string> & all_bipart){
  ifstream fin;
  string str, type;
  int pos;
  unsigned int nbiparts;
  fin.open(infilename.c_str());
  if (!fin){
    cerr << "cannot open file " << infilename << " for reading!" << endl;
    exit(2);
  }

  for (unsigned int i  =0; i < 4; i++)
    getline(fin, str); //skip first three lines, and read the fourth
  parse_and_get(str, "NBIPARTITIONS", nbiparts);

  for (unsigned int i =0; i < nbiparts; i++){
    getline(fin, str);
    pos = str.find_first_of(" ");
    type  = str.substr(0,pos);
    all_bipart.push_back(type);
  }
  fin.close(); //close the first file
  //repeat the first part of the above procedure with the second file
  fin.open(mergefile.c_str());
  if (!fin){
    cerr << "cannot open file " << mergefile << " for reading!" << endl;
    exit(2);
  }
  getline(fin,str); //skip taxa line
  getline(fin, str); //skip num_trees line
  getline(fin, str); //skip unique line
  getline(fin, str); //read in bipartition line
  pos = str.find_first_of(" ");
  type = str.substr(0,pos);
  if (type != "NBIPARTITIONS"){
    cerr << "cannot find bipartition line in " << infilename << "!" <<  endl;
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  nbiparts = atoi(str.c_str());
  for (unsigned int i =0; i < nbiparts; i++){
    getline(fin, str);
    pos = str.find_first_of(" ");
    type  = str.substr(0, pos);
    unsigned int j =0;
    while (j != all_bipart.size()){
      if (all_bipart[j] == type)
	break;
      j++;
    }
    if (j == all_bipart.size())
      all_bipart.push_back(type);
  }

  //let's print out the master list of bipartitions as a way to check!
  //fprintf(stderr, "Found %u total bipartitions!\n", (unsigned int)all_bipart.size());
  NUMBIPART = all_bipart.size();
  fin.close();
  //if (HETERO)
  // sort(all_bipart.begin(), all_bipart.end(), sort_by_ones);
}

string convert_id(string temp_branch, unsigned int x, unsigned int encode_size){
  unsigned int strpad, tmpad;
  string newtemp, myprefix, tmplen, result;
  stringstream stm;
  int pos;
  strpad = 0;

  pos = temp_branch.find_first_of(" ");
  result = temp_branch.substr(0, pos);
  if (result == "")
    return temp_branch;

  stm << x;
  stm >> myprefix;
  temp_branch = temp_branch.substr(encode_size); //get rid of encoded_id part
  if (myprefix.size() < (encode_size*2))
    strpad = (encode_size*2)-myprefix.size();
  while (strpad > 0){
    myprefix = "0" + myprefix; 
    strpad--;
  }
  stringstream tempstr;
  if (encode_size > 1){
    strpad = 0;
    while ((strpad < encode_size)){
      tmplen = myprefix.substr(0,2);
      tmpad = atoi(tmplen.c_str());
      tmpad+=33;
      tempstr << char(tmpad);
      strpad++;
      myprefix = myprefix.substr(2);
    }
  }
  else{
    tmpad = atoi(myprefix.c_str());
    tmpad+=33;
    tempstr << char(tmpad);
  }  
  tempstr << temp_branch;
  result = tempstr.str();
  return result;
}

bool null_case(string infilename, string mergefile, string outfile, int option){
  ifstream fin;
  ifstream fin2;
  fin.open(infilename.c_str());
  if (!fin){
    cerr << "cannot open file for reading!";
    exit(2);
  }
  fin2.open(mergefile.c_str());
  if (!fin2){
    cerr << "cannot open merge file for reading!";
    exit(2);
  }
  unsigned int size1, size2, count;
  fin.seekg(0, ios::end);
  size1 = fin.tellg();
  fin2.seekg(0, ios::end);
  size2 = fin2.tellg();
  fin.close();
  fin2.close();
  if (size1 && size2){
    ifstream fin, fin2;
    string line, line2;
    unsigned int t1, t2;
    fin.open(infilename.c_str());
    fin2.open(mergefile.c_str());
    for (unsigned int i = 0; i < 2; i++)
      getline(fin, line);
    for (unsigned int i  = 0; i < 2; i++)
      getline(fin2, line2);
    t1 = get_ntrees(line);
    t2 = get_ntrees(line2);
    fin.close();
    fin2.close();
    if (t1 && t2)
      return false;
    if (!t1 && t2){
      ofstream fout;
      string myline;
      count = 0;
      fin2.open(mergefile.c_str());
      if (outfile == "")
	outfile = "merged.trz";
      fout.open(outfile.c_str());
      getline(fin2, myline);
      fout << myline << endl;
      if (option == 1){
	while (!fin2.eof()){
	  getline(fin2, myline);
	  if (fin2.eof())
	    break;
	  fout << myline << endl;
	  count++;
	  if (count == 1){
	    unsigned int nt;
	    nt = get_ntrees(myline);
	    cout << "= Total Number of Trees: " << nt << endl; 
	  }
	}
      }
      else{
	cout << "= Total Number of Trees: " << t2 << endl; 
	fout << "NTREES 0";
	if (HETERO)
	  fout << " H";
	fout << endl;
	fout << "UNIQUE_T 0" << endl;
	fout << "NBIPARTITIONS 0" << endl;
	fout << "DUPLICATES 0" << endl;
      } 
      fout.close();
      fin2.close();
      return true;
    }
    if (!t2 && t1){
      ofstream fout;
      string myline;
      fin.open(infilename.c_str());
      if (outfile == "")
	outfile = "merged.trz";
      fout.open(outfile.c_str());
      count = 0;
      getline (fin, myline);
      fout << myline << endl;
      if ((option == 1)||(option == 3)){
	while (!fin.eof()){
	  getline(fin, myline);
	  if (fin.eof())
	    break;
	  fout << myline << endl;
	  count++;
	  if (count == 1){
	    unsigned int nt = get_ntrees(myline);
	    cout << "= Total Number of Trees: " << nt << endl; 
	  }
	}
      }
      else{
	cout << "= Total Number of Trees: " << t1 << endl; 
	fout << "NTREES 0";
	if (HETERO)
	  fout << " H";
	fout << endl;
	fout << "UNIQUE_T 0" << endl;
	fout << "NBIPARTITIONS 0" << endl;
	fout << "DUPLICATES 0" << endl;
      }
      fout.close();
      fin.close();
      return true;
    }
    return true;
  }
  else{
    if (!size1){ //if first file is empty
      ofstream fout;
      count = 0;
      fin2.open(mergefile.c_str());
      if (outfile == "")
	outfile = "merged.trz";
      fout.open(outfile.c_str());
      if (option == 1){
	string myline;
	while (!fin2.eof()){
	  getline(fin2, myline);
	  if (fin2.eof())
	    break;
	  fout << myline << endl;
	  count++;
	  if (count == 2){
	    unsigned int nt = get_ntrees(myline);
	    cout << "= Total Number of Trees: " << nt << endl; 
	  }
	}
      } 
      fout.close();
      fin2.close();
    }
    if (!size2){
      ofstream fout;
      fin.open(infilename.c_str());
      if (outfile == "")
	outfile = "merged.trz";
      fout.open(outfile.c_str());
      count = 0;
      if ((option == 1)||(option == 3)){
	string myline;
	while (!fin.eof()){
	  getline(fin, myline);
	  if (fin.eof())
	    break;
	  fout << myline << endl;
	  count++;
	  if (count == 2){
	    unsigned int nt = get_ntrees(myline);
	    cout << "= Total Number of Trees: " << nt << endl; 
	  }
	}
      }
      fout.close();
    }	
    fin.close();
    fin2.close();
    return true;
  }
}
void setup_tmerge(string infilename, string mergefile, string outfile, int option, bool quiet){
  unsigned int offset;
  LabelMap lm;
  string tempfile = infilename;
  //cout << "tempfile is: " << tempfile << endl;
  //make sure we can validly merge these two trz files
  unsigned int total_trees = 0;
  offset = validate_merge(infilename, mergefile, lm, option, total_trees);
  vector<string> all_bipart;
  vector < vector<string> > all_branch;
  all_branch.resize(total_trees); //this is for the weighted case
  if (!quiet)
    cout << "= Total Number of Trees: " << total_trees << endl; 
  unsigned int maxbranch = NUM_TAXA+(NUM_TAXA-2);
  for (unsigned int i = 0; i < total_trees; i++)
   all_branch[i].resize(maxbranch);
  vector<unsigned int> branch_sizes(total_trees, 0);
   unsigned int max_bipart = 2*NUM_TAXA - 1;
  unsigned int ** inverted_index = (unsigned int **)malloc(NUM_TREES*sizeof(unsigned int *));
  unsigned int * inverted_index_sizes = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  unsigned int t1, t2; //number of trees in file 1 and 2 respectively
  bool * is_dup = (bool *)malloc(total_trees*sizeof(bool));
  for (unsigned int i = 0; i < total_trees; i++)
    is_dup[i] = 0;

  unsigned int * found = (unsigned int *)malloc(total_trees*sizeof(unsigned int));  
  vector<string> dups;
  for (unsigned int i = 0; i  < NUM_TREES; i++){
    inverted_index[i] = (unsigned int *)malloc(max_bipart * sizeof(unsigned int));
    inverted_index_sizes[i] = 0;
  }

 
  int * true_ids = (int *)malloc(NUM_TREES*sizeof(int));
  if (!quiet)
    fprintf(stderr, "step 1: merging bipartitions.\n");
  merge_bipartitions(infilename, mergefile, all_bipart); //merge bipartitions
  getdups(infilename, t1, dups); //now, get the duplicates
  dups.push_back("next");
  getdups(mergefile, t2, dups);

  //process duplicates
  if (!quiet)
    fprintf(stderr, "step 2: processing duplicates.\n");
  string mydup;
  bool useoffset = false;
  unsigned int decode_size, dec_loc, dec_val;
  //need to edit this
  vector< vector<unsigned int> > alldups;
  alldups.resize(total_trees);
 
  for (unsigned int i = 0; i < dups.size(); i++){
    mydup = dups[i];
    if (mydup == "next"){
      useoffset = true;
      continue;
    }
    decode_size = decode(mydup, found); //decode the line                       
    dec_loc = found[0];
    if (useoffset)
      dec_loc += t1;
    for (unsigned int j = 1; j < decode_size; j++){                           
      dec_val = found[j];
      if (useoffset)
	dec_val+= t1;
      alldups[dec_loc].push_back(dec_val);
      is_dup[dec_val] = 1; //also set those locations in our bool structure to a 1                                
    }  
  }   
  infilename = tempfile;
  //cout << "filename is: " << infilename << endl;
 
  
  /* for (unsigned int i  =0; i < total_trees; i++){
    cout << i << ": ";
    for(unsigned int j = 0; j < alldups[i].size(); j++){
      cout << alldups[i][j] << " ";
    }
    cout << endl;
 }
 exit(1);*/
  int tempj = 0;                                              
  for (unsigned int i = 0; i < total_trees; i++){  
    if (tempj >= total_trees)
      break;
    if (is_dup[tempj]){
      while (is_dup[tempj])                                                     
        tempj++;                                                         
    }           
    if (tempj < total_trees)
      true_ids[i] = tempj;                                                        
    tempj++;                                                                    
  }      
  free(is_dup);
  
  //debug
  /* cout << "is_dup: [ ";
  for (unsigned int i = 0; i < total_trees; i++)
    cout <<  is_dup[i] << " ";
    cout << "]" << endl;
  cout << "true_ids: " << endl;                                               
  for (unsigned int i  = 0; i < NUM_TREES; i++)                           
    cout << i << "-->" << true_ids[i] << endl;                                  
  cout << endl;                                                                
                                                                                
  cout << "duplicates: " << endl;                                                     
  for (unsigned int i = 0; i < total_trees; i++){                                 
    cout << i << ": ";                                                          
    for (unsigned int j = 0; j < alldups[i].size(); j++){                          
      cout << alldups[i][j] << " ";                                                
    }                                                                           
    cout << endl;                                                               
    }
  exit(1);
  */

  if (!quiet)
    fprintf(stderr, "step 3: parsing and loading each file.\n");
  unsigned int encode_old, encode_new, numtrees_strsize, true_offset;
  stringstream tmpss;
  string numtrees_str;
  unsigned int total_ntrees;
  total_ntrees = t1+t2;
  //cout << "total_ntrees is: " << total_ntrees << endl;
  encode_old = 0;
  encode_new = 0; 
  tmpss << t1;
  tmpss >> numtrees_str;
  numtrees_strsize = numtrees_str.size();
  encode_old = (numtrees_strsize+1)/2;    
  //cout << "encode_old is: " << encode_old << endl;
  tmpss.clear();
  //cout << "t1 is: " << t1 << " and t2 is: " << t2 << " (total is: " << total_ntrees << ")." << endl;
  tmpss << total_ntrees;
  tmpss >> numtrees_str;
  //cout << "numtrees_str is: " << numtrees_str << endl;
  numtrees_strsize = numtrees_str.size();
  encode_new = (numtrees_strsize+1)/2;
  //cout << "encode_new is: " << encode_new << endl;
  true_offset = 0;
  tmerge_parse_and_load(infilename, offset, 0, inverted_index, inverted_index_sizes, all_bipart, t1, encode_old, encode_new, all_branch, branch_sizes, true_ids, alldups, total_trees, true_offset);
  tmpss.clear();
  //exit(0);
  //now, do the second file
  tmpss << t2;
  tmpss >> numtrees_str;
  numtrees_strsize = numtrees_str.size();
  encode_old = (numtrees_strsize+1)/2;    
  tmpss.clear();
  //encode_new = (numtrees_strsize+1)/2;
  //cout << "encode_new is: " << encode_new << endl;
  unsigned int secsize = NUM_TREES - offset;
  true_offset = t1;
  tmerge_parse_and_load(mergefile, secsize, offset, inverted_index, inverted_index_sizes, all_bipart, t2, encode_old, encode_new, all_branch, branch_sizes, true_ids, alldups, total_trees, true_offset);

  //exit(0);
  //debug
  /* 
  cout << "bipartition list: " << endl;
  for (unsigned int i  = 0; i < all_bipart.size(); i++)
    cout << i << ": " << all_bipart[i] << endl;

  cout << "inverted index: " << endl;
  for (unsigned int i  =0; i < NUM_TREES; i++){
    cout << "T" << i << ": ";
    for (unsigned int j  = 0; j < inverted_index_sizes[i]; j++){
      cout << inverted_index[i][j] << " ";
    }
    cout << endl;
  }

  if (WEIGHTED){
    cout << "inverted index (branches):" << endl;
    for (unsigned int i = 0; i < total_trees; i++){
      cout << "T" << i << ": ";
      for (unsigned int j = 0; j < all_branch[i].size(); j++){
	cout << all_branch[i][j] << " ";
      }
      cout << endl;
    }
  }
*/
  int * chaindup = (int *)malloc(total_trees*sizeof(int));
  if (chaindup == NULL){
    cerr << "could not allocate!" << endl;
    exit(2);
  }
  for (unsigned int i  =0; i < total_trees; i++)
    chaindup[i] = -1;

  unsigned int pl;
  unsigned int alldup_total, alldup_curr;
  alldup_total = alldups.size();
  for (unsigned int i  = 0; i < alldup_total; i++){
    pl = i;
    alldup_curr = alldups[i].size();
    for (unsigned int j = 0; j < alldup_curr; j++){
      chaindup[pl] = alldups[i][j];
      pl = alldups[i][j];
    }
  }

  cout << "done." << endl;
  /*
  cout << "chaindups: " << endl;
  for (unsigned int i = 0; i < total_trees; i++){
    cout << i << "-->" << chaindup[i] << endl;
  }
  cout << endl;
  
  exit(0);
  for (unsigned int i  = 0; i < total_trees; i++){
    cout << "T[" << i << "] = ";
    for (unsigned int j = 0; j < all_branch[i].size(); j++){
      cout << all_branch[i][j] << " ";
    }
    cout << endl;
  }
  */
  unsigned int * foundun = (unsigned int *)malloc(total_trees*sizeof(unsigned int));
  if (foundun == NULL){
    cerr << "cannot allocate array!" << endl;
    exit(2);
  }
 for (unsigned int i  =0; i < total_trees; i++)
    foundun[i] = 0;

  useoffset = false;
  

 //seed new hash functions
  if (!quiet)
    fprintf(stderr, "step 4: seed new hash functions, and rehash trees\n");
  HashRFMap tmerge_hash;
  tmerge_hash.uhashfunc_init_unique(NUM_TREES, NUM_TAXA, all_bipart.size(), C, NEWSEED);
  unsigned long long sm1 = tmerge_hash._HF.getM1();
  unsigned long long sm2 = tmerge_hash._HF.getM2();
  
  tmerge_hash._hashtab2.resize(sm1); //resize the hash table
 
  //rehash trees
  unsigned int loc, trueid;
  trueid = 0;
  for (unsigned int i = 0; i < NUM_TREES; i++){ 
    unsigned int pos = 0;
    loc = inverted_index_sizes[i];
    unsigned long long temp1 = 0;
    unsigned long long temp2 = 0;
    bool * cbs = new bool[NUMBIPART];
    for (unsigned int j = 0; j < NUMBIPART; j++) {
      cbs[j] = 0;
    }
    for (unsigned int j = 0; j < loc; j++){
      pos = inverted_index[i][j];
      cbs[pos] = 1;
      temp1+= tmerge_hash._HF.getA1(pos);
      temp2+= tmerge_hash._HF.getA2(pos);
    }
    temp1 = temp1 % sm1;
    temp2 = temp2 % sm2;
    //while(true_ids[trueid] == -1)
    //  trueid++;
    //cout << "the next true id is: " << trueid << endl;
    tmerge_hash.hashing_bs_unique(i, NUMBIPART, temp1, temp2, cbs);

    delete [] cbs;
  }
 
  //step 4: determine which tree ids matter
  vector<unsigned int> want;
  switch(option){
  case 1: 
    tmerge_union(tmerge_hash, want, true_ids, chaindup, quiet);
    break;
  case 2: tmerge_intersection(tmerge_hash, want, offset, quiet);
    break;
  case 3: tmerge_setDifference(tmerge_hash, want, offset, quiet);
    break;
  default: 
    cerr << "Invalid option!" << endl; exit(4);
  }


  /*cout << "chaindups is now: " << endl;
  for (unsigned int i = 0; i < total_trees; i++){
    cout << i << "-->" << chaindup[i] << endl;
  }
  cout << endl;*/

  vector<string> newdups;
  stringstream sscon;
  string wantstring;
  unsigned int encode_size, wantstringsize;
  sscon << want.size();
  sscon >> wantstring;
  wantstringsize = wantstring.size();
  encode_size = (wantstringsize+1)/2;    
  unsigned int max = want.size();
  if (option == 5){
    unsigned int x = 0;
    x = 0;
    //now, using chained dups, figure out the chains and encode them
    unsigned int xpos = 0;
    unsigned int xsize = 0;
    unsigned int pl2;
    string encoded_dup;
    for (unsigned int i = 0; i < total_trees; i++){
      pl = i;
      foundun[xpos] = i;
      xpos++;
      xsize++;
      //cout << i;
      while (chaindup[pl] != -1){
	//cout << "-->" << chaindup[pl];
	foundun[xpos] = chaindup[pl];
	pl2 = chaindup[pl];
	chaindup[pl] = -1;
	pl = pl2;
	xpos++;
	xsize++;
      }
      chaindup[pl] = -1;
      //cout << endl;
      if (xsize > 1){ //encode
	encoded_dup = encode(foundun, xsize);
	//cout << "the encoded duplicate line is: " << encoded_dup << endl;
	newdups.push_back(encoded_dup);
      }
      xpos = 0;
      xsize = 0;
    }
  }
  //allocate the bipartition table

  vector<vector<unsigned int> > bipart_table;
  bipart_table.resize(NUMBIPART);
  for (unsigned int i  =0; i < NUMBIPART; i++)
    bipart_table[i].resize(0);

  vector< vector<string> > branch_table;
  vector<string> encode_branches;
  if (WEIGHTED){
    branch_table.resize(NUMBIPART);
    for (unsigned int i  = 0; i < NUMBIPART; i++)
      branch_table[i].resize(0);

    encode_branches.resize(NUMBIPART);
    for (unsigned int i = 0; i < NUMBIPART; i++)
      encode_branches[i] = "";
  }



  //given a list of tree ids that matter, rebuild bipartition table, discarding inverted index in the process
  //cout << "the trees I want are: [";
  //for (unsigned int i = 0; i < want.size(); i++)
  //  cout << want[i] << " ";
  //cout << "]" << endl;
  vector < unsigned int > branch_positions(total_trees, 0);

  unsigned int x = 0;
  unsigned int temp, mytrue_id, mypos;
  string temp_branch;
  //cout << "NUM_TREES is: " << NUM_TREES << endl;
  //cout << "want.size() is: " << want.size() << endl;
  //cout << "NUMBIPART is: " << NUMBIPART << endl;
  //cout << "encode_size is: " << encode_size << endl;
  //cout << "encode_new is: " << encode_new << endl;
  if (want.size() >  0){
    for (unsigned int i =0; i < NUM_TREES; i++){
      //cout << "i is: " << i << endl;
      if (i == want[x]){
	//cout << "i is: " << i << endl;
	mytrue_id = true_ids[i];
	//cout << "mytrue_id is: " << mytrue_id << endl;
	for (unsigned int j = 0; j < inverted_index_sizes[i]; j++){
	  temp = inverted_index[i][j];
	  bipart_table[temp].push_back(x);
	  if (WEIGHTED){
	    mypos = branch_positions[mytrue_id];
	    branch_positions[mytrue_id]++;
	    //cout << "mypos is: " << mypos << endl;
	    temp_branch = all_branch[mytrue_id][mypos];
	    //cout << "temp_branch is: " << temp_branch << endl;
	    if ((mytrue_id != x) || (encode_size != encode_new)){  //if id is not correct, convert temp branch's id to x
	      //temp_branch = convert_id(temp_branch, x, encode_size);
	      unsigned int strpad, tmpad;
	      strpad = 0;
	      string newtemp, myprefix, tmplen, result;
	      stringstream stm;
	      int pos;	      
	      pos = temp_branch.find_first_of(" ");
	      result = temp_branch.substr(0, pos);	      
	      if (result != ""){ //check if we have integral values
		stm << x;
		stm >> myprefix; //convert the integral values into a string
		stm.str("");
		//cout << "temp_branch is: " << temp_branch << endl;
		//cout << "myprefix is: " << myprefix << endl;
		temp_branch = temp_branch.substr(encode_new); //get rid of encoded_id part
		if (myprefix.size() < (encode_size*2))
		  strpad = (encode_size*2)-myprefix.size();
		else 
		  strpad = 0;
		while (strpad > 0){
		  myprefix = "0" + myprefix; 
		  strpad--;
		}
		stringstream tempstr;
		if (encode_size > 1){ 
		  strpad = 0;
		  while ((strpad < encode_size)){ 
		    tmplen = myprefix.substr(0,2);
		    tmpad = atoi(tmplen.c_str());
		    tmpad+=33;
		    tempstr << char(tmpad);
		    strpad++;
		    myprefix = myprefix.substr(2);
		  }
		}
		else{
		  tmpad = atoi(myprefix.c_str());
		  tmpad+=33;
		  tempstr << char(tmpad);
		}  
		tempstr << temp_branch;
		temp_branch = tempstr.str();
		tempstr.str("");
	      }
	    } // end new optimization

	    branch_table[temp].push_back(temp_branch);
	  }
	}

	if (x < max)
	  x++;
      } 
      free(inverted_index[i]);
    }
    free(inverted_index);
    free(inverted_index_sizes);
  }
 
  //exit(0);    

  if (WEIGHTED){
    for (unsigned int i = 0; i < NUMBIPART; i++){
      string subs, subs2, tempstr;
      stringstream intss, fracss;
      int pos;
      for (unsigned int j = 0; j < branch_table[i].size(); j++){
	tempstr = branch_table[i][j];
	pos = tempstr.find_first_of(" ");
	subs = tempstr.substr(0, pos);
	if (subs != ""){
	  intss << subs;
	}
	tempstr = tempstr.substr(pos+1);
	fracss << tempstr;
      }
      subs2 = intss.str() + " " + fracss.str();
      encode_branches[i] = subs2;
      intss.str("");
      fracss.str("");
    }
  }

  //output file contents
  if (outfile == "") //if no outfile is specified, use merged.trz, else use specified file
    outfile = "merged.trz";
  ofstream fout;
  fout.open(outfile.c_str()); 
  if (!fout){
     cerr << "cannot open file for writing!" << endl;
    exit(3);
  }
  fout << "TAXA " << lm.name(0);
  for (unsigned int i  =1; i < NUM_TAXA; i++){
    fout << ":" << lm.name(i);
  }
  fout << endl;

  ifstream fin, fin2;
  string newLine, ltype;
  unsigned int truesize = 0;
  if (option == 5){
    fin.open(infilename.c_str());
    fin2.open(mergefile.c_str());
    getline(fin, newLine);
    getline(fin, newLine);
    getline(fin2, newLine); //repeat with second file
    getline(fin2, newLine);
    fout << "NTREES " << t1+t2;
    if (HETERO)
      fout << " H";
    fout << endl;
    string unq1, unq2;
    unsigned int unc1, unc2;
    getline(fin, newLine);
    unc1 = get_unique(newLine, t1);
    getline(fin2, newLine);
    unc2 = get_unique(newLine, t2);
     //unc1 will be correct -> in a merge, the the number of unique trees in the first file is the same
    //what can change is the number of unique trees in the second file
    truesize = want.size();
    if (truesize < (t1+t2)){ //need to encode tree ids 
      fout << "UNIQUE_T " << truesize << endl;
    }
    else
      fout << "UNIQUE_T all"<< endl; 
  }
  else{ //if option is not 5
    fout << "NTREES " << want.size();
    if (HETERO)
      fout << " H";
    fout << endl;
    fout << "UNIQUE_T all"<< endl;
  }
  //count up the new number of expressed bipartitions
  unsigned int new_biparts = 0;
  for (unsigned int i = 0; i < NUMBIPART; i++){
    if (bipart_table[i].size()){
      new_biparts++;
    }
  }
  //cout << "new_biparts is: " << new_biparts << endl;
  //if (want.size() > 0)
  //  assert(new_biparts == NUMBIPART);
 
  fout << "NBIPARTITIONS " << new_biparts << endl;
  //next sort the bipartitions by the number of ones
  vector< pair<unsigned int, pair<string, unsigned int> > > ordered_bitstrings;
  unsigned nOnes;
  string mybitstr;
  for (unsigned int i = 0; i < NUMBIPART; i++){
    nOnes = count_ones(all_bipart[i]);
    convert_bitstring(all_bipart[i], mybitstr);
    ordered_bitstrings.push_back(make_pair(nOnes, make_pair(mybitstr, i)));
  }
  sort(ordered_bitstrings.begin(), ordered_bitstrings.end(), sort_b);
  
  //for every line in the bipart table,
  string compressed; 
  vector< pair<unsigned int, pair<string, unsigned int> > >::iterator b_itr = ordered_bitstrings.begin();
  unsigned int pos = 0;
  //debug
  /*while (b_itr != ordered_bitstrings.end()){
    pair<string, unsigned int> temp_pair = b_itr -> second;
    cout << b_itr->first << ": " << temp_pair.first << endl;
    b_itr++;
  }
  exit(0);*/
  while (b_itr != ordered_bitstrings.end()){
    pair<string, unsigned int> temp_pair = b_itr -> second;
    pos = temp_pair.second;
    if (bipart_table[pos].size()){
      fout << all_bipart[pos] << " ";
      compressed = compress_line(bipart_table[pos], want.size());
      fout << compressed;
      if (WEIGHTED)
	fout << " " << encode_branches[pos]; 
      fout << endl;
    }
    b_itr++;
  }

  if (option == 5){
    cout << "TrueSize is!! " << truesize << endl;
    if (truesize < total_trees){
      fout << "DUPLICATES " << newdups.size() << endl;
      for (unsigned int i = 0; i < newdups.size(); i++)
	fout << newdups[i] << endl;
    }
    else
      fout << "DUPLICATES 0" << endl;
  }
  else{
    fout << "DUPLICATES 0" << endl;
  }

}


void parse_integral(string temp, vector<string> &hold_integral, unsigned int encode_old, unsigned int encode_new, unsigned int offset){

  unsigned int m = 0;
  unsigned int intpos = 0;
  unsigned char v;
  unsigned int strpad, tmpad;
  unsigned int treeloc = 0;
  string newtemp, myprefix, tmplen;
  stringstream tempstr;
  tempstr.str("");

  while (intpos < temp.size()){
    m = 0;
    treeloc = 0;
    while (m < encode_old){ //while we're less than the encoded size   
      v = temp[intpos];
      treeloc+= v;
      treeloc-=33;
      m++;
      if (m < encode_old)
	treeloc *= 100;
      intpos++; 
    }
    //cout << "treeloc is: " << treeloc << endl;
    treeloc += offset; //calculate new tree id
    //cout << "treeloc is now: " << treeloc << endl;
    //stm.str("");    
    stringstream stm;              
    stm << treeloc;
    stm >> myprefix; //covert tree id into string
    //cout << "myprefix is first: " << myprefix << endl;
    if (myprefix.size() < (encode_new*2)){
      strpad = (encode_new*2)-myprefix.size();  
      while (strpad > 0){     
	myprefix = "0" + myprefix; 
	strpad--;
      } 
    }
    //cout << "myprefix is: " << myprefix << endl;
    //cout << "encode_new is: " << encode_new << endl;
    if (encode_new > 1){
      strpad = 0;
      while ((strpad < encode_new)){
	tmplen = myprefix.substr(0,2);
	tmpad = atoi(tmplen.c_str());
	tmpad+=33;
	//fprintf(fout, "%c", tmpad); 
	//printf("%c", tmpad);
	tempstr << char(tmpad);
	strpad++;
	myprefix = myprefix.substr(2);
      }
    }
    else{
      tmpad = atoi(myprefix.c_str());
      tmpad+=33;
      //fprintf(fout, "%c", tmpad);
      //printf("%c", tmpad);
      tempstr << char(tmpad);
    }  

    //cout << "tempstr is now: " << tempstr.str();
    //cout << "intpos is now: " << intpos << endl;
    v = temp[intpos];
    //cout << "v is: " << v << endl;
    tempstr << v;
    intpos++;
    if (v == '+'){
      v = temp[intpos];
      while (v != '+'){
	v = temp[intpos];
	tempstr << v;
	intpos++;
      }
    }
    //cout << "inserting " << tempstr.str() << " into position " << treeloc << endl;
    hold_integral[treeloc] = tempstr.str();
    tempstr.str("");
    //exit(0);
  }
  //exit(0);
  //temp = tempstr.str(); //put new value of temp into temp
  //tempstr.str("");
  //cout << "temp is finally: " << temp << endl;  
}

void parse_and_add_line(string treeline, unsigned int ** inverted_index, unsigned int * inverted_index_sizes, unsigned int bp, unsigned int uTrees, unsigned int offset,  vector< vector<string> > &all_branch, vector<unsigned int> & branch_sizes, unsigned int encode_old, unsigned int encode_new, unsigned int nTrees, int * true_ids, vector< vector<unsigned int> > alldups, unsigned int totalTrees, unsigned int trueOffset){
  //parse the set of trees associated with the line
  //cout << "nTrees is: " << nTrees << endl;
  //cout << "treeline is: " << treeline << endl;
  int pos;
  unsigned int count, beg;
  string line_type, mycount, ids, branches, integral, fractional;
  vector<string> hold_integral;
  stringstream ss;
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
      fprintf(stderr, "WEIGHTED trees detected!\n");
      WEIGHTED= true;
    }
    ids = treeline.substr(0, pos);
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of("\n");
    branches = treeline.substr(0, pos);
    pos = branches.find_first_of(" ");
    integral = branches.substr(0,pos);
    fractional = branches.substr(pos+1);
    //cout << "integral is: " << integral << endl;
    //cout << "fractional is: " << fractional << endl;
    //use integral portion to fill hold_integral
    hold_integral.resize(totalTrees);
    for (unsigned int i  = 0; i < totalTrees; i++)
      hold_integral[i] = "";
    parse_integral(integral, hold_integral, encode_old, encode_new, trueOffset);
  }
  else{
    pos = treeline.find_first_of("\n");
    ids = treeline.substr(0, pos);
  }
  
  unsigned int loc_size, myloc;
  if (count == 0) { 
    unsigned int true_loc;
    for (unsigned int i = 0; i < uTrees; i++) {
      true_loc = i + offset; 
      loc_size = inverted_index_sizes[true_loc];
      inverted_index[true_loc][loc_size] = bp;
      inverted_index_sizes[true_loc]++;
    }
    if (WEIGHTED){
      string mybranch;
      beg = 0;
      for (unsigned int i = 0; i < nTrees; i++){
	true_loc = i+trueOffset;
	mybranch = hold_integral[true_loc];
	string tempstr(fractional, beg, 3);
	mybranch.append(" ").append(tempstr);
	myloc = branch_sizes[true_loc];
	branch_sizes[true_loc]++;
	all_branch[true_loc][myloc] = mybranch;
	//all_branch[true_loc].push_back(mybranch);
	beg+=3;
      }
    }
  }
  else{
    unsigned int dec_size = 0;
    unsigned int * found = (unsigned int *)malloc(uTrees*sizeof(unsigned int));
    if (found == NULL){
      cerr << "cannot allocate found!" << endl;
      exit(2);
    }
    dec_size = decode(ids, found);
    vector<unsigned int> all_id;

    //add to inverted index
    if (line_type == "-") { //compressed line
      bool * check = (bool *)malloc(nTrees*sizeof(bool));
      for (unsigned int i = 0; i < uTrees; i++)
	check[i] = false;

      for (unsigned int i = 0; i < dec_size; ++i) {
	unsigned int temp = found[i];
	check[temp] = true;
      }
      unsigned int true_loc;
      unsigned int my_trueid;
      for (unsigned int i = 0; i < uTrees; ++i) {
	if (check[i] == false) {
	  true_loc = i+offset;
	  my_trueid = true_ids[true_loc];
	  all_id.push_back(my_trueid); 
	  loc_size = inverted_index_sizes[true_loc];
	  inverted_index[true_loc][loc_size] = bp;
	  inverted_index_sizes[true_loc]++;
	  //cout << "adding " << bp << " to inverted_index[" << true_loc << "][" << loc_size << "]." << endl;
	}
      }
      free(check);
      if (WEIGHTED){
	unsigned int size = all_id.size();
	unsigned int temploc, secloc;
	for (unsigned int i = 0; i < size; i++){
	  temploc = all_id[i];
	  if (alldups[temploc].size()){
	      for (unsigned int j = 0; j < alldups[temploc].size(); j++){
		secloc = alldups[temploc][j];
		all_id.push_back(secloc);
	      }
	  }
	}
	sort(all_id.begin(), all_id.end());
	//now, using the all_ids vector, please do the following:
	
	string mybranch;
	unsigned int true_loc;
	beg = 0;
	for (unsigned int i = 0; i < all_id.size(); i++){
	  true_loc = all_id[i];
	  mybranch = hold_integral[true_loc];
	  string tempstr(fractional, beg, 3);
	  mybranch.append(" ").append(tempstr);
	  myloc = branch_sizes[true_loc];
	  branch_sizes[true_loc]++;
	  all_branch[true_loc][myloc] = mybranch;
	  //all_branch[true_loc].push_back(mybranch);
	  beg+=3;
	}
		
      } //end weighted
    } //end bout line
    else{
      for (unsigned int i = 0; i < count; ++i) { 
	unsigned int temp = found[i];
	temp+=offset;
	loc_size = inverted_index_sizes[temp];
	inverted_index[temp][loc_size] = bp;
	inverted_index_sizes[temp]++;
      } 
      if (WEIGHTED){
	unsigned int my_trueid, sec_id;
	for (unsigned int i = 0; i < count; i++){
	  unsigned int temp = found[i];
	  temp += offset;
	  my_trueid = true_ids[temp];
	  all_id.push_back(my_trueid);
	  if (alldups[my_trueid].size()){
	    for (unsigned int j = 0; j < alldups[my_trueid].size(); j++)
	      sec_id = alldups[my_trueid][j];
	    all_id.push_back(sec_id);
	  }
	}
	sort(all_id.begin(), all_id.end());
	
	string mybranch;
	unsigned int true_loc;
	beg = 0;
	for (unsigned int i = 0; i < all_id.size(); i++){
	  true_loc = all_id[i];
	  mybranch = hold_integral[true_loc];
	  string tempstr(fractional, beg, 3);
	  mybranch.append(" ").append(tempstr);
	  myloc = branch_sizes[true_loc];
	  branch_sizes[true_loc]++;
	  all_branch[true_loc][myloc] = mybranch;
	  //all_branch[true_loc].push_back(mybranch);
	  beg+=3;
	}
		
      } //end if weighted
    } //end else uncompressed line (bin)
    free(found);
  } //end else count > 0
  
}

void getdups(string filename, unsigned int &nt, vector<string> & dups){
  ifstream fin;
  string str, abipart, ltype;
  int pos;
  unsigned int ndup, bipart;
  fin.open(filename.c_str());
  if (!fin){
    cerr << "cannot open file " << filename << " for reading!" << endl;
    exit(1);
  }
  getline(fin, str); //skip taxa
  getline(fin, str); //get num trees
  nt = get_ntrees(str);
  getline(fin, str); //skip unique
  getline(fin, str);  //get bipart
  parse_and_get(str, "NBIPARTITIONS", bipart);
  //cout << "bipartitions: " << bipart << endl;
  for (unsigned int i = 0; i < bipart; i++)
    getline(fin, str); //skip all the bipartitions
  getline(fin, str); //read duplicates line;
  parse_and_get(str, "DUPLICATES", ndup);
  for (unsigned int i  = 0; i < ndup; i++){
    getline(fin, str);
    pos = str.find_first_of("\n");
    ltype = str.substr(0, pos);
    dups.push_back(ltype);
  }  
  fin.close();
}

void tmerge_parse_and_load(string filename, unsigned int uTrees, unsigned int offset, unsigned int ** inverted_index, unsigned int * inverted_index_sizes, vector<string>& all_bipart, unsigned int nTrees, unsigned int encode_old, unsigned int encode_new, vector< vector<string> > &all_branch, vector<unsigned int> & branch_sizes, int * true_ids, vector< vector<unsigned int> > alldups, unsigned int totalTrees, unsigned int trueOffset){

  //next, parse the TRZ files -- build inverted index
  ifstream fin;
  string str, abipart, ltype;
  int pos;
  unsigned int bipart;
  fin.open(filename.c_str());

  if (!fin){
    cerr << "cannot open file " << filename << " for reading!" << endl;
    exit(1);
  }

  for (unsigned int i = 0; i < 4; i++)
    getline(fin, str); //skip first four lines

  parse_and_get(str, "NBIPARTITIONS", bipart);
  //cout << "bipartitions is: " << bipart << endl;
  for (unsigned int i = 0; i  < bipart; i++){ //should be bipart!!!!
    getline(fin, str);
    pos = str.find_first_of(" ");
    abipart = str.substr(0, pos);
    //cout << "our bipartition is: " << abipart << endl;
    unsigned int j = 0; 
    while (j != NUMBIPART){
      if (abipart == all_bipart[j])
	break;
      j++;
    }
    if (j == NUMBIPART){
      cerr << "could not find desired bipartition!"; 
      exit(2);
    }
    //cout << "the id of this bipartition is: " << j << endl;
    //j is now the id that we will add to our inverted index
    str = str.substr(pos+1);

    //------------ parse_and_add_line starts here!
    int pos;
    unsigned int count, beg;
    string line_type, mycount, ids, branches, integral, fractional;
    vector<string> hold_integral;
    stringstream ss;
    string treeline = str;
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
	fprintf(stderr, "WEIGHTED trees detected!\n");
	WEIGHTED= true;
      }
      ids = treeline.substr(0, pos);
      treeline = treeline.substr(pos+1);
      pos = treeline.find_first_of("\n");
      branches = treeline.substr(0, pos);
      pos = branches.find_first_of(" ");
      integral = branches.substr(0,pos);
      fractional = branches.substr(pos+1);
      //cout << "integral is: " << integral << endl;
      //cout << "fractional is: " << fractional << endl;
      //use integral portion to fill hold_integral
      hold_integral.resize(totalTrees);
      for (unsigned int i  = 0; i < totalTrees; i++)
	hold_integral[i] = "";
      parse_integral(integral, hold_integral, encode_old, encode_new, trueOffset);
    }
    else{
      pos = treeline.find_first_of("\n");
      ids = treeline.substr(0, pos);
    }
    
    unsigned int loc_size, myloc;
    if (count == 0) { 
      unsigned int true_loc;
      for (unsigned int i = 0; i < uTrees; i++) {
	true_loc = i + offset; 
	loc_size = inverted_index_sizes[true_loc];
	inverted_index[true_loc][loc_size] = j;
	inverted_index_sizes[true_loc]++;
      }
      if (WEIGHTED){
	string mybranch;
	beg = 0;
	for (unsigned int i = 0; i < nTrees; i++){
	  true_loc = i+trueOffset;
	  mybranch = hold_integral[true_loc];
	  string tempstr(fractional, beg, 3);
	  mybranch.append(" ").append(tempstr);
	  myloc = branch_sizes[true_loc];
	  branch_sizes[true_loc]++;
	  all_branch[true_loc][myloc] = mybranch;
	  //all_branch[true_loc].push_back(mybranch);
	  beg+=3;
	}
      }
    }
    else{
      unsigned int dec_size = 0;
      unsigned int * found = (unsigned int *)malloc(uTrees*sizeof(unsigned int));
      if (found == NULL){
	cerr << "cannot allocate found!" << endl;
	exit(2);
      }
      dec_size = decode(ids, found);
      vector<unsigned int> all_id;
      
      //add to inverted index
      if (line_type == "-") { //compressed line
	bool * check = (bool *)malloc(nTrees*sizeof(bool));
	for (unsigned int i = 0; i < uTrees; i++)
	  check[i] = false;
	
	for (unsigned int i = 0; i < dec_size; ++i) {
	  unsigned int temp = found[i];
	  check[temp] = true;
	}
	unsigned int true_loc;
	unsigned int my_trueid;
	for (unsigned int i = 0; i < uTrees; ++i) {
	  if (check[i] == false) {
	    true_loc = i+offset;
	    my_trueid = true_ids[true_loc];
	    all_id.push_back(my_trueid); 
	    loc_size = inverted_index_sizes[true_loc];
	    inverted_index[true_loc][loc_size] = j;
	    inverted_index_sizes[true_loc]++;
	    //cout << "adding " << bp << " to inverted_index[" << true_loc << "][" << loc_size << "]." << endl;
	  }
	}
	free(check);
	if (WEIGHTED){
	  unsigned int size = all_id.size();
	  unsigned int temploc, secloc;
	  for (unsigned int i = 0; i < size; i++){
	    temploc = all_id[i];
	    if (alldups[temploc].size()){
	      for (unsigned int j = 0; j < alldups[temploc].size(); j++){
		secloc = alldups[temploc][j];
		all_id.push_back(secloc);
	      }
	    }
	  }
	  sort(all_id.begin(), all_id.end());
	  //now, using the all_ids vector, please do the following:
	  
	  string mybranch;
	  unsigned int true_loc;
	  beg = 0;
	  for (unsigned int i = 0; i < all_id.size(); i++){
	    true_loc = all_id[i];
	    mybranch = hold_integral[true_loc];
	    string tempstr(fractional, beg, 3);
	    mybranch.append(" ").append(tempstr);
	    myloc = branch_sizes[true_loc];
	    branch_sizes[true_loc]++;
	    all_branch[true_loc][myloc] = mybranch;
	    //all_branch[true_loc].push_back(mybranch);
	    beg+=3;
	  }
	  
	} //end weighted
      } //end bout line
      else{
	for (unsigned int i = 0; i < count; ++i) { 
	  unsigned int temp = found[i];
	  temp+=offset;
	  loc_size = inverted_index_sizes[temp];
	  inverted_index[temp][loc_size] = j;
	  inverted_index_sizes[temp]++;
	} 
	if (WEIGHTED){
	  unsigned int my_trueid, sec_id;
	  for (unsigned int i = 0; i < count; i++){
	    unsigned int temp = found[i];
	    temp += offset;
	    my_trueid = true_ids[temp];
	    all_id.push_back(my_trueid);
	    if (alldups[my_trueid].size()){
	      for (unsigned int j = 0; j < alldups[my_trueid].size(); j++){
		sec_id = alldups[my_trueid][j];
		all_id.push_back(sec_id);
	      }
	    }
	  }
	  sort(all_id.begin(), all_id.end());
	  string mybranch;
	  unsigned int true_loc;
	  beg = 0;
	  for (unsigned int i = 0; i < all_id.size(); i++){
	    true_loc = all_id[i];
	    mybranch = hold_integral[true_loc];
	    string tempstr(fractional, beg, 3);
	    mybranch.append(" ").append(tempstr);
	    myloc = branch_sizes[true_loc];
	    branch_sizes[true_loc]++;
	    all_branch[true_loc][myloc] = mybranch;
	    //all_branch[true_loc].push_back(mybranch);
	    beg+=3;
	  }
	  
	} //end if weighted
      } //end else uncompressed line (bin)
      free(found);
    } //end else count > 0
  }
  fin.close();   
}

void tmerge_union(HashRFMap tmerge_hash, vector<unsigned int> & want, int * true_ids, int * chaindup, bool quiet){
  if (!quiet)  
    cout << "Selected option: UNION" << endl;
  unsigned int pl, tmpid, tmpid2;
  for (unsigned int hti=0; hti<tmerge_hash._hashtab2.size(); ++hti) {
    unsigned int sizeVec = tmerge_hash._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; ++i) {
	tmpid = tmerge_hash._hashtab2[hti][i]._vec_treeidx[0];
	//cout << tmpid << " ";
        want.push_back(tmpid);
	pl = true_ids[tmpid];
	unsigned int debug = tmerge_hash._hashtab2[hti][i]._vec_treeidx.size();
	for (unsigned int j = 1; j < debug; j++){
	  tmpid = tmerge_hash._hashtab2[hti][i]._vec_treeidx[j];
	  tmpid2 = true_ids[tmpid];
	  while (chaindup[pl]!= -1){
	    tmpid = chaindup[pl];
	    pl = tmpid;
	  }
	  chaindup[pl] = tmpid2;
	  //cout << tmerge_hash._hashtab2[hti][i]._vec_treeidx[j] << " ";
	}
	//cout << endl;
      }
    } //end if sizevec
    
  } //end hashtable traversal

  sort(want.begin(), want.end());
  if (!quiet)
    cout << "found: " << want.size() << " unique trees!" << endl;
}


void tmerge_intersection(HashRFMap tmerge_hash, vector<unsigned int> & want, unsigned int offset, bool quiet){
  bool first, second;
  unsigned int temp_id;
  if (!quiet)
    cout << "Selected option: INTERSECTION" << endl;
  for (unsigned int hti=0; hti<tmerge_hash._hashtab2.size(); hti++) {
    unsigned int sizeVec = tmerge_hash._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; i++) {
        unsigned int sizeTreeIdx = tmerge_hash._hashtab2[hti][i]._vec_treeidx.size();
        first = false;
        second = false;
        for (unsigned int j = 0; j < sizeTreeIdx; j++){  
          temp_id = tmerge_hash._hashtab2[hti][i]._vec_treeidx[j];
          if (temp_id < offset)
            first = true;
          if (temp_id >= offset)
            second = true;
        }
        if (first && second)
          want.push_back(tmerge_hash._hashtab2[hti][i]._vec_treeidx[0]);
      }
    } //end if sizevec
  } //end hashtable traversal
  sort(want.begin(), want.end());
  if (!quiet)
    cout << "found: " << want.size() << " unique trees in the intersection!" << endl;
}
void tmerge_setDifference(HashRFMap tmerge_hash, vector<unsigned int> & want, unsigned int offset, bool quiet){
  bool first, second;
  unsigned int temp_id;
  if (!quiet)
    cout << "Selected option: SET DIFFERENCE" << endl;
  for (unsigned int hti=0; hti<tmerge_hash._hashtab2.size(); hti++) {
    unsigned int sizeVec = tmerge_hash._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; i++) {
        unsigned int sizeTreeIdx = tmerge_hash._hashtab2[hti][i]._vec_treeidx.size();
        first = false;
        second = false;
        for (unsigned int j = 0; j < sizeTreeIdx; j++){  
          temp_id = tmerge_hash._hashtab2[hti][i]._vec_treeidx[j];
          if (temp_id < offset)
            first = true;
          if (temp_id >= offset)
            second = true;
        }
        if (first && !second)
          want.push_back(tmerge_hash._hashtab2[hti][i]._vec_treeidx[0]);
      }
    } //end if sizevec
  } //end hashtable traversal
  sort(want.begin(), want.end());
  if (!quiet)
    cout << "found: " << want.size() << " unique trees in the set difference!" << endl;
  
}


