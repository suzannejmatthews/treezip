#include <iostream>
#include <stdlib.h>
#include "global.h"
#include "bmerge.h"
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <limits>
#include <stdio.h>

//char tocode[11] = "ABCDEFGHIJ";
string tocode = "ABCDEFGHIJ";

string rle(bool * temp, unsigned int length){
    unsigned int curr = temp[0];
    unsigned int currcnt = 0;
    unsigned int newvalb = 0;
    stringstream ss;
    for (unsigned int j =  0; j < length; j++) { 
      newvalb = temp[j];      
      if (newvalb != curr) { 
        if( currcnt > 1) { 
          if (curr == 1) { 
            //fprintf(fout, "K%d", currcnt);
            ss << "K" << currcnt;
          }
          else{
            //fprintf(fout, "L%d", currcnt);
            ss << "L" << currcnt;
          }
        }
        else {
          if (curr == 1) {  
            //fprintf(fout, "B");
            ss << "B";
          }
          else{
            //fprintf(fout, "A");
            ss << "A";
          }
        }
        curr = newvalb;
        currcnt = 1;
      }
      else{
        currcnt++;
      }
    }
    
    if (currcnt > 1) { 
      if (curr == 1) { 
        //fprintf(fout, "K%d", currcnt);
        ss << "K" << currcnt;
      }
      else{
        //fprintf(fout, "L%d", currcnt);
        ss << "L" << currcnt;
      }
}
    else{
      if (curr == 1) { 
        //fprintf(fout, "B");
        ss << "B";
      }
      else{
        //fprintf(fout, "A");
        ss << "A";
      }
    }
    return ss.str();
}

void decode_bitstring(string bitstring, bool * bs, unsigned int maxLength){
  short bitstring_size = bitstring.size();
  unsigned int j = 0;
  int x = 0;
  bool typeb = false;
  int ypos = 0;
  int y = 0;
  while ( x < bitstring_size ) {
    int val = bitstring[x];
    val-= 65;
    if (val < 2) { 
      bs[j] = (bool)val;
      ++j;
    }
    else{
      char type = bitstring[x];
      if (type == 'K')
	typeb = true;
      if (type == 'L')
	typeb = false;
      
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
	bs[j] = typeb;
	++j;
      }
    }
    ++x;
  }
  assert(j == maxLength);
} 

unsigned int get_bitstring_length(string bitstring){ //determines the size of the bitstring 
  short bitstring_size = bitstring.size();
  int x = 0;
  int ypos = 0;
  int y = 0;
  unsigned int count = 0;
  int val;
  while ( x < bitstring_size ) {
    char type = bitstring[x];
    if ((type == 'A') || (type == 'B')) //if it is either A or B, just increment by one
      count++;
    else if ((type == 'K') || (type == 'L')){ //
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
      count += iamt;
    }
    x++;
  } //end go through bitstring
  return count;
} 

//parses the ntrees line of the TRZ file
unsigned int get_ntrees(string str){
  unsigned int trees = 0;
  int pos = str.find_first_of(" ");
  string line_type = str.substr(0, pos);
  if (line_type != "NTREES"){
    cerr << "Error! Number of trees field not found. Exiting..\n";
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  pos = str.find_first_of(" ");
  if (pos != -1){
    string ntrees = str.substr(0, pos);
    str = str.substr(pos+1);
    if ( (str == "H") && (!HETERO) ){
      cerr << "Warning! Detected Heterogeneous collection of trees!" << endl;
      HETERO = true;
    }
    trees = atoi(ntrees.c_str());
  }
  else
    trees = atoi(str.c_str()); //this is the total number of trees
  return trees;
} 

//parses the unique line of the TRZ file
unsigned int get_unique(string str, unsigned int ntrees){
  int pos = str.find_first_of(" ");
  string line_type = str.substr(0, pos);
  unsigned int num_unique = ntrees;
  assert(num_unique != 0);
  if (line_type != "UNIQUE_T"){
    cerr << "Error! Unique trees line not found! Exiting...\n";
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  if (str != "all")
    num_unique = atoi(str.c_str());
  return num_unique;
}

//used for parsing nbipartitions and duplicates lines
void parse_and_get(string str, string check, unsigned int & var){
  int pos;
  pos = str.find_first_of(" ");
  string type =str.substr(0, pos);
  if (type != check){
    cerr << "Error! " << check << " field not found!" << endl;
    cout << str << endl;
    exit(1);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  var = atoi(str.c_str());
}

bool does_taxa_match(string infile, string mergefile, LabelMap & lm){
  ifstream fin, fin2;
  string str, str2, line_type;
  fin.open(infile.c_str());
  fin2.open(mergefile.c_str());
  if (!fin){
    cerr << "Cannot open file " << infile << " for reading!" << endl;
    exit(2);
  }
  if (!fin2){
    cerr << "Cannot open file " << mergefile << " for reading!" << endl;
  }
  getline(fin, str);
  getline(fin2, str2);
  int pos = str.find_first_of(" ");
  line_type = str.substr(0, pos);
  if (line_type != "TAXA"){
    cerr << "Error! No taxa labels identified for file! Exiting..." << endl;
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str2.find_first_of(" ");
  line_type = str2.substr(0, pos);
  if (line_type != "TAXA"){
    cerr << "Error! No taxa labels identified for file! Exiting..." << endl;
    exit(2);
  }
  str2 = str2.substr(pos+1);
  
  //now, check the labels
  stringstream ss(str);
  unsigned int ntaxa, ntaxa2, count;
  ntaxa = 0;
  ntaxa2 = 0;
  string taxa;
  while(getline(ss, taxa, ':')){
    lm.push(taxa);
    ntaxa++;
  }
  ss.clear();
  count =  0;
  ss.str(str2);
  while (getline(ss, taxa, ':')){
    if (taxa != lm.name(count)){
      return false;
    }
    count++;
    ntaxa2++;
  }
    
  if (ntaxa!= ntaxa2){
    cout << "number of taxa do not match!" << endl;
    return false;
  }
  NUM_TAXA = ntaxa;
  fin.close();
  fin2.close();
  return true;
}

void  rewrite_file(string file, vector<string> taxa){
  //create a temporary file that will store rewrite information
  ifstream fin;
  ofstream fout;
  string str;

  fin.open(file.c_str());
  if (!fin){
    cerr << "Cannot open file " << file << " for reading!"  << endl;
    exit(1);
  }
  fout.open(".tmptaxafile");
  if (!fout){
    cerr << "Error creating temporary file! Exiting..." << endl;
    exit(1);
  }
  fout << "TAXA " << taxa[0];
  //write out global set of taxa to temporary file
  for (unsigned int i = 1; i < taxa.size(); i++)
    fout << ":" << taxa[i];
  fout << endl;
  //allocate master bit-string of size global_taxa
  unsigned int taxasize = taxa.size();
  bool * global_bs = new bool[taxasize];   
  for (unsigned int i = 0; i < taxasize; i++)
    global_bs[i] = 0;
  getline(fin, str);
  int pos = str.find_first_of(" ");
  string l_type = str.substr(0, pos);
  if (l_type != "TAXA"){
    cerr << "Cannot find taxa!" << endl;
    exit(1);
  }
  str = str.substr(pos+1);
  stringstream ss(str);
  string mytaxon;
  vector<string> temptaxa;
  while (getline(ss, mytaxon, ':')){
    temptaxa.push_back(mytaxon);
   }
  //create helper array that is the size of current set of taxa
  unsigned int *helper = new unsigned int[temptaxa.size()];
  //shows what position the taxa is in the global_taxa set
  for (unsigned int i = 0; i < temptaxa.size(); i++){
    mytaxon = temptaxa[i];
    unsigned int j = 0;
    //cout << "mytaxon: " << mytaxon << endl;
    while (j != taxa.size()){
      if (taxa[j] == mytaxon)
	break;
      j++;
    }
    if (j == taxa.size()){
      cerr << "Error! Taxon not found!" << endl;
      exit(3);
    }
    helper[i] = j;
  }
  //(DEBUG) print out helper and check to make sure it matches!
  /*for (unsigned int i = 0; i < temptaxa.size(); i++)
    cout << helper[i] << " ";
    cout << endl;*/

  getline(fin, str);
  fout << str << endl;  //write out NTREES
  getline(fin, str);
  fout << str << endl;  //write out UNIQUE_T
  getline(fin, str);   //write out NBIPARTITIONS
  fout << str << endl;
  pos = str.find_first_of(" ");
  string bipart = str.substr(pos+1);
  unsigned int nbipart = atoi(bipart.c_str());
  unsigned int count  = 0;
  unsigned int maxLength = 0;
  unsigned int truepos = 0;
  unsigned int maxpos = 0;
  string newbipart;
  cout << "NBPIART is: " << nbipart << endl;
  while (count < nbipart){   //for every bipart:
    getline(fin, str);
    pos = str.find_first_of(" ");
    bipart = str.substr(0, pos);   //find first space, grab bipart string
    str = str.substr(pos+1); // remove the original bipart from str
    maxLength = get_bitstring_length(bipart);//convert string to bitstring
    bool * bs = new bool[maxLength];
    //cout << "bipart is: " << bipart << endl;
    decode_bitstring(bipart, bs, maxLength);
    /*cout << "decoded: ";
    for (unsigned int i = 0; i < maxLength; i++)
      cout << bs[i];
      cout << endl;*/
    maxpos = 0;
    for (unsigned int i = 0; i < maxLength; i++){
      //cout << "i is: " << i << endl;
      if (bs[i]){ //if bit is 1
	truepos = helper[i];
	//cout << "truepos is: " << truepos << endl;
	if (truepos > maxpos)
	  maxpos = truepos;     //keep track of the highest position
	global_bs[truepos] = 1;     //copy to appropriate place in master bit-string
      }
    }
    maxpos++;
    /* cout << "new: ";
    for (unsigned int i = 0; i < maxpos; i++)
      cout << global_bs[i];
      cout << endl;*/
    newbipart = rle(global_bs, maxpos);     //create a new bitstring based on range 
    //cout << "newbipart is: " << newbipart << endl;
    fout << newbipart << " " << str << endl;      //output bitstring + rest of line to temp file
    for (unsigned int i = 0; i < maxpos; i++)
      global_bs[i] = 0;   //reset master bit-string from 0 to highest
    delete [] bs;
    count++;
  }  //end for every bipart
  //output rest of file
  while (!fin.eof()){
    if (fin.eof())
      break;
    getline(fin, str);
    fout << str << endl;
  }
  fin.close();
  fout.close();
  //cout << "we get here!" << endl;
  //cout << "file is: " << file.c_str() << endl;
  rename(".tmptaxafile", file.c_str()); //rename temp file to file
  delete [] global_bs; //deallocate global
  delete [] helper; //deallocate helper
}

void  merge_taxa(string infile, string mergefile, LabelMap &lm){
  ifstream fin, fin2;
  string str, str2, line_type;
  unsigned int ntaxa = 0;
  fin.open(infile.c_str());
  fin2.open(mergefile.c_str());
  if (!fin){
    cerr << "Cannot open file " << infile << " for reading!" << endl;
    exit(2);
  }
  if (!fin2){
    cerr << "Cannot open file " << mergefile << " for reading!" << endl;
  }
  getline(fin, str);
  getline(fin2, str2);
  int pos = str.find_first_of(" ");
  line_type = str.substr(0, pos);
  if (line_type != "TAXA"){
    cerr << "Error! No taxa labels identified for file! Exiting..." << endl;
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str2.find_first_of(" ");
  line_type = str2.substr(0, pos);
  if (line_type != "TAXA"){
    cerr << "Error! No taxa labels identified for file! Exiting..." << endl;
    exit(2);
  }
  str2 = str2.substr(pos+1);
  vector<string> taxa1;
  vector<string> taxa2;   //load the taxa into the vectors
  stringstream ss(str);
  stringstream ss2(str2);
  string taxa;
  while(getline(ss, taxa, ':')){
    taxa1.push_back(taxa);
  }
  while(getline(ss2, taxa, ':')){
    taxa2.push_back(taxa);
  }
  fin.close();
  fin2.close();

  /*cout << "printing out contents of the arrays..." << endl;
  cout << "T1 = { ";
  for (unsigned int i  =0; i < taxa1.size(); i++)
    cout << taxa1[i] << " ";
  cout << "}" << endl;
  cout << "T2 = { ";
  for (unsigned int i  =0; i < taxa2.size(); i++)
    cout << taxa2[i] << " ";
  cout << "}" << endl;
  vector<string>::reverse_iterator t_itr = taxa2.rbegin();
  cout << "[ ";
  while (t_itr != taxa2.rend()){
    cout << (*t_itr) << " ";
    t_itr++;
  }
  cout << "]" << endl;*/
  if (taxa1.size() < taxa2.size()){
    vector<string>::reverse_iterator t_itr = taxa1.rbegin();
    vector<string>::iterator t2_itr = taxa2.begin();
    string taxon, taxon2;
    while (t_itr != taxa1.rend()){
      taxon = (*t_itr);
      t2_itr = taxa2.begin();
      while(t2_itr != taxa2.end()){
	taxon2 = (*t2_itr);
	if (taxon == taxon2)
	  break;
	t2_itr++;
      }
      if (t2_itr != taxa2.end()){ //we found it in the other list
	taxa2.erase(t2_itr);
      }
      t2_itr = taxa2.begin();
      taxa2.insert(t2_itr, taxon);
      t_itr++;
    }
    //display new taxa ordering:
    /* for (unsigned int i = 0; i < taxa2.size(); i++)
      cout << "t[" << i << "]: " << taxa2[i] << " ";
      cout << endl;*/
    for (unsigned int  i = 0; i < taxa2.size(); i++){
      lm.push(taxa2[i]);
      ntaxa++;
    }
    rewrite_file(mergefile, taxa2);
  }
  else{
    vector<string>::reverse_iterator t_itr = taxa2.rbegin();
    vector<string>::iterator t2_itr = taxa1.begin();
    string taxon, taxon2;
    while (t_itr != taxa2.rend()){
      taxon = (*t_itr);
      t2_itr = taxa1.begin();
      while(t2_itr != taxa1.end()){
	taxon2 = (*t2_itr);
	if (taxon == taxon2)
	  break;
	t2_itr++;
      }
      if (t2_itr != taxa1.end()){ //we found it in the other list
	taxa1.erase(t2_itr);
      }
      t2_itr = taxa1.begin();
      taxa1.insert(t2_itr, taxon);
      t_itr++;
    }
    //display new taxa ordering:
    /*for (unsigned int i = 0; i < taxa1.size(); i++)
      cout << "t[" << i << "]: " << taxa1[i] << " ";
      cout << endl;*/
    rewrite_file(infile, taxa1);
    for (unsigned int  i = 0; i < taxa1.size(); i++){
      lm.push(taxa1[i]);
      ntaxa++;
    }
  }
  NUM_TAXA = ntaxa; //this is the total number of taxa
}

unsigned int validate_merge(string & infile, string & mergefile, LabelMap &lm, unsigned int option, unsigned int & total_trees){
  /*opens up the two files and checks to make sure that:
     1. They are valid TRZ files
     2. The taxa match
     3. which file is smaller (based on number of bipartitions) --> switch merge and infile
         accordingly
     4. Note the offset. 
  */
  //first check to make sure both files are valid TRZ files (right now, check extension)
  total_trees = 0;
  string ext = infile.substr(infile.size()-4);
  string ext2 = mergefile.substr(mergefile.size()-4);
  if (ext != ".trz") {
    cerr << "Invalid extension. Please input a valid .trz for file 1\n";
    exit(2);
  }
  else if (ext2 != ".trz"){
   cerr << "Invalid extension. Please input a valid .trz for file 2\n";
   exit(2);
  }

  if (!HETERO){
    //next, populate taxa labels and ensure the taxa match
    bool match = does_taxa_match(infile, mergefile, lm);
    if (!match){
      cerr << "Taxa do not match! Assuming heterogenous collection..." << endl;
      HETERO = true;    
      lm.clear();
    }
  }
  if (HETERO){
    merge_taxa(infile, mergefile, lm);
  }
  ifstream fin, fin2;
  string str, str2, line_type;
  fin.open(infile.c_str());
  fin2.open(mergefile.c_str());
  if (!fin){
    cerr << "Cannot open file " << infile << " for reading!" << endl;
    exit(2);
  }
  if (!fin2){
    cerr << "Cannot open file " << mergefile << " for reading!" << endl;
  }
  getline(fin, str); //skip taxa line in file1
  getline(fin2, str2); //skip taxa line in file2

  //now that we know that the taxa match, we need note the offset 
  //(aka, the number of trees in the larger file).
  unsigned int ntrees1, ntrees2;
  getline(fin, str);
  ntrees1 = get_ntrees(str);
  total_trees += ntrees1;
  getline(fin2, str2);
  ntrees2 = get_ntrees(str2);
  total_trees += ntrees2; //at this point, totaltrees is populated
  getline(fin,str); //read in unique trees
  ntrees1 = get_unique(str, ntrees1); //update if it all trees not unique
  getline(fin2,str2); //read in unique trees
  ntrees2 = get_unique(str2, ntrees2); //also update if all trees not unique  
  NUM_TREES = ntrees1 + ntrees2;
  getline(fin, str);
  unsigned int nbiparts, nbiparts2;
  parse_and_get(str, "NBIPARTITIONS", nbiparts);
  getline(fin2, str2); //now do it for file2
  parse_and_get(str2, "NBIPARTITIONS", nbiparts2);
  
  //  if (option == 3)
  // NUM_TREES = ntrees1;
  
  return ntrees1;
}

void setup_bmerge(string infilename, string mergefile, int option){
  unsigned int offset = 0;
  LabelMap lm;
  FILE * fout;

  cout << "setting up merging.." << endl;
  //makes sure that everything is valid for merging, and sets up our varibles etc.
  unsigned int total = 0;
  offset = validate_merge(infilename, mergefile, lm, option, total);
  fout = fopen("merged.trz", "w");
  fprintf(fout, "TAXA ");
  fprintf(fout, "%s", lm.name(0).c_str());
  for (unsigned int i = 1; i < NUM_TAXA; ++i) { 
    fprintf(fout, ":%s", lm.name(i).c_str());
  }
  fprintf(fout, "\n");
  fprintf(fout, "NTREES %d", NUM_TREES);
  if (HETERO)
    fprintf(fout, " H");
  fprintf(fout, "\n");
  switch (option){
  case 5: bmerge_union(infilename, mergefile, offset, fout);
    break;
  case 6: bmerge_intersection(infilename, mergefile, offset, fout);
    break;
  case 7: bmerge_setDifference(infilename, mergefile, offset, fout);
    break;
  default: 
    cerr << "Invalid option! Exiting..." << endl;  exit(2);
  }
  fclose(fout);
}

unsigned int count_ones(string bipart){
  unsigned int count = 0;
  unsigned int size = bipart.size();
  unsigned int pos = 0;
  unsigned int temp, x, dec;
  //cout << "bipart is: " << bipart << endl;
  //cout << "size is: " << size << endl;
  while (pos < size){
    //cout << "pos: " << pos << ", bipart[" << pos << "]: " << bipart[pos] << endl;
    if (bipart[pos] == 'B'){
      count++;
      pos++;
    }
    else if (bipart[pos] == 'K'){ //special case!
      pos++;
      x = int(bipart[pos]);
      temp = 0;
      dec = 1;
      //cout << "special case!" << endl;
      while ((x < 58) && (x > 47)){
	//cout << "x is: " << x << endl;
	x -= 48;
	//cout << "x is now: " << x << endl;
	//cout << "dec is: " << dec << endl;

	temp = temp*dec + x;
	dec = 10;
	pos++;
	//cout << "temp is: " << temp << endl;
	x = int(bipart[pos]);
	//cout << "next x is: " << x << endl;
      }
      //cout << "out of loop." << endl;
      count+= temp;
      //cout << "count is now: " << count << endl;
    }
    else
      pos++;
  }
  //cout << "count is: " << count << endl;
  //exit(0);
  return count;
}


unsigned int decode(string encoded, unsigned int * found){
  //transform string into another string that is "intermediary"
  //cout << "string to decode is: " << encoded << endl;
  unsigned int ids_length = encoded.size();
  char *tmp = new char[std::numeric_limits <int>::digits10 + 2];
  string myids;

  for (unsigned int a = 0; a < encoded.size()-1; ++a) {
    int trueval = encoded[a];
    int trueval2 = encoded[a+1];
    trueval -= 65;
    trueval2 -= 65;
    if (trueval >= 0) { //we have a number!
      if (trueval < 10) { 
	sprintf(tmp, "%d", trueval);
	myids+=tmp;         
      }
      else{
	myids+="1:";
      }
    }
    else{
      myids+=encoded[a];
    }
    if (trueval2 >=0){ 
      myids+=" ";
    }
  }
  int trueval = encoded[ids_length-1];
  trueval-=65;
  if (trueval >= 0){
    sprintf(tmp, "%d", trueval);
    myids+=tmp;
  }
  else{
    myids+=encoded[ids_length-1];
  }
  delete[] tmp;


  //decode the string my ids:
  stringstream splitarray(myids);
  int pos;

  unsigned place = 0;
  string myid;
  unsigned int current = 0;
  while(splitarray >> myid) { 
    pos = myid.find_first_of(":");
    if (pos == -1){ 
      unsigned int temp = atoi(myid.c_str()); 
      temp += current;            
      found[place] = temp;
      place++;
      current = temp;
    }
    else{
      string id_str = myid.substr(0,pos);
      unsigned int id = atoi(id_str.c_str());
      myid = myid.substr(pos+1);
      unsigned int times = atoi(myid.c_str());
      for (unsigned int i = 0; i < times; ++i) {
	found[place] = current+id;
	current+=id;
	place++;
      }
    }
  }
  return place;
}

string encode(unsigned int * new_vec, unsigned int newlength){
  unsigned int val = new_vec[0];
  unsigned int curr_val  = new_vec[0];
  unsigned int curr_count = 1; 
  stringstream stm, newss;
  string tempval;
  short digit;
  string encoded_result;
  for (unsigned int j = 0; j < newlength-1; ++j)  {
    val = new_vec[j+1] - new_vec[j];
    stm << curr_val;
    stm >> tempval;
    stm.clear();
    digit = tempval[0];
    digit -= 48;
    tempval[0] = tocode[digit];
    if (val != curr_val) { 
      if(curr_count > 1){
	if (tempval.compare("B") == 0) { 
	  //fprintf(fout, "K%d", curr_count);
	  newss << "K" << curr_count;
	}
	else { 
	  //sprintf(fout, "%s:%d", tempval.c_str(), curr_count);
	  newss << tempval << ":" << curr_count;
	}
      }
      else{
	//fprintf(fout, "%s", tempval.c_str());  
	newss << tempval; 
      }  
      curr_val = val;
      curr_count = 1;
    }
    else { 
      curr_count++;
    }     
  } //end for
  stm.str("");
  stm << curr_val;
  stm >> tempval;
  digit = tempval[0];
  digit -= 48;
  tempval[0] = tocode[digit];
  if (curr_count > 1){
    if (tempval.compare("B")==0) { 
      //fprintf(fout, "K%d", curr_count);
      newss << "K" << curr_count;
    }
    else{
      //fprintf(fout, "%s:%d", tempval.c_str(), curr_count);
      newss << tempval << ":" << curr_count;
    }
  }
  else{
    //fprintf(fout, "%s", tempval.c_str());
    newss << tempval;
  }
  encoded_result = newss.str();
  newss.str("");
  return encoded_result;

}

string decode_and_reverse(string encoded, unsigned int nTrees, unsigned int line_size){
 
  //cout << "line_size is: " << line_size << endl;
 
  unsigned int * found, tempfound;
  bool * check = (bool *)malloc(nTrees*sizeof(bool));
  for (unsigned int i = 0; i < nTrees; i++)
    check[i] = 0;
  string new_encode;
  if (line_size != 0){
    found = (unsigned int *)malloc(line_size*sizeof(unsigned int));
    if (found == NULL){
      fprintf(stderr, "error! was not allocated!");
      exit(0);
    }
    
    decode(encoded, found);
 
    //reverse
    tempfound = 0;
    for (unsigned int i =0; i < line_size; i++){
      tempfound = found[i];
      check[tempfound] = true;
    }
  }
 
  unsigned int newlength = nTrees-line_size;
  unsigned int * new_vec = (unsigned int *)malloc(newlength * sizeof(unsigned int));
  unsigned int nx = 0;
  for (unsigned int i = 0; i < nTrees; i++){
    if (!check[i]){
      new_vec[nx] = i;
      nx++;
    }
  }
  free(check);
  //encode
  new_encode = encode(new_vec, newlength);
  //cout << "new_encode is: " << new_encode << endl;
  free(new_vec);
  return new_encode;
}

string correct_offset(string first, string second, unsigned int fsize, unsigned int secsize, unsigned int offset){
  string final_encode = "";
  stringstream ss;
  unsigned int final_size = 0;
  unsigned int * final;
  unsigned * found1, * found2;
  found1 = NULL; 
  found2 = NULL;
  if ( (fsize != 0) && (secsize != 0) ){
    found1 = (unsigned int *)malloc(fsize*sizeof(unsigned int));
    if (found1 == NULL){
      fprintf(stderr, "error! was not allocated!");
      exit(0);
    }
  }
  if (secsize != 0){ 
    found2 = (unsigned int *)malloc(secsize*sizeof(unsigned int));
    if (found2 == NULL){
      fprintf(stderr, "error! was not allocated!");
      exit(0);
    }
  }
  if (!fsize && !secsize)
    return final_encode; //if they are both blank, there is nothing to merge
  else if (fsize && !secsize) // if we have an element merged with a blank, just return the first
    return first;
  else if (!fsize && secsize){
    assert(found2 != NULL);
    decode(second, found2);
    final_size = secsize;
    final = (unsigned int *)malloc(final_size*sizeof(unsigned int));
    if (final == NULL){
      fprintf(stderr, "error! was not allocated!");
      exit(0);
    }
    for (unsigned int i  =0; i < final_size; i++){
      final[i] = found2[i]+offset;
    }
    free(found2);
  }
  else{
    assert(found1 != NULL);
    assert(found2 != NULL);
    //decode first and second into their respective tree id vectors
    decode(first, found1);
    decode(second, found2);

    //merge arrays
    final_size = fsize+secsize;
    final = (unsigned int *)malloc(final_size*sizeof(unsigned int));
    if (final == NULL){
      fprintf(stderr, "error! was not allocated!");
      exit(0);
    }
    for (unsigned int i  =0; i < fsize; i++){
      final[i] = found1[i];
    }
    for (unsigned int i = fsize; i < final_size; i++){
      final[i] = found2[i-fsize]+offset;
    }
    free(found1);
    free(found2);
  }
  //encode
  final_encode = encode(final, final_size);
  //remove after debugging is done
  //ss << first << "|" << second;
  //final_encode = ss.str();
  return final_encode;
}
string merge_biparts(string first, string second, unsigned int offset){
  string merged_string, temp1, temp2, branches, integral; 
  stringstream fs; //final stream
  unsigned int size1, size2;
  //first, need to determine the line types of each
  int pos, pos2;
  //cout << "first is: " << first << endl;
  //cout << "second is: " << second << endl;
  pos  = first.find_first_of(" ");
  if (pos != -1){ //if it is found!
    if (!WEIGHTED){
      WEIGHTED = true;
      cout << "Branch lengths detected!!" << endl;
    }
    pos++;
    pos2 = first.find_first_of("\n");
    branches = first.substr(pos, pos2);
    pos--;
    first = first.substr(0,pos);
    pos2 = branches.find_first_of(" ");
    integral = branches.substr(0, pos2);
    branches = branches.substr(pos2+1);

    //now do it for the second one
    pos = second.find_first_of(" ");
    if (pos == -1){
     cerr << "error! second file does not have branch lengths! exiting..." << endl;
     exit(1);
    }
    pos++;
    pos2 = second.find_first_of("\n");
    temp1 = second.substr(pos, pos2); //this is just temporary
    pos--;
    second = second.substr(0,pos);
    //parse out the integral
    pos2 = temp1.find_first_of(" ");
    integral += temp1.substr(0, pos2);
    temp1 = temp1.substr(pos2+1);
    branches += temp1;
  }
  //cout << "first is: " << first << endl;
  //cout << "second is: " << second << endl;
  //cout << "integral+branches is: " << integral << branches << endl;
  //exit(0);
 stringstream ss(first);
 stringstream ss2(second);
  getline(ss, temp1, ':');
  getline(ss2, temp2, ':');
  //cout << "temp1 is: " << temp1 <<  ", and temp2 is: " << temp2 << endl;
  if (temp1 == temp2){
    string type = temp1;
    //cout << "compact merge!" << endl;
    getline(ss, temp1, ':');
    getline(ss2, temp2, ':');
    size1 = atoi(temp1.c_str());
    size2 = atoi(temp2.c_str());
    size1 += size2; //new size
    getline(ss, temp1, ':');
    getline(ss2, temp2, ':');
    temp1 = correct_offset(temp1, temp2, (size1-size2), size2, offset);
    if (type == "-")
      fs << "-:" << size1 << ":" << temp1;
    else 
      fs << "+:" << size1 << ":" << temp1;
  }
  else {
    unsigned int second_size = NUM_TREES - offset;
    unsigned int new_size = 0;
    if (temp1 == "+"){ //temp1 is + , temp2 is -
      getline(ss2, temp2, ':');
      size2 = atoi(temp2.c_str());
      size2 = second_size - size2;
      getline(ss, temp1, ':');
      size1 = atoi(temp1.c_str());
      new_size = size1 + size2;
      //cout << "new_size is: " << new_size << endl;
      //cout << "NUM_TREES is: " << NUM_TREES << endl;
      if (new_size > NUM_TREES/2){ //temp1 should be -
	new_size = NUM_TREES - new_size;
	getline(ss, temp1, ':');
	temp1 = decode_and_reverse(temp1, offset, size1);
	getline(ss2, temp2, ':');
	temp1 = correct_offset(temp1, temp2, (offset-size1), (second_size-size2), offset);
	fs << "-:" << new_size << ":" << temp1;
      }
      else{//temp2 should be +
	//cout << "do we get here?" << endl;
	getline(ss2, temp2, ':');
	temp2 = decode_and_reverse(temp2, second_size, (second_size-size2));
	getline(ss, temp1, ':');
	temp1 = correct_offset(temp1, temp2, size1, size2, offset);
	fs << "+:" << new_size << ":" << temp1;
      }
    }
    else { //temp1 is -, temp2 is +
      getline(ss, temp1, ':');
      size1 = atoi(temp1.c_str());
      size1 = offset - size1;
      getline(ss2, temp2, ':');
      size2 = atoi(temp2.c_str());
      new_size = size1 + size2;
      if (new_size > NUM_TREES/2){ //temp2 should be -
	new_size = NUM_TREES - new_size;
	getline(ss2, temp2, ':');
	temp2 = decode_and_reverse(temp2, second_size, size2);
	getline(ss, temp1, ':');
	temp1 = correct_offset(temp1, temp2, (offset-size1), (second_size-size2), offset);
	fs << "-:" << new_size << ":" << temp1;
      }
      else{
	getline(ss, temp1, ':'); //temp1 should be +
	temp1 = decode_and_reverse(temp1, offset, (offset-size1));
	getline(ss2, temp2, ':');
	temp1 = correct_offset(temp1, temp2, size1, size2, offset);
	fs << "+:" << new_size << ":" << temp1;
      }
    }
  }
  if (WEIGHTED)
    fs << " " << integral << " " << branches;
  merged_string = fs.str();
  //cout << "merged string is: " << merged_string << endl;
  //exit(0);
  return merged_string;
}

void get_num_biparts(ifstream & fin, ifstream & fin2, unsigned int & bipart, unsigned int & bipart2){
  string str, str2;
  //skip taxa line and trees lines
  getline(fin, str);
  getline(fin, str);
  getline(fin2, str2);
  getline(fin2, str2);

  //get number of bipartitions
  getline(fin, str);
  int pos = str.find_first_of(" ");
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  bipart = atoi(str.c_str());
  getline(fin2, str2);
  pos = str2.find_first_of(" ");
  str2 = str2.substr(pos+1);
  pos = str2.find_first_of("\n");
  str2 = str2.substr(0, pos);
  bipart2 = atoi(str2.c_str());
}

void load_file_into_vectors(ifstream & fin, vector<string> &bipart, vector<string> &file, vector<unsigned int>&count, unsigned int b){
  unsigned int x = 0;
  string mybipart, str;
  unsigned int pos;
  for (unsigned int i  = 0; i < b; i++){
    getline(fin, str);
    pos = str.find_first_of(" ");
    mybipart = str.substr(0,pos);
    x = count_ones(mybipart);
    count.push_back(x);
    bipart.push_back(mybipart);
    str = str.substr(pos+1);
    pos = str.find_first_of("\n");
    file.push_back(str.substr(0, pos));
  }
  fin.close();
}

//to fix: when merging back (adding biparts from second file to merged.trz),
//we need to update the set of tids! the id need to be corrected by offset
void bmerge_union(string infilename, string mergefile, unsigned int offset, FILE * fout){
  cout << "Beginning union!" << endl;
  //count  up bipartitions as we print it to a string
  ifstream fin, fin2;
  unsigned int bipart, bipart2;
  fin.open(infilename.c_str());
  if (!fin){
    cerr << "cannot open " << infilename <<"!" << endl;
    exit(2);
  }
  fin2.open(mergefile.c_str());
  if (!fin){
    cerr << "cannot open " << mergefile << "!" << endl;
    exit(2);
  }
  string str, str2;
  cout << "offset is: " << offset << endl;
  get_num_biparts(fin, fin2, bipart, bipart2);
  cout << "Biparts are: " << bipart << " and " << bipart2 << endl;

  //with union, we include it no matter what
  //however, we need to print it out in the right order of the number of ones.
 
  vector<string> bipart_list, bipart_list2, file1, file2;
  vector<unsigned int> counts1, counts2;
  //load contents of infile and mergefile into vectors
  load_file_into_vectors(fin, bipart_list, file1, counts1, bipart);
  load_file_into_vectors(fin2, bipart_list2, file2, counts2, bipart2);
  
  assert(file1.size() == bipart_list.size());
  assert(file2.size() == bipart_list2.size());

  /*cout << "printing out bipartitions!" << endl;
  cout << "biparts in first file:" << endl;
  for (unsigned int i = 0; i < bipart; i++){
    cout << bipart_list[i] << " " << counts1[i] << endl;
  }
  cout << "biparts in second file:" << endl;
  for (unsigned int i = 0; i < bipart2; i++){
    cout << bipart_list2[i] << " "  << counts2[i] << endl;
    }*/

  //now update the number of bipartitions by comparing the mergefile to infile
  unsigned int new_biparts=  bipart;
  unsigned int a;
  bool found, check;
  found = false;
  check = false;

  unsigned int minpos = 0;
  unsigned int temp;
  bool added;
  bool * printed = (bool*)malloc(bipart2*sizeof(bool));
  if (printed == NULL){
    cout << "cannot allocate!" << endl;
  }
  for (unsigned int i = 0; i < bipart2; i++)
    printed[i] = 0;

  for (unsigned int i  = 0; i < bipart; i++){
    a = counts1[i];
    temp = 0;
    check = false;
    found = false;
    added = false;
    for (unsigned int j = minpos; j < bipart2; j++){
      if (a < counts2[j]){
	if (printed[j] == 0){
	  new_biparts++;
	  printed[j] = 1;
	}
	else{
	  temp += (j+1);
	  continue;
	}
      }
      else if (a == counts2[j]){ 	//actually check
	check = true;
	if (bipart_list[i] == bipart_list2[j]){
	  found = true;
	  added = true;
	  printed[j] = 1;
	  break;
	}
      }
      else{ //a is > counts2[j]
	added = true;
	break;
      }
    }
    //minpos += temp;
  }
  for (unsigned int j = minpos; j < bipart2; j++){
    if (printed[j] == 0){
      new_biparts++;
    }
    else
      printed[j] = 0;
  }
  fprintf(fout, "NBIPARTITIONS %u\n", new_biparts);
  
 //now, begin merge
  string merged_bipart;
  string mytemp;
  int pos2;
  for (unsigned int i  = 0; i < bipart; i++){
    a = counts1[i];
    //cout << "a is: " << a << endl;
    temp = 0;
    check = false;
    found = false;
    added = false;
    //cout << "bipartition we are looking for.." << bipart_list[i] << endl;
    for (unsigned int j = minpos; j < bipart2; j++){
      //cout << "counts2[" <<  j << "] is: " << counts2[j] << endl;
      if (a < counts2[j]){
	//cout << "time to check the other way!" << endl;
	//cout << "minpos is: " << minpos << ", and j is: " << j << endl;
	//cout << "for " << bipart_list2[j] << ": ";
	if (printed[j] == 0){
	  //cout << "printing!" << endl;
	  //need to fix this:
	  mytemp = file2[j];
	  pos2 = mytemp.find_first_of(" ");
	  if (pos2 != -1)
	    merged_bipart = "+:0:  ";
	  else 
	    merged_bipart = "+:0:";
	  merged_bipart = merge_biparts(merged_bipart, file2[j], offset);
	  //cout << "merged_bipart is: " << merged_bipart << endl;
	  fprintf(fout, "%s %s\n", bipart_list2[j].c_str(), merged_bipart.c_str());
	  printed[j] = 1;
	}
	else{
	  //cout << "already found!" << endl;	
	  temp += (j+1);
	  continue;
	}
      }
      else if (a == counts2[j]){ 	//actually check
	check = true;
	if (bipart_list[i] == bipart_list2[j]){
	  //cout << "Found!" << endl;
	  //cout << "need to merge bipartitions " << bipart_list[i] << " and ";
	  //cout << bipart_list2[j] << endl;
	  merged_bipart = merge_biparts(file1[i], file2[j], offset);
	  found = true;
	  fprintf(fout, "%s %s\n",bipart_list[i].c_str(),  merged_bipart.c_str());
	  added = true;
	  printed[j] = 1;
	  break;
	}
      }
      else{ //a is > counts2[j]
	//if ((found == false) && (check == true))
	//cout << "could not find " << bipart_list[i] << ". printing out." << endl;
	fprintf(fout, "%s %s\n", bipart_list[i].c_str(), file1[i].c_str());
	added = true;
	break;
      }
    }
    //minpos += temp;
    if (!added){
      //cout << "could not find.. adding." << endl;
      fprintf(fout, "%s %s\n", bipart_list[i].c_str(), file1[i].c_str()); 
    }   
    //cout << "next element in loop.." << endl;
  }
  for (unsigned int j = minpos; j < bipart2; j++){
    if (printed[j] == 0){
      //cout << "did not find match for " << bipart_list2[j] << endl;
      mytemp = file2[j];
      pos2 = mytemp.find_first_of(" ");
      if (pos2 != -1)
	merged_bipart = "+:0:  ";
      else 
	merged_bipart = "+:0:";
      merged_bipart = merge_biparts(merged_bipart, file2[j], offset);
      fprintf(fout, "%s %s\n", bipart_list2[j].c_str(), merged_bipart.c_str()); 
    }
  }
}



void bmerge_intersection(string infilename, string mergefile, unsigned int offset, FILE * fout){
  cout << "Beginning intersection!" << endl;
  ifstream fin, fin2;
  unsigned int bipart, bipart2;
  fin.open(infilename.c_str());
  if (!fin){
    cerr << "cannot open " << infilename <<"!" << endl;
    exit(2);
  }
  fin2.open(mergefile.c_str());
  if (!fin){
    cerr << "cannot open " << mergefile << "!" << endl;
    exit(2);
  }
  string str, str2;
  cout << "offset is: " << offset << endl;
  get_num_biparts(fin, fin2, bipart, bipart2);
  cout << "Biparts are: " << bipart << " and " << bipart2 << endl;

  //with intersection, we include it only if the bipartition is in both files
  //however, we need to print it out in the right order of the number of ones.
  vector<string> bipart_list, bipart_list2, file1, file2;
  vector<unsigned int> counts1, counts2;
  //load contents of infile and mergefile into vectors
  load_file_into_vectors(fin, bipart_list, file1, counts1, bipart);
  load_file_into_vectors(fin2, bipart_list2, file2, counts2, bipart2);

  assert(file1.size() == bipart_list.size());
  assert(file2.size() == bipart_list2.size());
 
  /*cout << "printing out bipartitions!" << endl;
  cout << "biparts in first file:" << endl;
  for (unsigned int i = 0; i < bipart; i++){
    cout << bipart_list[i] << " " << counts1[i] << endl;
  }
  cout << "biparts in second file:" << endl;
  for (unsigned int i = 0; i < bipart2; i++){
    cout << bipart_list2[i] << " " << counts2[i] << endl;
    }*/

  //now update the number of bipartitions by comparing the mergefile to infile
  unsigned int new_biparts=  0; //it can never be more than the number in the smaller file
  unsigned int a;
  bool found, check;
  for (unsigned int i  =0; i < bipart; i++){
    a = counts1[i];
    found = false;
    check = false;
    //cout << "a is: " << a << " (" << bipart_list[i] << ")" <<   endl;
    for (unsigned int j = 0; j < bipart2; j++){
      //cout << "counts2[" << j << "] is " << counts2[j] << " (" << bipart_list2[j] << ")" << endl;
      if (a < counts2[j])
	continue; //go to next element
      else if (a == counts2[j]){
	//cout << "Checking..." << endl;
	check = true;
	if (bipart_list[i] == bipart_list2[j]){
	  found = true;
	  new_biparts++;
	  break;
	}
      }
      else
	break;
    }
  }

  fprintf(fout, "NBIPARTITIONS %u\n", new_biparts);
  //exit(0);
  //now, begin merge
  string merged_bipart;
  for (unsigned int i  =0; i < bipart; i++){
    a = counts1[i];
    found = false;
    check = false;
    //cout << "a is: " << a << " (" << bipart_list[i] << ")" <<   endl;
    for (unsigned int j = 0; j < bipart2; j++){
      //cout << "counts2[" << j << "] is " << counts2[j] << " (" << bipart_list2[j] << ")" << endl;
      if (a < counts2[j])
	continue; //go to next element
      else if (a == counts2[j]){
	//cout << "Checking..." << endl;
	check = true;
	if (bipart_list[i] == bipart_list2[j]){
	  found = true;
	  //cout << "Found!" << endl;
	  //cout << "need to merge bipartitions " << bipart_list[i] << " and ";
	  //cout << bipart_list2[j] << endl;
	  merged_bipart = merge_biparts(file1[i], file2[j], offset);
	  fprintf(fout, "%s %s\n",bipart_list[i].c_str(),  merged_bipart.c_str());
	  break;
	}
      }
      else
	break;
    }
  }
}
void bmerge_setDifference(string infilename, string mergefile, unsigned int offset, FILE * fout){
  cout << "Beginning set difference!" << endl;
  ifstream fin, fin2;
  unsigned int bipart, bipart2;
  fin.open(infilename.c_str());
  if (!fin){
    cerr << "cannot open " << infilename <<"!" << endl;
    exit(2);
  }
  fin2.open(mergefile.c_str());
  if (!fin){
    cerr << "cannot open " << mergefile << "!" << endl;
    exit(2);
  }
  string str, str2;
  cout << "offset is: " << offset << endl;
  get_num_biparts(fin, fin2, bipart, bipart2);
  cout << "Biparts are: " << bipart << " and " << bipart2 << endl;

  //with set difference, we include it only if the bipartition is in the first, but not in the second
  //however, we need to print it out in the right order of the number of ones.
  vector<string> bipart_list, bipart_list2, file1, file2;
  vector<unsigned int> counts1, counts2;
  //load contents of infile and mergefile into vectors
  load_file_into_vectors(fin, bipart_list, file1, counts1, bipart);
  load_file_into_vectors(fin2, bipart_list2, file2, counts2, bipart2);

  assert(file1.size() == bipart_list.size());
  assert(file2.size() == bipart_list2.size());
 
  /*cout << "printing out bipartitions!" << endl;
  cout << "biparts in first file:" << endl;
  for (unsigned int i = 0; i < bipart; i++){
    cout << bipart_list[i] << " " << counts1[i] << endl;
  }
  cout << "biparts in second file:" << endl;
  for (unsigned int i = 0; i < bipart2; i++){
    cout << bipart_list2[i] << " " << counts2[i] << endl;
    }*/

  //now update the number of bipartitions by comparing the mergefile to infile
  unsigned int new_biparts=  bipart; //it can never be more than the number in the bigger file
  unsigned int a;
  bool found, check,added;
  found = false;
  check = false;
  //basically if the particular bipartition is not found, we decrement new_biparts
  for (unsigned int i  =0; i < bipart; i++){
    a = counts1[i];
    found = false;
    check = false;
    //cout << "a is: " << a << " (" << bipart_list[i] << ")" <<   endl;
    for (unsigned int j = 0; j < bipart2; j++){
      //cout << "counts2[" << j << "] is " << counts2[j] << " (" << bipart_list2[j] << ")" << endl;
      if (a < counts2[j])
	continue; //go to next element
      else if (a == counts2[j]){
	//cout << "Checking..." << endl;
	check = true;
	if (bipart_list[i] == bipart_list2[j]){
	  found = true;
	  new_biparts--;
	  //cout << "Found!" << endl;
	  break;
	}
      }
      else{
	break;
      }
    }
  }

  cout << "new biparts is: " << new_biparts << endl;
  fprintf(fout, "NBIPARTITIONS %u\n", new_biparts);

  //now, begin merge
  string merged_bipart;
  found = false;
  check = false;
  added = false;
  //basically if the particular bipartition is not found, we decrement new_biparts
  for (unsigned int i = 0; i < bipart; i++){
    a = counts1[i];
    found = false;
    check = false;
    added = false;
    cout << "a is: " << a << " (" << bipart_list[i] << ")" <<   endl;
    for (unsigned int j = 0; j < bipart2; j++){
      cout << "counts2[" << j << "] is " << counts2[j] << " (" << bipart_list2[j] << ")" << endl;
      if (a < counts2[j])
	continue; //go to next element
      else if (a == counts2[j]){
	//cout << "Checking..." << endl;
	check = true;
	if (bipart_list[i] == bipart_list2[j]){
	  found = true;
	  //cout << "Found!" << endl;
	  break;
	}
      }
      else{
	cout << "could not find. adding.." << endl;
	fprintf(fout, "%s %s\n", bipart_list[i].c_str(), file1[i].c_str());
	added = true; 
	break;
      }
    }
    if ((check == true) && (found == false) && (added ==false)){
      cout << "could not find. adding.." << endl;
      fprintf(fout, "%s %s\n", bipart_list[i].c_str(), file1[i].c_str()); 
    }
  }
}
