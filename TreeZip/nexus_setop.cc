#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include "nexus_setop.hh"
#include "tmerge.h"

using namespace std;

bool lex (string a, string b){
  return (a > b);
}

void get_taxa(ifstream & fin, vector<string> & taxa, unsigned int & count){
  string line, type;
  stringstream ss;
  type = "";
  while (type != "TRANSLATE"){
    getline(fin, line);
    ss.str(line);
    ss >>type;
    count++;
  }
  ss.clear();
  getline(fin, line); //get taxa for this block
  count++;
  ss.str(line);
  ss >> type;
  if (type != "TAXAL2"){
    cerr << "Taxa labels not found!" << endl;
    exit(1);
  }
  string mytaxa, taxon;
  ss >> mytaxa;

  stringstream ss2(mytaxa);
  while(getline(ss2, taxon, ':')) 
    taxa.push_back(taxon);
}

bool subset_check(vector<string> first, vector<string> sec){
  unsigned int num_correspond = 0;
  if (first.size() <= sec.size()){
    for (unsigned int i = 0; i < first.size(); i++){
      if ( binary_search(sec.begin(), sec.end(), first[i]) )
	num_correspond++;
      if (num_correspond > 3)
	break;
    }
    if (num_correspond > 3)
      return true;
    else
      return false;
  }
  else{
    for (unsigned int i = 0; i < sec.size(); i++){
      if ( binary_search(first.begin(), first.end(), sec[i]) )
	num_correspond++;
      if (num_correspond > 3)
	break;
    }
    if (num_correspond > 3)
      return true;
    else
      return false;
  }
}

void find_corresponding_blocks(string infile, string mergefile, vector< pair<unsigned int, unsigned int> > & blocks){
  //returns a list of blocks (based on ID of the block, or line number) that correspond
  //open up infile
  ifstream fin; 
  fin.open(infile.c_str());
  if (!fin){
    cerr << "cannot open file for reading!" << endl;
    exit(1);
  }
  string line, type;
  stringstream ss;
  vector<string> first_taxa, sec_taxa;
  unsigned int count1, count2, dummy, dummy2;
  count1 = 0;
  dummy = 0;
  dummy2= 0;
  while (!fin.eof()){
    getline(fin, line);
    if (fin.eof())
      break;
    count1++;
    ss.str(line);
    ss >> type;
    ss.clear();
    if (type == "TREES"){   //read until you hit a tree block
      get_taxa(fin, first_taxa, dummy);
      ifstream min;
      min.open(mergefile.c_str());  //open up mergefile
      if (!min){
	cerr << "cannot open mergefile for reading!" << endl;
	exit(1);
      }
      string line2, type2;
      stringstream ms;
      count2 = 0;
      while (!min.eof()){
	getline(min, line2);
	if (min.eof())
	  break;
	count2++;
	ms.str(line2);
	ms >> type2;
	ms.clear();
	//cout << "type2 is: " << type2 << endl;
	if (type2 == "TREES"){ //read until you hit a tree block
	  get_taxa(min, sec_taxa, dummy2);  //get the taxa
	  bool correspond = subset_check(first_taxa, sec_taxa);
	  sec_taxa.clear();
	  if (correspond){
	    //cout << "blocks at: " << count1 << " and " << count2 << " correspond." << endl;
	    blocks.push_back(make_pair(count1,count2)); //create pair and add to vector
	    count2 += dummy2;
	    dummy2 = 0;
	  }
	}
      }
      min.close();
    }
    first_taxa.clear();  //clear taxa1 vector
    count1 += dummy;
    dummy = 0;
  }
  fin.close();
}

 void extract_trz(ifstream & fin, string outfile, vector<string> taxa_labels){
  ofstream fout; 
  fout.open(outfile.c_str());
  if (!fout){
    cerr << "Cannot open temporary file for writing!" << endl;
    exit(1);
  }
  string line, type;
  type = "";
  while (type != "TAXA"){
    getline(fin, line);
    stringstream ss(line);
    ss >> type;
  }
  fout << "TAXA " << taxa_labels[0];
  for (unsigned int i  = 1; i < taxa_labels.size(); i++)
    fout << ":" << taxa_labels[i];
  fout << endl;
  //fout << line << endl; //output taxa line
  for (unsigned int i  = 0; i < 3; i++){
    getline(fin, line);
    fout << line << endl;
  }
  stringstream ss(line);
  string bipart;
  ss >> type;
  if (type != "NBIPARTITIONS"){
    cerr << "MALFORMED TRZ FILE!" << endl;
    exit(1);
  }
  ss >> bipart;
  unsigned int nbipart = atoi(bipart.c_str());
  for (unsigned int i = 0; i < nbipart; i++){
    getline(fin, line);
    fout << line << endl;
  }
  getline(fin, line); //get duplicates line
  ss.clear();
  ss.str(line);
  ss >> type;
  if (type != "DUPLICATES"){
    cerr << "MALFORMED TRZ FILE!" << endl;
    exit(1);
  }
  fout << line << endl;
  string dups;
  ss >> dups;
  unsigned int ndups = atoi(dups.c_str());
  for (unsigned int i = 0; i < ndups; i++){
    getline(fin, line);
    fout << line << endl;
  }
  fout.close(); 
}
  
void perform_op_on_blocks(string infile, string mergefile,  pair<unsigned int, unsigned int> & block, unsigned int option, ofstream & fout){
  ifstream fin, fin2;
  fin.open(infile.c_str());   //open infile, mergefile
  fin2.open(mergefile.c_str());
  if (!fin || !fin2){
    cerr << "cannot open files for reading!" << endl;
    exit(1);
  }
  string line, line2;
  unsigned int count, count2;
  count = block.first;
  count2 = block.second;
  for (unsigned int i = 0; i < count; i++)  //go to corresponding blocks
    getline(fin, line);
  for (unsigned int i = 0; i < count2; i++)
    getline(fin2, line2);
  vector<string> first_taxa, sec_taxa;
  unsigned int dummy = 0;
  get_taxa(fin, first_taxa, dummy);
  get_taxa(fin2, sec_taxa, dummy);
  dummy = 0;
  extract_trz(fin, ".tmpfirst.trz", first_taxa);
  extract_trz(fin2, ".tmpsec.trz", sec_taxa);
  //do setop on it 
  if (!null_case(".tmpfirst.trz", ".tmpsec.trz", ".tmpmerge.trz", option)){
    if (option < 4){
      setup_tmerge(".tmpfirst.trz", ".tmpsec.trz", ".tmpmerge.trz",  option, true);
      //fprintf(stderr, "out of setup_tmerge function\n");
    }
    else{
      cerr << "Invalid option! Exiting..." << endl;
      exit(1);
    }
  }
  fin.close();
  fin2.close();
  fin.open(".tmpmerge.trz");
  if (!fin){
    cerr << "cannot open temporary file for reading!" << endl;
    exit(1);
  }
  getline(fin, line); //get taxa line
  string type;
  stringstream ss(line);
  ss >> type;
  if (type != "TAXA"){
    cerr << "Cannot find taxa line!" << endl;
    exit(1);
  }
  fout << "TREES" << endl;

  string taxa,taxon;
  ss >> taxa;
  ss.str("");
  switch(option){
  case 1:
    fout << "    TITLE Union_TREES_Block;" << endl;
    fout << "    TRANSLATE" << endl;
    fout << "TAXAL2 " << taxa << endl;
    fout << "TID:    TREE union" << endl;
    break;
  case 2:
    fout << "    TITLE Intersection_TREES_Block;" << endl;
    fout << "    TRANSLATE" << endl;
    fout << "TAXAL2 " << taxa << endl;
    fout << "TID:    TREE intersection" << endl;
    break;
  case 3: 
    fout << "    TITLE SetDifference_TREES_Block;" << endl;
    fout << "    TRANSLATE" << endl;
    fout << "TAXAL2 " << taxa << endl;
    fout << "TID:    TREE setDifference" << endl;
    break;
  default:
    cerr << "Invalid option!" << endl;
    exit(1);
  }
  fout << "TYPE: " << endl;
  stringstream ss2(taxa);
  unsigned int taxaCount = 0;
  while(getline(ss2, taxon, ':'))
    taxaCount++;
  fout << "TAXA 1";
  for (unsigned int i  =1; i < taxaCount; i++)
    fout << ":" << i+1;
  fout << endl;
  while (!fin.eof()){
    getline(fin, line); 
    if (fin.eof())
      break;
    fout << line << endl;  //write the contents to fout
  }
//close both files and remove temporary files
  fin.close();
  unlink(".tmpfirst.trz");
  unlink(".tmpsec.trz");
  unlink(".tmpmerge.trz");
}
 
