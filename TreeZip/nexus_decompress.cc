#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include "nexus_decompress.hh"
#include "compressfunc.h"

using namespace std;

void taxa_decompress(ifstream & fin, ofstream & fout, string myline){
  int pos = myline.find_first_of(" ");
  string title = myline.substr(0, pos);
  myline = myline.substr(pos+1);
  if (title != "")
    fout << "TITLE " << title << endl; 
  stringstream ss(myline);
  string taxon;
  vector<string> taxa;
  while (getline(ss, taxon, ':')){
    taxa.push_back(taxon);
  }
  fout << "DIMENSIONS NTAX=" << taxa.size() << ";" << endl;
  fout << "TAXLABELS" << endl;
  for (unsigned int x = 0; x < taxa.size(); x++){
    fout <<"    " <<  taxa[x] << endl;
  }
}

void trees_decompress(ifstream & fin, ofstream & fout){
  string myline = "";
  int pos;
  string type = "";

  //create temporary file 
  while (type != "TRANSLATE"){
    getline(fin, myline);
    stringstream ss(myline);
    ss >> type;
    if (type == "TRANSLATE")
      break;
    if (fin.eof()){
      cerr << "reached end of file!" << endl;
      return;
    }
    fout << myline << endl;
  }
  if (type == "TRANSLATE"){
    getline(fin, myline);//taxal2
    pos = myline.find_first_of(" ");
    type = myline.substr(0, pos);
    if (type != "TAXAL2"){
      cerr << "TAXAL2 field not found!" << endl;
      exit(1);
    }
    myline = myline.substr(pos+1);
    stringstream ss(myline);
    if (myline != ""){
      string taxon;
      vector<string> mytaxa;
      fout << "    TRANSLATE" << endl; //output translate only if there are taxa
      while (getline(ss, taxon, ':')){
	//getline(s2, id, ':'); //parse each of the taxa names and the numbers
	mytaxa.push_back(taxon);
	//fout << "        " << id << "= " << taxon << endl;    //output to file
      } 
      for (unsigned int i = 0; i < mytaxa.size()-1; i++)
	fout << "        " << i+1 << " " << mytaxa[i] << "," << endl;    //output to file
      fout << "        " << mytaxa.size() << " " <<  mytaxa[mytaxa.size()-1] << ";" << endl;    //output to file
    }
    string tempt;
    getline(fin, myline); //tid
    pos = myline.find_first_of(":");
    tempt = myline.substr(0, pos);
    ss.clear();
    ss.str(tempt);
    ss >> type;
    if (type != "TID"){
      cerr << "TID field not found!" << endl;
      exit(1);
    }
    myline = myline.substr(pos+1);
    string tid = myline;
    getline(fin, myline); //type
    pos = myline.find_first_of(":");
    tempt = myline.substr(0, pos);
    ss.clear();
    ss.str(tempt);
    ss >> type;
    if (type != "TYPE"){
      cerr << "TYPE field not found!" << endl;
      exit(1);
    }
    myline = myline.substr(pos+1);
    string my_type = myline;
    ofstream ftmp;
    ftmp.open(".mytemp.trz");     //open temporary file
    if (!ftmp){
      cerr << "cannot open file!" << endl;
      exit(1);
    }
    getline(fin, myline); //taxa line
    int pos = myline.find_first_of(" ");
    string type = myline.substr(0, pos);
    if (type != "TAXA"){
      cerr << "MALFORMED TRZ FILE!" << endl;
      exit(1);
    }
    ftmp << myline << endl; //write out taxa
    getline(fin, myline); 
    pos = myline.find_first_of(" ");
    type = myline.substr(0, pos);
    if (type != "NTREES"){
      cerr << "MALFORMED TRZ FILE!" << endl;
      exit(1);
    }
    ftmp << myline << endl;
    myline = myline.substr(pos+1); //write out ntrees
    unsigned int ntrees= atoi(myline.c_str()); //store ntrees
    getline(fin, myline);
    pos = myline.find_first_of(" ");
    type = myline.substr(0, pos);
    if (type != "UNIQUE_T"){
      cerr << "MALFORMED TRZ FILE!" << endl;
      exit(1);
    }
    ftmp << myline << endl; //write out unique
    getline(fin, myline);   
    //write out nibiparts
    pos = myline.find_first_of(" ");
    type = myline.substr(0, pos);
    if (type != "NBIPARTITIONS"){
      cerr << "MALFORMED TRZ FILE!" << endl;
      exit(1);
    }
    ftmp << myline << endl;
    myline = myline.substr(pos+1);
    unsigned int nbiparts = atoi(myline.c_str());
    for (unsigned int x = 0; x < nbiparts; x++){
      getline(fin, myline);     //write out the next nbipart lines
      ftmp << myline << endl;
    }
    getline(fin, myline);   
    pos = myline.find_first_of(" ");
    type = myline.substr(0, pos);
    if (type != "DUPLICATES"){
      cerr << "MALFORMED TRZ FILE!" << endl;
      exit(1);
    }
    ftmp << myline << endl;   //write out duplicates
    myline = myline.substr(pos+1);
    unsigned int ndups = atoi(myline.c_str());
    for (unsigned x = 0; x < ndups; x++){
      getline(fin, myline);   //write out n duplicates
      ftmp << myline << endl;
    }
    ftmp.close();  //close temporary file
    decompress(".mytemp.trz", 0, true, 0, 'n', true); //decompress file
    ifstream f2tmp;
    f2tmp.open(".mytemp.e");
    if (!f2tmp){
      cerr << "Cannot open file for reading/writing trees!" << endl;
      exit(1);
    }
    for (unsigned int x  = 0; x < ntrees; x++){ //for every tree
      getline(f2tmp, myline);
      fout << tid << "=" << my_type << myline << endl; //write out tid and tree
    }   
    //get rid of temporary files
    unlink(".mytemp.e");
    unlink(".mytemp.trz");    
  }
  else{
    cerr << "not trees found!" << endl;
    return;
  }
}

void char_decompress(ifstream & fin, ofstream & fout, unsigned int lines){
 string myline;
  for (unsigned int i = 0; i < lines; i++){
    getline(fin, myline);
    fout << myline << endl;
  }
}

void dist_decompress(ifstream & fin, ofstream & fout, unsigned int lines){
 string myline;
  for (unsigned int i = 0; i < lines; i++){
    getline(fin, myline);
    fout << myline << endl;
  }
}

void common_decompress(ifstream & fin, ofstream & fout, unsigned int lines){
  string myline;
  for (unsigned int i = 0; i < lines; i++){
    getline(fin, myline);
    fout << myline << endl;
  }
}

