#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <unistd.h>
#include "nexus_compress.hh" 
#include "compressfunc.h" 

using namespace std;
bool sortTaxa(const string & a, const string & b){
  return a < b;
}

void gather_lines(ifstream & fin, vector<string> & lines){
  string line = "";
  string myline;
  int pos = -1;
  getline(fin, line);
  while (pos == -1){
    lines.push_back(line);
    getline(fin, line);
    myline = line;
    transform(myline.begin(), myline.end(), myline.begin(), ::toupper);
    pos = myline.find("END;");
    if (pos != -1)
      break;
    pos = myline.find("ENDBLOCK;");
    if (pos != -1)
      break;
    if (fin.eof()){
      cerr << "never found end statement!" << endl;
      exit(3);
    }
  }
}

void taxa_compress(ifstream & fin, bool & etaxa, short & ntaxa){
  string line, taxon, delim, title, type;
  bool found_taxa = false;
  type = "";
  vector<string> taxa;
  while ( (type != "END;") || (type != "ENDBLOCK;") ){
    getline(fin, line);

    stringstream mytype(line);
    mytype >> type;
    transform(type.begin(), type.end(), type.begin(), ::toupper);
    if ( (type == "END;") || (type == "ENDBLOCK;"))
      break;
    if (type == "TITLE"){
      mytype >> title;
    }
    else if (type == "DIMENSIONS"){ //ignored, because we can figure this out later
      continue;
    }
    else if (type == "TAXLABELS"){
      found_taxa = true;
      string taxon;
      while(mytype >> taxon) {
	taxa.push_back(taxon);
      }
    }
    else{
      if (found_taxa){
	taxa.push_back(type);
	while(mytype >> taxon) {
	  taxa.push_back(taxon);
	}
      }
      else{
	cerr << "line is " << line << ". unknown type!" << endl;
	exit(2);
      }
    }
  }
  if (found_taxa != true){
    cerr << "cannot find taxa!" << endl;
    exit(2);
  }
  ofstream fout;
  fout.open(".nexus-taxa", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(2);
  }
  fout << "TAXAL " << title << " ";
  fout << taxa[0];
  for (unsigned int i = 1; i < taxa.size(); i++)
    fout << ":" << taxa[i];
  fout << endl;
  fout.close();
  etaxa = true;
  ntaxa++;
}

void trees_compress(ifstream & fin, bool & etree, short & ntree){
  ofstream fout, fout2;
  vector<string> lines;
  string myline, templine;
  gather_lines(fin, lines); //get the lines
  fout.open(".nexus-trees", ios::app);
  if (!fout){
    cerr << "Cannot open temporary file for writing!" << endl;
    exit(1);
  }
  fout << "TREES " << endl;
  fout2.open(".mytemptre");
  if (!fout2){
    cerr<< "Cannot open temporary tree file for writing!" << endl;
    exit(1);
  }
  int pos = 0;
  int pos2 = 0;
  bool found = false;
  bool ready = false;
  string taxa, taxal, taxon, label, tid, type;
  vector<string> mytaxa;
  taxa = "";
  taxal = "";
  for (unsigned int i = 0; i < lines.size(); i++){ //write lines to temporary file
    myline = lines[i];
    templine = myline;
    pos = templine.find('=');
    if (pos != -1 && !found){
      fout << "TRANSLATE" << endl;
      found = true;
    }
    if (!found){ //until you hit translate, print out everything to the main temporary file
      transform(templine.begin(), templine.end(), templine.begin(), ::toupper);
      pos = templine.find("TRANSLATE");
      if ( pos == -1){
	fout << lines[i] << endl;
      }
      else {
	found = true;
	fout << templine << endl;
	continue;
      }
    }
    else if ( (found) && (!ready) ){ //read the taxa in next
       pos = myline.find_first_of('=');
       if (pos == -1){ //we are still dealing with taxa
	 stringstream ss(myline);
	 ss >> taxon >>  label;
	 pos2 = label.find_first_of(',');
	 if (pos2 == -1){
	   pos2 = label.find_first_of(';');
	   if (pos2 != -1)
	     label = label.substr(0, pos2);
	 }
	 else
	   label = label.substr(0, pos2);
	 mytaxa.push_back(label);
	//ss.str("");
	//if (taxa == ""){
	//taxa += taxon;
	//taxal += label;
	//}
	//else{
	//taxa += ":";
	//  taxa += taxon;
	//taxal += ":";
	//  taxal += label;
	//}
       } //end if pos = -1
       else{
	 ready = true;
	 tid = myline.substr(0, pos);
	 myline = myline.substr(pos+1);
	 pos = myline.find_first_of(']');
	 type = myline.substr(0, pos+1);
	 myline = myline.substr(pos+1);
	 fout2 << myline << endl;
       }	
    }
    else if (ready){
      pos = myline.find_first_of('=');
      myline = myline.substr(pos+1);
      pos = myline.find_first_of(']');
      if (pos != -1)
	myline = myline.substr(pos+1);
      fout2 << myline << endl;
    }
  }
  fout2.close();
  //run compress on the file
  compress(".mytemptre", false, true, mytaxa, false);
  ifstream myfin;
  string aline;
  if (mytaxa.size()){
    sort(mytaxa.begin(), mytaxa.end(), sortTaxa);
    fout << "TAXAL2 " << mytaxa[0];
    for (unsigned int i = 1; i < mytaxa.size(); i++)
      fout << ":" << mytaxa[i]; 
    fout << endl;
  }
  else
    fout << "TAXAL2 " << endl;
  fout << "TID:" << tid << endl;
  fout << "TYPE:" << type << endl;
  myfin.open(".mytemptre.trz"); //now open up the trz file
  if (!myfin){
    cerr << "cannot find the tree file for reading!" << endl;
    exit(1);
  }
  while (!myfin.eof()){
    getline(myfin, aline);
    if (myfin.eof())
      break;
    fout << aline << endl;
  }
  fout.close(); //close the file
  unlink(".mytemptre.trz"); //get rid of temporary files
  unlink(".mytemptre");
  etree = true;
  ntree++;
}

void char_compress(ifstream & fin, bool & echar, short & nchar){
  ofstream fout;
  fout.open(".nexus-chars", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "CHARACTERS " << lines.size() << " " << " 0" << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  echar = true;
  nchar++;
}

void dist_compress(ifstream & fin, bool & edist, short & ndist){
  ofstream fout;
  fout.open(".nexus-dist", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "DIST " << lines.size() << " " << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  edist = true;
  ndist++;
}

void unalign_compress(ifstream & fin, bool & eunalign, short & nunalign){
  ofstream fout;
  fout.open(".nexus-unalign", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "UNALIGNED " << lines.size() << " "  << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  eunalign = true;
  nunalign++;
}

void data_compress(ifstream & fin, bool & edata, short & ndata){
  ofstream fout;
  fout.open(".nexus-data", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "CHARACTERS " << lines.size() << " " << " 1" << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  edata = true;
  ndata++;
}

void set_compress(ifstream & fin, bool & eset, short & nset){
  ofstream fout;
  fout.open(".nexus-sets", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "SETS " << lines.size() << " " << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  eset = true;
  nset++;
}

void codon_compress(ifstream & fin, bool & ecodon, short & ncodon){
  ofstream fout;
  fout.open(".nexus-codons", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "CODONS " << lines.size() << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  ecodon = true;
  ncodon++;
}

void note_compress(ifstream & fin, bool & enote, short & nnote){
  ofstream fout;
  fout.open(".nexus-notes", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "NOTES " << lines.size() << " " << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  enote = true;
  nnote++;
}

void assume_compress(ifstream & fin, bool & eassume, short & nassume){
  ofstream fout;
  fout.open(".nexus-assume", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "ASSUME " << lines.size() << " " << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  eassume = true;
  nassume++;
}

void paup_compress(ifstream & fin, bool & epaup, short & npaup){
  ofstream fout;
  fout.open(".nexus-paup", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "PAUP " << lines.size() << " " << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  epaup = true;
  npaup++;
}

void mb_compress(ifstream & fin, bool & emb, short & nmb){
  ofstream fout;
  fout.open(".nexus-mb", ios::app);
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  vector<string> lines;
  gather_lines(fin, lines);
  fout << "MRBAYES " << lines.size() << " " << endl;
  for (unsigned int i = 0; i < lines.size(); i++)
    fout << lines[i] << endl;
  fout.close();
  emb = true;
}
