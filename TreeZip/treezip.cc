//*****************************************************************/
/*
This is TreeZip, a compression software for phylogenetic trees. 
It is based on HashRF and HashCS, developed by SeungJin Sul.

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
//#include <bitset>

#include "label-map.hh"
#include "hashfunc.hh"
#include "hash.hh"
#include <tclap/CmdLine.h>
				    //#include <tclap/CmdLine.h>
#include "SCTree.h"
#include "parsing.h"
#include "global.h"
#include "compressfunc.h"
#include "bmerge.h"
#include "tmerge.h"
#include "nexus_setop.hh"

// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;

void copyFile(string from, string to){
  //from and to are file names with paths
  ifstream fromstr(from.c_str(), fstream::binary);
  ofstream tostr(to.c_str(), fstream::trunc|fstream::binary);
  tostr << fromstr.rdbuf ();
}

int main(int argc, char** argv) {
  string infilename;
  string mergefile;
  //int OFFSET = 0;
  bool DECOMPRESS = false;
  bool ALTERNATE = false;
  bool MERGING = false;
  bool QUIET = false;
  bool bUbid = false; // for counting the number of unique bipartitions  
  int random_decompress = 0;
  bool branch = true;
  bool calc_unique_trees = false; 
  unsigned char consensus = 'n';
  unsigned int option = 0;
  // TCLAP
  try {

    // Define the command line object.
    string 	helpMsg  = "treezip [-d] <input file>\n";

    helpMsg += "Input file: \n";
    helpMsg += "   The current version of TreeZip only supports the Newick format.\n";

    helpMsg += "Example of Newick tree: \n";
    helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
    helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n\n";

    helpMsg += "Examples: \n";
    helpMsg += "  compression: treezip foo.tre\n";
    helpMsg += "  decompression: treezip -d foo.tre.trz\n";
    helpMsg += "  decompression (unique trees): treezip -du foo.tre.trz\n";
    helpMsg += "  decompression (strict consensus): treezip -d -c s foo.tre.trz\n";
    helpMsg += "  decompression (majority consensus): treezip -d -c m foo.tre.trz\n"; 
    helpMsg += "  set op (union): treezip foo.tre.trz -m bar.tre.trz 1\n";    
    helpMsg += "  set op (intersection): treezip foo.tre.trz -m bar.tre.trz 2\n";    
    helpMsg += "  set op (set difference): treezip foo.tre.trz -m bar.tre.trz 3\n";    
    TCLAP::CmdLine cmd(helpMsg, ' ', "3.0.0");

    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "Input tree file name"  );
    cmd.add( fnameArg );

    TCLAP::UnlabeledValueArg<int>  numtreeArg( "numtree", "number of trees", false, 0, "Number of trees"  );
    cmd.add( numtreeArg );

    TCLAP::ValueArg<char> cArg("c", "consensus", "consensus", false, 'n', "consensu(s");
    cmd.add( cArg );

    TCLAP::ValueArg<unsigned int> pArg("p", "pvalue", "percent to compress for hash table", false, 50, "percent to compress");
    cmd.add( pArg );

    TCLAP::SwitchArg weightedSwitch("w", "weighted", "weighted trees", false);
    cmd.add( weightedSwitch );    

    TCLAP::SwitchArg nArg("n", "nobranch", "decompress with no branches", false);
    cmd.add( nArg );    
    TCLAP::SwitchArg dArg("d", "decompress", "decompression", false);
    cmd.add( dArg );
    TCLAP::SwitchArg aArg("a", "alternative", "alternate approach", false);
    cmd.add( aArg );
    TCLAP::SwitchArg uArg("u", "unique", "unique trees", false);
    cmd.add( uArg );
    TCLAP::SwitchArg qArg("q", "quiet", "quiet treezip", false);
    cmd.add( qArg );
    TCLAP::ValueArg<int> seedArg("s", "seedvalue", "user specified seed value", false, 1000, "user specified seed value");
    cmd.add( seedArg );
    TCLAP::ValueArg<int> rArg("r", "random", "random percentage of trees to decompress", false, 0, "random number to decompress");
    cmd.add( rArg );
    TCLAP::ValueArg<string>  mergeArg( "m", "mfile", "m file", false, "", "optional second trz file for merging"  );
    cmd.add( mergeArg );
    TCLAP::ValueArg<unsigned int> oArg("o", "operation", "set operation for merging", false, 0, "operation choice");
    cmd.add( oArg );

    cmd.parse( argc, argv );

    NUM_TREES = numtreeArg.getValue();
    DECOMPRESS = dArg.getValue();
    ALTERNATE = aArg.getValue();
    QUIET = qArg.getValue();
    infilename = fnameArg.getValue();

    if (cArg.getValue())
      consensus = cArg.getValue();
    
    if (pArg.getValue()){
      P = pArg.getValue();
      if (P > 100) { 
	cerr << "error: P value is greater than 100!" << endl;
	return 1;
      }
    }

    P = 50; //enforced thresholding

    if (weightedSwitch.getValue())
      WEIGHTED = true;

    if (seedArg.getValue()) { 
      NEWSEED = seedArg.getValue();
    }

    if (rArg.getValue())
      random_decompress = rArg.getValue();
   
    if (nArg.getValue())
      branch=false;
 
    if (uArg.getValue())
      calc_unique_trees = true;

    //merging options
    mergefile = mergeArg.getValue();
    if (mergefile != ""){
      MERGING = true;
      fprintf(stderr, "Compare mode on!\n");
      fprintf(stderr, "TRZ file to compare: %s\n", mergefile.c_str());
      if (oArg.getValue()){
	option = oArg.getValue();
	if (option == 4){ //equivalence option
	  test_equivalence(infilename, mergefile);
	  return 0;
	}
	else if (option < 4)
	  fprintf(stderr, "selected option is: %d\n", option);
	else{
	  cerr << "Selected option must be from [0..3]!" << endl;
	  exit(1);
	}
      }
      else{
	cerr << "no option specified! quitting!" << endl;
	exit(1);
      }
      if (MERGING && DECOMPRESS){
	cerr << "cannot merge and decompress at the same time! Please choose on or the other.." << endl;
	exit(1);
      }
    }
    if (oArg.getValue() &&!MERGING){
      cerr << "cannot specify a set operation without specifying a file to merge! quitting!" << endl;
      exit(1);
    }
  } catch (TCLAP::ArgException &e) 
  { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }


  //cout << "random_decompress is: " << random_decompress << endl;
  //cout << "branch is: " << branch << endl;
  // cout << "calculate unique trees is: " << calc_unique_trees << endl;
  //cout << "consensus is: " << consensus << endl;
  //##################################################3
  if (is_nexus(infilename) && !DECOMPRESS && !MERGING){
    cout << "The file is a NEXUS file! Following NEXUS parsing procedure..." << endl;
    parse_nexus(infilename);
    cout << "File compressed successfully." << endl;
    return 0;
  } 
  if (WEIGHTED){
    fprintf(stderr, "Results are Weighted!\n");
  }
  if (MERGING){
    fprintf(stderr, "Beginning merging procedure!\n");
    if (is_nexus(infilename) && is_nexus(mergefile)){
      // new procedure goes here
      vector< pair<unsigned int, unsigned int> > my_blocks;
      //determine if blocks correspond
      find_corresponding_blocks(infilename, mergefile, my_blocks); 
      cout << "Found " << my_blocks.size() << " corresponding block(s)." << endl;
      if (my_blocks.size() == 0){
	cerr << "Error! No corresponding tree blocks in input files!" << endl;
	cerr<< "Exiting.." << endl;
	exit(1);
      }
      ofstream fout; 
      fout.open("merged.trz"); //create file to store results
      if (!fout){
	cerr << "Cannot open file for outputting set ops!" << endl;
	exit(1);
      }
      fout << "#NEXUS" << endl;
      for (unsigned int i  =0; i < my_blocks.size(); i++){
	pair< unsigned int, unsigned int> block = my_blocks[i];
	perform_op_on_blocks(infilename, mergefile, block, option, fout);
      }
      fout.close();
    }
    else if (is_nexus(infilename) || is_nexus(mergefile)){
      cerr << "Error! For now, both files must either be PHYLIP or NEXUS files!" << endl;
    }
    else{  //perform merge regularly
      string outfile = "";
      if (null_case(infilename, mergefile, outfile, option))
	return 0;
      if (option < 4){
	setup_tmerge(infilename, mergefile, outfile, option, false);
	//fprintf(stderr, "out of setup_tmerge function\n");
      }
      else{
	if (option < 8){
	  //setup_bmerge(infilename, mergefile, option);
	  fprintf(stderr, "out of setup_bmerge function\n");
	}
	else{
	  cerr << "Invalid option! Exiting..." << endl;
	  exit(1);
	}
      }
    }
  }
  else if (DECOMPRESS) { 
    if (!QUIET)
      fprintf(stderr, "Decompressing file...\n");
    if (is_nexus(infilename))
      decompress_nexus(infilename, random_decompress, branch, calc_unique_trees, consensus, QUIET);
    else
      decompress(infilename, random_decompress, branch, calc_unique_trees, consensus, QUIET);
    if (!QUIET)
      fprintf(stderr, "File decompressed successfully\n");
  }
  else { //compress
    if (!QUIET)
      fprintf(stderr, "Compressing file...\n");
    vector<string> dummy;
    compress(infilename, bUbid, QUIET, dummy, false);
    if (!QUIET)
      fprintf(stderr, "File Compressed successfully\n");
  }
  
  return 0;
}
