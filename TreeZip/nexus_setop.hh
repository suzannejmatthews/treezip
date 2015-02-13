#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unistd.h>

using namespace std;

#ifndef NEXUS_SETOP_HH
#define NEXUS_SETOP_HH

void get_taxa(ifstream & fin, vector<string> & taxa);
void find_corresponding_blocks(string infile, string mergefile, vector< pair<unsigned int, unsigned int> > & blocks);
void perform_op_on_blocks(string infile, string mergefile,  pair<unsigned int, unsigned int> & block, unsigned int option, ofstream & fout);

#endif
