#include <iostream>
#include <fstream>

using namespace std;
#ifndef _NEXUS_DECOMPRESS_HH
#define _NEXUS_DECOMPRESS_HH

void taxa_decompress(ifstream & fin, ofstream & fout, string line);
void trees_decompress(ifstream & fin, ofstream & fout);
void char_decompress(ifstream & fin, ofstream & fout, unsigned int lines);
void dist_decompress(ifstream & fin, ofstream & fout, unsigned int lines);
void common_decompress(ifstream & fin, ofstream & fout, unsigned int lines);

#endif

