
#include <iostream>
#include <fstream>

#ifndef _NEXUS_COMPRESS_HH
#define _NEXUS_COMPRESS_HH

using namespace std;

void taxa_compress(ifstream & fin, bool & etaxa, short & ntaxa);
void trees_compress(ifstream & fin, bool & etree, short & ntree);
void char_compress(ifstream & fin, bool & echar, short & nchar);
void unalign_compress(ifstream & fin, bool & eunalign, short & nunalign);
void data_compress(ifstream & fin, bool & edata, short & ndata);
void set_compress(ifstream & fin, bool & eset, short & nset);
void codon_compress(ifstream & fin, bool & ecodon, short & ncodon);
void note_compress(ifstream & fin, bool & enote, short & nnote);
void assume_compress(ifstream & fin, bool & eassume, short & nassume);
void paup_compress(ifstream & fin, bool & epaup, short & npaup);
void mb_compress(ifstream & fin, bool & emb, short & nmb);
void dist_compress(ifstream & fin, bool & edist, short & ndist);

#endif
