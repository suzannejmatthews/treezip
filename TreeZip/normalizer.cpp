//	Copyright (C) 2007-2008 Mark T. Holder
//
//	This file is part of NCL (Nexus Class Library).
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc.,
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

/*******************************************************************************
 * This file contains code for 4 executables:
 *		* NEXUSnormalizer - writes a NEXUS version of the file with consistent
 *			ordering of blocks and commands. Ideally 2 equivalent files will
 *			produce the same normalized output. This version of tthe program is
 *			less ambitious. The goal is to be able to run (for any valid NEXUS
 *			in.nex file):
 *				$ NEXUSnormalizer in.nex > outOrig.nex
 *				$ NEXUSnormalizer outOrig.nex > outSecond.nex
 *				$ diff outOrig.nex outSecond.nex
 *			and find no differences.
 * All other code has to do with reading command line arguments and other
 * 	user-interface concerns.
 */
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"
#include <cassert>
using namespace std;

bool gQuietMode = false;
std::ofstream gCommonFileStream;
std::ostream * gCommonOstream = 0L;

bool gAltNexus = false;

void writeAsNexus(PublicNexusReader & nexusReader, ostream & os);

long gStrictLevel = 2;
bool gUnderscoresToSpaces = false;
bool gValidateInternals = true;
bool gTreesViaInMemoryStruct = true;
long gInterleaveLen = -1;
bool blocksReadInValidation = false;
bool gSuppressingNameTranslationFile = false;

enum ProcessActionsEnum {
  REPORT_BLOCKS,
  OUTPUT_NORMALIZED_NEXUS,
  OUTPUT_ANY_FORMAT,
  OUTPUT_NEXML,
  VALIDATE_ONLY
};

void processContent(PublicNexusReader & nexusReader, ostream *os, ProcessActionsEnum currentAction);
MultiFormatReader * instantiateReader();

#	if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
MultiFormatReader * gNexusReader = NULL;
#	endif


void writeAsNexus(PublicNexusReader & nexusReader, ostream & os) {
  BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();
  os << "#NEXUS\n";
  for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt) {
    NxsBlock * b = *bIt;
    if (b)
      b->WriteAsNexus(os);
  }
}



////////////////////////////////////////////////////////////////////////////////
// Takes NxsReader that has successfully read a file, and processes the
//	information stored in the reader.
//
// The caller is responsibel for calling DeleteBlocksFromFactories() to clean
//	up (if the reader uses the factory API).
////////////////////////////////////////////////////////////////////////////////
void processContent(PublicNexusReader & nexusReader, ostream *os, ProcessActionsEnum currentAction) {
  BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();
  writeAsNexus(nexusReader, *os);
}

MultiFormatReader * instantiateReader() {
  MultiFormatReader * nexusReader = new MultiFormatReader(-1, NxsReader::WARNINGS_TO_STDERR);
  if (gQuietMode)
    nexusReader->SetWarningOutputLevel(NxsReader::SKIPPING_CONTENT_WARNING);
  if (gStrictLevel != 2)
    nexusReader->SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
  if (gUnderscoresToSpaces) 
    nexusReader->SetCoerceUnderscoresToSpaces(true);
  NxsCharactersBlock * charsB = nexusReader->GetCharactersBlockTemplate();
  NxsDataBlock * dataB = nexusReader->GetDataBlockTemplate();
  charsB->SetAllowAugmentingOfSequenceSymbols(true);
  dataB->SetAllowAugmentingOfSequenceSymbols(true);
  if (gInterleaveLen > 0) {
    assert(charsB);
    charsB->SetWriteInterleaveLen(gInterleaveLen);
    dataB->SetWriteInterleaveLen(gInterleaveLen);
  }
  
  NxsTreesBlock * treesB = nexusReader->GetTreesBlockTemplate();
  assert(treesB);
  if (gStrictLevel < 2)
    treesB->SetAllowImplicitNames(true);
  treesB->SetWriteFromNodeEdgeDataStructure(gTreesViaInMemoryStruct);
  treesB->setValidateInternalNodeLabels(gValidateInternals);
  if (gAltNexus)
    treesB->setWriteTranslateTable(false);
  if (gStrictLevel < 2) {
    NxsStoreTokensBlockReader *storerB =  nexusReader->GetUnknownBlockTemplate();
    assert(storerB);
    storerB->SetTolerateEOFInBlock(true);
  }
  nexusReader->conversionOutputRecord.addNumbersToDisambiguateNames = true;
  
  if (gSuppressingNameTranslationFile)
    nexusReader->conversionOutputRecord.writeNameTranslationFile = false;
  return nexusReader;
}

////////////////////////////////////////////////////////////////////////////////
// Creates a NxsReader, and tries to read the file `filename`.  If the
//	read succeeds, then processContent will be called.
//	\returns 0 on success
////////////////////////////////////////////////////////////////////////////////
#if !defined(MULTIFILE_NEXUS_READER) ||  !MULTIFILE_NEXUS_READER
int processFilepath(
		    const char * filename, // name of the file to be read
		    ostream * os, // output stream to use (NULL for no output). Not that cerr is used to report errors.
		    MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
		    ProcessActionsEnum currentAction) // enum that is passed on to processContent to indicate what should be done with the content of the file.
#else
  int processFilepath(
		      const char * filename, // name of the file to be read
		      ostream * , // output stream to use (NULL for no output). Not that cerr is used to report errors.
		      MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
		      ProcessActionsEnum ) // enum that is passed on to processContent to indicate what should be done with the content of the file.
#endif 
{
  assert(filename);
  try
    {
      MultiFormatReader * nexusReader;
#		if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
      nexusReader = gNexusReader;
#		else
      nexusReader = instantiateReader();
#		endif
	    
      if (!gQuietMode)
	cerr << "Executing" << endl;
      try {
	nexusReader->DemoteBlocks();
	nexusReader->ReadFilepath(filename, fmt);
#			if !defined(MULTIFILE_NEXUS_READER) ||  !MULTIFILE_NEXUS_READER
	                processContent(*nexusReader, os, currentAction);
#			endif
      }
      catch(...) {
	nexusReader->DeleteBlocksFromFactories();
#	if ! defined(MULTIFILE_NEXUS_READER) || !MULTIFILE_NEXUS_READER
	delete nexusReader;
#	endif
	throw;
      }
#     if !defined(MULTIFILE_NEXUS_READER) ||  !MULTIFILE_NEXUS_READER
      nexusReader->DeleteBlocksFromFactories();
      delete nexusReader;
#     endif
      return 0;
    }
  catch (const NxsException &x)
    {
      cerr << "Error:\n " << x.msg << endl;
      if (x.line > 0 || x.pos > 0)
	cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
      return 2;
    }
}

/*! \returns 0 on success*/
int readFilepathAsNEXUS(const char *filename, const char *outfilename) {
  MultiFormatReader::DataFormatType fmt(MultiFormatReader::NEXUS_FORMAT);
  ProcessActionsEnum currentAction = OUTPUT_NORMALIZED_NEXUS;
  ofstream fout;
  fout.open(outfilename);
  if (!fout){
    cerr << "Cannot open output file for writing!" << endl;
    return 1;
  } 
  try {
    ostream * outStream;
    outStream = &fout;
    return processFilepath(filename, outStream, fmt, currentAction);
    fout.close();
  }
  catch (...) {
    cerr << "Normalizing of " << filename << " failed (with an exception)" << endl;
    return 1;
  }
}

/*! \returns 0 on success*/


const char * gExeName = "NEXUSnormalizer";
	

int do_main(int argc, char *argv[]) {
  NxsReader::setNCLCatchesSignals(true);
# if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
  gNexusReader = instantiateReader();
  gNexusReader->cullIdenticalTaxaBlocks(true);
# endif
  bool readfile = false;
  if (argc!= 3){
    cerr << "Usage: ./ncl <input> < output>" << endl;
    cerr << "Normalizes NEXUS files." << endl;
    return 1;
  }
  const char * filepath = argv[1];
  const char * outfile = argv[2];
  readfile = true;
  int rc = readFilepathAsNEXUS(filepath, outfile);
  if (rc != 0)
    return rc;    

# if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
  if (gNexusReader) {
    processContent(*gNexusReader, &std::cout, OUTPUT_NORMALIZED_NEXUS);
    gNexusReader->DeleteBlocksFromFactories();
    delete gNexusReader;
  }
# endif
  return 0;
}

/*int main(int argc, char *argv[]) {
  int rc = do_main(argc, argv);
  if (gCommonOstream != 0L && gCommonOstream == &gCommonFileStream)
    gCommonFileStream.close();
  return rc;
  }*/
