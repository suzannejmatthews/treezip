comdir = ${top_srcdir}/Common
randir = ${comdir}/RandomLib
tclap  = ${top_srcdir}/tclap

AM_CXXFLAGS = -Wall -O2 -funroll-loops -fexpensive-optimizations

AM_CXXFLAGS += -I${top_srcdir}
AM_CXXFLAGS += -I${comdir}

# real program stuff
bin_PROGRAMS = treezip

#common sources
treezip_SOURCES =  
treezip_SOURCES += ${comdir}/hashfunc.hh ${comdir}/randomc.h ${comdir}/hashfunc.cc
treezip_SOURCES += ${comdir}/label-map.hh ${comdir}/label-map.cc
treezip_SOURCES += ${comdir}/hash.hh ${comdir}/hash.cc

# common sources: Newick parser
treezip_SOURCES += ${comdir}/newick.h ${comdir}/newick.c

# common sources: Mersenne Twister
treezip_SOURCES += ${comdir}/Random.cpp
treezip_SOURCES += ${randir}/Random.hpp
treezip_SOURCES += ${randir}/RandomAlgorithm.hpp
treezip_SOURCES += ${randir}/RandomCanonical.hpp
treezip_SOURCES += ${randir}/RandomEngine.hpp
treezip_SOURCES += ${randir}/RandomMixer.hpp
treezip_SOURCES += ${randir}/RandomPower2.hpp
treezip_SOURCES += ${randir}/RandomSeed.hpp
treezip_SOURCES += ${randir}/RandomType.hpp

# common sources: Command Line parser
treezip_SOURCES += ${tclap}/Arg.h
treezip_SOURCES += ${tclap}/ArgException.h
treezip_SOURCES += ${tclap}/CmdLine.h
treezip_SOURCES += ${tclap}/CmdLineInterface.h
treezip_SOURCES += ${tclap}/CmdLineOutput.h
treezip_SOURCES += ${tclap}/Constraint.h
treezip_SOURCES += ${tclap}/HelpVisitor.h
treezip_SOURCES += ${tclap}/IgnoreRestVisitor.h
treezip_SOURCES += ${tclap}/MultiSwitchArg.h
treezip_SOURCES += ${tclap}/MultiArg.h
treezip_SOURCES += ${tclap}/OptionalUnlabeledTracker.h
treezip_SOURCES += ${tclap}/StdOutput.h
treezip_SOURCES += ${tclap}/SwitchArg.h
treezip_SOURCES += ${tclap}/UnlabeledMultiArg.h
treezip_SOURCES += ${tclap}/UnlabeledValueArg.h
treezip_SOURCES += ${tclap}/ValueArg.h
treezip_SOURCES += ${tclap}/ValuesConstraint.h
treezip_SOURCES += ${tclap}/VersionVisitor.h
treezip_SOURCES += ${tclap}/Visitor.h
treezip_SOURCES += ${tclap}/XorHandler.h

# common sources: Tree Building
treezip_SOURCES += ${comdir}/SCNode.h ${comdir}/SCNode.cc
treezip_SOURCES += ${comdir}/SCTree.h ${comdir}/SCTree.cc
treezip_SOURCES += ${comdir}/buildtree.h ${comdir}/buildtree.cc

treezip_SOURCES += ${comdir}/global.h ${comdir}/global.cc
treezip_SOURCES += ${comdir}/parsing.h ${comdir}/parsing.cc

# TreeZip files
treezip_SOURCES += nexus_compress.hh nexus_compress.cc
treezip_SOURCES += nexus_decompress.hh nexus_decompress.cc
treezip_SOURCES += nexus_setop.hh nexus_setop.cc
treezip_SOURCES += compressfunc.h compressfunc.cc
treezip_SOURCES += bmerge.h bmerge.cc
treezip_SOURCES += tmerge.h tmerge.cc
treezip_SOURCES += treezip.cc
