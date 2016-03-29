TreeZip 3.0

Phylogenetic trees are tree structures that depict relationships between organisms. Popular analysis techniques often produce large collections of hypothetical trees, which can be expensive to store. This can also make the exchange of such data also difficult. TreeZip compresses phylogenetic trees based on the shared evolutionary relationships. In our experiments, TreeZip has been shown to be very effective, typically compressing a tree file to less than 2% of its original size. When coupled with standard compression methods such as 7zip, TreeZip can compress a file to less than 1% of its original size. 

We have extended our TreeZip algorithm to include the compression of heterogeneous trees. To the best of our knowledge, TreeZip is the only algorithm capable of identifying relationships in large collections of heterogeneous trees. The TreeZip format is meant to serve as an alternative format to the Newick representation. While Newick is still preferable for single trees, TreeZip is preferable for large collections of trees. First, all trees and their underlying evolutionary relationships are stored uniquely in the TreeZip Compressed (TRZ) format, resulting in a robustness to branch rotations and a smaller file size.  Second, TreeZip enables operations such as the consensus to be computer in mere seconds.

TreeZip's previous home was: treezip.googlecode.com. Please visit that page for previous versions of TreeZip. This site is new, and we will be updating the information here to include usage information.

For detailed command usage, compile the TreeZip program and type: 
./treezip --help

If you use TreeZip 3.0, please cite our paper:

Matthews, SJ. Heterogeneous Compression of Large Collections of Evolutionary Trees. IEEE/ACM Transactions on Computational Biology and Bioinformatics:Special Issues on Software and Databases, to appear. November 2014. DOI: 10.1109/TCBB.2014.2366756 
