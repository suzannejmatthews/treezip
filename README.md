
## TreeZip 3.0


Phylogenetic trees are tree structures that depict relationships between organisms. Popular analysis techniques often produce large collections of hypothetical trees, which can be expensive to store. This can also make the exchange of such data also difficult. TreeZip compresses phylogenetic trees based on the shared evolutionary relationships. In our experiments, TreeZip has been shown to be very effective, typically compressing a tree file to less than 2% of its original size. When coupled with standard compression methods such as 7zip, TreeZip can compress a file to less than 1% of its original size. 

We have extended our TreeZip algorithm to include the compression of heterogeneous trees. To the best of our knowledge, TreeZip is the only algorithm capable of identifying relationships in large collections of heterogeneous trees. The TreeZip format is meant to serve as an alternative format to the Newick representation. While Newick is still preferable for single trees, TreeZip is preferable for large collections of trees. First, all trees and their underlying evolutionary relationships are stored uniquely in the TreeZip Compressed (TRZ) format, resulting in a robustness to branch rotations and a smaller file size.  Second, TreeZip enables operations such as the consensus to be computer in mere seconds.

TreeZip's previous home was: [treezip.googlecode.com](https://code.google.com/archive/p/treezip/). Please visit that page for previous versions of TreeZip. This site is new, and we will be updating the information here to include usage information.

## Usage
TreeZip is very easy to use. To compress tree files, simply type
```bash
./treezip file.tre
```
where `file.tre` is the tree file (in Newick format) that you wish to compress. This will produce the TRZ file `file.tre.trz`.

### Popular operations
Once in the TRZ format, you can perform several different operations:
* Get the unique trees in a collection: `treezip -du file.tre.trz`
* Get the strict consensus tree associated with a collection: `treezip -d -c s file.tre.trz`
* Get the majority consensus tree associated with a collection: `treezip -d -c m foo.tre.trz`

### Comparing tree collections
The TRZ file is a 1:1 representation of a collection of trees. As such, TreeZip enables you to compare your large tree collections:
* Compute the union between two tree collections `treezip file1.tre.trz -m file2.tre.trz 1`
* Compute the intersection between two tree collections `treezip file1.tre.trz -m file2.tre.trz 2`
* Compute the set difference between two tree collections `treezip foo.tre.trz -m bar.tre.trz 3`

### Other commands
For detailed command usage, compile the TreeZip program and type: 
```bash
./treezip --help
```

## Have an idea for a feature that TreeZip doesn't yet have?
Please feel free to contribute to the TreeZip project! You are always welcome and encouraged to request new features. TreeZip is designed to serve the phylogenetics community! Please also reach out to me if you run into any issues.

## To Cite:
If you use TreeZip 3.0, please cite our paper:

Matthews SJ. Heterogeneous Compression of Large Collections of Evolutionary Trees". IEEE/ACM Transactions on Computational Biology and Bioinformatics, 12(4), pp. 807-814, July-Aug. 2015. DOI: 10.1109/TCBB.2014.2366756 
