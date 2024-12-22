## Comparison of Dengue Variants via Phylogenetic Organization

This Python repo processes FASTA files, aligns them using ClustalW, builds phylogenetic trees using IQ-TREE, converts the resulting Newick-formatted tree files into ASCII trees, and saves the ASCII representations to text files. The script uses the Biopython library for tree parsing and rendering.




### Dependencies:

1. Python 3.6 or higher

2. Biopython library

3. ClustalW software

4. IQ-TREE software



   


### Installation:

**Biopython**

  To install Biopython, run:

    pip install biopython

**IQ-TREE:**
  Download and install from the [official site](http://www.iqtree.org/). In order to use iqtree2 command in your command-line, the IQ-Tree bin file must be added to your path.
  
**ClustalW:** (Optional) 
  Download MegaX software and run ClustalW alignment. MegaX's ClustalW alignment tool was used to align the sequences, all fasta files in this repository were aligned  already for our experiment. 
  


### Usage:

**Running IQ-Tree**

  The script processes FASTA files and runs the IQTree software with the FASTA file. The provided FASTA files include:

    all_serotypes_aligned.fasta

    serotype_1.fasta
  
    serotype_3.fasta
  
    serotype_3_2017.fasta
  
    serotype_3_2018.fasta

  These are all the fasta files that were run in our experiment.

  To run IQTree, run:

    python iq_tree_run.py fasta_name.fasta


  Replace `fasta_name.fasta` with any of the FASTA files above. This will generate several reports and output files from IQ-Tree. The file with the .iqtree extension is the full report plus text representation of the tree. the file with the .bionj is the Newick format of the generated tree.

**Visualizing Trees**

  We utilized multiple tools to visualize our resulting trees. We used an online software [iTOL](https://itol.embl.de/upload.cgi) to generate our large trees. This was especially useful for our serotype 1 dataset. We also used BioPhython to generate readable trees for our smaller datasets.

  The script is executed from the command line and requires a single argument: the path to the input tree file.
  
  Run:
    ```
    python tree_gen.py nwk_file.txt
    ```

  To replicate our trees, use `tree_3.txt`, `tree_3_2017.txt` or `tree_3_2018.txt` instead of `nwk_file.txt`. Or, use the .bionj file output from IQTree. This should generate a .png file of the same name that visualizes the tree.



Acknowledgments

  The Biopython library (https://biopython.org/) is used for tree parsing and visualization.
  
  IQ-TREE is an essential tools for tree construction.

