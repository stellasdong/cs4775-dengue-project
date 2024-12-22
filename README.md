Comparison of Dengue Variants via Phylogenetic Organization

This Python repo processes FASTA files, aligns them using ClustalW, builds phylogenetic trees using IQ-TREE, converts the resulting Newick-formatted tree files into ASCII trees, and saves the ASCII representations to text files. The script uses the Biopython library for tree parsing and rendering.




Dependencies:

1. Python 3.6 or higher

2. Biopython library

3. ClustalW software

4. IQ-TREE software



   


Installation:

Biopython

  To install Biopython, run:

    pip install biopython
  
ClustalW and IQ-TREE
  
  ClustalW: Download MegaX software and run ClustalW alignment.
  
  IQ-TREE: Download and install from the official site. In order to use iqtree2 command in your command-line, the IQ-Tree bin file must be added to your path.
  





Usage:

Input Files

  The script processes FASTA files and expects them to be properly formatted. Example FASTA files include:
  
    serotype_1.fasta
  
    serotype_3.fasta
  
    serotype_3_2017.fasta
  
    serotype_3_2018.fasta






Workflow:

Alignment with ClustalW
  The FASTA files are aligned using ClustalW to create multiple sequence alignments.

Tree Construction with IQ-TREE
  The aligned files are processed with IQ-TREE using the following command:

    iqtree2 -s <aligned_file.fasta> -m MFP

  This generates Newick-formatted tree files. Example output files:

    tree_1.txt
  
    tree_3.txt






Newick to ASCII Tree Conversion:
The Newick tree files are read, visualized as ASCII trees using BioPython, and saved to text files.

Running the Script

  The script is executed from the command line and requires a single argument: the path to the input FASTA file.
  
  Example:
  
    python script.py input_file.txt

  To replicate our trees, use tree_1.txt or tree_3.txt instead of input_file.txt.

Output

Terminal Output: Displays the parsed tree as an ASCII representation.

File Output: Saves the ASCII tree to a file named <base_name>_output.txt.


Error Handling

  If the specified input file does not exist, the script will print:
  
  Error: File '<input_file>' not found.


Acknowledgments

  The Biopython library (https://biopython.org/) is used for tree parsing and visualization.
  
  ClustalW and IQ-TREE are essential tools for alignment and tree construction.

