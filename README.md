Tree Generator and Visualizer

This Python script processes FASTA files, aligns them using ClustalW, builds phylogenetic trees using IQ-TREE, converts the resulting Newick-formatted tree files into ASCII trees, and saves the ASCII representations to text files. The script uses the Biopython library for tree parsing and rendering.




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
  
  ClustalW: Download and install from the official site.
  
  IQ-TREE: Download and install from the official site.
  





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
The Newick tree files are read, visualized as ASCII trees, and saved to text files.

Running the Script

  The script is executed from the command line and requires a single argument: the path to the input FASTA file.
  
  Example:
  
    python script.py input_file.fasta

Output

Terminal Output: Displays the parsed tree as an ASCII representation.

File Output: Saves the ASCII tree to a file named <base_name>_output.txt.

Example

  Input
  
  A FASTA file, serotype_1.fasta, is processed as described above. The Newick file tree_1.txt contains:
  
  ((A,B),(C,D));
  
  Output
  
  Terminal Output:
  
  Tree successfully parsed:
    _______ A
   |
  _|_______ B
   |
   |_______ C
           |
           D

  File Output:
  The file serotype_1_output.txt will contain the same ASCII representation.


Error Handling

  If the specified input file does not exist, the script will print:
  
  Error: File '<input_file>' not found.

Notes

  Ensure ClustalW and IQ-TREE are installed and accessible from your command line.
  
  Only properly formatted FASTA files are supported.


Acknowledgments

  The Biopython library (https://biopython.org/) is used for tree parsing and visualization.
  
  ClustalW and IQ-TREE are essential tools for alignment and tree construction.

