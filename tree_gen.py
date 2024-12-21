import argparse
import os
from Bio import Phylo
from io import StringIO

def convert(input_file):
    try:
        tree = Phylo.read(input_file, "newick")
        
        print("Tree successfully parsed:")
        Phylo.draw_ascii(tree)
        
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_output.txt"
        
        ascii_tree = StringIO()
        Phylo.draw_ascii(tree, file=ascii_tree)
        ascii_tree_output = ascii_tree.getvalue()
        
        with open(output_file, 'w') as file:
            file.write(ascii_tree_output)
        
        return tree
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a Newick file and convert it to a tree.")
    parser.add_argument("input_file", help="Path to the input Newick file.")
    
    args = parser.parse_args()
    
    convert(args.input_file)
