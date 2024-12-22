import argparse
import os
from Bio import Phylo
import matplotlib.pyplot as plt

def convert(input_file):
    try:
        tree = Phylo.read(input_file, "newick")
        print("Tree successfully parsed.")
        
        base_name = os.path.splitext(input_file)[0]
        output_image = f"{base_name}_tree.png"
        
        plt.figure(figsize=(10, 8))
        Phylo.draw(tree, do_show=False)
        plt.savefig(output_image, format="png")
        plt.close()
        
        print(f"Tree image saved to '{output_image}'")
        return tree
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a Newick file and create an image of the tree.")
    parser.add_argument("input_file", help="Path to the input Newick file.")
    
    args = parser.parse_args()
    
    convert(args.input_file)
