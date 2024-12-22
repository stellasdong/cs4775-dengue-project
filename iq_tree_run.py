import subprocess
import argparse

parser = argparse.ArgumentParser(description="Run IQ-TREE with a specified FASTA file.")
parser.add_argument(
    "fasta_file", 
    type=str, 
)

args = parser.parse_args()

command = [
    "iqtree2", 
    "-s", args.fasta_file, 
    "-m", "MFP", 
]

try:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    print("IQ-TREE finished successfully.")
    print("Standard Output:", result.stdout)
    print("Standard Error:", result.stderr)
except subprocess.CalledProcessError as e:
    print("Error occurred:", e)
    print("Standard Output:", e.stdout)
    print("Standard Error:", e.stderr)
