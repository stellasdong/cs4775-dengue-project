import subprocess

# IQ Tree Inputs
command = [
    "iqtree2", 
    "-s", "all_serotypes.fasta", 
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
