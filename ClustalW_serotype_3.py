import subprocess

# CLustalW requirements
command = [
    "clustalw2", 
    "-INFILE=serotype_3.fasta"
]

try:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    print("ClustalW finished successfully.")
    print("Standard Output:", result.stdout)
    print("Standard Error:", result.stderr)
except subprocess.CalledProcessError as e:
    print("Error occurred:", e)
    print("Standard Output:", e.stdout)
    print("Standard Error:", e.stderr)
