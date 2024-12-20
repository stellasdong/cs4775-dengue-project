def process_genomic_file(input_file, output_file):
    try:
        with open(input_file, 'r') as infile:
            lines = infile.readlines()

        with open(output_file, 'w') as outfile:
            for line in lines:
                line = line.strip()  # Remove leading/trailing whitespaces and line breaks
                if line.startswith('>'):
                    # Write header lines as is
                    outfile.write(line + '\n')
                else:
                    # Write sequence lines without line breaks and quotation marks
                    outfile.write(line.replace('"', ''))

        print(f"Processed file saved to {output_file}")

    except FileNotFoundError:
        print(f"Error: The file {input_file} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
# Replace 'input.txt' with the path to your input file and 'output.fasta' with your desired output path
process_genomic_file('sequence-3.txt', 'output.txt')
