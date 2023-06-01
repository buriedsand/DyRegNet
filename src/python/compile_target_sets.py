import sys

input_files = sys.argv[1:-2]
output_file = sys.argv[-2]
refseq_file = sys.argv[-1]

# List of transcription regulators (TR) options
with open(refseq_file, "r") as f:
    gene_list = {line.split("\t", 1)[0].split(":", 1)[1] for line in f}

try:
    with open(output_file, "w") as output:
        for input_file in input_files:
            tr_name = input_file.split("/")[-1].split(".")[0]
            try:
                with open(input_file, "r") as f:
                    try:
                        # Skip header line
                        next(f)
                    except StopIteration:
                        # Empty file, skip processing
                        print(f"Empty file: {input_file}")
                        continue
                    for line in f:
                        target_gene = line.strip().split("\t")[0]
                        if target_gene and (target_gene in gene_list):
                            output.write(f"{tr_name},{target_gene}\n")
            except FileNotFoundError:
                print(f"File not found: {input_file}")
except Exception as e:
    print(f"Error writing to output file: {e}")
