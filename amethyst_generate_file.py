import numpy as np
import pandas as pd
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import shutil
import os
import sys
import platform
import csv
import seaborn

#CONFIGURATION

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def get_user_config():
    print("\n--- Configuration Setup ---")
    fixed_seq = input("Enter Fixed Sequence (nucleotide sequence before ORF, e.g., gggatcc): ").strip()
    
    while True:
        try:
            aa_length = int(input("Enter Amino Acid Length (e.g., 155): ").strip())
            break
        except ValueError:
            print("Please enter a valid integer.")
            
    reference_seq = input("Enter Reference Sequence (one letter amino acid sequence): ").strip()
    
    return fixed_seq, aa_length, reference_seq

def generate_filenames(base_name):
    if base_name.lower().endswith('.fastq.gz'):
        base = base_name[:-9]
    elif base_name.lower().endswith('.fastq'):
        base = base_name[:-6]
    elif base_name.lower().endswith('.txt'):
        base = base_name[:-4]
    else:
        base = base_name
        
    return {
        'translated': f"{base}_translated.txt",
        'one_mut': f"{base}_1_one_mut.txt",
        'counts': f"{base}_1_counts.txt",
        'matrix': f"{base}_1_mutation_count_matrix.csv"
    }

#STEP 0: NGmerge (Paired-End Merging)
def run_ngmerge(r1, r2, output_file):
    print(f"--- Merging Paired-End Reads ---")
    print(f"R1: {r1}")
    print(f"R2: {r2}")
    print(f"Output: {output_file}")

    #Check for NGmerge executable
    ngmerge_cmd = None
    
    #Check PyInstaller Bundle Path (Onefile or Onedir)
    if getattr(sys, 'frozen', False):
        # If onefile, _MEIPASS is the temp dir.
        # If onedir, dirname(sys.executable) is the app dir (Contents/MacOS).
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(sys.executable))

        possible_paths = [
            os.path.join(base_path, "NGmerge", "NGmerge"),
            os.path.join(base_path, "_internal", "NGmerge", "NGmerge"), # Sometimes _internal is used
        ]

        if platform.system() == 'Darwin':
             possible_paths.append(os.path.join(base_path, "..", "Frameworks", "NGmerge", "NGmerge"))

        for p in possible_paths:
            p = os.path.abspath(p) # Normalize path
            if os.path.isfile(p) and os.access(p, os.X_OK):
                ngmerge_cmd = p
                break

    #Check Input PATH or CWD
    if ngmerge_cmd is None:
        if shutil.which("NGmerge"):
            ngmerge_cmd = "NGmerge"
        elif os.path.isfile("./NGmerge") and os.access("./NGmerge", os.X_OK):
            ngmerge_cmd = "./NGmerge"
        elif os.path.isfile("./NGmerge/NGmerge") and os.access("./NGmerge/NGmerge", os.X_OK):
            ngmerge_cmd = "./NGmerge/NGmerge"
    
    if ngmerge_cmd is None:
         debug_info = f"Base: {getattr(sys, '_MEIPASS', os.path.dirname(sys.executable) if getattr(sys, 'frozen', False) else 'Local')}"
         raise Exception(f"NGmerge executable not found. Searched standard paths and bundle locations. {debug_info}")

    #Construct command: NGmerge -1 <R1> -2 <R2> -o <Output> -y
    cmd = [ngmerge_cmd, "-1", r1, "-2", r2, "-o", output_file, "-y"]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("Merge successful.")
        return True
    except subprocess.CalledProcessError as e:
        error_msg = f"NGmerge failed with exit code {e.returncode}.\nStderr: {e.stderr}\nStdout: {e.stdout}"
        print(error_msg)
        raise Exception(error_msg)
    except Exception as e:
        print(f"Unexpected error: {e}")
        raise e

#STEP 1: Process FASTQ to Translated Text
def step1_fastq_to_translated(input_fastq, output_fasta, fixed_seq, aa_length):
    print(f"--- Processing {input_fastq} ---")
    nt_length = aa_length * 3
    count = 0
    
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fastq, "fastq"):
            seq = str(record.seq).lower()

            pos = seq.find(fixed_seq)
            if pos == -1:
                continue

            start = pos + len(fixed_seq)
            if start + nt_length > len(seq):
                continue

            coding_seq = seq[start : start + nt_length]
            protein_seq = Seq(coding_seq).translate(to_stop=False)
            out_f.write(f"{protein_seq}\n")
            count += 1
            
    print(f"Extracted {count} sequences.")
    return True

#STEP 2: Filter and Count Mutations
def step2_filter_and_count(input_file, output_seqs, output_counts_file, reference_seq):
    print(f"--- Counting mutations ---")
    
    mutation_counts = Counter()
    total_reads = 0
    
    with open(input_file, "r") as f, open(output_seqs, "w") as fout:
        for line_num, line in enumerate(f, 1):
            read_seq = line.strip()
            total_reads += 1
            
            if len(read_seq) != len(reference_seq):
                continue

            mutations = sum(1 for a, b in zip(read_seq, reference_seq) if a != b)
            mutation_counts[mutations] += 1

            if mutations == 1:
                fout.write(read_seq + "\n")

    with open(output_counts_file, "w") as f:
        for mut_count in sorted(mutation_counts.keys()):
            f.write(f"{mutation_counts[mut_count]}\n")

    print(f"Total reads processed: {total_reads}")
    return True

# STEP 3: Generate Mutation Matrix
def step3_generate_matrix(input_file, output_csv, reference_seq):
    print(f"--- Generating matrix ---")
    
    mutation_matrix = pd.DataFrame(
        np.nan,
        index=range(len(reference_seq)),
        columns=AA_LIST
    )
    
    flag = 0
    with open(input_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            seq = line.strip()
            if len(seq) != len(reference_seq):
                continue

            for i, (ref_aa, mut_aa) in enumerate(zip(reference_seq, seq)):
                if ref_aa != mut_aa:
                    if mut_aa not in AA_LIST:
                        flag += 1
                        continue

                    if pd.isna(mutation_matrix.at[i, mut_aa]):
                        mutation_matrix.at[i, mut_aa] = 1
                    else:
                        mutation_matrix.at[i, mut_aa] += 1
                    break

    mutation_matrix.to_csv(output_csv)
    return True

#MAIN
if __name__ == "__main__":
    print("Welcome to Amethyst (Paired-End Edition with NGmerge).")
    fixed_seq, aa_length, reference_seq = get_user_config()
    
    print("\n--- Input Files ---")
    print("Please provide Paired-End FASTQ files (Read 1 and Read 2).")
    print("Supports .fastq and .fastq.gz files.")

    r1 = input("Enter path to Control Read 1 FASTQ/FASTQ.GZ: ").strip()
    r2 = input("Enter path to Control Read 2 FASTQ/FASTQ.GZ: ").strip()

    merged_fastq = "merged_control.fastq"
    
    if not run_ngmerge(r1, r2, merged_fastq):
        print("Aborting due to merge failure.")
        exit(1)
        
    m_files = generate_filenames(merged_fastq)
    
    if step1_fastq_to_translated(merged_fastq, m_files['translated'], fixed_seq, aa_length):
        if step2_filter_and_count(m_files['translated'], m_files['one_mut'], m_files['counts'], reference_seq):
            step3_generate_matrix(m_files['one_mut'], m_files['matrix'], reference_seq)

    files_to_delete = [merged_fastq, m_files['translated'], m_files['one_mut']]
    for f_path in files_to_delete:
        try:
            if os.path.exists(f_path):
                os.remove(f_path)
        except Exception as e:
            print(f"Could not delete {f_path}: {e}")
