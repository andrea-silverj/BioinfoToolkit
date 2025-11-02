#!/usr/bin/env python3

# Splitting a FASTA file into multiple files with n sequences and k folders

__author__ = 'Andrea Silverj'
__version__='2.0'
__date__='1 November 2025'

# This script contains content that was partially generated using a Generative AI model.
# Tools used: Claude Sonnet v4.5, Google Gemini 2.5 Pro, Grok 4 Expert
# Date of generation/assistance: December 31, 2025
# Purpose of AI usage: Generating boilerplate code, debugging assistance, creating initial function structure.
# All AI-generated content has been reviewed and edited by a human to ensure quality and relevance.

from Bio import SeqIO
from Bio.SeqIO import FastaIO
import argparse as ap
import os
import sys
import math


def read_args(args):

    parser = ap.ArgumentParser(description='# Splitting a FASTA file into multiple files with n sequences and k folders #\n')

    required = parser.add_argument_group('required arguments')

    required.add_argument('-f',
                          required=True,
                          metavar='fasta_file',
                          nargs='?',
                          help="FASTA file",
                          type=str)

    required.add_argument('-o',
                          required=True,
                          metavar='output_folder',
                          nargs='?',
                          help='name of the output folder',
                          type=str)

    optional = parser.add_argument_group('optional arguments')

    optional.add_argument('-n',
                          required=False,
                          metavar='NseqsPerFile',
                          default=1,
                          help='number of sequences per output file (default: 1)',
                          type=int)

    optional.add_argument('-k',
                          required=False,
                          metavar='Nfolders',
                          default=1,
                          help='number of folders to distribute output files (default: 1)',
                          type=int)

    optional.add_argument('-w',
                          '--wrap',
                          required=False,
                          metavar='wrap_length',
                          default=None,
                          help='wrap FASTA sequences at specified length (default: single-line, no wrapping)',
                          type=int)

    return vars(parser.parse_args())


def check_missing_files(args):

    if not os.path.isfile(args['f']):
        print("Error: file '"+args['f']+"' is not accessible!")
        sys.exit(1)


def create_output_structure(outname, k_folders):
    """Create output directory and subdirectories if needed"""

    # Create main output directory
    if not os.path.exists(outname):
        os.makedirs(outname)

    # Create subdirectories if k > 1
    folder_paths = []
    if k_folders > 1:
        for i in range(k_folders):
            folder_path = os.path.join(outname, f"folder_{i+1}")
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            folder_paths.append(folder_path)
    else:
        folder_paths.append(outname)

    return folder_paths


def write_fasta_record(file_handle, seq_id, sequence, wrap_length=None):
    """Write a FASTA record with optional wrapping"""

    file_handle.write(f">{seq_id}\n")

    if wrap_length is None:
        # Single line output
        file_handle.write(f"{str(sequence)}\n")
    else:
        # Wrapped output
        seq_str = str(sequence)
        for i in range(0, len(seq_str), wrap_length):
            file_handle.write(f"{seq_str[i:i+wrap_length]}\n")


def split_fasta(fasta_file, n_seqs, k_folders, outname, wrap_length):
    """Main function to split FASTA file"""

    # 1. Read all sequences FIRST
    print(f"Reading sequences from {fasta_file}...\n")
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    total_seqs = len(sequences)
    print(f"Total sequences found: {total_seqs}")

    if total_seqs == 0:
        print("Warning: No sequences found in the input file!")
        # Create only the main output directory, but no subfolders
        if not os.path.exists(outname):
            os.makedirs(outname)
        print(f"Created empty output directory: {outname}")
        return # Stop processing

    # 2. Calculate number of output files
    n_files = math.ceil(total_seqs / n_seqs)

    # 3. Perform validation checks *before* creating folders
    if n_seqs > total_seqs:
        print(f"Warning: Requested {n_seqs} sequences per file, but only {total_seqs} sequence(s) available.")
        print(f"Will create 1 file with all {total_seqs} sequence(s).")

    if k_folders > n_files:
        print(f"\nError: Requested {k_folders} folder(s), but only {n_files} file(s) will be created.")
        print(f"This would result in {k_folders - n_files} empty folder(s).")
        print(f"Maximum folders allowed: {n_files}")
        print(f"\nSuggestions:")
        print(f"  - Reduce number of folders to {n_files} or fewer")
        print(f"  - Reduce sequences per file (current: {n_seqs})")
        sys.exit(1) # Exit *before* any folders are created

    # 4. Create output directory structure *after* validation has passed
    folder_paths = create_output_structure(outname, k_folders)

    print(f"Creating {n_files} output file(s) with {n_seqs} sequence(s) per file")
    print(f"Distributing files across {k_folders} folder(s)\n")

    if wrap_length:
        print(f"Wrapping sequences at {wrap_length} characters per line")
    else:
        print("Using single-line format (no wrapping)")

    # Calculate files per folder for sequential even distribution
    base = n_files // k_folders
    extra = n_files % k_folders
    files_per_folder = [base + 1 if i < extra else base for i in range(k_folders)]

    # Split sequences into files
    global_file_counter = 1  # Global counter for continuous numbering
    seq_counter = 0
    current_batch = []
    current_folder = 0
    current_folder_file_count = 0

    for i, record in enumerate(sequences):
        current_batch.append(record)
        seq_counter += 1

        # Write file when batch is complete or last sequence is reached
        if seq_counter == n_seqs or i == total_seqs - 1:
            # Get current folder
            output_folder = folder_paths[current_folder]

            # Create filename
            if n_seqs == 1:
                rec_for_name = current_batch[0]
                # Use first token of ID, sanitize for filename
                base_id = rec_for_name.id.split()[0]
                base_id = base_id.replace("/", "_").replace("\\", "_").replace(" ", "_").replace("|", "_")
                output_file = os.path.join(output_folder, f"{base_id}.fa")
                # Avoid overwrite if duplicate IDs
                counter = 1
                original_file = output_file
                while os.path.exists(output_file):
                    output_file = os.path.join(output_folder, f"{base_id}_{counter}.fa")
                    counter += 1
            else:
                output_file = os.path.join(output_folder, f"sequences_{global_file_counter}.fa")

            # Write sequences to file
            with open(output_file, "w") as f_out:
                for rec in current_batch:
                    write_fasta_record(f_out, rec.description, rec.seq, wrap_length)

            print(f"  Written: {output_file} ({len(current_batch)} sequence(s))")

            # Update folder counters
            current_folder_file_count += 1
            if current_folder_file_count == files_per_folder[current_folder]:
                current_folder += 1
                current_folder_file_count = 0

            # Reset batch counters
            global_file_counter += 1
            seq_counter = 0
            current_batch = []

    print(f"\nDone! Created {global_file_counter-1} file(s) in {outname}/")


if __name__ == '__main__':

    args = read_args(sys.argv)

    check_missing_files(args)

    fasta = args['f']
    n_seqs = args['n']
    k_folders = args['k']
    outname = args['o']
    wrap_length = args['wrap']

    # Validate arguments
    if n_seqs < 1:
        print("Error: Number of sequences per file must be at least 1!")
        sys.exit(1)

    if k_folders < 1:
        print("Error: Number of folders must be at least 1!")
        sys.exit(1)

    if wrap_length is not None and wrap_length < 1:
        print("Error: Wrap length must be at least 1!")
        sys.exit(1)

    # Run the splitting function
    try:
        split_fasta(fasta, n_seqs, k_folders, outname, wrap_length)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
