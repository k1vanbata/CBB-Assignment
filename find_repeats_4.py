import re

import gzip

from multiprocessing import Pool

from Bio import SeqIO

from Bio.Seq import Seq


# Function to read DNA sequence from a gzipped FASTA file using Biopython

def read_dna_sequence_from_gz(file_path, num_genes):

    sequence = ""
    count = 0

    with gzip.open(file_path, 'rt') as file:

        for record in SeqIO.parse(file, "fasta"):
            if count < num_genes:
                sequence += str(record.seq)
                count += 1
            else:
                break

    return sequence



# Function to get the reverse complement of a DNA sequence using Biopython

def reverse_complement(seq):

    return str(Seq(seq).reverse_complement())



# Optimized function to find inverted repeats with no spacer, counting with dictionary

def find_inverted_repeats_no_spacer(sequence, min_half_length=4):

    repeat_counts = {}

    seq_len = len(sequence)

    seen_substrings = {}  # Stores seen substrings and their reverse complements



    # Iterate over possible half lengths

    for half_len in range(min_half_length, seq_len // 2 + 1):

        for i in range(seq_len - 2 * half_len + 1):

            first_half = sequence[i:i + half_len]

            second_half = sequence[i + half_len:i + 2 * half_len]



            # Calculate reverse complement only once and store in a dictionary

            if first_half not in seen_substrings:

                seen_substrings[first_half] = reverse_complement(first_half)



            # Check if second_half matches the reverse complement of first_half

            if second_half == seen_substrings[first_half]:

                full_repeat = f"{first_half} | {second_half}"

                repeat_counts[full_repeat] = repeat_counts.get(full_repeat, 0) + 1



    return repeat_counts



# Optimized function to find inverted repeats with spacers, counting with dictionary

def find_inverted_repeats_with_spacer(sequence, min_half_length=4, min_spacer=1, max_spacer=10):

    repeat_counts = {}

    seq_len = len(sequence)

    seen_substrings = {}



    # Iterate over possible spacer lengths

    for spacer_len in range(min_spacer, max_spacer + 1):

        for half_len in range(min_half_length, (seq_len - spacer_len) // 2 + 1):

            for i in range(seq_len - 2 * half_len - spacer_len + 1):

                first_half = sequence[i:i + half_len]

                second_half = sequence[i + half_len + spacer_len:i + half_len + spacer_len + half_len]



                # Precompute reverse complement for each unique first_half

                if first_half not in seen_substrings:

                    seen_substrings[first_half] = reverse_complement(first_half)



                # Check if second_half matches the reverse complement of first_half

                if second_half == seen_substrings[first_half]:

                    full_repeat = f"{first_half} | {second_half}"

                    repeat_counts[full_repeat] = repeat_counts.get(full_repeat, 0) + 1



    return repeat_counts



# Function to analyze a chunk of DNA sequence for both types of repeats

def analyze_chunk(sequence_chunk):

    no_spacer_repeats = find_inverted_repeats_no_spacer(sequence_chunk)

    spacer_repeats = find_inverted_repeats_with_spacer(sequence_chunk)

    return no_spacer_repeats, spacer_repeats



# Use multiprocessing to analyze the DNA sequence in parallel

def process_dna_file_parallel(file_path, output_file, num_genes, num_processes=4):

    # Step 1: Read the DNA sequence from the gzipped FASTA file

    dna_sequence = read_dna_sequence_from_gz(file_path, num_genes)

    print(f"DNA sequence length: {len(dna_sequence)}")



    # Step 2: Divide sequence into chunks for parallel processing

    chunk_size = len(dna_sequence) // num_processes

    sequence_chunks = [dna_sequence[i:i + chunk_size] for i in range(0, len(dna_sequence), chunk_size)]



    # Step 3: Process each chunk in parallel

    with Pool(num_processes) as pool:

        results = pool.map(analyze_chunk, sequence_chunks)



    # Combine results from all processes

    no_spacer_counts = {}

    spacer_counts = {}

    for no_spacer, spacer in results:

        # Combine no spacer repeats

        for repeat, count in no_spacer.items():

            no_spacer_counts[repeat] = no_spacer_counts.get(repeat, 0) + count

        # Combine spacer repeats

        for repeat, count in spacer.items():

            spacer_counts[repeat] = spacer_counts.get(repeat, 0) + count



    # Write results to a text file

    with open(output_file, 'w') as f:

        f.write("Inverted repeats with no spacer:\n")

        for repeat, count in no_spacer_counts.items():

            f.write(f"Repeat: {repeat}, Count: {count}\n")



        f.write("\nInverted repeats with spacer:\n")

        for repeat, count in spacer_counts.items():

            f.write(f"Repeat: {repeat}, Count: {count}\n")



    print(f"Results written to {output_file}")


# Main function to run the parallel processing

if __name__ == '__main__':

    file_path = 'C:\Users\k1van\Box\bme580\coding\coding\coding_shuffled_001.gz'  # Update this with the actual file path

    output_file = 'C:\Users\k1van\Box\bme580\coding\coding\repeat_analysis_output.txt'  # Update with desired output file path

    num_genes = 5

    process_dna_file_parallel(file_path, output_file, num_genes, num_processes=4)