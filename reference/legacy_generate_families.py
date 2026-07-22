#!/usr/bin/env python3

import argparse
import os
import logging
import time
import re
import typing
import shutil

import pandas as pd
import pyfamsa
import pyhmmer
import pytrimal
import numpy as np

# --- Constants ----------------------------------------------------------------

# use a singleton protein alphabet
ALPHABET = pyhmmer.easel.Alphabet.amino()

# logger
logging.basicConfig(
    filename="execution.log",
    filemode="a",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# --- Data handling utilities --------------------------------------------------

class Sequence(typing.NamedTuple):
    id: str
    seq: str

class SizedIterator:
    def __init__(self, it, length):
        self.it = iter(it)
        self.length = length

    def __len__(self):
        return self.length

    def __iter__(self):
        return self

    def __next__(self):
        item = next(self.it)
        self.length -= 1
        return item

def pyfamsa_to_pyhmmer(ali: pyfamsa.Alignment) -> pyhmmer.easel.DigitalMSA:
    return pyhmmer.easel.TextMSA(
        sequences=[
            pyhmmer.easel.TextSequence(name=seq.id, sequence=seq.sequence.decode())
            for seq in ali
        ]
    ).digitize(ALPHABET)

def pytrimal_to_pyhmmer(ali: pytrimal.Alignment) -> pyhmmer.easel.TextMSA:
    return pyhmmer.easel.TextMSA(
        sequences=[
            pyhmmer.easel.TextSequence(name=name, sequence=sequence) 
            for name, sequence in zip(ali.names, ali.sequences)
        ]
    )

# ---- Process -----------------------------------------------------------------

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Process clustering data and extract sequences.")
    
    parser.add_argument("-c", "--clusters_chunk", required=True, type=str, help="Path to the clusters chunk file.")
    parser.add_argument("-f", "--fasta_file", required=True, type=str, help="Path to the MGnifams input FASTA file.")
    parser.add_argument("-p", "--cpus", required=True, type=int, help="Number of CPUs to use.")
    parser.add_argument("-n", "--chunk_num", required=True, type=str, help="Chunk number used for naming output files.")
    parser.add_argument("--discard_min_rep_length", required=True, type=int, help="Minimum allowed representative sequence length, to keep family.")
    parser.add_argument("--discard_max_rep_length", required=True, type=int, help="Maximum allowed representative sequence length, to keep family.")
    parser.add_argument("--discard_min_starting_membership", required=True, type=float, help="Minimum allowed recruited initial cluster sequence membership percentage, to keep family.")
    parser.add_argument("--max_seq_identity", required=True, type=float, help="Maximum sequence identity, to filter out sequences from an alignment.")
    parser.add_argument("--max_seed_seqs", required=True, type=int, help="Maximum number of allowed seed MSA sequences.")
    parser.add_argument("--max_gap_occupancy", required=True, type=float, help="Maximum allowed gap occupancy for trimming.")
    parser.add_argument("--recruit_evalue_cutoff", required=True, type=float, help="E-value cutoff for sequence recruitment.")
    parser.add_argument("--recruit_hit_length_percentage", required=True, type=float, help="Minimum allowed recruited sequence env hit length percentage with family HMM.")

    args = parser.parse_args(args)
    return args

def define_globals(args):
    global log_file, refined_families_tsv_file, \
        discarded_clusters_file, successful_clusters_file, \
        converged_families_file, family_metadata_file, family_reps_file, \
        tmp_folder, seed_msa_folder, full_msa_folder, hmm_folder, rf_folder, \
        tmp_family_sequences_path, tmp_seed_msa_path, tmp_full_msa_path, \
        tmp_hmm_path, tmp_sequences_to_remove_path, tmp_rf_path
    logs_folder                = "logs"
    refined_families_folder    = "refined_families"
    discarded_clusters_folder  = "discarded_clusters"
    successful_clusters_folder = "successful_clusters"
    converged_families_folder  = "converged_families"
    family_metadata_folder     = "family_metadata"
    family_reps_folder         = "family_reps"
    tmp_folder                 = "tmp"
    seed_msa_folder            = "seed_msa_sto"
    full_msa_folder            = "full_msa_sto"
    hmm_folder                 = "hmm"
    rf_folder                  = "rf"

    for folder in [tmp_folder, seed_msa_folder, full_msa_folder, hmm_folder, rf_folder, \
        logs_folder, refined_families_folder, discarded_clusters_folder, \
        successful_clusters_folder, converged_families_folder, family_metadata_folder, family_reps_folder
    ]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    tmp_seed_msa_path = os.path.join(tmp_folder, 'seed_msa.sto')
    tmp_full_msa_path = os.path.join(tmp_folder, 'full_msa.sto')

    log_file                  = os.path.join(logs_folder               , f'{args.chunk_num}.txt')
    refined_families_tsv_file = os.path.join(refined_families_folder   , f'{args.chunk_num}.tsv')
    discarded_clusters_file   = os.path.join(discarded_clusters_folder , f'{args.chunk_num}.csv')
    successful_clusters_file  = os.path.join(successful_clusters_folder, f'{args.chunk_num}.txt')
    converged_families_file   = os.path.join(converged_families_folder , f'{args.chunk_num}.txt')
    family_metadata_file      = os.path.join(family_metadata_folder    , f'{args.chunk_num}.csv')
    family_reps_file          = os.path.join(family_reps_folder        , f'{args.chunk_num}.fasta')

def create_empty_output_files():
    open(refined_families_tsv_file, 'w').close()
    open(discarded_clusters_file  , 'w').close()
    open(successful_clusters_file , 'w').close()
    open(converged_families_file  , 'w').close()
    open(family_metadata_file     , 'w').close()
    open(family_reps_file         , 'w').close()

def load_clusters_df(clusters_chunk: str) -> pd.DataFrame:
    return pd.read_csv(clusters_chunk, sep='\t', header=None, names=['representative', 'member'], dtype=str)

def log_time(start_time, text: str):
    duration = time.time() - start_time
    logging.info(f"{text} took {duration:.6f} seconds")

def get_next_family(clusters_df: pd.DataFrame) -> typing.Tuple[typing.Optional[str], typing.Optional[typing.List[str]]]:
    start_time = time.time()

    if clusters_df.empty:
        return None, None

    # Select the next family (first row in the dataFrame)
    next_family_rep = clusters_df.iloc[0]['representative']
    next_family_members = clusters_df.loc[clusters_df['representative'] == next_family_rep, 'member'].tolist()

    # Remove the selected family from the dataframe
    clusters_df.drop(clusters_df[clusters_df['representative'] == next_family_rep].index, inplace=True)

    # Logging execution time and results
    logging.info(f"get_next_family: {time.time() - start_time}")
    logging.info(f"S: {next_family_rep}, s: {len(next_family_members)}")
    
    return next_family_rep, next_family_members

def read_pyhmmer_seqs(mgnifams_input_fasta_file: str) -> pyhmmer.easel.DigitalSequenceBlock:
    start_time = time.time()

    # pre-fetching targets - fast, but needs the whole target database in memory
    with pyhmmer.easel.SequenceFile(mgnifams_input_fasta_file, digital=True) as seq_file:
        seqs = seq_file.read_block()

    log_time(start_time, "read_pyhmmer_seqs (pyhmmer): ")

    return seqs

def build_sequence_dict(seqs: pyhmmer.easel.DigitalSequenceBlock) -> dict[str, pyhmmer.easel.DigitalSequence]:
    """Safely build a name → sequence mapping without crashing like .indexed"""
    start_time = time.time()

    seq_dict = {seq.name.decode(): seq for seq in seqs}
    
    log_time(start_time, "build_sequence_dict: ")

    return seq_dict

def get_fasta_sequences(
    seq_dict: dict[str, pyhmmer.easel.DigitalSequence],
    headers: typing.Iterable[str],
) -> typing.List[Sequence]:
    """Retrieve sequences from a prebuilt name → sequence dictionary."""
    results = []

    for header in headers:
        seq = seq_dict.get(header)
        if seq is not None:
            decoded_seq = ALPHABET.decode(seq.sequence)
            results.append(Sequence(header, decoded_seq))

    return results

def run_initial_msa(
    members: typing.Iterable[str],
    seqs: pyhmmer.easel.DigitalSequenceBlock,
    cpus: int = 0,
) -> pyhmmer.easel.DigitalMSA:
    """Aligns sequences using pyfamsa and writes the result to a FASTA file."""
    start_time = time.time()
    
    # Convert sequences to pyfamsa.Sequence objects
    seq_dict = {seq.name: seq for seq in seqs if seq.name in map(str.encode, members)}

    sequences = [
        pyfamsa.Sequence(member, seq_dict[member].textize().sequence.encode())
        for member in map(str.encode, members)
        if member in seq_dict
    ]
    if not sequences:
        raise ValueError("No valid sequences found for alignment.")

    # Create an Aligner object
    aligner = pyfamsa.Aligner()
    
    # Perform the multiple sequence alignment
    alignment = aligner.align(sequences)

    log_time(start_time, "run_initial_msa (pyfamsa): ")
    return pyfamsa_to_pyhmmer(alignment)

def run_hmmbuild(chunk: int, iteration: int, seed_msa: pyhmmer.easel.DigitalMSA, hand: bool = False) -> pyhmmer.plan7.HMM:
    """Runs HMMER's hmmbuild using pyhmmer."""
    start_time = time.time()

    seed_msa.name = f"{chunk}_{iteration}".encode()

    architecture: typing.Literal["hand", "fast"] = "hand" if hand else "fast"
    builder = pyhmmer.plan7.Builder(ALPHABET, architecture=architecture)
    background = pyhmmer.plan7.Background(ALPHABET)
    hmm, _, _ = builder.build_msa(seed_msa, background)

    log_time(start_time, "run_hmmbuild (pyhmmer): ")

    return hmm

def mask_sequence(sequence: Sequence, env_from: int, env_to: int) -> Sequence:
    # Unpack the tuple (header, seq)
    header, seq = sequence

    # Modify header and sequence
    header = f"{header}/{env_from}_{env_to}"
    seq = seq[env_from - 1:env_to]  # -1 for 0-based indexing

    # Replace the old tuple with the new one
    return Sequence(header, seq)

def run_hmmsearch(
    hmm: pyhmmer.plan7.HMM, 
    pyhmmer_seqs: pyhmmer.easel.DigitalSequenceBlock, 
    seq_dict: dict[str, pyhmmer.easel.DigitalSequence],
    exit_flag: bool,
    recruit_evalue_cutoff: float,
    recruit_hit_length_percentage: float,
    cpus: int = 0,
) -> typing.List[Sequence]: # hmm obj changes, pyhmmer_seqs is the same, whole protein library
    """Runs HMMER's hmmsearch using pyhmmer."""
    start_time = time.time()

    filtered_sequences = []

    for top_hits in pyhmmer.hmmer.hmmsearch(hmm, pyhmmer_seqs, cpus=cpus, E=recruit_evalue_cutoff):
        qlen = top_hits.query.M
        for hit in top_hits:
            sequence_name = hit.name.decode()
            tlen = hit.length
            for subhit in hit.domains:
                env_length = subhit.env_to - subhit.env_from + 1
                if (exit_flag or env_length >= recruit_hit_length_percentage * qlen):
                    sequence = get_fasta_sequences(seq_dict, [sequence_name])[0]
                    if (env_length < tlen):
                        sequence = mask_sequence(sequence, subhit.env_from, subhit.env_to)
                    filtered_sequences.append(sequence)

    log_time(start_time, "run_hmmsearch (pyhmmer): ")

    return filtered_sequences

def run_hmmalign(hmm: pyhmmer.plan7.HMM, family_sequences: typing.Iterable[Sequence]) -> typing.Tuple[pyhmmer.easel.TextMSA, int, int]:
    """Runs HMMER's hmmalign using pyhmmer."""
    start_time = time.time()

    num_seqs_result = 0
    non_gap_seq_length = 0

    seqs = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=seq.id.encode(), sequence=seq.seq)
        for seq in family_sequences
    ).digitize(ALPHABET)   
    hmmalign_res = typing.cast(pyhmmer.easel.TextMSA, pyhmmer.hmmer.hmmalign(hmm, seqs, trim=False))

    num_seqs_result = len(hmmalign_res.names)
    non_gap_seq_length = len(re.sub(r"[.\-~]", "", hmmalign_res.alignment[0]))

    log_time(start_time, "run_hmmalign (pyhmmer): ")

    return hmmalign_res, num_seqs_result, non_gap_seq_length

def check_rep_length(non_gap_seq_length, iteration, args):
    """Check if representative length is within limits.

    Returns:
        discard_flag (bool), discard_reason (str or None), discard_value (int or None)
    """
    if non_gap_seq_length < args.discard_min_rep_length:
        logging.info(f"Discard-Warning: {iteration} rep length is only {non_gap_seq_length}")
        return True, "family representative length too small", non_gap_seq_length
    elif non_gap_seq_length > args.discard_max_rep_length:
        logging.info(f"Discard-Warning: {iteration} rep length is over {non_gap_seq_length}")
        return True, "family representative length too large", non_gap_seq_length
    return False, None, None

def extract_first_part(sequence_name: str) -> str:
    return sequence_name.split('/')[0]

def unmask_sequence_names(sequences: typing.Iterable[Sequence]) -> typing.List[str]:
    return [extract_first_part(name) for name, _ in sequences]

def run_pytrimal_reps(full_msa: pyhmmer.easel.TextMSA, threshold: float, max_seed_seqs: int) -> pyhmmer.easel.TextMSA:
    start_time = time.time()

    # keep reference line
    assert full_msa.reference is not None
    rf = np.array(list(full_msa.reference.decode()))

    # Load the MSA - Stream sequences instead of materializing all at once
    sequences = SizedIterator(
        (seq.upper().replace('.', '-').replace('~', '-') for seq in full_msa.alignment),
        len(full_msa.alignment)
    )
    ali = pytrimal.Alignment(full_msa.names, sequences)

    # and remove redundant sequences
    repTrimmer = pytrimal.RepresentativeTrimmer(identity_threshold=threshold)
    ali = repTrimmer.trim(ali)
    rf = rf[ali.residues_mask]

    logging.info("run_pytrimal_reps finished.")
    number_of_remaining_sequences = len(list(ali.names))
    logging.info(f"Remaining sequences: {str(number_of_remaining_sequences)}")

    if (number_of_remaining_sequences > max_seed_seqs):
        logging.info(f"Capping to {str(max_seed_seqs)} most representative")
        repTrimmer = pytrimal.RepresentativeTrimmer(clusters=max_seed_seqs)
        ali = repTrimmer.trim(ali)
        rf = rf[ali.residues_mask]

    log_time(start_time, "run_pytrimal_reps: ")

    out_msa = pytrimal_to_pyhmmer(ali)
    out_msa.reference = ''.join(rf).encode()
    return out_msa

def calculate_trim_positions(sequence_matrix: np.ndarray, occupancy_threshold: float):
    numeric_matrix = np.where(sequence_matrix == '-', 0, 1)
    num_rows = numeric_matrix.shape[0]
    column_sums = np.sum(numeric_matrix, axis=0)
    column_sums_percentage = column_sums / num_rows
    start_position = np.argmax(column_sums_percentage > occupancy_threshold)
    end_position = len(column_sums_percentage) - np.argmax(column_sums_percentage[::-1] > occupancy_threshold) - 1

    return start_position, end_position

def clip_ends(msa: pyhmmer.easel.TextMSA, occupancy_threshold: float) -> pyhmmer.easel.TextMSA:
    sequence_matrix = np.array([list(row) for row in msa.alignment])
    start_position, end_position = calculate_trim_positions(sequence_matrix, occupancy_threshold)
    
    return msa.select(columns=range(start_position, end_position))

def check_seed_membership(original_sequence_names, filtered_seq_names):
    original_first_parts = set(map(extract_first_part, original_sequence_names))
    filtered_first_parts = set(map(extract_first_part, filtered_seq_names))
    common_first_parts = original_first_parts & filtered_first_parts
    common_count = len(common_first_parts)
    percentage_membership = common_count / len(original_sequence_names)

    return percentage_membership

def parse_protein_name(seq_name, seq_length, seq_whole_name, original_length, start, end):
    splits = seq_name.split('_')

    if (end - start) == original_length:
        if len(splits) == 3: # example 199394972_125_238
            seq_name = f"{splits[0]}/{splits[1]}-{splits[2]}"
        return seq_name
    else:
        old_length = len(seq_whole_name)
        if len(splits) == 3: # example 199394972_125_238
            new_start = start + int(splits[1])
            new_end = new_start + seq_length - 1
            new_name = f"{splits[0]}/{new_start}-{new_end}"
        else:
            new_name = f"{splits[0]}/{start + 1}-{end}"

        # Pad or trim new_name to match old_length exactly
        if len(new_name) < old_length:
            new_name = new_name.ljust(old_length)
        elif len(new_name) > old_length:
            # If longer, trim spaces at the end to fit old_length
            # Trim only trailing spaces, but if still too long, hard cut
            new_name = new_name.rstrip()
            if len(new_name) > old_length:
                new_name = new_name[:old_length]
            else:
                new_name = new_name.ljust(old_length)

        return new_name

def renumber_sto_msa(
    in_sto_file,
    out_sto_file,
    seq_dict,
    write_metadata=False,
    iteration=None,
    chunk_num=None,
    full_msa_num_seqs=None,
    consensus=None,
    converged=None
):
    if write_metadata:
        assert iteration is not None and chunk_num is not None

        familyfile = open(refined_families_tsv_file, 'a')
        metadatafile = open(family_metadata_file, 'a')
        repsfile = open(family_reps_file, 'a')

    with open(in_sto_file, 'r') as infile, open(out_sto_file, 'w') as outfile:
        seq_name = ""
        rep_flag = True
        seen_seq_names = set()  # Found edge cases that contain duplicates

        for line in infile:
            split_line = line.split()

            if split_line: # If not empty
                if split_line[0] == '//' or split_line[1] == 'RF':
                    outfile.write(line)
                elif split_line[1] == 'STOCKHOLM':
                    outfile.write(f'{line}\n')
                elif split_line[0] not in ['#=GF', '#=GS', '#=GC', '#=GR']:
                    raw_seq_name = split_line[0].split("/")[0]
                    seq = re.sub(r"[.\-~]", "", split_line[1]).upper()
                    original_seq = get_fasta_sequences(seq_dict, [raw_seq_name])[0][1]
                    start = original_seq.find(seq)
                    end = start + len(seq)
                    seq_name = parse_protein_name(raw_seq_name, len(seq), split_line[0], len(original_seq), start, end)

                    if seq_name.strip() in seen_seq_names:
                        continue  # Skip duplicates
                    seen_seq_names.add(seq_name.strip())  # Mark as seen

                    line = line.replace(split_line[0], seq_name, 1)
                    outfile.write(line)

                    if write_metadata:
                        familyfile.write(f"{iteration}\t{seq_name}\n")
                        if rep_flag:
                            splits = seq_name.split("/")
                            region = splits[1].strip() if "/" in seq_name else "-"
                            metadatafile.write(f"{iteration},{full_msa_num_seqs},\"{splits[0]}\",{region},{len(seq)},{seq},{consensus},{converged}\n")
                            repsfile.write(f">{seq_name.strip()}\t{chunk_num}_{iteration}\n{seq}\n")
                            rep_flag = False

    if write_metadata:
        familyfile.close()
        metadatafile.close()
        repsfile.close()

def remove_tmp_files() -> None:
    for item in os.listdir(tmp_folder):
        item_path = os.path.join(tmp_folder, item)
        if os.path.isfile(item_path):
            os.remove(item_path)

def main():
    args = parse_args()
    define_globals(args)
    create_empty_output_files()
    clusters_df = load_clusters_df(args.clusters_chunk)
    pyhmmer_seqs = read_pyhmmer_seqs(args.fasta_file)
    seq_dict = build_sequence_dict(pyhmmer_seqs)

    iteration = 0
    while True:
        iteration += 1
        family_rep, family_members = get_next_family(clusters_df)
        if not family_members:
            logging.info("Exiting all...")
            break
        
        original_sequence_names = family_members

        seed_msa = run_initial_msa(family_members, pyhmmer_seqs, cpus=args.cpus)

        total_checked_sequences = []
        full_msa_num_seqs = 0
        discard_flag = False
        discard_reason = ""
        discard_value = 0.0
        exit_flag = False
        converged = False
        hmm = ""
        family_iteration = 0
        while True:
            family_iteration += 1
            if (family_iteration > 3):
                exit_flag = True

            if (exit_flag):
                logging.info("Exiting-3 loops.")
            logging.info(f"{str(family_iteration)}")

            if not exit_flag: # main strategy branch
                hmm = run_hmmbuild(args.chunk_num, iteration, seed_msa)

                filtered_seqs = run_hmmsearch(
                    hmm, 
                    pyhmmer_seqs,
                    seq_dict,
                    exit_flag=exit_flag, 
                    cpus=args.cpus, 
                    recruit_evalue_cutoff=args.recruit_evalue_cutoff,
                    recruit_hit_length_percentage=args.recruit_hit_length_percentage,
                )
                if (len(filtered_seqs) == 0): # low complexity sequence, confounding cluster, discard and move on to the next
                    discard_flag = True
                    discard_reason = "low complexity model - confounding cluster"
                    discard_value = 0.0
                    break

                new_recruited_sequences = set(unmask_sequence_names(filtered_seqs)) - set(total_checked_sequences) # new_recruited_sequences always has something at first turn since total_checked_sequences starts empty []
                total_checked_sequences += list(new_recruited_sequences)

                if not new_recruited_sequences:
                    exit_flag = True
                    converged = True
                    logging.info("Exiting-CONVERGED: no new sequences recruited.")
                    with open(converged_families_file, 'a') as file:
                        file.write(f"{iteration}\n")

            if exit_flag: # exit strategy branch
                logging.info("Exiting branch strategy:")

                final_hmm = run_hmmbuild(args.chunk_num, iteration, seed_msa, hand=True)

                filtered_seqs = run_hmmsearch(
                    hmm, 
                    pyhmmer_seqs,
                    seq_dict,
                    exit_flag, 
                    args.recruit_evalue_cutoff, 
                    args.recruit_hit_length_percentage, 
                    args.cpus
                )
                if (len(filtered_seqs) == 0): # low complexity sequence, confounding cluster, discard and move on to the next
                    discard_flag = True
                    discard_reason = "low complexity model - confounding cluster"
                    discard_value = 0.0
                    break

                membership_percentage = check_seed_membership(original_sequence_names, unmask_sequence_names(filtered_seqs))
                if (membership_percentage < args.discard_min_starting_membership):
                    discard_flag = True
                    discard_reason = "few seed sequences remained"
                    discard_value = membership_percentage
                    logging.info(f"Discard-Warning: {iteration} seed percentage in MSA is {membership_percentage}")
                    break
                elif (membership_percentage < 1):
                    logging.info(f"Warning: {iteration} seed percentage in MSA is {membership_percentage}")

                full_msa, full_msa_num_seqs, non_gap_seq_length = run_hmmalign(hmm, filtered_seqs) # final full MSA, including smaller sequences
                discard_flag, discard_reason, discard_value = check_rep_length(non_gap_seq_length, iteration, args)
                if discard_flag:
                    break

                break # break from main strategy

            # main strategy continue, if not converged
            hmmalign_res, _, non_gap_seq_length = run_hmmalign(hmm, filtered_seqs)
            discard_flag, discard_reason, discard_value = check_rep_length(non_gap_seq_length, iteration, args)
            if discard_flag:
                break

            full_msa = run_pytrimal_reps(hmmalign_res, args.max_seq_identity, args.max_seed_seqs) # removes redundant sequences
            full_msa = clip_ends(full_msa, args.max_gap_occupancy) # removes gaps above threshold at ends
            seed_msa = full_msa.digitize(ALPHABET)

        # Exiting family loop
        if (discard_flag): # unsuccessfully
            logging.info(f"Discarding cluster {family_rep}")
            with open(discarded_clusters_file, 'a') as outfile:
                outfile.write(str(family_rep) + "," + discard_reason + "," + str(discard_value) + "\n")
            iteration -= 1 # keep proper track of family ids
        else: # successfully
            # write successful cluster
            with open(successful_clusters_file, 'a') as outfile:
                outfile.write(str(family_rep) + "\n")
            # write reference line
            with open(os.path.join(rf_folder, f'{args.chunk_num}_{iteration}.txt'), "w") as f:
                f.write(seed_msa.reference.decode())
            # write HMM
            with open(os.path.join(hmm_folder, f'{args.chunk_num}_{iteration}.hmm'), "wb") as f:
                final_hmm.write(f)
            # write alignments
            with open(tmp_seed_msa_path, "wb") as dst:
                seed_msa.write(dst, format="pfam")
            with open(tmp_full_msa_path, "wb") as dst:
                full_msa.write(dst, format="pfam")

            renumber_sto_msa(
                in_sto_file=tmp_seed_msa_path,
                out_sto_file=os.path.join(seed_msa_folder, f'{args.chunk_num}_{iteration}.sto'),
                seq_dict=seq_dict
            )
            renumber_sto_msa(
                in_sto_file=tmp_full_msa_path,
                out_sto_file=os.path.join(full_msa_folder, f'{args.chunk_num}_{iteration}.sto'),
                seq_dict=seq_dict,
                write_metadata=True,
                iteration=iteration,
                chunk_num=args.chunk_num,
                full_msa_num_seqs=full_msa_num_seqs,
                consensus=final_hmm.consensus,
                converged=converged
            )

        remove_tmp_files()

    # # End of all families
    logging.info("DONE.")
    logging.shutdown()
    shutil.move("execution.log", log_file)

if __name__ == "__main__":
    main()
