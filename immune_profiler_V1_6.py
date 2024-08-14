#!/usr/bin/env python3
# Author: Sean Sun

script_name = 'Cogent_NGS_immune_profiler'
version = 'CogentIP v1.6'
desc = """
Description:
This script is a part of the Cogent NGS immune profiling pipeline (CogentIP v1.6). 
It leverages tools like ImmuneProfiler, MIGEC, and MiXCR to perform immune repertoire profiling.
The script includes functionalities for:
- Generating input files for VDJviz.
- Adding custom UMI cutoffs via command-line arguments.
- Enabling human and mouse species profiling.
- Error handling and reporting for various steps in the pipeline.
- Managing directories and files during execution.

Done:
- Integrated ImmuneProfiler v1.6 with MIGEC 1.2.9 and MiXCR 3.
- Added BCRv1, BCRv2, TCRv1, and TCRv2 as selectable options.
- Included a custom UMI cutoff via arguments.
- Enhanced error reporting and handling mechanisms.
- Adjusted stats filename generation.

Todo:
- Fully enable species selection (human, mouse) with proper error handling.
- Remove analysis steps not involving MIGEC.
- Refactor functions for running MiXCR and MIGEC.
- Improve user-friendly error messages.

"""

# ---------- Import Required Modules ---------- #

import argparse
import csv
import logging
import re
from collections import OrderedDict
import datetime
import openpyxl
import subprocess
from contextlib import ExitStack
import os
import fnmatch
import sys
import copy
import gzip
import platform
import time
from shutil import copyfile
from itertools import chain
import shutil
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Align
#import pairwise2

# ---------- Function Definitions ---------- #

def required_python_module_check():
    """
    Checks for the presence of required Python packages and the correct Java version.

    Outputs:
    - A message indicating whether all required packages are installed.
    - Verifies that Java version is compatible with CogentIP.
    - Exits the script if any package is missing or Java version is incompatible.
    """
    # List of required packages
    required_pkgs = ["argparse", "csv", "logging", "re", "collections",
                     "datetime", "openpyxl", "subprocess", "contextlib", "os",
                     "fnmatch", "sys", "copy", "gzip", "platform", "numpy",
                     "time", "shutil", "matplotlib.pyplot", "Bio", "itertools"]

    # Track missing packages
    not_installed = []

    sys.stdout.write('\nChecking for Python packages required by Immune Profiler...\n')

    # Check if required packages are installed
    for pkg in required_pkgs:
        try:
            __import__(pkg)
        except ImportError:
            not_installed.append(pkg)

    # Report missing packages
    if not_installed:
        sys.stdout.write('\n')
        for package in not_installed:
            sys.stdout.write(f'\t-- {package} is not installed.\n')
        sys.stdout.write('\nSee https://docs.python.org/3/installing/ for help with package installation.\n')
        sys.exit(0)
    else:
        sys.stdout.write('\nAll required packages are installed.\n\n')

    # Check Java version compatibility
    vjava = subprocess.check_output(['java', '-version'], stderr=subprocess.STDOUT)
    vjava = vjava.decode("utf-8").split("\n")[0]
    vjava_version = int(vjava.split('"')[1].split('.')[0])
    
    if int(re.match(r'([a-z]+)([0-9]+)', p_args.cogentIP_version, re.I).groups()[1]) < 2 and vjava_version > 15:
        print(f'your Java version: {vjava}, required v1.8 or v8-15')
        sys.exit(0)

    # Ensure installation path has no spaces
    if ' ' in p_args.meta_file:
        print('The address for CogentIP installation should be a string without spaces')
        sys.exit(0)

    return vjava

def check_dir(d, exp):
    """
    Checks if a directory exists. Exits with an error if not.

    Args:
    - d (str): Directory path.
    - exp (str): Error message to display if the directory does not exist.
    
    Returns:
    - True if directory exists.
    """
    if not os.path.isdir(d):
        print(f'{get_time()} [ERROR] {d} {exp}', flush=True)
        sys.exit(1)
    return True

def check_file(f, exp):
    """
    Checks if a file exists. Exits with an error if not.

    Args:
    - f (str): File path.
    - exp (str): Error message to display if the file does not exist.
    
    Returns:
    - True if file exists.
    """
    if not os.path.isfile(f):
        print(f'{get_time()} [ERROR] {f} {exp}', flush=True)
        sys.exit(1)
    return True

def get_time():
    """
    Returns the current timestamp in the format 'YYYY-MM-DD HH:MM:SS'.

    Returns:
    - str: Current timestamp.
    """
    ts = time.time()
    return datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

def timer(start_time):
    """
    Calculates the time elapsed since a given start time.

    Args:
    - start_time (datetime): The time to start timing from.
    
    Returns:
    - str: Time difference in 'Hh Mm Ss' format.
    """
    end_time = datetime.datetime.now()
    elapsed_time = end_time - start_time
    secs = round(elapsed_time.total_seconds() % 60, 2)
    mins = int(elapsed_time.total_seconds() // 60) % 60
    hours = int(elapsed_time.total_seconds() // 3600)
    diff = ':'.join([str(hours) + 'h', str(mins) + 'm', str(secs) + 's'])
    return diff

def load_meta(meta_file, fastq_dir):
    """
    Load metadata from a CSV file into an ordered dictionary.

    Args:
        meta_file (str): Path to the metadata CSV file.
        fastq_dir (str): Directory containing FASTQ files.

    Returns:
        OrderedDict: A dictionary where keys are sample names and values are lists of FASTQ file names.
    
    Raises:
        SystemExit: If there are errors in the metadata file or if required files are not found.
    """
    meta_dict = OrderedDict()
    try:
        infile = open(meta_file, 'r')
        next(infile)  # Skip header line
        for line in infile:
            words = [word.strip() for word in line.rstrip().split(',')]
            sname = ''
            r1_fq = ''
            r2_fq = ''
            if len(words) == 3:
                sname, r1_fq, r2_fq = words[0], words[1], words[2]
            elif len(words) < 3:
                print(get_time() + '[ERROR] sample ' + sname + ' does not have paired-fastqs')
            if sname in meta_dict:
                print(get_time() + '[ERROR] Found duplicate sample ID: ' + sname +
                      ', please rename it in meta data file ' + meta_file, flush=True)
                sys.exit(1)
            elif '_' in sname:
                print(get_time() + '[ERROR] Found _ in sample ID: ' + sname +
                      ', please rename it in meta data file ' + meta_file, flush=True)
                sys.exit(1)
            elif len(sname) > 20:
                print(get_time() + '[ERROR] There are more than 20 characters in sample ID: '
                      + sname + ', please rename', flush=True)
                sys.exit(1)
            meta_dict[sname] = [r1_fq, r2_fq]
            for fq in meta_dict[sname]:
                if find_file(fq, fastq_dir) == []:
                    print(get_time() + '[ERROR] The input FASTQ: ' + fq + ' is not found. \n' +
                        'Please make sure the FASTQ directory is correct & ' +
                        'FASTQ name in metadata file is correct (case sensitive)', flush=True)
                    sys.exit(1)
        infile.close()
        if '' in meta_dict.keys():
            del meta_dict['']
            print(get_time() + ' [WARN] Delete empty sample ID in meta file: ' + meta_file, flush=True)
            print(get_time() + ' [WARN] After deletion, ' + str(len(meta_dict.keys())) + ' samples to process',
                  flush=True)
            print(get_time() + ' [WARN] Please check and modify meta file if the sample number is incorrect',
                  flush=True)
            print(get_time() + ' [WARN] Check if there is any extra empty line in meta file', flush=True)
        return meta_dict
    except EnvironmentError:
        print(get_time() + ' [ERROR] Unable to open meta file: ' + meta_file, flush=True)
        sys.exit(1)


def load_bcs(bc_file):
    """
    Load barcode sequences from a CSV file.

    Args:
        bc_file (str): Path to the barcode CSV file.

    Returns:
        list: A list of barcode sequences.
    """
    bcs = []
    with open(bc_file) as f:
        for line in f:
            l = line.rstrip().split(',')[0]
            bcs.append(l)
    return bcs


def load_pseq(ps_file):
    """
    Load a mapping of current sequences to original sequences from a file.

    Args:
        ps_file (str): Path to the sequence mapping file.

    Returns:
        OrderedDict: A dictionary mapping current sequences to original sequences.
    """
    with open(ps_file) as f:
        ps_dict = OrderedDict()
        for line in f:
            (cur_seq, ori_seq) = re.split(',', line.strip())
            ps_dict[cur_seq] = ori_seq
    return ps_dict


def update_meta(bcs, meta_dict, nchain):
    """
    Update the metadata dictionary by appending barcode sequences to each sample entry.

    Args:
        bcs (list): List of barcode sequences.
        meta_dict (OrderedDict): Metadata dictionary to be updated.
        nchain (int): Number of barcode chains to assign per sample.

    Returns:
        None
    """
    idx = 0
    for k in meta_dict.keys():
        meta_dict[k].append(bcs[idx:(idx + nchain)])
        idx += nchain
    return


def load_sd(sd_fname):
    """
    Load sample data from a tab-separated file.

    Args:
        sd_fname (str): Path to the sample data file.

    Returns:
        OrderedDict: A dictionary where keys are sample names and values are lists of tuples with sample data.
    """
    sd = OrderedDict()
    with open(sd_fname) as f:
        for line in f:
            segs = line.rstrip().split('\t')
            key = segs[0][:-1]
            values = []
            for i in range(1, len(segs)):
                items = segs[i].split(',')
                item0 = items[0]
                item1 = int(items[1]) if items[1].isnumeric() else None
                item2 = float(items[2]) if items[2].replace(".", "").isnumeric() else None
                mytuple = (item0, item1, item2)
                values.append(mytuple)
            sd[key] = values
    return sd


def change_dir(dirname, distname):
    """
    Rename a directory.

    Args:
        dirname (str): Current directory name.
        distname (str): New directory name.

    Returns:
        None

    Raises:
        SystemExit: If the directory cannot be renamed or does not exist.
    """
    if os.path.isdir(dirname):
        try:
            os.rename(dirname, distname)
        except OSError as err:
            print('ERROR: Unable to change directory from:' + dirname + " to " + distname, flush=True)
            sys.exit(1)
    else:
        print('ERROR: The folder ' + dirname + ' does not exist', flush=True)
        sys.exit(1)


def create_dir(dirname):
    """
    Create a new directory.

    Args:
        dirname (str): Directory name to create.

    Returns:
        None

    Raises:
        SystemExit: If the directory cannot be created or already exists.
    """
    if os.path.isdir(dirname):
        print('ERROR: The folder ' + dirname + ' exists', flush=True)
        sys.exit(1)
    else:
        try:
            os.makedirs(dirname)
        except OSError as err:
            print('ERROR: Unable to create directory:' + dirname, flush=True)
            sys.exit(1)


def load_fastq(f):
    """
    Open a FASTQ file, supporting gzip-compressed files.

    Args:
        f (str): Path to the FASTQ file.

    Returns:
        file object: The opened FASTQ file.
    """
    if f.endswith('.gz'):
        f = gzip.open(f, mode='rt')
    else:
        f = open(f)
    return f


def get_sample_bcs(sample, meta_dict):
    """
    Retrieve barcode sequences for a given sample.

    Args:
        sample (str): Sample ID.
        meta_dict (OrderedDict): Metadata dictionary.

    Returns:
        list: List of barcode sequences associated with the sample.
    """
    bcs = meta_dict[sample][2]
    return bcs


def ham_check(hh, bc):
    """
    Check if a barcode sequence exists in the provided dictionary.

    Args:
        hh (dict): Dictionary of barcode sequences.
        bc (str): Barcode sequence to check.

    Returns:
        str: The barcode sequence if found, otherwise an empty string.
    """
    try:
        return hh[bc]
    except KeyError:
        return ''


def get_errorbs_bytwoseqs(seq1, seq2, errbases):
    """
    Calculate the number of errors between two sequences up to a maximum number of allowed errors.

    Args:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
        errbases (int): Maximum number of allowed errors.

    Returns:
        tuple: (sequence with errors, number of errors)
    """
    errbcount = 0
    j = 0
    seqL = len(seq1)
    while errbcount <= errbases and j < seqL:
        if seq1[j] != seq2[j]:
            errbcount += 1
        j += 1
    if errbcount <= errbases:
        return seq2, errbcount
    return '', errbcount


def get_alia_seqs(seq):
    """
    Generate all possible alias sequences for a sequence containing ambiguous characters.

    Args:
        seq (str): Sequence with ambiguous characters.

    Returns:
        list: List of all possible alias sequences.
    """
    seqs = [seq]
    while len(seqs) > 0 and any(char in seqs[0] for char in 'RSIYKM'):
        del_num = len(seqs)
        for i in range(len(seqs)):
            for char, replacements in {'R': 'AG', 'S': 'CG', 'I': 'ACT', 'Y': 'CT', 'K': 'GT', 'M': 'AC'}.items():
                if char in seqs[i]:
                    for replacement in replacements:
                        myseq = seqs[i].replace(char, replacement, 1)
                        seqs.append(myseq)
        if len(seqs) > del_num:
            del seqs[:del_num]
    return seqs


def ham_check_flex(prs, bc):
    """
    Flexible check to match a barcode sequence with sequences in a dictionary, considering alias sequences.

    Args:
        prs (dict): Dictionary of possible sequences and their original forms.
        bc (str): Barcode sequence to match.

    Returns:
        str: Matching sequence if found, otherwise an empty string.
    """
    errbs = int(len(bc) // 4)
    seqs = list(prs.keys())
    for seq in seqs:
        mappedseq, mappederr = get_errorbs_bytwoseqs(bc, seq, errbs)
        if len(mappedseq) > 0:
            return mappedseq
        if any(char in seq for char in 'RSIYKM'):
            seqalias = get_alia_seqs(seq)
            for alias in seqalias:
                mappedseq, mappederr = get_errorbs_bytwoseqs(bc, alias, errbs)
                if len(mappedseq) > 0:
                    return alias
    return ''


def get_dyn_align_score(bc, seq, match, mismatch, opengap, extgap):
    """
    Calculate the dynamic alignment score between two sequences.

    Args:
        bc (str): First sequence.
        seq (str): Second sequence.
        match (int): Score for matching characters.
        mismatch (int): Penalty for mismatched characters.
        opengap (int): Penalty for opening a gap.
        extgap (int): Penalty for extending a gap.

    Returns:
        tuple: (score with penalty, maximum score)
    """
    maxscore = aligner.align(bc, seq)[0].score
    return maxscore - 2, maxscore


def best_dynamic_score(bc, seq, match, mismatch, opengap, extgap):
    """
    Find the best dynamic alignment score between a barcode sequence and a list of possible sequences.

    Args:
        bc (str): Barcode sequence.
        seq (str): Sequence to match.
        match (int): Score for matching characters.
        mismatch (int): Penalty for mismatched characters.
        opengap (int): Penalty for opening a gap.
        extgap (int): Penalty for extending a gap.

    Returns:
        int: Best alignment score found.
    """
    score = -1
    seqs = get_alia_seqs(seq)
    for myseq in seqs:
        myscore, _ = get_dyn_align_score(bc, myseq, match, mismatch, opengap, extgap)
        if myscore > len(bc):
            return myscore
        elif myscore > score:
            score = myscore
    return score


def save_to_pass_ambiguous_chain_map(bc, seq, score, p_args, sample, prs):
    """
    Save successful matches with ambiguous sequences to a CSV file.

    Args:
        bc (str): Barcode sequence.
        seq (str): Matching sequence.
        score (int): Alignment score.
        p_args (argparse.Namespace): Argument parser namespace containing output directory.
        sample (str): Sample ID.
        prs (dict): Dictionary of possible sequences and their original forms.

    Returns:
        None
    """
    if p_args.internal_test and score == len(bc) + 1:
        fname = os.path.join(p_args.out_dir, 'pass_ambiguous_chain_map_' + sample + '.csv')
        with open(fname, 'a') as f:
            if not os.path.exists(fname) or os.stat(fname).st_size == 0:
                f.write("seq,chain\n")
            f.write(bc + ',' + seq + '(' + prs[seq] + ')\n')

def align_prs_bestscore(sample, prs, bc):
    """
    Aligns a barcode sequence with sequences in a dictionary and selects the best matching sequence based on alignment scores.

    Args:
        sample (str): Sample identifier.
        prs (dict): Dictionary containing possible sequences and their original forms.
        bc (str): Barcode sequence to align.

    Returns:
        str: The best matching sequence if alignment score criteria are met, otherwise an empty string.
    """
    # Define scoring parameters
    match = 2
    mismatch = -2
    opengap = -2
    extgap = -1

    # Skip alignment if the barcode contains 'GGGGGG'
    if bc.find('GGGGGG') >= 0:
        return ''

    seqs = list(prs.keys())
    best_score = 0
    best_id = -1

    for i in range(len(seqs)):
        seq = seqs[i]
        score, maxscore = get_dyn_align_score(bc, seq, match, mismatch, opengap, extgap)

        # If the score exceeds the length of the barcode, save the match and return the sequence
        if score > len(bc):
            save_to_pass_ambiguous_chain_map(bc, seq, score, p_args, sample, prs)
            return seq

        # Check for ambiguous characters and use flexible alignment if needed
        if any(char in seq for char in 'RSIYKM'):
            score = best_dynamic_score(bc, seq, match, mismatch, opengap, extgap)
            if score > len(bc):
                save_to_pass_ambiguous_chain_map(bc, seq, score, p_args, sample, prs)
                return seq

        if score > best_score:
            best_score = score
            best_id = i

    # Save failed ambiguous matches if criteria are met
    if p_args.internal_test and best_score == len(bc):
        fname = os.path.join(p_args.out_dir, 'fail_ambiguous_chain_map_' + sample + '.csv')
        with open(fname, 'a') as f:
            if not os.path.exists(fname) or os.stat(fname).st_size == 0:
                f.write("seq,chain\n")
            f.write(bc + ',' + seqs[best_id] + '(' + prs[seqs[best_id]] + ')\n')

    return ''

def get_partial_primer_start_end(p_args, prs):
    """
    Determines the start and end positions for partial primer sequences based on receptor kit and species.

    Args:
        p_args (argparse.Namespace): Argument parser namespace containing receptor kit and species information.
        prs (dict): Dictionary of possible sequences and their original forms.

    Returns:
        tuple: Start and end positions for the primer sequences.
    """
    # Set start and end positions based on receptor kit and species
    if p_args.receptor_kit == 'TCRv1':
        start, end = (7, 13) if p_args.species == 'mouse' else (18, 24)
    elif p_args.receptor_kit == 'TCRv2':
        start, end = (7, 13) if p_args.species == 'mouse' else (5, 11)
    elif p_args.receptor_kit == 'BCRv1':
        start, end = (16, 23) if p_args.species == 'mouse' else (16, 22)
    elif p_args.receptor_kit in ['BCRv2', 'BCRv2sub', 'BCRv2hy']:
        if p_args.species == 'human' and len(list(prs.keys())[0]) > 9:
            start, end = 14, 26
        elif p_args.species == 'mouse':
            start, end = 13, 27
        else:
            start, end = 17, 24

    # Adjust start and end positions if necessary
    if end - start != len(list(prs.keys())[0]):
        if len(list(prs.keys())[0]) > 11 and p_args.receptor_kit == 'BCRv1' and p_args.species == 'mouse':
            start -= 2
        end = start + len(list(prs.keys())[0])

    return start, end

def get_subseq(sample, read, repo_type, species, pseq_dict, prs, len_cutoff=30, mismatch=1):
    """
    Extracts and processes a subsequence from a read based on specified parameters.

    Args:
        sample (str): Sample identifier.
        read (tuple): Tuple containing the read information (header, sequence, etc.).
        repo_type (str): Type of repository (not used in function).
        species (str): Species type (not used in function).
        pseq_dict (dict): Dictionary of possible sequences (not used in function).
        prs (dict): Dictionary of possible sequences and their original forms.
        len_cutoff (int): Minimum length cutoff for subsequences.
        mismatch (int): Allowed number of mismatches.

    Returns:
        str: Processed subsequence or 'short' if sequence is too short.
    """
    start, end = get_partial_primer_start_end(p_args, prs)

    # Extract subsequence based on length cutoff
    subseq = read[1][start:end] if len(read[1]) >= len_cutoff else 'short'

    # Check for exact or flexible matches
    if mismatch == 0 or subseq in prs.keys() or subseq == 'short':
        pass
    else:
        rawseq = ham_check_flex(prs, subseq)
        if len(rawseq) == 0 and p_args.species == 'mouse' and p_args.receptor_kit.find('BCR') >= 0:
            rawseq = align_prs_bestscore(sample, prs, subseq)
        if len(rawseq) > 0:
            subseq = rawseq

    return subseq

def get_umi(read, umi_start=0, len_umi=12):
    """
    Extracts a Unique Molecular Identifier (UMI) from a read.

    Args:
        read (tuple): Tuple containing the read information (header, sequence, etc.).
        umi_start (int): Start position of the UMI in the sequence.
        len_umi (int): Length of the UMI.

    Returns:
        str: Extracted UMI or 'na' if the UMI length is insufficient.
    """
    return read[1][umi_start:umi_start + len_umi] if len(read[1]) >= len_umi else 'na'

def get_index(read):
    """
    Extracts index information from a read header.

    Args:
        read (tuple): Tuple containing the read information (header, sequence, etc.).

    Returns:
        tuple: Extracted index values (three integers).
    """
    line = read[0].split()[0]
    segs = line.split(":")
    return int(segs[-3]), int(segs[-2]), int(segs[-1])

def process_read(r1, r2, bcs, bc_idx, umi, p_args, workflow2_link):
    """
    Processes read information and modifies headers based on the receptor kit and UMI.

    Args:
        r1 (tuple): Tuple containing the R1 read information (header, sequence, etc.).
        r2 (tuple): Tuple containing the R2 read information (header, sequence, etc.).
        bcs (list): List of barcodes.
        bc_idx (int): Index of the barcode to use.
        umi (str): Unique Molecular Identifier.
        p_args (argparse.Namespace): Argument parser namespace with receptor kit and library kit information.
        workflow2_link (str): Path to workflow JAR file to determine processing.

    Returns:
        tuple: Modified R1 and R2 read tuples.
    """
    if p_args.receptor_kit in ['TCRv2', 'BCRv2', 'BCRv2sub', 'BCRv2hy'] or (p_args.receptor_kit == "BCRv1" and p_args.species == "human"):
        if 'ImmuneProfiler.jar' in workflow2_link and p_args.library_kit == 'takara_smartseq':
            bc = bcs[bc_idx].replace('"', '')
            r1[0] = ' '.join(r1[0].split() + ['_bc_' + bc + '_umi_' + umi])
            r2[0] = ' '.join(r2[0].split() + ['_bc_' + bc + '_umi_' + umi])
            r2[1] = bc + r2[1]
            r2[3] = 'I' * len(bc) + r2[3]
        elif "ImmuneProfilerv2.jar" in workflow2_link:
            tmp = r1[0].split()
            segs = tmp[0].split(':')
            segs[2] += '_umi_' + umi
            r1[0] = ' '.join(':'.join(segs))
            r2[0] = ' '.join(':'.join(segs))
        elif p_args.library_kit == 'takara_smartseq':
            r1[0] = ' '.join(r1[0].split() + ['_umi_' + umi])
            r2[0] = ' '.join(r2[0].split() + ['_umi_' + umi])

    return r1, r2

def correct_linker(read, start=12, linker='GTAC'):
    """
    Checks if the linker sequence in a read matches the expected linker.

    Args:
        read (tuple): Tuple containing the read information (header, sequence, etc.).
        start (int): Start position of the linker sequence in the read.
        linker (str): Expected linker sequence.

    Returns:
        bool: True if the linker sequence does not match, otherwise False.
    """
    return read[1][start:start + len(linker)] != linker

def get_subigindex(lines, p_s, p_e, umi):
    """
    Performs a binary search to find the position of a UMI in a list of lines.

    Args:
        lines (list): List of lines containing UMI information.
        p_s (int): Start index for the search.
        p_e (int): End index for the search.
        umi (str): Unique Molecular Identifier to find.

    Returns:
        int: Position of the UMI in the list of lines, or -1 if not found.
    """
    pos = -1
    mid = (p_s + p_e) // 2
    mid = (mid // 4) * 4

    if p_e > p_s:
        myumi = lines[mid + 1][:12]
        if umi > myumi:
            pos = get_subigindex(lines, mid, p_e, umi)
        elif umi < myumi:
            pos = get_subigindex(lines, p_s, mid, umi)
        else:
            pos = mid

    return pos

def get_umi_index(lines, pos):
    """
    Extracts UMI and index information from a list of lines based on a position.

    Args:
        lines (list): List of lines containing UMI and index information.
        pos (int): Position in the list to extract information from.

    Returns:
        list: List containing UMI and index values (three integers).
    """
    umi = lines[pos + 1][:12]
    line = lines[pos].split()[0]
    segs = line.split(":")
    return [umi, int(segs[-3]), int(segs[-2]), int(segs[-1])]

def do_merge2(read2, r2):
    """
    Merges two R2 read records (dummy implementation).

    Args:
        read2 (tuple): Tuple containing the R2 read information (header, sequence, etc.).
        r2 (list): List of lines to merge with R2 read.

    Returns:
        tuple: Merged R2 read tuple.
    """
    return read2

def do_merge1(read1, ref_subr1, pos):
    """
    Merges an R1 read record with a reference R1 read based on a position.

    Args:
        read1 (tuple): Tuple containing the R1 read information (header, sequence, etc.).
        ref_subr1 (str): Path to the reference R1 file.
        pos (int): Position in the reference file to merge with.

    Returns:
        tuple: Merged R1 read tuple.
    """
    with open(ref_subr1, "r") as fr:
        lines = fr.readlines()
    r1 = lines[pos:pos + 4]
    read1 = do_merge2(read1, r1)
    return read1

def merge_subIG_head(numl, fr1, fr2, sample, subseq, read1, read2, repo_type, species, pseq_dict, prs, len_cutoff=30, mismatch=1):
    """
    Merges two FASTQ files based on UMI information and sequence subsequence.

    Args:
        numl (int): Number of lines to read (not used in function).
        fr1 (file object): File object for the first FASTQ file.
        fr2 (file object): File object for the second FASTQ file.
        sample (str): Sample identifier.
        subseq (str): Subsequence to check for merging criteria.
        read1 (tuple): Tuple containing R1 read information (header, sequence, etc.).
        read2 (tuple): Tuple containing R2 read information (header, sequence, etc.).
        repo_type (str): Type of repository (not used in function).
        species (str): Species type (not used in function).
        pseq_dict (dict): Dictionary of possible sequences (not used in function).
        prs (dict): Dictionary of possible sequences and their original forms.
        len_cutoff (int): Minimum length cutoff for subsequences.
        mismatch (int): Allowed number of mismatches.

    Returns:
        tuple: Merged R1 and R2 read tuples.
    """
    infos = get_umi_index(read2, 0)
    lines2 = [fr2.readline() for _ in range(4)]
    lines1 = [fr1.readline() for _ in range(4)]
    myinfos = get_umi_index(lines2, 0)

    while myinfos[0] < infos[0]:
        lines2 = [fr2.readline() for _ in range(4)]
        lines1 = [fr1.readline() for _ in range(4)]
        myinfos = get_umi_index(lines2, 0)

    if myinfos[0] == infos[0]:
        do_merge2(read2, lines2)
        do_merge2(read1, lines1)
    else:
        print("get most similar r2")

    return read1, read2

def merge_subIG_mid(sample, subseq, read1, read2, repo_type, species, pseq_dict, prs, len_cutoff=30, mismatch=1):
    """
    Merges two FASTQ files based on UMI information and subsequence match.

    Args:
        sample (str): Sample identifier.
        subseq (str): Subsequence to check for merging criteria.
        read1 (tuple): Tuple containing R1 read information (header, sequence, etc.).
        read2 (tuple): Tuple containing R2 read information (header, sequence, etc.).
        repo_type (str): Type of repository (not used in function).
        species (str): Species type (not used in function).
        pseq_dict (dict): Dictionary of possible sequences (not used in function).
        prs (dict): Dictionary of possible sequences and their original forms.
        len_cutoff (int): Minimum length cutoff for subsequences.
        mismatch (int): Allowed number of mismatches.

    Returns:
        tuple: Merged R1 and R2 read tuples.
    """
    ref_subr1 = os.path.join(p_args.out_dir, "preprocess", sample + "_subIGG_R1_sorted.fastq")
    ref_subr2 = os.path.join(p_args.out_dir, "preprocess", sample + "_subIGG_R2_sorted.fastq")
    if subseq == "GGGAAGA":
        ref_subr1 = os.path.join(p_args.out_dir, "preprocess", sample + "_subIGA_R1_sorted.fastq")
        ref_subr2 = os.path.join(p_args.out_dir, "preprocess", sample + "_subIGA_R2_sorted.fastq")

    infos = get_umi_index(read2, 0)
    with open(ref_subr2, "r") as fr:
        lines = fr.readlines()

    pos = get_subigindex(lines, 0, len(lines) - 1, infos[0])
    if 0 <= pos < len(lines):
        while get_umi_index(lines, pos) < infos:
            pos += 4
        while get_umi_index(lines, pos) > infos:
            pos -= 4

        if get_umi_index(lines, pos) == infos:
            do_merge2(read2, lines[pos:pos + 4])
            do_merge1(read1, ref_subr1, pos)
        else:
            print("get most similar r2")

    else:
        print("have no the corresponding the sub r2; give up")

    return read1, read2

def save_sortedfastq(umis, r1_fq, r2_fq, p_args):
    """
    Saves sorted FASTQ files based on UMI information.

    Args:
        umis (list): List of UMIs and associated index values.
        r1_fq (str): Path to the R1 FASTQ file.
        r2_fq (str): Path to the R2 FASTQ file.
        p_args (argparse.Namespace): Argument parser namespace containing FASTQ directory information.

    Returns:
        None
    """
    for key in range(2):
        fpath = r2_fq if key == 0 else r1_fq
        fdname = os.path.dirname(fpath)
        fname = os.path.basename(fpath)
        fname_sorted = fname.replace(".fastq", "_sorted.fastq")

        with open(fpath, "r") as fr:
            lines = fr.readlines()

        with open(os.path.join(fdname, fname_sorted), 'w') as fw:
            for umi in umis:
                index = umi[-1]
                for j in range(4):
                    fw.write('%s\n' % lines[index * 4 + j].rstrip())

def write_sample_qc(p_args, sd):
    """
    Writes a sample QC statistics report to a CSV file.

    Args:
        p_args (argparse.Namespace): Argument parser namespace containing output directory and name information.
        sd (dict): Dictionary containing sample data for QC statistics.

    Returns:
        None
    """
    file_dir = p_args.out_dir
    file_name = p_args.out_name

    with open(os.path.join(file_dir, 'report', file_name + '_sample_QC_stats.csv'), 'w', newline='') as qc_file:
        writer = csv.writer(qc_file, delimiter=',', lineterminator="\n")
        col_names = ['Sample_ID']
        key1 = list(sd)[0]
        for item in sd[key1]:
            col_names.append(item[0].split('_')[0])
            col_names.append(item[0].split('_')[0] + '%')
        writer.writerow(col_names)
        
        for key, val in sd.items():
            sample_id = key
            l1 = [str(item[1]) for item in val]
            l2 = ['{:.1%}'.format(item[2]) for item in val]
            stats = list(chain.from_iterable(zip(l1, l2)))
            writer.writerow([sample_id] + stats)

def write_airr_sample_qc(p_args, sd):
    """
    Writes an AIRR sample QC statistics report to a CSV file.

    Args:
        p_args (argparse.Namespace): Argument parser namespace containing output directory and name information.
        sd (dict): Dictionary containing sample data for QC statistics.

    Returns:
        None
    """
    file_dir = p_args.out_dir
    file_name = p_args.out_name

    with open(os.path.join(file_dir, 'airr_report', file_name + '_sample_QC_stats.csv'), 'w', newline='') as qc_file:
        writer = csv.writer(qc_file, delimiter=',', lineterminator="\n")
        col_names = ['organism_id', 'sample_processing_id']
        key1 = list(sd)[0]
        for item in sd[key1]:
            col_names.append(item[0].split('_')[0])
            col_names.append(item[0].split('_')[0] + '%')
        writer.writerow(col_names)
        
        organism = p_args.species
        for key, val in sd.items():
            sample_id = key
            l1 = [str(item[1]) for item in val]
            l2 = ['{:.1%}'.format(item[2]) for item in val]
            stats = list(chain.from_iterable(zip(l1, l2)))
            writer.writerow([organism] + [sample_id] + stats)

def create_fd(sd):
    """
    Creates a file dictionary mapping sample names to their respective FASTQ files.

    Args:
        sd (dict): Dictionary containing sample data.

    Returns:
        OrderedDict: Mapping of sample names to their respective FASTQ files and barcodes.
    """
    sublist = ['short', 'undetermined', 'flc', 'total']
    if workflow2_link.find('ImmuneProfilerv2.jar') >= 0 and p_args.internal_test:
        sublist = ['short', 'total']

    pd = OrderedDict((k, sd[k]) for k in sd.keys())
    for k, v in pd.items():
        l1 = [i for i in range(len(v)) if v[i][1] > 0]
        l2 = [i for i in range(len(v)) if v[i][0] not in sublist]
        l = [value for value in l1 if value in l2]
        pd[k] = [v[idx] for idx in l]

    fd = OrderedDict()
    for k, v in pd.items():
        sample = k
        bc = ''
        for i in range(len(v)):
            chain = v[i][0].split('_')[0]
            if chain in ['undetermined', 'flc']:
                bc = chain
            else:
                bc = v[i][0].split('_')[1]
            name = sample + '_' + chain
            fd[name] = [name + '_R1.fastq', name + '_R2.fastq', bc]

    return fd

def dos2unix(fq_file):
    """
    Converts a file from DOS (Windows) format to Unix (Linux) format.

    Args:
        fq_file (str): Path to the file to be converted.

    Description:
        This function reads the file specified by `fq_file`, which may contain
        line endings in DOS format (`\r\n`). It writes the content back to the file
        with Unix line endings (`\n`), which is the format required by MiXCR for
        processing FASTQ files.

    Exceptions:
        EnvironmentError: If there is an issue opening or writing to the file.
    """
    try:
        text = open(fq_file).read()
        open(fq_file, 'w', newline='\n').write(text)
    except EnvironmentError as err:
        print('error', 'Unable to open file: ' + fq_file + '\n', flush=True)
    return

def prep_migec_input(p_args, fd):
    """
    Prepares the input file for MiGEC by creating a barcodes.txt file.

    Args:
        p_args (argparse.Namespace): Namespace containing output directory and other arguments.
        fd (dict): Dictionary where keys are sample names and values are lists containing
                   file names and barcodes.

    Description:
        This function generates a tab-delimited file named 'barcodes.txt' in the
        'run_migec' directory. Each line in this file includes sample information,
        UMI barcodes, and paths to the R1 and R2 FASTQ files.

    Returns:
        None
    """
    with open(os.path.join(p_args.out_dir, 'run_migec', 'barcodes.txt'), 'w', newline='') as bc_file:
        writer = csv.writer(bc_file, delimiter='\t')
        for key, val in fd.items():
            sample = key
            r1_fq = val[0]
            r2_fq = val[1]
            bc = val[2]
            writer.writerow(
                [sample] +
                [bc + 'N' * umi_length] + [''] +
                [os.path.join(p_args.out_dir, 'preprocess', r2_fq)] +
                [os.path.join(p_args.out_dir, 'preprocess', r1_fq)]
            )

def run_system_cmd(cmd, outfile, errfile, err_message, logname='NA'):
    """
    Executes a system command and logs the output.

    Args:
        cmd (str): The system command to run.
        outfile (str): Path to the file where standard output will be written.
        errfile (str): Path to the file where standard error will be written.
        err_message (str): Error message to log if the command fails.
        logname (str, optional): Logger name for error logging. Default is 'NA'.

    Description:
        This function runs a system command using `subprocess.run()`, directing
        the standard output and error to specified files. If the command fails,
        it logs the error message using either a logger or standard print, based
        on the `logname` parameter.

    Returns:
        None
    """
    with open(outfile, 'w') as out, open(errfile, 'w') as err:
        try:
            subprocess.run(cmd.split(), stdout=out, stderr=err, check=True)
        except subprocess.CalledProcessError:
            if logname == 'NA':
                print(err_message)
            else:
                logger = logging.getLogger(logname)
                logger.error(err_message + "\n" + str(err))
            sys.exit(1)
    return

def get_num_reads(fq):
    """
    Counts the number of reads in a FASTQ file.

    Args:
        fq (str): Path to the FASTQ file.

    Returns:
        int: Number of reads in the FASTQ file.

    Description:
        This function counts the number of reads in the FASTQ file by reading all
        lines and dividing the total number of lines by 4 (each read is represented
        by 4 lines in a FASTQ file).
    """
    with open(fq, 'r') as fp:
        num_reads = len(fp.readlines()) // 4
        return num_reads

def run_migec(p_args, migec_dir, log_name, errfile_name):
    """
    Executes the MiGEC workflow for processing FASTQ files, including checking out, generating histograms,
    and assembling results.

    Args:
        p_args (argparse.Namespace): Namespace containing the output directory, receptor kit, species, and UMI cutoff.
        migec_dir (str): Subdirectory within the output directory where MiGEC results will be stored.
        log_name (str): Name of the log file for MiGEC operations.
        errfile_name (str): Name of the error file for MiGEC operations.

    Description:
        This function orchestrates the MiGEC workflow by:
        1. Running the MiGEC checkout process.
        2. Generating histograms from the checkout results.
        3. Assembling results from the histogram data.
        It handles special cases for different receptor kits and species, and performs file operations specific
        to the operating system (Windows vs Unix).

        - **Checkout Phase**: 
          Constructs and executes a command to perform the MiGEC checkout. If the receptor kit is "TCRv1" or
          "BCRv1" and the species is "mouse", it copies files into the 'checkout_all' folder. Otherwise, it runs
          the checkout command using `run_system_cmd()`.

        - **Histogram Generation**:
          Constructs and executes a command to generate histograms from the checkout results. If the receptor kit is
          "TCRv1" or "BCRv1" and the species is "mouse", it creates an 'estimates.txt' file with sample information.
          Otherwise, it runs the histogram command using `run_system_cmd()`.

        - **Assembly Phase**:
          Constructs and executes a command to assemble the results. If a UMI cutoff is specified, it includes it in
          the command. If the receptor kit is "TCRv1" or "BCRv1" and the species is "mouse", it copies files into
          the 'assemble' folder. Otherwise, it runs the assemble command using `run_system_cmd()`.

        - **File Format Handling**:
          For Windows systems, it converts FASTQ files in the 'assemble' folder from DOS to Unix format using `dos2unix()`.
          For Unix systems, no conversion is performed.

        - **Cleanup**:
          If `p_args.keep_inter_file` is `False`, it deletes FASTQ files from the 'checkout_all' directory.

    Returns:
        None
    """
    out_log = os.path.join(p_args.out_dir, migec_dir, log_name)
    run_migec_err = os.path.join(p_args.out_dir, migec_dir, errfile_name)
    barcode_file = os.path.join(p_args.out_dir, migec_dir, 'barcodes.txt')
    checkout_folder = os.path.join(p_args.out_dir, migec_dir, 'checkout_all')
    checkout_cmd = workflow1_link + ' workflow1 CheckoutBatch -uto --skip-undef ' + barcode_file + ' ' + checkout_folder
    logger.info('Launching MIGEC check out using command:\n' + checkout_cmd)
    print(get_time() + ' [INFO] Processing fastq file check')
    
    if (p_args.receptor_kit == "TCRv1") or (p_args.receptor_kit == "BCRv1" and p_args.species == "mouse"):
        print("copy file into folder checkout_all")
        for key, value in fd.items():
            src = os.path.join(p_args.out_dir, "preprocess", value[0])
            dst = os.path.join(p_args.out_dir, "run_migec", "checkout_all", value[0])
            if not os.path.isdir(os.path.join(p_args.out_dir, "run_migec", "checkout_all")):
                os.makedirs(os.path.join(p_args.out_dir, "run_migec", "checkout_all"))
            copyfile(src, dst)
    else:
        run_system_cmd(checkout_cmd, out_log, run_migec_err, 'checkout failed', LOG_NAME)
    
    histogram_folder = os.path.join(p_args.out_dir, migec_dir, 'histogram')
    histogram_cmd = workflow1_link + ' workflow1 Histogram ' + checkout_folder + ' ' + histogram_folder
    logger.info('Launching MIGEC histogram generation using command:\n' + histogram_cmd)
    print(get_time() + ' [INFO] Processing UMI check')
    
    if (p_args.receptor_kit == "TCRv1") or (p_args.receptor_kit == "BCRv1" and p_args.species == "mouse"):
        print("create estimates.txt file in folder histogram")
        if not os.path.isdir(os.path.join(p_args.out_dir, "run_migec", "histogram")):
            os.makedirs(os.path.join(p_args.out_dir, "run_migec", "histogram"))
        dst = os.path.join(p_args.out_dir, "run_migec", "histogram", "estimates.txt")
        line = "#SAMPLE_ID\tSAMPLE_TYPE\tTOTAL_READS\tTOATL_MIG\tOVERSEQ_THRESHOLD\tCOLLISION_THRESHOLD\tUMI_QUAL_THRESHOLD\tUMI_LEN"
        with open(dst, 'w') as fr:
            fr.write("%s\n" % line)
            for key, value in fd.items():
                total_reads = get_num_reads(os.path.join(p_args.out_dir, "run_migec", "checkout_all", value[0]))
                line = key + "\tpaired\t" + str(total_reads) + "\t" + str(total_reads) + "\t1\t1\t1\t0"
                fr.write("%s\n" % line)
    else:
        run_system_cmd(histogram_cmd, out_log, run_migec_err, 'histogram failed', LOG_NAME)
    
    assemble_folder = os.path.join(p_args.out_dir, migec_dir, 'assemble')
    if p_args.umi_cutoff == '':
        assemble_cmd = workflow1_link + ' workflow1 AssembleBatch ' + checkout_folder + ' ' + histogram_folder + ' ' + assemble_folder
    elif p_args.umi_cutoff.isdigit():
        assemble_cmd = workflow1_link + ' workflow1 AssembleBatch --force-overseq ' + str(p_args.umi_cutoff) + \
                       ' ' + checkout_folder + ' ' + histogram_folder + ' ' + assemble_folder
    logger.info('Launching MIGEC assembling using command:\n' + assemble_cmd)
    print(get_time() + ' [INFO] Processing consensus read assembly')
    
    if (p_args.receptor_kit == "TCRv1") or (p_args.receptor_kit == "BCRv1" and p_args.species == "mouse"):
        print("create assemble.txt file in folder assemble")
        if not os.path.isdir(os.path.join(p_args.out_dir, "run_migec", "assemble")):
            os.makedirs(os.path.join(p_args.out_dir, "run_migec", "assemble"))
        dst = os.path.join(p_args.out_dir, "run_migec", "assemble", "assemble.txt")
        line = "#SAMPLE_ID\tUMI_CUT_OFF\tTCR_THRESHOLD\tBCR_THRESHOLD"
        with open(dst, 'w') as fr:
            fr.write("%s\n" % line)
            for key, value in fd.items():
                line = key + "\t1\t0\t0"
                fr.write("%s\n" % line)
    else:
        run_system_cmd(assemble_cmd, out_log, run_migec_err, 'assemble failed', LOG_NAME)
    
    if platform.system() == 'Windows':
        migec_fqs = []
        for f in os.listdir(os.path.join(p_args.out_dir, 'run_migec', 'assemble')):
            if f.endswith('.fastq'):
                migec_fqs.append(f)
        for fq in migec_fqs:
            dos2unix(os.path.join(p_args.out_dir, 'run_migec', 'assemble', fq))
    else:
        pass
    
    if not p_args.keep_inter_file:
        fq_list = find_file('*.fastq', os.path.join(p_args.out_dir, 'run_migec', 'checkout_all'))
        if len(fq_list) > 0:
            for item in fq_list:
                os.remove(item)

def summ_migec(p_args):
    """
    Summarizes MIGEC (Molecular Identifier Groups-based Error Correction) results by extracting and 
    appending relevant statistics from the 'estimates.txt' file into a dictionary.

    Args:
        p_args (argparse.Namespace): An object containing command-line arguments, including the output directory 
                                      and UMI (Unique Molecular Identifier) cutoff.

    Returns:
        None: The function updates the 'fd' dictionary with data from the 'estimates.txt' file.
    """
    if os.path.isfile(os.path.join(p_args.out_dir, 'run_migec', 'histogram', 'estimates.txt')):
        infile = open(os.path.join(p_args.out_dir, 'run_migec', 'histogram', 'estimates.txt'), 'r')
        next(infile)  # Skip the header line
        for line in infile:
            items = line.rstrip().split('\t')
            sample = items[0]
            total_reads = items[2]
            total_migs = items[3]
            overseq_threshold = items[4]
            collision_threshold = items[5]
            umi_qual_threshold = items[6]
            if sample in fd.keys():
                fd[sample].append(total_reads)
                fd[sample].append(total_migs)
                if p_args.umi_cutoff == '':
                    fd[sample].append(overseq_threshold)
                else:
                    fd[sample].append(str(p_args.umi_cutoff))
                fd[sample].append(collision_threshold)
                fd[sample].append(umi_qual_threshold)
        infile.close()

def find_file(pattern, path):
    """
    Recursively searches for files matching a specified pattern within a directory.

    Args:
        pattern (str): The filename pattern to match.
        path (str): The directory path to search within.

    Returns:
        list: A list of file paths that match the given pattern.
    """
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def getsubdirsandfiles(distfd):
    """
    Retrieves the list of subdirectories and files within a specified directory.

    Args:
        distfd (str): The directory path to list contents from.

    Returns:
        tuple: A tuple containing two lists - one for subdirectories and one for files.
    """
    complist = os.listdir(distfd)
    complist.sort()
    dirslist = []
    fileslist = []
    for i in range(len(complist)):
        fullpath = os.path.join(distfd, complist[i])
        if os.path.isfile(fullpath):
            fileslist.append(complist[i])
        elif os.path.isdir(fullpath):
            dirslist.append(complist[i])
    return dirslist, fileslist

def get_umi_cutoff(refinefile):
    """
    Extracts the UMI (Unique Molecular Identifier) cutoff value from a refine report file.

    Args:
        refinefile (str): The path to the refine report file.

    Returns:
        str: The extracted UMI cutoff value.
    """
    umi_cutoff = 1
    f = open(refinefile, 'r')
    lines = f.readlines()
    f.close()
    line = lines[len(lines)-2]
    if line.find("Effective threshold:") >= 0:
        umi_cutoff = str(int(float(line.split(": ")[1].split('\n')[0])))
    return umi_cutoff

def summ_results(analysis_dir, folder, inputDict, align_stats_dict, clone_stats_dict):
    """
    Summarizes the results of sequence alignment and clonotype assembly by extracting and updating relevant statistics.

    Args:
        analysis_dir (str): The directory containing the analysis results.
        folder (str): The specific folder within the analysis directory to process.
        inputDict (dict): A dictionary containing input data to be updated with summary results.
        align_stats_dict (dict): A dictionary template for storing alignment statistics.
        clone_stats_dict (dict): A dictionary template for storing clonotype statistics.

    Returns:
        dict: The updated inputDict containing alignment and clonotype statistics.
    """
    full_path = os.path.join(analysis_dir, 'run_mixcr', folder)
    alignkey = '.align.report.txt'
    cloneskey = '.assemble.report.txt'
    clonesall = '.clones_*.tsv'
    refinekey = '.refine.report.txt'
    if workflow2_link.find("ImmuneProfilerv1.6.jar") >= 0:
        alignkey = '_align_report.txt'
        cloneskey = '_clones_report.txt'
        clonesall = '_clones_all.txt'
    for key, val in inputDict.items():
        align_report = find_file(key + '_' + folder + alignkey, full_path)[0]
        infile1 = open(align_report, 'r')
        cur_align_dict = copy.deepcopy(align_stats_dict)
        for line in infile1:
            line_content = line.rstrip().split(':')
            if line_content[0] in cur_align_dict.keys():
                cur_align_dict[line_content[0]] += ':'
                cur_align_dict[line_content[0]] += line_content[1]
        for subkey, subvalue in cur_align_dict.items():
            inputDict[key].append(subvalue)
            if subkey == "Total sequencing reads" and val[3] == None:
                val[3] = subvalue.split(": ")[1]
            elif subkey == "Successfully aligned reads" and val[4] == None:
                val[4] = subvalue.split(": ")[1].split(" ")[0]
        if val[5] == None:
            if p_args.umi_cutoff == '':
                files = find_file(key + '_' + folder + refinekey, full_path)
                if len(files) > 0:
                    val[5] = val[6] = get_umi_cutoff(files[0])
            else:
                val[5] = val[6] = p_args.umi_cutoff
            val[7] = "15"
        infile1.close()

        clone_report = find_file(key + '_' + folder + cloneskey, full_path)[0]
        infile2 = open(clone_report, 'r')
        cur_clone_dict = copy.deepcopy(clone_stats_dict)
        for line in infile2:
            line_content = line.rstrip().split(':')
            if line_content[0] in cur_clone_dict.keys():
                cur_clone_dict[line_content[0]] += ':'
                cur_clone_dict[line_content[0]] += line_content[1]
        for subkey, subvalue in cur_clone_dict.items():
            inputDict[key].append(subvalue)
        infile2.close()

        tmp = find_file(key + '_' + folder + clonesall, full_path)
        if len(tmp) == 0:
            logger.info('No clonos_all file generated stop analysis \n')
            next
        else:
            for i in range(len(tmp)):
                clones_infile = tmp[i]
                clones_outfile = clones_infile[:len(clones_infile)-4] + '.csv'  # Output as CSV format
                if p_args.cogentIP_version.find('v2') < 0:
                    clones_outfile = clones_infile.replace('_clones_all.txt', '_clones_result.csv')
                parse_clones_report(clones_infile, clones_outfile, key, inputDict)
    return inputDict

def report_stats(r_type, file_dir, file_name, metaDict):
    """
    Generates a CSV report of mapping statistics based on summarized data.

    Args:
        r_type (str): The receptor type (e.g., 'BCRv1', 'TCRv2').
        file_dir (str): The directory where the report should be saved.
        file_name (str): The base name of the report file.
        metaDict (dict): A dictionary containing metadata and statistics for each sample.

    Returns:
        None: The function writes the report to a CSV file.
    """
    subheads = ['total MIG', 'UMI threshold', 'number of reads after MIG collapse']
    if list(metaDict.values())[0][5] == None or p_args.receptor_kit == 'TCRv1' or (p_args.receptor_kit == 'BCRv1' and p_args.species == 'mouse'):
        subheads = ['aligned reads', 'without UMI', 'number of reads available in file']
    elif workflow2_link.find('ImmuneProfilerv2.jar') >= 0:
        subheads = ['aligned reads', 'UMI threshold', 'number of reads available in file']
    with open(os.path.join(file_dir, file_name + '_mapping_stats.csv'), 'w', newline='') as stats_file:
        writer = csv.writer(stats_file, delimiter=',', lineterminator="\n")  # Write as CSV format
        heads = ['sample type', 'total read'] + subheads + ['aligned', 'pair-read overlap', 'overlapped and aligned']
        if r_type in ['BCRv1', 'BCRv2', 'BCRv2sub', 'BCRv2hy']:
            if p_args.internal_test:
                heads += ['reads per clonotype', 'reads in clonotypes', 'mapped low quality reads', 'pcr error correction', 
                          'failed mapping', 'clonotype count', 'IgG', 'IgM', 'IgK', 'IgL', 'IgA', 'IgD', 'IgE', 'IgH(lack constant region)']
            else:
                heads += ['clonotype count', 'IgG', 'IgM', 'IgK', 'IgL', 'IgA', 'IgD', 'IgE', 'IgH(lack constant region)']
        elif r_type in ['TCRv1', 'TCRv2']:
            if p_args.internal_test:
                heads += ['reads per clonotype', 'reads in clonotypes', 'mapped low quality reads', 'pcr error correction', 
                          'failed mapping', 'clonotype count', 'TRA', 'TRB']
            else:
                heads += ['clonotype count', 'TRA', 'TRB']
        writer.writerow(heads)
        lheads = len(heads)
        for key, val in metaDict.items():
            sample_id = key
            stats2report = val[3:6] + [item.split(':')[1].split('(')[0].replace(' ', '') for item in val[8:]]
            if len(stats2report) != lheads-1:
                for i in range(len(stats2report)):
                    if stats2report[i] == None or not stats2report[i].isnumeric():
                        stats2report[i] = '0'
                stats2report = stats2report + [0] * (lheads-1-len(stats2report))
            writer.writerow([sample_id] + stats2report)

def stats_csv_to_airr(in_csv, airr_csv):
    """
    Converts a CSV file to AIRR format by adding headers and relevant data.
    
    Parameters:
    in_csv (str): Path to the input CSV file.
    airr_csv (str): Path to the output AIRR-formatted CSV file.
    """
    rows = []
    with open(in_csv, 'rt', encoding='utf8') as f:
        reader = csv.reader(f, delimiter=',')
        count = 0
        for row in reader:
            if count == 0:
                # Modify the header row
                row = ["organism_id", "sample_processing_id"] + airrhead(row)
            else:
                # Extract sample ID and prepend relevant data
                sid = row[0].split("_")[0]
                row = [p_args.species, sid] + row
            count += 1
            rows.append(row)
    save_csv(airr_csv, rows)

def save_csv(file_path, datarows):
    """
    Saves data to a CSV file.
    
    Parameters:
    file_path (str): Path to the output CSV file.
    datarows (list): List of rows to write to the CSV file.
    """
    with open(file_path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',', lineterminator="\n")
        for row in datarows:
            writer.writerow(row)

def airrhead(items):
    """
    Processes the header row for AIRR formatting by modifying and standardizing header names.
    
    Parameters:
    items (list): List of header names.
    
    Returns:
    list: Processed header names.
    """
    items = [sub.replace(' ', '_') for sub in items]
    items = [sub.lower() for sub in items]
    for i in range(len(items)):
        # Specific header name replacements
        if items[i] == "amino_acid_length":
            items[i] = "junction_length_aa"
        elif items[i] == "igh(lack_constant_region)":
            items[i] = "IGH(lack_constant_region)"
        if items[i] in ["igh", "igi", "igk", "igl"]:
            items[i] = items[i].upper()
        elif items[i] in ["iga", "igd", "ige", "igg", "igm"]:
            items[i] = items[i].upper()
        elif items[i] in ["tra", "trb", "trd", "trg"]:
            items[i] = items[i].upper()
    return items

def merge_airr_csv(sample_id, process_type, report_name, file_dir, output_dir, airr_dir):
    """
    Merges multiple CSV files into an AIRR-formatted CSV and recalculates fractions.
    
    Parameters:
    sample_id (str): ID of the sample.
    process_type (str): Type of processing.
    report_name (str): Name of the report to generate.
    file_dir (str): Directory containing input files.
    output_dir (str): Directory for output files.
    airr_dir (str): Directory for AIRR files.
    """
    mergedflist = []
    stats_file = find_file('*_mapping_stats.csv', file_dir)[0]
    with open(stats_file, 'rt', encoding='utf8') as f:
        reader = csv.reader(f, delimiter=',')
        file_prefix_list = []
        ct = 0
        for row in reader:
            if ct == 0:
                pass  # Skip the header row
            elif row[0].split('_')[0] == sample_id:
                file_prefix_list.append(row[0])
            ct += 1
    
    for key in file_prefix_list:
        sid = key.split('_')[0]
        airr_report_path = os.path.join(p_args.out_dir, 'airr_report', sid)
        if not os.path.isdir(airr_report_path):
            os.makedirs(airr_report_path)
        
        infile = find_file(key + '_' + process_type + '_clones_*.csv', output_dir)[0]
        frowlist = []
        with open(infile, 'rt', encoding='utf8') as f:
            reader = csv.reader(f, delimiter=',')
            i = 0
            for row in reader:
                if i == 0:
                    row = ["organism_id", "sample_processing_id"] + airrhead(row)
                else:
                    row = [p_args.species, sid] + row
                frowlist.append(row)
                i += 1
        save_csv(os.path.join(airr_report_path, key + '_' + process_type + '_clones_result.csv'), frowlist)
        
        if not mergedflist:
            mergedflist = frowlist
        else:
            mergedflist += frowlist[1:]  # Append data excluding the header row
    
    # Recalculate the fractions
    total_reads = sum(float(row[2]) for row in mergedflist[1:])
    for i in range(1, len(mergedflist)):
        mergedflist[i][3] = float(mergedflist[i][2]) / total_reads
    
    output_file = os.path.join(p_args.out_dir, 'airr_report', p_args.out_name + '_' + sid + '_' + process_type + '_report.csv')
    save_csv(output_file, mergedflist)

def merge_csv(sample_id, process_type, report_name, file_dir, output_dir):
    """
    Merges multiple CSV files into an Excel workbook with multiple sheets.
    
    Parameters:
    sample_id (str): ID of the sample.
    process_type (str): Type of processing.
    report_name (str): Name of the report to generate.
    file_dir (str): Directory containing input files.
    output_dir (str): Directory for output files.
    
    Returns:
    str: Path to the statistics file.
    """
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = 'stats'
    
    stats_file = find_file('*_mapping_stats.csv', file_dir)[0]
    with open(stats_file, 'rt', encoding='utf8') as f:
        reader = csv.reader(f, delimiter=',')
        file_prefix_list = []
        ct = 0
        for row in reader:
            if ct == 0:
                ws.append(row)  # Keep header row
            elif row[0].split('_')[0] == sample_id:
                ws.append(row)
                file_prefix_list.append(row[0])
            ct += 1
    
    for key in file_prefix_list:
        sid = key.split('_')[0]
        report_path = os.path.join(p_args.out_dir, 'report', sid)
        if not os.path.isdir(report_path):
            os.makedirs(report_path)
        
        pattern = '_clones_*.csv'
        if workflow2_link.find("ImmuneProfilerv2.jar") >= 0:
            pattern = '.clones_*.csv'
        
        infiles = find_file(key + '_' + process_type + pattern, file_dir)
        for infile in infiles:
            sheet_name = key + '_clone'
            if key.find("_TR") < 0 and key.find("_IG") < 0:
                sheet_name = key + infile[-9:-4] + '_clone'
            
            ws = wb.create_sheet(sheet_name)
            with open(infile, 'rt', encoding='utf8') as f:
                reader = csv.reader(f, delimiter=',')
                for row in reader:
                    ws.append(row)
            
            os.rename(infile, os.path.join(report_path, os.path.basename(infile)))
    
    output_file = os.path.join(output_dir, report_name + '_report.xlsx')
    wb.save(output_file)
    
    return stats_file

def getsubIg_dict(input_dict, subIg):
    """
    Filters a dictionary to retain only keys that contain a specific substring.
    
    Parameters:
    input_dict (dict): Dictionary to filter.
    subIg (str): Substring to search for in the keys.
    
    Returns:
    dict: Filtered dictionary containing only the matching keys.
    """
    return {key: val for key, val in input_dict.items() if subIg.upper() in key}

def findmatchedpos(seg, subs, fixpos):
    """
    Finds the position of a substring in a string, with some flexibility.
    
    Parameters:
    seg (str): The string to search within.
    subs (str): The substring to search for.
    fixpos (int): Position within the substring to keep fixed.
    
    Returns:
    int: The last position where the substring (or a variation) is found, or -1 if not found.
    """
    poslist = []
    pos = seg.find(subs)
    while pos >= 0:
        poslist.append(pos)
        pos = seg.find(subs, pos + 1)
    
    if not poslist:
        # Generate variations of the substring
        l_subs = len(subs)
        atcg = ["A", "T", "C", "G"]
        for count in range(l_subs):
            if count != fixpos:
                for base in atcg:
                    temp = subs[:count] + base + subs[count+1:]
                    pos = seg.find(temp)
                    while pos >= 0:
                        if pos not in poslist:
                            poslist.append(pos)
                        pos = seg.find(temp, pos + 1)
    
    poslist.sort()
    return poslist[-1] if poslist else -1

def inbin(v0, v1, s, l):
    """
    Checks if the difference between two values (v1 - v0) lies within a specified range [s, l].

    Parameters:
    v0 (int): The first value.
    v1 (int): The second value.
    s (int): The minimum difference allowed.
    l (int): The maximum difference allowed.

    Returns:
    bool: True if the difference between v1 and v0 is within the range [s, l], otherwise False.
    """
    if (v0 >= 0 and v1 >= 0 and v1 - v0 >= s and v1 - v0 <= l):
        return True
    else:
        return False

def adjustPs(p0, p1, p2, p3, p4):
    """
    Adjusts positions p0 and p1 based on their relations to other positions p2, p3, and p4,
    ensuring that the pairs satisfy certain conditions defined by the inbin function.

    Parameters:
    p0 (int): First position.
    p1 (int): Second position.
    p2 (int): Third position.
    p3 (int): Fourth position.
    p4 (int): Fifth position.

    Returns:
    tuple: The adjusted positions (p0, p1, p2, p3, p4).
    """
    if p0 >= 0 and p1 >= 0:
        if inbin(p1, p2, 44, 46) or inbin(p1, p3, 44, 46) or inbin(p1, p4, 44, 46):
            p0 = -1
        elif inbin(p0, p2, 44, 46) or inbin(p0, p3, 44, 46) or inbin(p0, p4, 44, 46):
            p1 = -1
    elif p0 >= 0:  # p1 < 0
        if not inbin(p0, p2, 44, 46) and not inbin(p0, p3, 44, 46) and not inbin(p0, p4, 44, 46):
            p0 = -1
    elif p1 >= 0:  # p0 < 0
        if not inbin(p1, p2, 44, 46) and not inbin(p1, p3, 44, 46) and not inbin(p1, p4, 44, 46):
            p1 = -1
    
    return p0, p1, p2, p3, p4

def get_subIgG(seg, count):
    """
    Identifies the subtype of IgG (Immunoglobulin G) based on the segment (seg) by matching 
    specific patterns and adjusting positions.

    Parameters:
    seg (str): The segment string to analyze.
    count (int): An index or counter used for processing.

    Returns:
    str: The identified IgG subtype (e.g., 'IGG3/4', 'IGG1/2', 'IGG4', 'IGG2', 'IGG1', 'IGG3').
    """
    subIgG = 'IGG'
    labs = ['GCTTCCAC','GCCTCCAC','AGGAGCACCTCT','AAGAGCACCTCT','AGGAGCACCTCC']
    p0 = findmatchedpos(seg, labs[0], 2)
    p1 = findmatchedpos(seg, labs[1], 2)
    p2 = findmatchedpos(seg, labs[2], 1)
    p3 = findmatchedpos(seg, labs[3], 1)
    p4 = findmatchedpos(seg, labs[4], 11)
    p0, p1, p2, p3, p4 = adjustPs(p0, p1, p2, p3, p4)
    
    if p0 >= 0:
        subIgG = 'IGG3/4'
    elif p1 >= 0:
        subIgG = 'IGG1/2'
    
    if (p0 >= 0 or p1 >= 0) and (p2 >= 0 or p3 >= 0 or p4 >= 0):
        if inbin(p0, p4, 44, 46):
            subIgG = 'IGG4'
        elif inbin(p1, p4, 44, 46):
            subIgG = 'IGG2'
        elif inbin(p0, p3, 44, 46) or inbin(p1, p3, 44, 46):
            subIgG = 'IGG1'
        elif inbin(p0, p2, 44, 46) or inbin(p1, p2, 44, 46):
            subIgG = 'IGG3'
    
    return subIgG

def analyzeallCAlignments(p_args, segs):
    """
    Analyzes all C alignments from segment strings to identify the IgA subtype based on specific patterns.

    Parameters:
    p_args (argparse.Namespace): Argument parser containing parameters and settings.
    segs (list): List of segment strings to analyze.

    Returns:
    str: The identified IgA subtype.
    """
    subIgAtype = []
    subIgAvalue = []
    pos = 5
    if workflow2_link.find('immune_profilerv2.java') >= 0:
        pos += 2
    allcalignments = segs[pos]
    segs = allcalignments.split(',')
    for i in range(len(segs)):
        items = segs[i].split('*00(')
        subIgAtype.append(items[0].replace('IGHA','IGA'))
        subIgAvalue.append(int(items[1][:len(items[1])-1]))
    subIgA = subIgAtype[0]
    for i in range(1, len(subIgAtype)):
        if subIgAvalue[i] == subIgAvalue[i-1]:
            subIgA += ',' + subIgAtype[i]
    
    return subIgA

def get_subIg(p_args, aligns, key, read_type, target_region, mixcr_dir, clonesall):
    """
    Processes alignment data to identify the immunoglobulin subtype (IgG or IgA) and writes the results to a file.

    Parameters:
    p_args (argparse.Namespace): Argument parser containing parameters and settings.
    aligns (str): Path to the alignments file.
    key (str): Key or identifier related to the alignment.
    read_type (str): Type of read being processed.
    target_region (str): Target region of interest.
    mixcr_dir (str): Directory containing MiXCR results.
    clonesall (str): Additional data related to clones (not used directly in this function).

    Returns:
    None: The function writes the identified subtypes to an output file.
    """
    process_type = f"{read_type}_{target_region}"
    subIgfile = os.path.join(p_args.out_dir, mixcr_dir, process_type,
                             f"{key}_{process_type}_subIg{key[-1]}.tsv")
    subIG = 'subIGG'
    if key.find('IGA') >= 0:
        subIG = 'subIGA'
    
    if not os.path.exists(aligns) and p_args.cogentIP_version == 'v2':
        aligns = aligns.replace('_alignments.tsv', '_readIds.tsv')
    
    with open(aligns, 'r') as fr:
        line = fr.readline().replace('targetSequences', subIG).replace('targetQualities', 'targetSequences')
        with open(subIgfile, 'w') as fw:
            fw.write(line)
            count = 0
            while True:
                line = fr.readline()
                if line == "":
                    break
                segs = line.split('\t')
                if key.find('IGG') >= 0:
                    subIG = get_subIgG(segs[0], count)
                elif key.find('IGA') >= 0:
                    subIG = analyzeallCAlignments(p_args, segs)
                line = subIG + '\t' + segs[0]
                for i in range(2, len(segs)):
                    line += '\t' + segs[i]
                fw.write(line)
                count += 1

def get_rightindex(miglines, mig_pos, asslines, ass_pos):
    """
    Aligns the positions in the migration lines and assembly lines, adjusting indices as needed.

    Parameters:
    miglines (list): List of migration lines.
    mig_pos (int): Current position in the migration lines.
    asslines (list): List of assembly lines.
    ass_pos (int): Current position in the assembly lines.

    Returns:
    tuple: Adjusted positions (mig_pos, ass_pos).
    """
    print("False:")
    print(f"{miglines[mig_pos+1][:30]}")
    print(f"{asslines[ass_pos][:30]}")
    print(f"mig_pos={mig_pos}")
    print(f"ass_pos={ass_pos}")
    
    if asslines[ass_pos][:8] != "GTACGGGG":
        return mig_pos + 4, ass_pos + 1
    
    mig_p = mig_pos
    ass_p = ass_pos
    
    while miglines[mig_p + 1][:30] != asslines[ass_p][:30] and ass_p < ass_pos + 30 and ass_p < len(asslines):
        ass_p += 1
    
    if ass_p < ass_pos + 30 and ass_p < len(asslines):
        ass_pos = ass_p
    else:
        while miglines[mig_p + 1][:30] != asslines[ass_p][:30] and mig_p < mig_pos + 30 and mig_p < len(miglines):
            mig_p += 4
        if mig_p < mig_pos + 30 and mig_p < len(miglines):
            mig_pos = mig_p
        else:
            ass_pos += 1
            mig_pos += 4
    
    print(f"new mig_pos={mig_pos}")
    print(f"new ass_pos={ass_pos}")
    return mig_pos, ass_pos

def scanindex(mig1line, mig2line, asslines, p_thread):
    """
    Scans for the correct index in the assembly lines that matches with the migration lines.

    Parameters:
    mig1line (str): First migration line.
    mig2line (str): Second migration line.
    asslines (list): List of assembly lines.
    p_thread (int): Thread index or position.

    Returns:
    int: The position in assembly lines that matches with migration lines, or -1 if not found.
    """
    ass_pos = -1
    m2line = mig2line.split('\n')[0][-16:]
    m1line = mig1line.split('\n')[0][:16]
    
    for i in range(1, len(asslines)):
        assline = asslines[i].split('\t')[0]
        pe = findmatchedpos(m2line, assline[-16:], -1)
        ps = findmatchedpos(m1line, assline[:16], -1)
        if pe >= 0 and ps >= 0:
            ass_pos = i
            break
    
    return ass_pos

def get_umis_by_readIds_oldright(p_args, readIds_fpath, sample, subig):
    """
    Extracts UMI (Unique Molecular Identifier) sequences based on read IDs from older formats.

    Parameters:
    p_args (argparse.Namespace): Argument parser containing parameters and settings.
    readIds_fpath (str): Path to the read IDs file.
    sample (str): Sample name.
    subig (str): Sub IG type (e.g., IGG, IGA).

    Returns:
    list: List of UMIs.
    """
    readIdslines = []
    umis = []
    
    if not os.path.exists(readIds_fpath):
        return umis
    
    with open(readIds_fpath, 'r') as readIds_fr:
        for readIdsline in readIds_fr:
            segs = readIdsline.split('\t')
            readIdslines.append(segs[-1])
    
    prefd = os.path.join(p_args.out_dir, 'run_migec', 'assemble')
    prelab = '_R2.t*'
    
    if p_args.cogentIP_version == 'v2':
        prefd = os.path.join(p_args.out_dir, 'preprocess')
        prelab = '_R1.fastq'
    
    ig_fpath = find_file(f"{sample}_{subig}{prelab}", prefd)[0]
    
    with open(ig_fpath, 'r') as ig_fr:
        igfqlines = ig_fr.readlines()
    
    for i in range(1, len(readIdslines)):
        index_fq = int(readIdslines[i].rstrip())
        if p_args.cogentIP_version == 'v2':
            umi_fq = igfqlines[4*index_fq].split(':')[2].split('_umi_')[1]
        else:
            umi_fq = igfqlines[4*index_fq].split(':')[1]
        umis.append(umi_fq)
    
    return umis

def get_umis_by_readIds(p_args, readIds_fpath, sample, subig):
    """
    Extracts UMI (Unique Molecular Identifier) sequences based on read IDs from the latest formats.

    Parameters:
    p_args (argparse.Namespace): Argument parser containing parameters and settings.
    readIds_fpath (str): Path to the read IDs file.
    sample (str): Sample name.
    subig (str): Sub IG type (e.g., IGG, IGA).

    Returns:
    list: List of UMIs.
    """
    readIdslines = []
    umis = []
    
    if not os.path.exists(readIds_fpath):
        return umis
    
    with open(readIds_fpath, 'r') as readIds_fr:
        for readIdsline in readIds_fr:
            segs = readIdsline.split('\t')
            readIdslines.append(segs[-1])
    
    prefd = os.path.join(p_args.out_dir, 'run_migec', 'assemble')
    prelab = '_R2.t*'
    
    if p_args.cogentIP_version == 'v2':
        prefd = os.path.join(p_args.out_dir, 'preprocess')
        prelab = '_R1.fastq'
    
    ig_fpath = find_file(f"{sample}_{subig}{prelab}", prefd)[0]
    count_fqline = -1
    igfqline = ''
    
    with open(ig_fpath, 'r') as ig_fr:
        for i in range(1, len(readIdslines)):
            index_fq = int(readIdslines[i].rstrip())
            while count_fqline < 4*index_fq:
                igfqline = ig_fr.readline()
                count_fqline += 1
            if p_args.cogentIP_version == 'v2':
                umi_fq = igfqline.split(':')[2].split('_umi_')[1]
            else:
                umi_fq = igfqline.split(':')[1]
            umis.append(umi_fq)
    
    return umis

def get_umis_seqs_by_readIds(p_args, readIds_fpath, sample, subig):
    """
    Extracts UMI (Unique Molecular Identifier) sequences and assembled sequences based on read IDs.

    Parameters:
    - p_args: Arguments passed to the function, including output directory and version info.
    - readIds_fpath (str): File path to the read IDs file.
    - sample (str): Sample name.
    - subig (str): Sub-group identifier.

    Returns:
    - umis (list): List of extracted UMIs.
    - assemseqs (list): List of assembled sequences.
    """
    readIdslines, assemseqs, umis = [], [], []

    if not os.path.exists(readIds_fpath):
        return umis, assemseqs

    with open(readIds_fpath, 'r') as readIds_fr:
        for readIdsline in readIds_fr:
            segs = readIdsline.split('\t')
            if segs[0] != 'targetSequences':
                assemseqs.append(segs[0])
                readIdslines.append(segs[-1])

    prefd = os.path.join(p_args.out_dir, 'run_migec', 'assemble')
    prelab = '_R2.t*'
    if p_args.cogentIP_version == 'v2':
        prefd = os.path.join(p_args.out_dir, 'preprocess')
        prelab = '_R1.fastq'

    ig_fpath = find_file(sample + '_' + subig + prelab, prefd)[0]

    count_fqline, igfqline = -1, ''
    ig_fr = open(ig_fpath, 'r')

    # Loop through read ID lines and extract corresponding UMI from the fastq file
    for i in range(len(readIdslines)):
        index_fq = int(readIdslines[i].rstrip())
        while count_fqline < 4 * index_fq:
            igfqline = ig_fr.readline()
            count_fqline += 1
        umi_fq = igfqline.split(':')[2].split('_umi_')[1] if p_args.cogentIP_version == 'v2' else igfqline.split(':')[1]
        umis.append(umi_fq)

    ig_fr.close()
    return umis, assemseqs


def alignmentsPrettySeq(assemPretty_f, index, migIgAGr1_f, migIgAGr2_f, count):
    """
    Placeholder function to pretty print alignment sequences.

    Parameters:
    - assemPretty_f: Pretty alignment file.
    - index: Index of the sequence.
    - migIgAGr1_f: R1 file for MIGEC alignment.
    - migIgAGr2_f: R2 file for MIGEC alignment.
    - count: Sequence count.

    Returns:
    - seq (str): The formatted sequence.
    """
    seq = ''
    return seq


def replace_mig_subIG_R1_by_assembledseq(p_args, sample, subig, process_type, mixcr_dir):
    """
    Replaces sequences in MIGEC R1 file with assembled sequences from MiXCR analysis.

    Parameters:
    - p_args: Arguments passed to the function, including output directory and version info.
    - sample (str): Sample name.
    - subig (str): Sub-group identifier.
    - process_type (str): Type of process (e.g., alignments, readIds).
    - mixcr_dir (str): Directory for MiXCR output files.

    Returns:
    - None
    """
    mydir = os.path.join(p_args.out_dir, mixcr_dir, process_type)
    readIds_fpath = find_file(sample + '_' + subig + '_' + process_type + '_readIds.tsv', mydir)[0]
    assem_fpath = find_file(sample + '_' + subig + '_' + process_type + '_alignments.tsv', mydir)[0]
    assemPretty_fpath = find_file(sample + '_' + subig + '_' + process_type + '_alignmentsPretty.tsv', mydir)[0]
    umis_readids = get_umis_by_readIds(p_args, readIds_fpath, sample, subig)
    umis_unique = sorted(set(umis_readids))

    umis_fpath = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{sample}_{subig}_{process_type}_umis.txt")
    if not os.path.exists(umis_fpath):
        save_list_to_file(umis_fpath, umis_readids, True)

    with open(assem_fpath, 'r') as assem_fr:
        assemlines = assem_fr.readlines()

    mydir = os.path.join(p_args.out_dir, 'run_migec', 'assemble')
    migsubr1_fpath = find_file(sample + '_sub' + subig + '_R1*.fastq', mydir)[0]
    migsubr2_fpath = find_file(sample + '_sub' + subig + '_R2*.fastq', mydir)[0]

    with open(migsubr1_fpath, 'r') as migsub_fr1, open(migsubr2_fpath, 'r') as migsub_fr2:
        migsub1lines = migsub_fr1.readlines()
        migsub2lines = migsub_fr2.readlines()

    subr1, subr2 = [], []
    for count in range(0, len(migsub1lines), 4):
        umi = migsub1lines[count].split(":")[1]
        if umi in umis_readids:
            index = umis_readids.index(umi)
            seq = assemlines[index + 1].split('\t')[0]
            if ',' not in seq:
                subr1.append(migsub1lines[count])
                subr1.append(seq)
                subr1.append(migsub1lines[count + 2])
                subr1.append('I' * len(seq))
                subr2.append(migsub2lines[count])
                subr2.append(migsub2lines[count + 1])
                subr2.append(migsub2lines[count + 2])
                subr2.append(migsub2lines[count + 3])

    print('New MIGEC sub R1 is ready to save')
    save_list_to_file(migsubr1_fpath, subr1, True)
    save_list_to_file(migsubr2_fpath, subr2, True)

def replace_pre_subIG_R2_by_assembledseq(p_args, sample, subig, process_type, mixcr_dir):
    """
    Replaces sequences in subIG R1 and R2 files with assembled sequences from a MiXCR output.

    Parameters:
    - p_args (Namespace): Argument parser object containing paths and other settings.
    - sample (str): The sample name.
    - subig (str): Sub-IG type (e.g., 'IGG', 'IGA').
    - process_type (str): The type of processing (e.g., 'align', 'assemble').
    - mixcr_dir (str): Directory containing MiXCR output files.

    This function performs the following tasks:
    1. Locates necessary MiXCR output files, including `readIds.tsv` and `alignmentsPretty.tsv`.
    2. Extracts UMIs and sequences corresponding to read IDs.
    3. Reads from sub-IG R1 and R2 files, and matches UMIs with assembled sequences.
    4. Replaces the sequences in the sub-IG R1/R2 files with corresponding assembled sequences.
    5. Renames the original files after processing to avoid overwriting.
    """

    mydir = os.path.join(p_args.out_dir, mixcr_dir, process_type)
    readIds_fpath = find_file(sample + '_' + subig + '_' + process_type + '_readIds.tsv', mydir)[0]
    readIds_subfpath = readIds_fpath.replace(sample + '_' + subig, sample + '_sub' + subig)
    assemPretty_fpath = find_file(sample + '_' + subig + '_' + process_type + '_alignmentsPretty.tsv', mydir)[0]

    umis_igreadids, assemlines = get_umis_seqs_by_readIds(p_args, readIds_fpath, sample, subig)

    # Load sub-read IDs if available
    subreadids = []
    if os.path.exists(readIds_subfpath):
        with open(readIds_subfpath, 'r') as subreadIds_fr:
            subreadIdsline = subreadIds_fr.readline().rstrip()
            while subreadIdsline:
                segs = subreadIdsline.split('\t')
                subreadids.append(int(segs[-1]))
                subreadIdsline = subreadIds_fr.readline().rstrip()

    # Process R1 and R2 files
    mydir = os.path.join(p_args.out_dir, 'preprocess')
    migsubr1_fpath = find_file(sample + '_sub' + subig + '_R1*.fastq', mydir)[0]
    migsubr2_fpath = find_file(sample + '_sub' + subig + '_R2*.fastq', mydir)[0]
    migsubr1_f = migsubr1_fpath.replace('_R1.fastq', '_matched_R1.fastq')
    migsubr2_f = migsubr2_fpath.replace('_R2.fastq', '_matched_R2.fastq')

    if not os.path.exists(migsubr1_fpath) or os.stat(migsubr1_fpath).st_size == 0:
        print("The file does not exist or is empty; the replacement cannot proceed.")
        return

    subr1, subr2 = [], []
    count = 0

    with open(migsubr1_fpath, 'r') as migsub_fr1, open(migsubr2_fpath, 'r') as migsub_fr2, ExitStack() as cm:
        fhs1 = cm.enter_context(open(migsubr1_f, 'a'))
        fhs2 = cm.enter_context(open(migsubr2_f, 'a'))

        for x, y in zip(migsub_fr1, migsub_fr2):
            subr1.append(x.rstrip())
            subr2.append(y.rstrip())
            count += 1

            if len(subr1) == 4:
                umi_subig = subr1[0].split(':')[2].split('_umi_')[1]

                if (count // 4 - 1 in subreadids) and (umi_subig in umis_igreadids):
                    index = umis_igreadids.index(umi_subig)
                    seq = assemlines[index + 1].split('\t')[0]

                    if ',' not in seq:
                        subr1[1] = seq
                        subr1[3] = 'I' * len(seq)

                        for item in subr1:
                            fhs1.write(f'{item}\n')
                        for item in subr2:
                            fhs2.write(f'{item}\n')

                subr1.clear()
                subr2.clear()

    # Rename original files to indicate processing is complete
    print('New sub R1/2 files are ready to rename.')
    migsubr1_ftemp = migsubr1_fpath.replace('_R1.fastq', '_temp_R1.fastq')
    migsubr2_ftemp = migsubr1_fpath.replace('_R2.fastq', '_temp_R2.fastq')
    os.rename(migsubr1_fpath, migsubr1_ftemp)
    os.rename(migsubr2_fpath, migsubr2_ftemp)
    os.rename(migsubr1_f, migsubr1_fpath)
    os.rename(migsubr2_f, migsubr2_fpath)

def repone_subIG_Rn_tx_cf_by_alignassem(Rn, p_args, sample, subig, process_type, mixcr_dir):
    """
    Replaces sequences in sub-IG R1 and R2 files with sequences from MiXCR alignments.

    Parameters:
    - Rn (str): Placeholder, not used in the function.
    - p_args (Namespace): Argument parser object containing paths and other settings.
    - sample (str): The sample name.
    - subig (str): Sub-IG type (e.g., 'IGG', 'IGA').
    - process_type (str): The type of processing (e.g., 'align', 'assemble').
    - mixcr_dir (str): Directory containing MiXCR output files.

    This function performs the following tasks:
    1. Identifies the MiXCR alignment and read IDs files.
    2. Extracts UMIs (Unique Molecular Identifiers) based on read IDs.
    3. Saves the extracted UMIs to a file if it does not already exist.
    4. Reads the MiXCR alignments file to get the sequences.
    5. Matches the UMIs in the sub-IG R1 and R2 files with those from the MiXCR alignments.
    6. Replaces the sequences in the sub-IG R1/R2 files with the corresponding sequences from the MiXCR alignments.
    7. Saves the modified sub-IG R1 and R2 files.

    Notes:
    - The function processes the sub-IG files in blocks of four lines at a time (FASTQ format).
    - If the sequence does not contain a comma, it is considered valid for replacement.
    """

    # Initialize variables to store file paths
    assem_f = "" 
    readIds_f = ""
    migsubr1_f = ""
    migsubr2_f = ""

    # Identify the MiXCR files for alignments and read IDs
    for mixf in os.listdir(os.path.join(p_args.out_dir, mixcr_dir, process_type)):
        if mixf.endswith(sample + '_' + subig + '_' + process_type + '_alignments.txt'):
            assem_f = mixf
        elif mixf.endswith(sample + '_' + subig + '_' + process_type + '_readIds.txt'):
            readIds_f = mixf

    # Extract UMIs based on read IDs
    umis = get_umis_by_readIds(p_args, mixcr_dir, process_type, readIds_f, sample, subig)
    print("UMIs extraction completed")

    # Save UMIs to a file if not already present
    umi_f = sample + "_" + subig + "_" + process_type + '_umis.txt'
    if not os.path.exists(os.path.join(p_args.out_dir, mixcr_dir, process_type, umi_f)):
        save_list_to_file(os.path.join(p_args.out_dir, mixcr_dir, process_type, umi_f), umis, True)

    # Read the MiXCR alignments file
    with open(os.path.join(p_args.out_dir, mixcr_dir, process_type, assem_f), 'r') as assem_fr:
        asslines = assem_fr.readlines()

    # Identify the sub-IG R1 and R2 files in the MIGEC directory
    nonsubig = subig.replace("sub", "")
    for migsubfs in os.listdir(os.path.join(p_args.out_dir, 'run_migec', 'assemble')):
        if migsubfs.find(sample + "_sub" + nonsubig + "_R1") >= 0:
            migsubr1_f = migsubfs
        elif migsubfs.find(sample + "_sub" + nonsubig + "_R2") >= 0:
            migsubr2_f = migsubfs

    # Read the sub-IG R1 and R2 files
    with open(os.path.join(p_args.out_dir, 'run_migec', 'assemble', migsubr1_f), 'r') as migsub_fr1:
        migsub1lines = migsub_fr1.readlines()
    with open(os.path.join(p_args.out_dir, 'run_migec', 'assemble', migsubr2_f), 'r') as migsub_fr2:
        migsub2lines = migsub_fr2.readlines()

    # Process the sub-IG R1 and R2 files, replacing sequences with those from MiXCR alignments
    subr1 = []
    subr2 = []
    for count in range(0, len(migsub1lines), 4):
        umi = migsub1lines[count].split(":")[1]
        if umi in umis:
            index = umis.index(umi)
            print("Processing read at count:", count)
            print("Original sequence:", migsub1lines[count + 1][:30])
            print("Aligned sequence:", asslines[index + 1][:30])
            seq = asslines[index + 1].split('\t')[0]
            if seq.find(",") < 0:
                if subig.find("sub") < 0:
                    subr1.extend([migsub1lines[count], seq, migsub1lines[count + 2], 'I' * len(seq)])
                    subr2.extend(migsub2lines[count:count + 4])
                else:
                    subr2.extend([migsub2lines[count], seq, migsub2lines[count + 2], 'I' * len(seq)])
                    subr1.extend(migsub1lines[count:count + 4])

    print("Processing completed. Saving the modified files.")

    # Save the modified sub-IG R1 and R2 files
    new_migsubr1_f = migsubr1_f
    new_migsubr2_f = migsubr2_f
    save_list_to_file(os.path.join(p_args.out_dir, "run_migec", "assemble", new_migsubr1_f), subr1, True)
    save_list_to_file(os.path.join(p_args.out_dir, "run_migec", "assemble", new_migsubr2_f), subr2, True)

def replace_mig_ig2(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    """
    Replaces sequences in IGG or IGA files with sequences from sub-IGG or sub-IGA files.

    Parameters:
    - p_args (Namespace): Argument parser object containing paths and other settings.
    - input_dict (dict): Dictionary containing sample names and associated file paths.
    - read_type (str): Type of read (e.g., 'R1', 'R2').
    - target_region (str): Target region for processing.
    - mixcr_dir (str): Directory containing MiXCR output files.
    - sub_dir (str): Directory containing sub-IG files.
    - log_name (str): Log file name.
    - errfile_name (str): Error file name.

    This function identifies whether the input files are IGG or IGA and calls the function
    `repone_subIG_Rn_tx_cf_by_alignassem` to replace sequences in sub-IG files with those
    from the IG files.
    """
    process_type = f"{read_type}_{target_region}"
    subig = "IGG"
    for key, value in input_dict.items():
        if "_subIGG" in key:
            subig = "subIGG"
        elif "_subIGA" in key:
            subig = "subIGA"
        sample = key[:-7]
        repone_subIG_Rn_tx_cf_by_alignassem("subR2", p_args, sample, subig, process_type, mixcr_dir)

def replace_subig_by_assembledseq(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    """
    Replaces sequences in sub-IG files with assembled sequences.

    Parameters:
    - p_args (Namespace): Argument parser object containing paths and other settings.
    - input_dict (dict): Dictionary containing sample names and associated file paths.
    - read_type (str): Type of read (e.g., 'R1', 'R2').
    - target_region (str): Target region for processing.
    - mixcr_dir (str): Directory containing MiXCR output files.
    - sub_dir (str): Directory containing sub-IG files.
    - log_name (str): Log file name.
    - errfile_name (str): Error file name.

    This function identifies whether the input files are IGG or IGA and calls the appropriate
    function to replace sequences in sub-IG files with assembled sequences.
    """
    logger.info('Launching function replace_subig_by_assembledseq to execute assembled seq to replace _subIg?_R1/2 \n')
    process_type = f"{read_type}_{target_region}"
    subig = "IGG"
    for key, value in input_dict.items():
        if "_IGG" in key or "_IGA" in key:
            subig = key[-3:]
            sample = key[:-4]
            if p_args.cogentIP_version == 'cenegtv1.5':
                replace_mig_subIG_R1_by_assembledseq(p_args, sample, subig, process_type, mixcr_dir)
            else:
                replace_pre_subIG_R2_by_assembledseq(p_args, sample, subig, process_type, mixcr_dir)

def save_list_to_file(filename, mylist, overwrite):
    """
    Saves a list to a file.

    Parameters:
    - filename (str): The name of the file.
    - mylist (list): The list to be saved.
    - overwrite (bool): If True, overwrites the file if it exists.

    The function writes the contents of the list to the specified file.
    """
    if overwrite or not os.path.exists(filename):
        with open(filename, 'w') as fr:
            for item in mylist:
                item = item.split('\n')[0]
                fr.write("%s\n" % item)
                #fr.write(f"{item.split('\n')[0]}\n")
        if platform.system() == 'Windows':
            dos2unix(filename)

def getumidict(reads_fq):
    """
    Creates a dictionary of UMIs from a FASTQ file.

    Parameters:
    - reads_fq (list): List of lines from a FASTQ file.

    Returns:
    - dict: Dictionary with UMIs as keys and their corresponding indices as values.

    The function parses the FASTQ file to extract UMIs and stores their positions.
    """
    umidect = {}
    for i in range(0, len(reads_fq), 4):
        umi = reads_fq[i].split(":")[1]
        if umi in umidect:
            print("error")
        umidect[umi] = i
    return umidect

def reverse_complement(seq):
    """
    Computes the reverse complement of a DNA sequence.

    Parameters:
    - seq (str): The DNA sequence.

    Returns:
    - str: The reverse complement of the sequence.

    The function replaces each nucleotide with its complement and reverses the string.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def match_R2umi_non_sub(fnonr2, fsubr2, p_args, sub_dir, assemble_folder):
    """
    Matches UMIs between non-sub and sub-IG R2 files and modifies sequences.

    Parameters:
    - fnonr2 (str): File path of the non-sub R2 FASTQ file.
    - fsubr2 (str): File path of the sub-IG R2 FASTQ file.
    - p_args (Namespace): Argument parser object containing paths and other settings.
    - sub_dir (str): Directory containing sub-IG files.
    - assemble_folder (str): Directory containing assembled files.

    Returns:
    - tuple: Paths of the modified non-sub and sub-IG R2 FASTQ files.

    The function compares UMIs between the non-sub and sub-IG R2 files and replaces sequences
    accordingly, saving the results in new FASTQ files.
    """
    uminonr2, umisubr2 = [], []
    with open(fnonr2, 'r') as fnonr2_fr:
        fnonr2lines = fnonr2_fr.readlines()
    for i in range(0, len(fnonr2lines), 4):
        uminonr2.append(fnonr2lines[i].split(":")[1])

    with open(fsubr2, 'r') as fsubr2_fr:
        fsubr2lines = fsubr2_fr.readlines()
    for i in range(0, len(fsubr2lines), 4):
        umisubr2.append(fsubr2lines[i].split(":")[1])

    non_r2lines, sub_r2lines = [], []
    for i, umi in enumerate(umisubr2):
        if umi in uminonr2:
            index = uminonr2.index(umi)
            for j in range(4):
                line = fsubr2lines[4*i+j] if j != 0 else fsubr2lines[4*i+j].replace("R2", "R1")
                non_r2lines.append(fnonr2lines[4*index+j] if j != 0 else line)
                sub_r2lines.append(fsubr2lines[4*i+j])
        else:
            for j in range(4):
                line = fsubr2lines[4*i+j] if j != 0 else fsubr2lines[4*i+j].replace("R2", "R1")
                non_r2lines.append(line)
                sub_r2lines.append(fsubr2lines[4*i+j])

    non_fq2 = fnonr2.replace("_R2", "_match_R1")
    sub_fq2 = fsubr2
    save_list_to_file(non_fq2, non_r2lines, True)
    save_list_to_file(sub_fq2, sub_r2lines, True)
    return non_fq2, sub_fq2

def run_mixcr_sub2(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    """
    Runs MiXCR on sub-IG files and replaces sequences in non-sub IG files.

    Parameters:
    - p_args (Namespace): Argument parser object containing paths and other settings.
    - input_dict (dict): Dictionary containing sample names and associated file paths.
    - read_type (str): Type of read (e.g., 'R1', 'R2').
    - target_region (str): Target region for processing.
    - mixcr_dir (str): Directory containing MiXCR output files.
    - sub_dir (str): Directory containing sub-IG files.
    - log_name (str): Log file name.
    - errfile_name (str): Error file name.

    This function:
    1. Creates directories and prepares sub-IG samples.
    2. Matches UMIs between non-sub and sub-IG files.
    3. Runs MiXCR and replaces sequences in the non-sub IG files with those from sub-IG files.
    """
    out_log = os.path.join(p_args.out_dir, mixcr_dir, log_name)
    run_mixcr_err = os.path.join(p_args.out_dir, mixcr_dir, errfile_name)
    process_type = f"{read_type}_{target_region}"
    sub_dict = OrderedDict()
    for key, value in input_dict.items():
        if "_subIGG" in key or "_subIGA" in key:
            sub_dict[key] = value
    for key, value in input_dict.items():
        if "_IGG" in key or "_IGA" in key:
            fq1 = find_file(f"{key}_R2*.fastq", os.path.join(p_args.out_dir, sub_dir, 'assemble'))[0]
            fq2 = find_file(f"{key.replace('_', '_sub')}_R2*.fastq", os.path.join(p_args.out_dir, sub_dir, 'assemble'))[0]
            match_fq1, match_fq2 = match_R2umi_non_sub(fq1, fq2, p_args, sub_dir, 'assemble')
            subkey = f"{key[:-3]}sub{key[-3:]}"
            sub_dict[subkey][0] = os.path.basename(match_fq1)
            sub_dict[subkey][1] = os.path.basename(match_fq2)
    run_mixcr_base(2, p_args, sub_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)
    replace_mig_ig2(p_args, sub_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)

def matched_subIg_fsize(key, fd):
    """
    Retrieves the file size of matched sub-IG FASTQ files.

    Parameters:
    - key (str): The sample key.
    - fd (str): Directory path.

    Returns:
    - int: Size of the matched sub-IG FASTQ file.

    The function finds the matched sub-IG FASTQ file in the specified directory and returns its size.
    """
    segs = key.split('_')
    key = f"{segs[0]}_sub{segs[1]}"
    fnames = find_file(f"{key}*.fastq", fd)
    return os.stat(fnames[0]).st_size if fnames else 0

def sort_fastqs_in_list(sub_dict, p_args, sub_dir):
    """
    Sorts the records in the FASTQ files based on sequence IDs.

    Parameters:
    - sub_dict (dict): Dictionary containing sample names and associated FASTQ file paths.
    - p_args (Namespace): Argument parser object containing paths and other settings.
    - sub_dir (str): Directory containing the FASTQ files.

    The function sorts the sequences in the FASTQ files by their IDs and saves them to new files.
    """
    for key, value in sub_dict.items():
        forward_R1fq = os.path.join(p_args.out_dir, sub_dir, value[0])
        reverse_R2fq = os.path.join(p_args.out_dir, sub_dir, value[1])
        forward_R1out = os.path.join(p_args.out_dir, sub_dir, value[0].replace('R1.fastq', 'sorted_R1.fastq'))
        reverse_R2out = os.path.join(p_args.out_dir, sub_dir, value[1].replace('R2.fastq', 'sorted_R2.fastq'))
        
        with open(forward_R1fq, "r") as forward_file, open(reverse_R2fq, "r") as reverse_file, \
             open(forward_R1out, "w") as forward_out, open(reverse_R2out, "w") as reverse_out:
            forward_records = sorted(SeqIO.parse(forward_file, "fastq"), key=lambda x: x.id)
            reverse_records = sorted(SeqIO.parse(reverse_file, "fastq"), key=lambda x: x.id)
            SeqIO.write(forward_records, forward_out, "fastq")
            SeqIO.write(reverse_records, reverse_out, "fastq")
        
        if platform.system() == 'Windows':
            dos2unix(forward_R1out)
            dos2unix(reverse_R2out)

def sort_fastqs_in_list_new(sub_dict, p_args, sub_dir):
    """
    Sort FASTQ files based on read IDs in the forward and reverse reads.

    Parameters:
    sub_dict (dict): Dictionary containing sample identifiers and corresponding FASTQ file names.
    p_args (Namespace): Parsed arguments including the output directory.
    sub_dir (str): Subdirectory within the output directory to locate the FASTQ files.

    Returns:
    None: Writes the sorted FASTQ files to disk.
    """
    for key, value in sub_dict.items():
        # Construct file paths for the forward and reverse reads
        readsR1f = os.path.join(p_args.out_dir, sub_dir, value[0])
        readsR2f = os.path.join(p_args.out_dir, sub_dir, value[1])
        
        # Parse the FASTQ files
        forward_reads = SeqIO.parse(readsR1f, "fastq")
        reverse_reads = SeqIO.parse(readsR2f, "fastq")
                
        # Sort the reads by their ID
        sorted_forward_reads = sorted(forward_reads, key=lambda read: read.id)
        sorted_reverse_reads = sorted(reverse_reads, key=lambda read: read.id)
        
        # Write the sorted reads back to the original files
        with open(readsR1f, "w") as forward_outfile:
            SeqIO.write(sorted_forward_reads, forward_outfile, "fastq")
        with open(readsR2f, "w") as reverse_outfile:
            SeqIO.write(sorted_reverse_reads, reverse_outfile, "fastq")
        
        # Convert files to Unix format if running on Windows
        if platform.system() == 'Windows':
            dos2unix(readsR1f)
            dos2unix(readsR2f)

def run_mixcr_preIgGA(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    """
    Run MiXCR processing for IgG/IgA regions.

    Parameters:
    p_args (Namespace): Parsed arguments including the output directory and species.
    input_dict (dict): Dictionary containing sample identifiers and corresponding FASTQ file names.
    read_type (str): Type of reads (e.g., paired-end, single-end).
    target_region (str): Target region for MiXCR analysis (e.g., 'cdr3', 'fl').
    mixcr_dir (str): Directory for MiXCR output.
    sub_dir (str): Subdirectory within the output directory for intermediate files.
    log_name (str): Log file name.
    errfile_name (str): Error file name.

    Returns:
    None: Processes input files, sorts them, and runs MiXCR, then updates the input dictionary.
    """
    sub_dict = OrderedDict()

    # Filter input_dict for IgG or IgA samples
    for key, value in input_dict.items():
        if key.find('_IGG') >= 0 or key.find('_IGA') >= 0:
            sub_dict[key] = value.copy()

    # Additional filtering if using ImmuneProfiler
    if workflow2_link.find('ImmuneProfilerv2.jar') >= 0:
        for key, value in input_dict.items():
            if key.find('_subIGG') >= 0 or key.find('_subIGA') >= 0:
                sub_dict[key] = value.copy()

    # Sort FASTQ files
    sort_fastqs_in_list_new(sub_dict, p_args, sub_dir)

    # Run MiXCR for the selected samples
    run_mixcr(p_args, sub_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)

    # Replace subIg by assembled sequences
    replace_subig_by_assembledseq(p_args, sub_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)

    # Remove processed keys from input_dict
    keysdel = []
    for key, value in input_dict.items():
        if key.find('_IGG') >= 0 or key.find('_IGA') >= 0:
            mydir = os.path.join(p_args.out_dir, sub_dir, 'assemble')
            if workflow2_link.find('ImmuneProfilerv2.jar') >= 0:
                mydir = os.path.join(p_args.out_dir, sub_dir)
            if matched_subIg_fsize(key, mydir) == 0:
                segs = key.split('_')
                key = segs[0] + '_sub' + segs[1]
            keysdel.append(key)

    for key in keysdel:
        del input_dict[key]

def run_mixcr_base(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    """
    Run the MiXCR analysis pipeline for IgG/IgA or other receptor kits. This includes alignment, assembly, and exporting results.
    
    Parameters:
    p_args (Namespace): Parsed command-line arguments, including output directory, species, receptor kit, and threads.
    input_dict (dict): Dictionary with sample identifiers as keys and FASTQ file names as values.
    read_type (str): Type of reads (e.g., 'rawUnd').
    target_region (str): Target region for analysis ('cdr3' or 'fl').
    mixcr_dir (str): Directory for MiXCR outputs.
    sub_dir (str): Subdirectory within the output directory where intermediate files are stored.
    log_name (str): Name of the log file.
    errfile_name (str): Name of the error file.

    Returns:
    None: Executes MiXCR commands and writes outputs to the specified directories.
    """
    # Define paths for log files
    out_log = os.path.join(p_args.out_dir, mixcr_dir, log_name)
    run_mixcr_err = os.path.join(p_args.out_dir, mixcr_dir, errfile_name)
    
    # Determine the process type and preset
    process_type = f"{read_type}_{target_region}"
    tmp_dict = copy.deepcopy(input_dict)
    is_initUnd = False
    preset = 'takara'
    
    if p_args.species == 'human':
        preset += '-human'
    elif p_args.species == 'mouse':
        preset += '-mouse'
    
    if p_args.receptor_kit == 'TCRv1':
        preset += '-tcr-V1' if p_args.species == 'human' else '-tcr'
    elif p_args.receptor_kit == 'TCRv2':
        preset += '-tcr-V2' if p_args.species == 'human' else '-tcr-V2'
    elif 'BCR' in p_args.receptor_kit:
        preset += '-bcr'
        if p_args.species == 'mouse' and p_args.receptor_kit == 'BCRv2':
            preset += '-V2'
    
    if 'cdr3' in target_region:
        preset += '-cdr3'
    else:
        preset += '-full-length'
    
    # Define output file paths
    samps_clns = os.path.join(p_args.out_dir, mixcr_dir, process_type, '*.clns')
    alignQc_pdf = os.path.join(p_args.out_dir, 'report', 'alignQc.pdf')
    chainUsage_pdf = os.path.join(p_args.out_dir, 'report', 'chainUsage.pdf')
    createPDF = False
    
    # Process each sample
    for key, val in input_dict.items():
        reg_str = 'CDR3' if target_region == 'cdr3' else 'Full_length'
        print(f"{get_time()} [INFO] Processing {reg_str} region of sample {key}...", flush=True)
        
        # Define file paths for FASTQ files
        readslayout = 'Opposite'
        fq1, fq2 = '', ''
        
        if sub_dir == 'rawUnd' or read_type == 'rawUnd':
            fq1 = os.path.join(p_args.fastq_dir, val[0])
            fq2 = os.path.join(p_args.fastq_dir, val[1])
        elif sub_dir == 'run_migec':
            fq1 = find_file(f"{key}_R2*.fastq", os.path.join(p_args.out_dir, sub_dir, 'assemble'))[0]
            fq2 = find_file(f"{key}_R1*.fastq", os.path.join(p_args.out_dir, sub_dir, 'assemble'))[0]
            readslayout = 'Collinear'
        else:
            fq1 = find_file(f"{key}_R1*.fastq", os.path.join(p_args.out_dir, 'preprocess'))[0]
            fq2 = find_file(f"{key}_R2*.fastq", os.path.join(p_args.out_dir, 'preprocess'))[0]
        
        # Define output file paths for MiXCR results
        vdjca = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}.vdjca")
        align_report = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_align_report.txt")
        clns = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}.clns")
        clone_report = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_clones_report.txt")
        clone_all = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_clones_all.txt")
        vdjviz_input = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_vdjviz.txt")
        airr_tsv = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_airr.tsv")
        output_id = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}")
        na_R1 = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_notAligned_R1.fastq")
        na_R2 = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_notAligned_R2.fastq")
        exportaligns = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_alignments.tsv")
        exportPretty = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_alignmentsPretty.tsv")
        exportreadIds = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_readIds.tsv")
        exportclonesall = os.path.join(p_args.out_dir, mixcr_dir, process_type, f"{key}_{process_type}_clones_all.tsv")
        
        # Determine chains based on receptor kit
        chains = 'TRA,TRB' if p_args.receptor_kit in ['TCRv1', 'TCRv2'] else 'IG'
        
        # Set species-specific parameters
        if p_args.species == 'rhesus_monkey':
            mjson = os.path.join(REPO_PATH, 'src', 'config', 'imgt.201946-3.sv6.json')
            copyfile(mjson, os.path.join(p_args.main_dir, 'imgt.201946-3.sv6.json'))
            sp_cmd = ' --library imgt.201946-3.sv6 -s rhesus_monkey '
        elif p_args.species == 'mouse':
            sp_cmd = ' -s mouse '
        else:
            sp_cmd = ' -s human '
        
        # Set thread count if specified
        threads = f' --threads {p_args.process_threads}' if len(p_args.process_threads) > 0 else ''
        
        # Construct and execute MiXCR commands based on workflow version
        if 'ImmuneProfilerv1.6.jar' in workflow2_link:
            align_cmd = f"{workflow2_link} workflow2 align {threads} -p kaligner2 -OreadsLayout={readslayout} " \
                        f"-OvParameters.geneFeatureToAlign=Vtranscript {sp_cmd} {fq1} {fq2} {vdjca} --report {align_report}"
            if p_args.internal_test:
                align_cmd = f"{workflow2_link} workflow2 align -f {threads} -p kaligner2 -OreadsLayout={readslayout} " \
                            f"-OvParameters.geneFeatureToAlign=Vtranscript -OsaveOriginalReads=True " \
                            f"--not-aligned-R1 {na_R1} --not-aligned-R2 {na_R2} {sp_cmd} {fq1} {fq2} {vdjca} --report {align_report}"
                if p_args.allow_partial:
                    align_cmd = f"{workflow2_link} workflow2 align -f {threads} -p kaligner2 -OreadsLayout={readslayout} " \
                                f"-OvParameters.geneFeatureToAlign=Vtranscript -OsaveOriginalReads=True " \
                                f"-OallowPartialAlignments=True -OallowNoCDR3PartAlignments=True " \
                                f"--not-aligned-R1 {na_R1} --not-aligned-R2 {na_R2} {sp_cmd} {fq1} {fq2} {vdjca} --report {align_report}"
            logger.info(f'Launching MiXCR alignment using command:\n{align_cmd}')
            run_system_cmd(align_cmd, out_log, run_mixcr_err, 'alignment failed', log_name)
            
            assem_cmd = f"{workflow2_link} workflow2 assemble {vdjca} {clns} --report {clone_report}"
            if target_region == 'fl':
                assem_cmd = f"{workflow2_link} workflow2 assemble -OassemblingFeatures=VDJRegion {vdjca} {clns} --report {clone_report}"
            logger.info(f'Launching MiXCR assembling using command:\n{assem_cmd}')
            run_system_cmd(assem_cmd, out_log, run_mixcr_err, 'assembling failed', log_name)
            
            expt_aligns_cmd = f"{workflow2_link} workflow2 exportAlignments {vdjca} {exportaligns}"
            logger.info(f'Launching MiXCR exportAlignments using command:\n{expt_aligns_cmd}')
            run_system_cmd(expt_aligns_cmd, out_log, run_mixcr_err, 'exportAlignments failed', log_name)
            
            expt_cmd = f"{workflow2_link} workflow2 exportClones --preset full -c {chains} {clns} {clone_all}"
            logger.info(f'Launching MiXCR clone exporting using command:\n{expt_cmd}')
            run_system_cmd(expt_cmd, out_log, run_mixcr_err, 'clone exporting failed', log_name)
        
        else: # Assuming 'ImmuneProfilerv2.jar' is being used
            cmd = f"{workflow2_link} workflow2 analyze {preset} "
            if p_args.umi_cutoff:
                cmd += f" -Massemble.consensusAssemblerParameters.assembler.minRecordsPerConsensus={p_args.umi_cutoff} " \
                       f"-MrefineTagsAndSort.parameters.postFilter=null "
            cmd += f"{fq1} {fq2} {output_id} {threads}"
            logger.info(f'Launching MiXCR alignment using command:\n{cmd}')
            run_system_cmd(cmd, out_log, run_mixcr_err, 'alignment failed', log_name)
            
            airr_cmd = f"{workflow2_link} workflow2 exportAirr {clns} {airr_tsv}"
            logger.info(f'Launching MiXCR exportAirr using command:\n{airr_cmd}')
            run_system_cmd(airr_cmd, out_log, run_mixcr_err, 'exportAirr failed', log_name)
        
        if (p_args.receptor_kit == 'BCRv2hy' and ('IGG' in key or 'IGA' in key)) or p_args.internal_test:
            expt_alignsPretty_cmd = f"{workflow2_link} workflow2 exportAlignmentsPretty {vdjca} {exportPretty}"
            logger.info(f'Launching MiXCR exportAlignmentsPretty using command:\n{expt_alignsPretty_cmd}')
            run_system_cmd(expt_alignsPretty_cmd, out_log, run_mixcr_err, 'exportAlignmentsPretty failed', log_name)
            
            expt_readIds_cmd = f"{workflow2_link} workflow2 exportAlignments -readIds {vdjca} {exportreadIds}"
            logger.info(f'Launching MiXCR exportAlignmentsreadIds using command:\n{expt_readIds_cmd}')
            run_system_cmd(expt_readIds_cmd, out_log, run_mixcr_err, 'exportAlignmentsreadIds failed', log_name)
            
            if '_subIGG' in key or '_subIGA' in key:
                get_subIg(p_args, exportaligns, key, read_type, target_region, mixcr_dir, exportclonesall)
        
        if p_args.vdjviz_inputs:
            vdjviz_convert(clone_all, vdjviz_input)
    
    if 'ImmuneProfilerv2.jar' in workflow2_link and createPDF:
        alignQc_cmd = f"{workflow2_link} workflow2 exportQc align {samps_clns} {alignQc_pdf}"
        logger.info(f'Launching MiXCR exportQc using command:\n{alignQc_cmd}')
        run_system_cmd(alignQc_cmd, out_log, run_mixcr_err, 'exportQc align failed', log_name)
        
        chainUsage_cmd = f"{workflow2_link} workflow2 exportQc chainUsage {samps_clns} {chainUsage_pdf}"
        logger.info(f'Launching MiXCR exportQc using command:\n{chainUsage_cmd}')
        run_system_cmd(chainUsage_cmd, out_log, run_mixcr_err, 'exportQc chainUsage failed', log_name)
    
    if not p_args.keep_inter_file:
        vdjca_list = find_file('*.vdjca', os.path.join(p_args.out_dir, mixcr_dir, process_type))
        clns_list = find_file('*.clns', os.path.join(p_args.out_dir, mixcr_dir, process_type))
        del_list = vdjca_list + clns_list
        if del_list:
            for item in del_list:
                os.remove(item)

def run_mixcr(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    """
    Runs the MiXCR pipeline for a specific read type and target region.

    Parameters:
    - p_args: Argument parser object containing various user-defined parameters.
    - input_dict: Dictionary with input data for processing.
    - read_type: Type of sequencing reads (e.g., paired, single).
    - target_region: Target region for analysis (e.g., CDR3, full length).
    - mixcr_dir: Directory where MiXCR results will be stored.
    - sub_dir: Sub-directory for additional processing.
    - log_name: Name of the log file.
    - errfile_name: Name of the error file.

    Steps:
    1. Create a directory for the specific process type (read type and target region).
    2. Run the base MiXCR process using the provided parameters.
    3. Summarize and report the results.
    4. Handle post-processing steps such as merging CSV files and managing AIRR output.
    """

    process_type = f"{read_type}_{target_region}"
    tmp_dict = copy.deepcopy(input_dict)
    
    create_dir(os.path.join(p_args.out_dir, mixcr_dir, process_type))
    
    run_mixcr_base(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)
    
    cur_dict = summ_results(p_args.out_dir, process_type, tmp_dict, align_stats_dict, clone_stats_dict)
    report_stats(p_args.receptor_kit, os.path.join(p_args.out_dir, mixcr_dir, process_type), f"{p_args.out_name}_{process_type}", cur_dict)
    
    sample_list = list(set([k.split("_")[0] for k in fd.keys()]))
    
    for sample_id in sample_list:
        stats_file = merge_csv(sample_id, process_type, f"{p_args.out_name}_{sample_id}_{process_type}",
                               os.path.join(p_args.out_dir, mixcr_dir, process_type),
                               os.path.join(p_args.out_dir, 'report'))
        
        if p_args.receptor_kit == 'BCRv2hy':
            subIgfs = find_file(f"{sample_id}_subIG*_subIg?.tsv", os.path.join(p_args.out_dir, mixcr_dir, process_type))
            for subIgf in subIgfs:
                os.rename(subIgf, os.path.join(p_args.out_dir, 'report', os.path.basename(subIgf)))
        
        if p_args.airr_create:
            merge_airr_csv(sample_id, process_type, f"{p_args.out_name}_{sample_id}_{process_type}",
                           os.path.join(p_args.out_dir, mixcr_dir, process_type),
                           os.path.join(p_args.out_dir, 'report'),
                           os.path.join(p_args.out_dir, 'airr_report'))
    
    os.rename(stats_file, os.path.join(p_args.out_dir, 'report', os.path.basename(stats_file)))
    
    reg_str = 'CDR3' if target_region != 'fl' else 'Full_length'
    print(f"{get_time()} [INFO] Finished {reg_str} region processing for all samples", flush=True)
    
    if p_args.airr_create:
        stats_csv_to_airr(os.path.join(p_args.out_dir, 'report', os.path.basename(stats_file)),
                          os.path.join(p_args.out_dir, 'airr_report', os.path.basename(stats_file)))

def run_raw_mixcr(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    """
    Runs the MiXCR pipeline for raw sequencing data.

    Parameters: (similar to `run_mixcr`)
    - p_args: Argument parser object containing various user-defined parameters.
    - input_dict: Dictionary with input data for processing.
    - read_type: Type of sequencing reads (e.g., paired, single).
    - target_region: Target region for analysis (e.g., CDR3, full length).
    - mixcr_dir: Directory where MiXCR results will be stored.
    - sub_dir: Sub-directory for additional processing.
    - log_name: Name of the log file.
    - errfile_name: Name of the error file.

    The function is similar to `run_mixcr` but is tailored for raw sequencing data, with some differences in file handling.
    """

    process_type = f"{read_type}_{target_region}"
    tmp_dict = copy.deepcopy(input_dict)
    
    create_dir(os.path.join(p_args.out_dir, mixcr_dir, process_type))
    
    run_mixcr_base(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)
    
    cur_dict = summ_results(p_args.out_dir, process_type, tmp_dict, align_stats_dict, clone_stats_dict)
    report_stats(p_args.receptor_kit, os.path.join(p_args.out_dir, mixcr_dir, process_type), f"{p_args.out_name}_{process_type}", cur_dict)
    
    sample_list = list(set([k.split("_")[0] for k in fd.keys()]))
    
    for sample_id in sample_list:
        stats_file = merge_csv(sample_id, process_type, f"{p_args.out_name}_{sample_id}_{process_type}",
                               os.path.join(p_args.out_dir, mixcr_dir, process_type),
                               os.path.join(p_args.out_dir, 'report'))
        
        if workflow2_link.find('ImmuneProfilerv2.jar') >= 0 and p_args.airr_create:
            airrfiles = find_file(f"{sample_id}_*_airr.tsv", os.path.join(p_args.out_dir, mixcr_dir, process_type))
            for airr_orig in airrfiles:
                airr_dist = os.path.join(p_args.out_dir, "report", sample_id, os.path.basename(airr_orig))
                os.rename(airr_orig, airr_dist)
    
    os.rename(stats_file, os.path.join(p_args.out_dir, 'report', os.path.basename(stats_file)))
    
    reg_str = 'CDR3' if target_region != 'fl' else 'Full_length'
    print(f"{get_time()} [INFO] Finished raw {reg_str} region processing for all samples", flush=True)

def change_run_mixcr_name(p_args, mytype):
    """
    Changes the name of the MiXCR processing directory and associated report directories.

    Parameters:
    - p_args: Argument parser object containing various user-defined parameters.
    - mytype: A string representing the type of processing being done (e.g., 'raw', 'filtered').

    The function handles directory renaming and cleanup after processing. It also manages the
    copying of specific files to the appropriate directories.
    """

    change_dir(os.path.join(p_args.out_dir, 'run_mixcr'), os.path.join(p_args.out_dir, f"{mytype}_mixcr"))
    create_dir(os.path.join(p_args.out_dir, 'run_mixcr'))
    
    change_dir(os.path.join(p_args.out_dir, 'report'), os.path.join(p_args.out_dir, f"{mytype}_report"))
    create_dir(os.path.join(p_args.out_dir, 'report'))
    
    src = os.path.join(p_args.out_dir, f"{mytype}_report", f"{p_args.out_name}_sample_QC_stats.csv")
    dst = os.path.join(p_args.out_dir, 'report', f"{p_args.out_name}_sample_QC_stats.csv")
    shutil.copy(src, dst)
    
    if p_args.airr_create and workflow2_link.find('ImmuneProfiler.jar') >= 0:
        change_dir(os.path.join(p_args.out_dir, 'airr_report'), os.path.join(p_args.out_dir, f"{mytype}_airr_report"))
        create_dir(os.path.join(p_args.out_dir, 'airr_report'))
        
        if 'raw' in mytype:
            src = os.path.join(p_args.out_dir, f"{mytype}_airr_report", f"{p_args.out_name}_sample_QC_stats.csv")
            dst = os.path.join(p_args.out_dir, 'airr_report', f"{p_args.out_name}_sample_QC_stats.csv")
            shutil.copy(src, dst)
    
    if p_args.receptor_kit == 'BCRv2hy' and not p_args.internal_test and not debug:
        shutil.rmtree(os.path.join(p_args.out_dir, f"{mytype}_mixcr"))
        shutil.rmtree(os.path.join(p_args.out_dir, f"{mytype}_report"))
        
        if p_args.airr_create:
            shutil.rmtree(os.path.join(p_args.out_dir, f"{mytype}_airr_report"))

def parse_clones_report(in_file, out_file, key, inputDict):
    """
    Parses a clones report file and processes the data to generate summary statistics.

    Args:
        in_file (str): Path to the input file containing the clone report data.
        out_file (str): Path to the output file where the processed results will be written.
        key (str): A key to access specific entries in the inputDict for storing processed data.
        inputDict (dict): A dictionary to store summary statistics for different chain types.

    Processing:
        - Reads the clone report file and skips the header.
        - Extracts relevant fields from each line of the report, including information about 
          clone count, sequence quality, and V/D/J/C segments.
        - Determines the chain type (IG or TR) and specific chain (e.g., G, M, K, L, A, D, E, A, B).
        - Aggregates clone counts based on chain types and stores the results in the inputDict.
        - Writes the processed data to the output CSV file.
        - Handles special cases for data formatting based on the species and receptor kit used.

    Special Handling:
        - The function differentiates between IG and TR chains and further categorizes IG chains 
          by their constant region (e.g., IgG, IgM).
        - Frame shifts and stop codons in the amino acid sequences are detected and noted.
        - The function also merges existing statistics with new data if they are already present in inputDict.

    Example:
        Given a file with clone report data, the function processes it and generates summary 
        statistics for different immune chains (e.g., IgG, TRA) which are then stored in inputDict.

    Notes:
        - The function is designed to work with output from different versions of the ImmuneProfiler tool.
        - It supports processing for both mouse and human species depending on the receptor kit used.

    Returns:
        None
    """

    with open(in_file, 'r') as infile, open(out_file, 'a') as outfile:
        next(infile)  # Skip the header line
        header = [
            'Read Count', 'Fraction', 'Clonal Sequence', 'Clonal Sequence Quality', 
            'CDR3 Min Quality', 'CDR3 Sequence', 'CDR3 Amino Acid Sequence', 
            'Clonal Type', 'Frame Shift', 'Stop Codon', 'Amino Acid Length', 
            'V segment', 'all V hits', 'D segment', 'all D hits', 
            'J segment', 'all J hits', 'C segment', 'all C hits'
        ]
        outfile.write(','.join(header) + '\n')

        # Initialize counters for various chain types
        chains, chain_types, clones = set(), set(), set()
        countNoC, countAll, countHnoC = 0, 0, 0
        countG, countM, countK, countL, countA, countD, countE = 0, 0, 0, 0, 0, 0, 0
        countTRA, countTRB = 0, 0
        countClnNoC, countClnAll, countClnHnoC = 0, 0, 0
        countClnG, countClnM, countClnK, countClnL, countClnA, countClnD, countClnE = 0, 0, 0, 0, 0, 0, 0
        countClnTRA, countClnTRB = 0, 0

        shift = 0
        if workflow2_link.find('ImmuneProfilerv2.jar') >= 0:
            shift = 2
            if p_args.species == 'mouse' and p_args.receptor_kit in ['BCRv1', 'TCRv1']:
                shift = 0

        for line in infile:
            line = line.rstrip().split('\t')
            cloneId, cloneCount, cloneFraction = line[0], line[1], line[2]
            clonalSequence = line[3 + shift]
            clonalSequenceQuality = line[4 + shift]

            nSeqCDR3 = line[13 + shift]
            minQualCDR3 = line[14 + shift]
            aaSeqCDR3 = line[15 + shift]
            if workflow2_link.find('ImmuneProfiler.jar') >= 0:
                nSeqCDR3 = line[23]
                minQualCDR3 = line[24]
                aaSeqCDR3 = line[32]

            aaframeshift = 'Frame Shift' if '_' in aaSeqCDR3 else ''
            aaStopCodon = 'Stop Codon' if '*' in aaSeqCDR3 else ''
            aaLength = str(len(aaSeqCDR3))

            vlist = line[5 + shift].split(',')
            Vsegments = [x.split('*')[0] for x in vlist][0]
            vmulti = ';'.join([x.split('*')[0] for x in vlist])

            dlist = line[6 + shift].split(',')
            Dsegments = [x.split('*')[0] for x in dlist][0]
            dmulti = ';'.join([x.split('*')[0] for x in dlist])

            jlist = line[7 + shift].split(',')
            Jsegments = [x.split('*')[0] for x in jlist][0]
            jmulti = ';'.join([x.split('*')[0] for x in jlist])

            c_list = line[8 + shift].split(',')
            Csegments = [x.split('*')[0] for x in c_list][0]
            cmulti = ';'.join([x.split('*')[0] for x in c_list])

            # Determine chain type and specific chain
            chain_type = Vsegments.split(',')[0][:2]
            if chain_type == 'TR':
                chain = Vsegments.split(',')[0][2]
                clonalType = Vsegments.split(',')[0][:3]
                if Csegments.split(',')[0][:2] == '':
                    countNoC += int(float(cloneCount))
                    countClnNoC += 1
            elif chain_type == 'IG':
                if Csegments.split(',')[0][:2] != '':
                    if Csegments.split(',')[0][2:3] == 'H':
                        chain = Csegments.split(',')[0][3]  # G/M
                        clonalType = Csegments.split(',')[0][:2] + Csegments.split(',')[0][3]
                    else:
                        chain = Csegments.split(',')[0][2]  # K/L
                        clonalType = Csegments.split(',')[0][:3]
                else:
                    chain = Vsegments.split(',')[0][2]  # K/L/H
                    clonalType = Vsegments.split(',')[0][:3]
                    countNoC += int(float(cloneCount))
                    countClnNoC += 1

            chain_types.add(chain_type)
            chains.add(chain)
            countAll += int(float(cloneCount))
            countClnAll += 1

            # Aggregate counts by chain type
            if 'IG' in chain_types:
                if chain == 'G':
                    countG += int(float(cloneCount))
                    countClnG += 1
                elif chain == 'M':
                    countM += int(float(cloneCount))
                    countClnM += 1
                elif chain == 'H':
                    countHnoC += int(float(cloneCount))
                    countClnHnoC += 1
                elif chain == 'K':
                    countK += int(float(cloneCount))
                    countClnK += 1
                elif chain == 'L':
                    countL += int(float(cloneCount))
                    countClnL += 1
                elif chain == 'A':
                    countA += int(float(cloneCount))
                    countClnA += 1
                elif chain == 'D':
                    countD += int(float(cloneCount))
                    countClnD += 1
                elif chain == 'E':
                    countE += int(float(cloneCount))
                    countClnE += 1

            if 'TR' in chain_types:
                if chain == 'A':
                    countTRA += int(float(cloneCount))
                    countClnTRA += 1
                elif chain == 'B':
                    countTRB += int(float(cloneCount))
                    countClnTRB += 1

            # Write the processed data to the output file
            result = [
                cloneCount, cloneFraction, clonalSequence, clonalSequenceQuality, 
                minQualCDR3, nSeqCDR3, aaSeqCDR3, clonalType, 
                aaframeshift, aaStopCodon, aaLength, 
                Vsegments, vmulti, Dsegments, dmulti, 
                Jsegments, jmulti, Csegments, cmulti
            ]
            outfile.write(','.join(result) + '\n')

        # Update inputDict with summary statistics
        length = len(inputDict[key])
        if 'IG' in chain_types:
            if inputDict[key][length-8] == None or inputDict[key][length-8].find('IgG:') < 0:
                inputDict[key].extend([
                    'IgG:' + str(countClnG), 'IgM:' + str(countClnM), 
                    'IgK:' + str(countClnK), 'IgL:' + str(countClnL), 
                    'IgA:' + str(countClnA), 'IgD:' + str(countClnD), 
                    'IgE:' + str(countClnE), 'IgH(lack constant region):' + str(countClnHnoC)
                ])
            else:
                merge_chain_value(inputDict, key, 'IgG:', countClnG, length-8)
                merge_chain_value(inputDict, key, 'IgM:', countClnM, length-7)
                merge_chain_value(inputDict, key, 'IgK:', countClnK, length-6)
                merge_chain_value(inputDict, key, 'IgL:', countClnL, length-5)
                merge_chain_value(inputDict, key, 'IgA:', countClnA, length-4)
                merge_chain_value(inputDict, key, 'IgD:', countClnD, length-3)
                merge_chain_value(inputDict, key, 'IgE:', countClnE, length-2)
                merge_chain_value(inputDict, key, 'IgH(lack constant region):', countClnHnoC, length-1)
        
        if 'TR' in chain_types:
            if inputDict[key][length-2].find('TRA:') < 0:
                inputDict[key].extend([
                    'TRA:' + str(countClnTRA), 'TRB:' + str(countClnTRB)
                ])
            else:
                merge_chain_value(inputDict, key, 'TRA:', countClnTRA, length-2)
                merge_chain_value(inputDict, key, 'TRB:', countClnTRB, length-1)

def merge_chain_value(inputDict, key, chain, value, index):
    """
    Merges a new value into an existing entry in inputDict by adding the new value to the existing one.

    Args:
        inputDict (dict): A dictionary containing chain information.
        key (str): The key used to access the specific dictionary entry.
        chain (str): The specific chain type (e.g., 'IgG:').
        value (int): The value to be added to the existing entry.
        index (int): The index within the entry where the chain value is stored.

    Returns:
        None
    """
    old_value = int(inputDict[key][index].split(chain)[1])
    inputDict[key][index] = chain + str(old_value + value)
    return

def vdjviz_convert(file_clones_all, chord_file):
    """
    Converts a MiXCR clones_all.txt file into a VDJviz input file format.

    Args:
        file_clones_all (str): Path to the MiXCR clones_all.txt file.
        chord_file (str): Path to the output file to be used by VDJviz.

    Processing:
        - Reads the input file and extracts relevant clone data.
        - Modifies and writes the data into a format compatible with VDJviz.
        - If the input file is empty or contains fewer than three lines, the output file is not generated.

    Returns:
        None
    """
    if os.path.getsize(file_clones_all) > 0:
        with open(chord_file, 'w') as out_file:
            headers = ['count', 'frequency', 'CDR3nt', 'CDR3aa', 'V', 'D', 'J']
            out_file.write('\t'.join(headers) + '\n')

        with open(file_clones_all) as all_clones:
            all_clones_line_num = 0
            for line in all_clones:
                all_clones_line_num += 1
                if all_clones_line_num != 1:
                    all_clones_line_content = line.rstrip().split("\t")
                    all_clones_line_content[5] = all_clones_line_content[5].split("*")[0]
                    all_clones_line_content[6] = all_clones_line_content[6].split("*")[0]
                    all_clones_line_content[7] = all_clones_line_content[7].split("*")[0]
                    for idx in [5, 6, 7]:
                        if len(all_clones_line_content[idx]) == 0:
                            all_clones_line_content[idx] = str(r".")
                    newline = [
                        str(all_clones_line_content[1]),  # count
                        str(all_clones_line_content[2]),  # frequency
                        str(all_clones_line_content[23]),  # CDR3nt
                        str(all_clones_line_content[32]),  # CDR3aa
                        str(all_clones_line_content[5]),  # V
                        str(all_clones_line_content[6]),  # D
                        str(all_clones_line_content[7])  # J
                    ]
                    with open(chord_file, "a") as out_file:
                        out_file.write("\t".join(newline) + "\n")

            if all_clones_line_num < 3:
                os.remove(chord_file)
    return

def get_mergedfqs(p_args, log_name, errfile_name):
    """
    Logs the completion of merging fastq files for a hybrid kit.

    Args:
        p_args (Namespace): Command-line arguments or parameters.
        log_name (str): Name of the log file.
        errfile_name (str): Name of the error file.

    Returns:
        None
    """
    out_log = os.path.join(p_args.out_dir, log_name)
    merge_err = os.path.join(p_args.out_dir, errfile_name)
    print("merged fastq for hybrid kit done")

def split_align_unalign_by_readIds(fq_fn, id_fn):
    """
    Splits a FASTQ file into aligned and unaligned reads based on provided read IDs.

    Args:
        fq_fn (str): Path to the input FASTQ file.
        id_fn (str): Path to the file containing read IDs to align.

    Processing:
        - Reads the input FASTQ file and read ID file.
        - Segregates reads into aligned and unaligned based on the provided IDs.
        - Writes aligned and unaligned reads into separate FASTQ files.

    Returns:
        None
    """
    readsaligned = []
    readsunaligned = []
    with open(fq_fn) as fq_f:
        fqlines = fq_f.readlines()
    with open(id_fn) as ids_f:
        idlines = ids_f.readlines()

    count = 0
    for i in range(1, len(idlines)):
        readId = int(idlines[i])
        for j in range(4):
            readsaligned.append(fqlines[readId*4+j])
        for diff in range(count, readId):
            for j in range(4):
                readsunaligned.append(fqlines[diff*4+j])
        count = readId + 1

    fq_fn_aligned = fq_fn[:-6] + "_aligned.fastq"
    fq_fn_unaligned = fq_fn[:-6] + "_unaligned.fastq"
    save_list_to_file(fq_fn_aligned, readsaligned, True)
    save_list_to_file(fq_fn_unaligned, readsunaligned, True)

def umi_polyGNs(umi):
    """
    Checks if a UMI (Unique Molecular Identifier) contains a poly-G or poly-N sequence.

    Args:
        umi (str): The UMI sequence.

    Returns:
        bool: True if the UMI contains poly-G or poly-N, False otherwise.
    """
    ispolyGs = False
    i = 0
    while i < len(umi) and not ispolyGs:
        if umi[i:i+5] == 'GGGGG' or umi[i:i+5] == 'NNNNN':
            ispolyGs = True
        i += 1
    return ispolyGs

def get_adapter_args(p_args):
    """
    Determines the start position and length of the adapter based on the library kit and species.

    Args:
        p_args (Namespace): Command-line arguments or parameters containing library kit and species.

    Returns:
        tuple: (adapter_start, adapter_len) - The start position and length of the adapter.
    """
    adapter_start = 0
    adapter_len = 0
    if p_args.library_kit == 'takara_smartseq':
        if p_args.species == 'human':
            if p_args.receptor_kit == 'TCRv1':
                adapter_start = 8
                adapter_len = 18  # +1
            elif p_args.receptor_kit == 'TCRv2':
                adapter_start = 5
                adapter_len = 21
            elif p_args.receptor_kit == 'BCRv1':
                adapter_start = 0
                adapter_len = 0
            elif p_args.receptor_kit == 'BCRv2':
                adapter_start = 7
                adapter_len = 19  # +5
        elif p_args.species == 'mouse':
            if p_args.receptor_kit == 'TCRv1':
                adapter_start = 7
                adapter_len = 20
            elif p_args.receptor_kit == 'TCRv2':
                adapter_start = 7
                adapter_len = 20
            elif p_args.receptor_kit == 'BCRv1':
                adapter_start = 0
                adapter_len = 0
            elif p_args.receptor_kit == 'BCRv2':
                adapter_start = 7
                adapter_len = 20  # +2
    return adapter_start, adapter_len

def remove_umi_rawUnd(p_args, valfqs):
    """
    Removes UMI (Unique Molecular Identifier) and optional adapter sequences from FASTQ files.

    Args:
    p_args (argparse.Namespace): Contains various arguments such as fastq_dir, library_kit, rm_adapter, etc.
    valfqs (list): List containing the names of the input FASTQ files (typically paired-end).

    Returns:
    tuple: 
        - umi_dict (OrderedDict): Dictionary containing UMIs as keys and their counts as values.
        - total (int): Total number of reads processed.
        - umiread_count (int): Total number of reads containing valid UMIs.
    """
    raw_noumi_fq = os.path.join(p_args.fastq_dir, valfqs[0])  # Path to R1 FASTQ file (without UMI)
    linker_range = 0
    umi_dict = OrderedDict()
    umi_start = 0
    len_umi = 12
    len_umi_link = len_umi + 7
    raw_umi_fq = os.path.join(p_args.fastq_dir, valfqs[1])  # Path to R2 FASTQ file (with UMI)
    
    # Adjust UMI length if using 'cellecta_air' kit
    if p_args.library_kit == 'cellecta_air':
        len_umi = 14
        len_umi_link = 14
    
    r1fq = []
    r2fq = []
    umiread_count = 0
    total = 0
    
    # Remove adapter from R1 if specified
    adapter_start = 0
    adapter_len = 0
    if p_args.rm_adapter:
        adapter_start, adapter_len = get_adapter_args(p_args)
    
    with load_fastq(raw_noumi_fq) as f1, load_fastq(raw_umi_fq) as f2:
        prefix1 = os.path.basename(raw_noumi_fq).split('.fastq')[0]
        out_noumi_fq = os.path.join(p_args.fastqRMumi_dir, prefix1 + '.fastq')
        prefix2 = os.path.basename(raw_umi_fq).split('.fastq')[0]
        out_rmumi_fq = os.path.join(p_args.fastqRMumi_dir, prefix2 + '.fastq')
        
        fw1 = open(out_noumi_fq, 'w')
        fw2 = open(out_rmumi_fq, 'w')
        
        for x, y in zip(f1, f2):
            r1fq.append(x.rstrip())
            r2fq.append(y.rstrip())
            
            if len(r1fq) == 4:  # Process one complete FASTQ entry (4 lines)
                total += 1
                umi = get_umi(r2fq, umi_start, len_umi)
                
                if not umi_polyGNs(umi):
                    # Validate UMI for 'takara_smartseq' kit or other kits
                    if (p_args.library_kit == 'takara_smartseq' and r2fq[1][len_umi:len_umi_link] == 'GTACGGG') or p_args.library_kit != 'takara_smartseq':
                        umiread_count += 1
                        if umi in umi_dict.keys():
                            umi_dict[umi] += 1
                        else:
                            umi_dict[umi] = 1
                
                # Trim UMI and linker from R2
                r2fq[1] = r2fq[1][umi_start+len_umi_link:]
                r2fq[3] = r2fq[3][umi_start+len_umi_link:]
                
                # Trim adapter from R1 if applicable
                r1fq[1] = r1fq[1][adapter_start + adapter_len:]
                r1fq[3] = r1fq[3][adapter_start + adapter_len:]
                
                # Write the processed reads to output files
                for item in r1fq:
                    fw1.write('%s\n' % item)
                for item in r2fq:
                    fw2.write('%s\n' % item)
                
                # Clear the read buffers for the next entry
                del r1fq[:]
                del r2fq[:]
        
        fw1.close()
        fw2.close()
    
    # Convert output files to Unix format if running on Windows
    if platform.system() == 'Windows':
        dos2unix(out_noumi_fq)
        dos2unix(out_rmumi_fq)
    
    print("Total reads=", total, "\nUMI reads=", umiread_count)
    return umi_dict, total, umiread_count

def insertmax_returnmin(list_key5max, list_val5max, key_seq, val):
    """
    Inserts a new UMI and its count into the top 5 maximum UMI list.

    Args:
    list_key5max (list): List of top 5 UMI sequences.
    list_val5max (list): List of corresponding counts for the top 5 UMIs.
    key_seq (str): UMI sequence to be inserted.
    val (int): Count of the UMI sequence.

    Returns:
    tuple:
        - (str): The UMI sequence with the minimum count in the top 5 list.
        - (int): The minimum count in the top 5 list.
    """
    len_5top = 5
    index = 0
    
    # Find the correct position for the new UMI based on its count
    while index < len_5top and val > list_val5max[index]:
        index += 1
    
    if index > 0:
        index -= 1
        for i in range(index):
            if i < len_5top:
                list_key5max[i] = list_key5max[i + 1]
                list_val5max[i] = list_val5max[i + 1]
        
        list_key5max[index] = key_seq
        list_val5max[index] = val
    
    return list_key5max[0], list_val5max[0]

def visualize_umi_distribution(p_args, target_range, umi_dict, total, umireads_count):
    """
    Visualizes the distribution of UMI counts and writes summary statistics to a file.

    Args:
    p_args (argparse.Namespace): Contains various arguments such as fastqRMumi_dir, etc.
    target_range (str): Name or range identifier used in file naming.
    umi_dict (OrderedDict): Dictionary containing UMIs as keys and their counts as values.
    total (int): Total number of reads processed.
    umireads_count (int): Total number of reads containing valid UMIs.

    Returns:
    None
    """
    numkeys = 5
    list_key5max = ['NNNNNNNNNNNN'] * numkeys
    list_val5max = [0] * numkeys
    
    # Find the top 5 UMIs by count
    keymax = list(umi_dict.keys())[0]
    valmax = umi_dict[keymax]
    keymax, valmax = insertmax_returnmin(list_key5max, list_val5max, keymax, valmax)
    
    for key_umi, val_umicount in umi_dict.items():
        if val_umicount > valmax:
            valmax = val_umicount
            keymax = key_umi
            keymax, valmax = insertmax_returnmin(list_key5max, list_val5max, keymax, valmax)
    
    umi_summary_f = os.path.join(p_args.fastqRMumi_dir, target_range + '_umi_summary.csv')
    with open(umi_summary_f, 'w') as fw:
        print('Total reads: %s' % str(total))
        fw.write('Total reads:' + str(total) + '\n')
        print('UMI reads: %s' % str(umireads_count))
        fw.write('UMI reads:' + str(umireads_count) + '\n')
        print('Number of UMIs: %s' % len(umi_dict.keys()))
        fw.write('Number of UMIs:' + str(len(umi_dict.keys())) + '\n')
        print('Top 5 UMIs=', list_key5max)
        fw.write('Top 5 UMIs=')
        fw.write(','.join(str(item) for item in list_key5max))
        fw.write('\n')
        print('Top 5 UMIs with count=', list_val5max)
        fw.write('Top 5 UMIs with count=')
        fw.write(','.join(str(item) for item in list_val5max))
        fw.write('\n')
    
    # Create distribution of UMI counts
    valmax = list_val5max[-1]
    umi_counts = [0] * (valmax + 1)
    
    for key_umi, val_umicount in umi_dict.items():
        umi_counts[val_umicount] += 1
    
    x = [*range(0, valmax + 1, 1)]
    y = umi_counts
    
    umi_range_f = os.path.join(p_args.fastqRMumi_dir, target_range + '_umi_distribution.csv')
    with open(umi_range_f, 'w') as fp:
        fp.write('umi_type_numbers,')
        fp.write(','.join(str(item) for item in x))
        fp.write('\n')
        fp.write('frequency,')
        fp.write(','.join(str(item) for item in y))
    
    drawplot(x, y, umi_range_f[:-4] + '.png')
    return

def drawplot(x, y, fname):
    """
    Plots the UMI type distribution and saves the plots in various formats.

    Args:
    x (list): List of UMI type numbers.
    y (list): List of frequencies corresponding to UMI types.
    fname (str): The filename to save the plot, without extension.

    Returns:
    None
    """
    len_exp = len(y)
    
    # Adjust length if the last frequency is 1
    if y[-1] == 1:
        len_exp -= 1
    
    # Trim the arrays to exclude trailing zeros
    while len_exp > 50 and (y[len_exp-1] == 0 or y[len_exp-2] == 0 or y[len_exp-3] == 0 or y[len_exp-4] == 0):
        len_exp -= 1
    
    x = x[:len_exp]
    y = y[:len_exp]
    
    basename = os.path.basename(fname)
    sampleID = basename.split('_')[0]
    
    # Determine setname based on filename
    setname = 'takara TCR '
    if 'outCellecta' in fname:
        setname = 'cellecta '
    else:
        setname = 'takara '
    
    if 'TCR' in fname:
        setname += 'TCR '
    else:
        setname += 'BCR '
    
    # Generate and save four different plots
    for i in range(4):
        plt.plot(x, y)
        plt.xlabel('UMI type numbers')
        plt.ylabel('Frequency of UMI types')
        plt.title(f'x=UMI types and y=Frequency of UMI types with {setname}{sampleID}')
        
        if i == 1:
            plt.xlim([0, 50])
        if i == 2:
            plt.xlim([0, 50])
            plt.ylim([0, 20000])
        if i == 3:
            plt.ylim([0, 2000])
        
        fwname = fname[:-4] + str(i) + '.png'
        if i == 0:
            fwname = fname[:-4] + '.png'
        
        plt.savefig(fwname)
        plt.close()

def drawUMIdistribution():
    """
    Draws UMI distributions from specified CSV files and saves the corresponding plots.

    Args:
    None

    Returns:
    None
    """
    for i in range(9, 11):
        frname = f'C:\\Users\\suns\\Downloads\\John_peak\\outmTCRv2ilku1\\John_peak_rmumi\\S{i}_umi_distribution.csv'
        with open(frname, 'r') as fr:
            a = fr.readline().replace('umi_count,', '').split(',')
            b = fr.readline().replace('frequency,', '').split(',')
        
        x = [int(i) for i in a]
        y = [int(i) for i in b]
        drawplot(x, y, frname[:-4] + '.png')

def get_meta_file(p_args):
    """
    Generates a metadata CSV file for FASTQ files in the specified directory.

    Args:
    p_args (argparse.Namespace): Contains various arguments including fastq_dir and receptor_kit.

    Returns:
    tuple: 
        - fastq_dir (str): Directory containing the FASTQ files.
        - meta_file (str): Path to the generated metadata CSV file.
    """
    basename = os.path.basename(p_args.fastq_dir)
    p_args.fastq_dir = os.path.dirname(p_args.fastq_dir)
    fastqs = find_file(basename, p_args.fastq_dir)
    fastqs.sort()
    
    p_args.meta_file = os.path.join(p_args.fastq_dir, p_args.receptor_kit + '_meta.csv')
    with open(p_args.meta_file, 'w') as f:
        f.write('sampleID,read1_file_name,read2_file_name\n')
        sampleID = ''
        for fqname in fastqs:
            fqname_base = os.path.basename(fqname)
            if fqname_base.split('_R')[0] != sampleID:
                sampleID = fqname_base.split('_R')[0]
                f.write(sampleID + ',' + fqname_base)
            else:
                f.write(',' + fqname_base + '\n')
    
    if platform.system() == 'Windows':
        dos2unix(p_args.meta_file)
    
    return p_args.fastq_dir, p_args.meta_file

def comp_mixcr_report_v3_v4_matrix(v3reportf, v4reportf, receptor):
    """
    Compares clonotype sequences between two MiXCR reports (v3 and v4).

    Args:
    v3reportf (str): Path to the MiXCR v3 report file.
    v4reportf (str): Path to the MiXCR v4 report file.
    receptor (str): The receptor type, either 'BCRv2' or other (assumed 'TCR').

    Returns:
    list: A list containing statistics on clonotype counts and overlaps.
    """
    v3matrix = []
    v4matrix = []
    results = []
    chains = ['TRA', 'TRB'] if receptor != 'BCRv2' else ['IGG', 'IGM', 'IGK', 'IGL', 'IGA', 'IGD', 'IGE', 'IGH']
    id_chain = 7
    id_seq = 5
    
    for _ in chains:
        v3matrix.append([])
        v4matrix.append([])
    
    for i, vf in enumerate([v3reportf, v4reportf]):
        with open(vf, 'r') as file:
            csvFile = csv.reader(file)
            for count, line in enumerate(csvFile):
                if count > 0:
                    chain_type = line[id_chain]
                    seq = line[id_seq]
                    for id, chain in enumerate(chains):
                        if chain_type == chain:
                            if i == 0:
                                v3matrix[id].append(seq)
                            else:
                                v4matrix[id].append(seq)
    
    for id, chain in enumerate(chains):
        v3set = set(v3matrix[id])
        v4set = set(v4matrix[id])
        interset = v3set.intersection(v4set)
        overlap = len(interset)
        results.append([len(v3matrix[id]), len(v3set), overlap, len(v4matrix[id]), len(v4set)])
        print(len(v3matrix[id]), len(v3set), overlap, len(v4matrix[id]), len(v4set))
    
    return results

def plot_multi_bar(results, fb):
    """
    Plots a multi-bar chart comparing clonotype counts from two MiXCR versions.

    Args:
    results (list): List of clonotype count statistics.
    fb (str): The filename to save the bar plot.

    Returns:
    None
    """
    N = 3 * len(results) - 1
    bases = (results[0][2], results[0][2])
    diffs = (results[0][1] - results[0][2], results[0][4] - results[0][2])
    
    for i in range(1, len(results)):
        bases += (0, results[i][2], results[i][2])
        diffs += (0, results[i][1] - results[i][2], results[i][4] - results[i][2])
    
    ind = np.arange(N)
    width = 0.35
    
    fig, _ = plt.subplots(figsize=(1 * N, 7))
    p1 = plt.bar(ind, bases, width)
    p2 = plt.bar(ind, diffs, width, bottom=bases)
    
    plt.ylabel('Clonotype count')
    tname = os.path.basename(fb).split('_bar.png')[0].split('_')[-1]
    
    if len(results) == 2:
        tname += ' TRA,TRB clonotype count MIXCR v3 vs v4'
        xlab = ('v3TRA', 'v4TRA', '', 'v3TRB', 'v4TRB')
    else:
        tname += ' IgG,IgM,IgK,IgL,IgA,IgD,IgE,IgH clonotype count MIXCR v3 vs v4'
        xlab = ('v3IgG', 'v4IgG', '', 'v3IgM', 'v4IgM', '', 'v3IgK', 'v4IgK', '', 'v3IgL', 'v4IgL', '', 
                'v3IgA', 'v4IgA', '', 'v3IgD', 'v4IgD', '', 'v3IgE', 'v4IgE', '', 'v3IgH', 'v4IgH')
    
    plt.title(tname)
    plt.xticks(ind, xlab)
    plt.legend((p1[0], p2[0]), ('bases', 'diffs'))
    
    plt.savefig(fb)
    plt.close()

def get_modify_fastq_head(ifd, ifr, barcode, ofd):
    """
    Modifies the header of FASTQ files to include a barcode.

    Args:
    ifd (str): Input directory containing the FASTQ files.
    ifr (str): Input FASTQ file name.
    barcode (str): Barcode to be included in the header.
    ofd (str): Output directory to save the modified FASTQ files.

    Returns:
    None
    """
    fwname = os.path.join(ofd, ifr)
    fr = os.path.join(ifd, ifr)
    
    with open(fr, 'r') as f, open(fwname, 'w') as fw:
        for count, line in enumerate(f, 1):
            if count % 4 == 1:
                segs = line.split(' ')
                newline = f'{segs[0]}_{barcode} {segs[1]}'
                fw.write(newline)
            else:
                fw.write(line)

def get_demux_counts_all_csv_fastq(ifpath, ofd):
    """
    Generates modified FASTQ files for all barcodes listed in a CSV file.

    Args:
    ifpath (str): Path to the input CSV file containing barcode and key information.
    ofd (str): Output directory to save the modified FASTQ files.

    Returns:
    None
    """
    ifd = os.path.dirname(ifpath)
    
    with open(ifpath, 'r') as f:
        for line in f:
            segs = line.split(',')
            if len(segs[0]) == 16:
                key = segs[1].replace('_', '')
                get_modify_fastq_head(ifd, f'demux_{segs[0]}_{key}_R1.fastq.gz', segs[0], ofd)
                get_modify_fastq_head(ifd, f'demux_{segs[0]}_{key}_R2.fastq.gz', segs[0], ofd)

def splitfastq_by_adapter(sample, fastq_dir, adapterseq):
    """
    Splits FASTQ files based on the presence of a specific adapter sequence.

    Parameters:
    - sample: The sample name (str).
    - fastq_dir: Directory containing FASTQ files (str).
    - adapterseq: Adapter sequence to search for in the reads (str).

    Creates two FASTQ files for each read pair:
    - One with reads containing the adapter sequence.
    - One without the adapter sequence.
    """
    r1, r2 = [], []
    raw_r1_fq = os.path.join(fastq_dir, sample + "_L001_R1_001.fastq.gz")
    raw_r2_fq = os.path.join(fastq_dir, sample + "_L001_R2_001.fastq.gz")
    out_r1_fqs, out_r2_fqs = [], []
    prs = ["with-adapter", "without-adapter"]
    
    # Create output directories if they don't exist
    for directory in prs:
        if not os.path.isdir(os.path.join(fastq_dir, directory)):
            create_dir(os.path.join(fastq_dir, directory))
        out_r1_fqs.append(os.path.join(fastq_dir, directory, sample + '_L001_R1_001.fastq'))
        out_r2_fqs.append(os.path.join(fastq_dir, directory, sample + '_L001_R2_001.fastq'))
    
    # Create empty FASTQ files for output
    for item in out_r1_fqs + out_r2_fqs:
        open(item, 'w+').close()
    
    fhs1, fhs2 = [], []
    # Process FASTQ files
    with load_fastq(raw_r1_fq) as f1, load_fastq(raw_r2_fq) as f2, ExitStack() as cm:
        for name in out_r1_fqs:
            fhs1.append(cm.enter_context(open(name, 'a')))
        for name in out_r2_fqs:
            fhs2.append(cm.enter_context(open(name, 'a')))
        
        # Split reads based on the presence of adapter sequence
        for x, y in zip(f1, f2):
            r1.append(x.rstrip())
            r2.append(y.rstrip())
            if len(r1) == 4:
                if r1[1].find(adapterseq) >= 0 or r2[1].find(adapterseq) >= 0:
                    for item in r1:
                        fhs1[0].write(f'{item}\n')
                    for item in r2:
                        fhs2[0].write(f'{item}\n')
                else:
                    for item in r1:
                        fhs1[1].write(f'{item}\n')
                    for item in r2:
                        fhs2[1].write(f'{item}\n')
                r1.clear()
                r2.clear()
    
    # Convert to Unix format if on Windows
    if platform.system() == 'Windows':
        for name in out_r1_fqs + out_r2_fqs:
            dos2unix(name)
    
    # Compress the output FASTQ files
    for directory in prs:
        for f in os.listdir(os.path.join(fastq_dir, directory)):
            if f.endswith('.fastq'):
                fq = os.path.join(fastq_dir, directory, f)
                cmd = f'tar -vczf {fq}.gz {fq}'
                run_system_cmd(cmd, 'run_gzip.log', 'run_gzip.error', 'gzip failed', 'NA')
                os.remove(fq)

def clonetypes_comparison(fp1, fp2):
    """
    Compares two CSV files containing clone types and merges them.

    Parameters:
    - fp1: Filepath to the first CSV file (str).
    - fp2: Filepath to the second CSV file (str).

    Returns:
    - count_same: Number of matching clone types (int).
    - non_matching_count: Number of non-matching clone types (int).
    """
    with open(fp1, 'r') as f1, open(fp2, 'r') as f2:
        flist1, flist2 = f1.readlines(), f2.readlines()
        fwlist = flist1.copy()
        count_same = 0 
        
        for j in range(1, len(flist2)):
            segs2 = flist2[j].split(",")
            i = 1
            is_same = False
            while not is_same and i < len(flist1):
                segs1 = flist1[i].split(",")
                if all(segs2[k] == segs1[k] for k in [2, 11, 13, 15, 17]):
                    count_same += 1
                    is_same = True
                i += 1
            if not is_same:
                fwlist.append(flist2[j])
    
    fw = fp1.replace('_IGG_', '_IGGMerged_')
    with open(fw, 'w') as f:
        f.writelines(fwlist)
    
    return count_same, len(flist2) - 1 - count_same

def dosformat(fpath):
    """
    Checks if a file is in DOS format.

    Parameters:
    - fpath: Filepath to check (str).

    Returns:
    - is_dosformat: Boolean indicating if the file is in DOS format (bool).
    """
    with open(fpath, 'rb') as f:
        line = f.readline()
        is_dosformat = any(sub in line for sub in [b'\r\n', b'\n\r', b'\r'])
    return is_dosformat

def plot_QC_chains(fname):
    """
    Generates pie charts from a CSV file for all QC chains.

    Parameters:
    - fname: Filepath to the CSV file (str).

    The function saves pie chart images for each sample in the CSV file.
    """
    with open(fname, 'r') as f:
        dirname = os.path.dirname(fname)
        flist = f.readlines()
        labels = []

        for i, line in enumerate(flist):
            segs = line.split(",")
            y, mylabels = [], []
            for j in range(1, len(segs) - 2):
                if i == 0:
                    if j % 2 == 1:
                        labels.append(segs[j])
                else:
                    if j % 2 == 1:
                        value = int(segs[j])
                        if value > 1:
                            y.append(value)
                    if j % 2 == 0 and value > 1:
                        mylabels.append(labels[j // 2 - 1] + '_' + segs[j])
            
            if i > 0:
                ofname = os.path.join(os.path.dirname(dirname), 'QC_pie_plot_' + segs[0] + ".png")
                y[1], y[2] = y[2], y[1]
                mylabels[1], mylabels[2] = mylabels[2], mylabels[1]
                if len(y) > 7:
                    y[5], y[6] = y[6], y[5]
                    mylabels[5], mylabels[6] = mylabels[6], mylabels[5]
                myexplode = [0] * len(y)
                myexplode[-1] = 0.2
                plt.pie(y, labels=mylabels, explode=myexplode)
                plt.title(f"QC plot sample_{segs[0]}")
                plt.savefig(ofname)
                plt.close()

def flc_analysis(fname):
    """
    Analyzes and generates pie charts based on a CSV file with specific data columns.

    Parameters:
    fname (str): The path to the CSV file to be analyzed.

    The function performs the following steps:
    1. Reads the CSV file line by line.
    2. Extracts and organizes data, skipping the first column and last two columns.
    3. Converts certain columns to appropriate data types for analysis.
    4. Creates and saves pie charts based on specific conditions in the data.

    Output:
    Saves pie charts as PNG files in the same directory as the input file.
    """
    # Extract directory name from the file path
    dirname = os.path.dirname(fname)
    
    count = 0
    heads = []
    data = []
    
    # Read the file line by line
    with open(fname, 'r') as fp:
        for line in fp:
            line = line.strip()
            segs = line.split(',')
            if count == 0:
                # Extract header row
                heads = segs.copy()
            else:
                # Extract data, skipping first and last two columns
                data.append(segs[1:len(segs)-2])
            count += 1
    
    # Create a DataFrame from the data
    df = pd.DataFrame(data=data, columns=heads[1:len(heads)-2])
    
    # Total number of rows in the DataFrame
    total = len(df)
    
    # Convert specific columns to appropriate data types
    df['pos_GTACGGG'] = df['pos_GTACGGG'].astype('int32')
    df['p0_40'] = df['p0_40'].astype('int32')
    df['sim(GTACGGG)'] = df['sim(GTACGGG)'].astype('float32')
    df['p0_40_sim_max'] = df['p0_40_sim_max'].astype('float32')
    
    # Loop through keys 'pos' and 'sim_max' to generate pie charts
    for key in ['pos', 'sim_max']:
        y = []
        mylabels = []
        
        # Customize behavior based on the key
        if key == 'pos':
            df = df.sort_values(by=['pos_GTACGGG', 'p0_40_sim_max'], ascending=False)
            labels = [
                'p-1(0_40sim<3):', 'p-1(0_40sim=3):', 'p-1(0_40sim=4):', 'p-1(0_40sim=5):', 
                'p-1(0_40sim=6):', 'p0-10(sim=7):', 'p11(sim=7):', 'p13(sim=7):', 
                'p14-40(sim=7):', 'p41-117(sim=7):'
            ]
            
            # Generate data for pie chart based on different conditions
            for i in range(len(labels)):
                if i == 0:
                    sub_df = df.loc[(df['pos_GTACGGG'] == -1) & (df['p0_40_sim_max'] < 3)]
                elif i == 1:
                    sub_df = df.loc[(df['pos_GTACGGG'] == -1) & (df['p0_40_sim_max'] == 3)]
                elif i == 2:
                    sub_df = df.loc[(df['pos_GTACGGG'] == -1) & (df['p0_40_sim_max'] == 4)]
                elif i == 3:
                    sub_df = df.loc[(df['pos_GTACGGG'] == -1) & (df['p0_40_sim_max'] == 5)]
                elif i == 4:
                    sub_df = df.loc[(df['pos_GTACGGG'] == -1) & (df['p0_40_sim_max'] == 6)]
                elif i == 5:    
                    sub_df = df.loc[(df['pos_GTACGGG'] > -1) & (df['pos_GTACGGG'] <= 10)]
                elif i == 6:
                    sub_df = df.loc[(df['pos_GTACGGG'] == 11)]
                elif i == 7:
                    sub_df = df.loc[(df['pos_GTACGGG'] == 13)]
                elif i == 8:
                    sub_df = df.loc[(df['pos_GTACGGG'] > 13) & (df['pos_GTACGGG'] <= 40)]
                elif i == 9:
                    sub_df = df.loc[(df['pos_GTACGGG'] > 40)]
                
                # Append counts and labels for pie chart
                y.append(sub_df.shape[0])
                mylabels.append(labels[i] + str(sub_df.shape[0]) + 
                                "(" + str(int(100 * len(sub_df) / total + 0.5)) + "%)")
        
        else:  # key == 'sim_max'
            df = df.sort_values(by=['p0_40_sim_max', 'pos_GTACGGG'], ascending=False)
            labels = ['0_40sim_max<3:', '0_40sim_max=3:', '0_40sim_max=4:', 
                      '0_40sim_max=5:', '0_40sim_max=6:', '0_40sim_max=7:']
            
            # Generate data for pie chart based on different conditions
            for i in range(len(labels)):
                if i == 0:
                    sub_df = df.loc[df['p0_40_sim_max'] < 3]
                elif i == 1:
                    sub_df = df.loc[df['p0_40_sim_max'] == 3]
                elif i == 2:    
                    sub_df = df.loc[df['p0_40_sim_max'] == 4]
                elif i == 3:
                    sub_df = df.loc[df['p0_40_sim_max'] == 5]
                elif i == 4:
                    sub_df = df.loc[df['p0_40_sim_max'] == 6]
                elif i == 5:
                    sub_df = df.loc[df['p0_40_sim_max'] == 7]
                
                # Append counts and labels for pie chart
                y.append(sub_df.shape[0])
                mylabels.append(labels[i] + str(sub_df.shape[0]))
        
        # Generate file name for the pie chart
        ofname = os.path.join(dirname, 'flc_pie_plot_' + key + '_' + segs[0] + ".png")
        
        # Create and save the pie chart
        plt.pie(y, labels=mylabels)
        plt.title("Pie plot sample_" + key + '_' + segs[0])
        plt.savefig(ofname, dpi=250)
        plt.close()
    
    return

# ---------- main ---------- #
if __name__ == '__main__':
    
    # ---------- setup | parser ---------- #
    # Initialize the argument parser with a description of the program.
    parser = argparse.ArgumentParser(description=desc)
    
    # Separate optional and required argument groups for better organization.
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)

    # ---------- setup | user options ---------- #
    
    # Add required arguments for the script.
    required.add_argument('-r', '--receptor_kit', dest='receptor_kit',
                          help='specify receptor kit: TCRv1, TCRv2, BCRv1, or BCRv2', 
                          choices=['TCRv1', 'TCRv2', 'BCRv1', 'BCRv2'],
                          required=True)
    required.add_argument('-f', '--fastq_dir', dest='fastq_dir',
                          help='Directory containing all input FASTQ files', 
                          required=True)
    required.add_argument('-m', '--meta_file', dest='meta_file',
                          help='File containing sample IDs and corresponding FASTQ pairs',
                          required=True)
    required.add_argument('-o', '--output_name', dest='out_name',
                          help='Name for the output directory and file prefix; should be less than 20 characters',
                          required=True)
    required.add_argument('-t', '--target_region', dest='target_region',
                          help='Specify the target regions reads should map to',
                          choices=['CDR3', 'Full_length', 'Both'],
                          required=True)

    # Add optional arguments for additional configurations.
    optional.add_argument('-k', '--keep_inter_file', dest='keep_inter_file',
                          help='Decide whether to keep intermediate files, including MiXCR files and preprocessed FASTQs [Default: False]',
                          action='store_true', default=False)
    optional.add_argument('-l', '--linker_correction', dest='linker_correction',
                          help='Decide whether to remove reads based on sequence match of linker [Default: False]',
                          action='store_true', default=False)
    optional.add_argument('-s', '--species', dest='species',
                          help='Specify the genome species: human, mouse [Default: human]',
                          default='human')
    optional.add_argument('-u', '--umi_cutoff', dest='umi_cutoff',
                          help='Specify an integer as the UMI cutoff [Default: \'\']',
                          default='')
    optional.add_argument('-e', '--memory_size', dest='memory_size',
                          help='Specify memory allocation for the project [Default: 32]',
                          default='32')

    # ---------- setup | parse/check user options ---------- #
    # Parse the arguments provided by the user.
    p_args = parser.parse_args()
    
    # Additional logic based on parsed arguments.
    if p_args.receptor_kit == 'TCRv1' or (p_args.receptor_kit == 'BCRv1' and p_args.species == 'mouse'):
        p_args.linker_correction = False
        p_args.umi_cutoff = ''
    
    # Ensure the FASTQ directory path does not end with a slash.
    if p_args.fastq_dir[-1] in ['/', '\\']:
        p_args.fastq_dir = p_args.fastq_dir[:-1]
    
    # Initialize additional parameters.
    p_args.internal_test = False
    p_args.cogentIP_version = 'v1.6'
    p_args.process_threads = ''
    p_args.library_kit = 'takara_smartseq'
    p_args.rm_adapter = False
    p_args.vdjviz_inputs = False
    p_args.allow_partial = False
    
    # Conditional settings based on the version of cogentIP.
    if p_args.cogentIP_version != 'v2':
        p_args.airr_create = True
    else:
        p_args.airr_create = False
    
    # Validate UMI cutoff value.
    if p_args.umi_cutoff and not p_args.umi_cutoff.isdigit():
        print(f"The umi_cutoff '{p_args.umi_cutoff}' specified by user is not an integer, please modify and relaunch")
        sys.exit(1)
    
    # Additional species and receptor kit checks.
    if p_args.species == 'mouse' and p_args.receptor_kit != 'TCRv2':
        print('Support for mouse species is only available for the TCRv2 receptor kit.')
        sys.exit(1)
    if p_args.receptor_kit == 'BCRv2hy' and p_args.species == 'mouse':
        print(f"The receptor_kit '{p_args.receptor_kit}' specified is for human, please modify and relaunch.")
        sys.exit(1)
    
    # ---------- setup | start ---------- #
    vjava = required_python_module_check()  # Check for required Python modules.
    start_time = datetime.datetime.now()  # Record the start time of the script.

    # Validate the output name.
    if '/' in p_args.out_name or '\\' in p_args.out_name:
        print('The output name should be a string, not a directory path; a new folder with this name will be created in the metadata file directory.')
        sys.exit(1)
    
    # --------- setup aligner | end ---------- #
    # Initialize the pairwise aligner with specific scoring parameters.
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    
    # Extract FASTQ files based on the directory name and meta file.
    basename = os.path.basename(p_args.fastq_dir)
    if '*' in basename or '?' in basename:
        p_args.fastq_dir, p_args.meta_file = get_meta_file(p_args)

    # Set the main and output directory paths.
    p_args.main_dir = os.path.split(os.path.abspath(p_args.meta_file))[0]
    p_args.out_dir = os.path.join(p_args.main_dir, p_args.out_name)

    # Check the existence of the FASTQ directory and meta file.
    check_dir(p_args.fastq_dir, 'does not exist!\n')
    if p_args.meta_file:
        check_file(p_args.meta_file, 'does not exist!\n')

    # ---------- setup | load meta, create outdir ---------- #
    # Load metadata from the provided file and create the output directory.
    meta_dict = load_meta(p_args.meta_file, p_args.fastq_dir)
    print(get_time() + ' [INFO] Loading meta data file specified by user', flush=True)

    # ---------- delete the folder only for internal_test or debug ----- #
    debug = False  # Set to True for debugging.
    logger = logging.getLogger('temp_recoder')
    
    # Remove the output directory if it's for internal testing or debugging.
    if os.path.isdir(p_args.out_dir) and (p_args.internal_test or debug):
        logging.shutdown()
        shutil.rmtree(p_args.out_dir)

    # Check if the output directory already exists or if the output name is invalid.
    if os.path.isdir(p_args.out_dir):
        print(get_time() + ' [INFO] Analysis dir already exists: ' + p_args.out_dir, flush=True)
        sys.exit(1)
    elif p_args.out_name == os.path.basename(os.path.normpath(p_args.main_dir)):
        print('The output folder name is identical to its upper folder, please rename', flush=True)
        sys.exit(1)
    else:
        try:
            os.makedirs(p_args.out_dir)
        except OSError as err:
            print(get_time() + ' [ERROR] Unable to create directory: ' + p_args.out_dir, flush=True)
            sys.exit(1)

    # ---------- setup | logger ---------- #
    # Initialize logging for the script.
    LOG_NAME = 'immune_profiler'
    LOG_FILENAME = os.path.join(p_args.out_dir, p_args.out_name + '_immune_profiler.log')
    LOG_FORMAT = '%(levelname)s - %(asctime)s - %(name)s - %(message)s'
    LOG_FORMAT_ARG_DATE = '%Y-%m-%d %H:%M:%S'

    logger = logging.getLogger(LOG_NAME)
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(LOG_FILENAME, mode='w')
    formatter = logging.Formatter(LOG_FORMAT, datefmt=LOG_FORMAT_ARG_DATE)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.info('Started script execution : ' + script_name + '_v' + version + ' with ' + vjava)
    logger.info('Original call : ' + ' '.join(sys.argv))

    # ----------- setup | define parameter & local resources----------- #
    # Define various parameters and paths for the workflow.
    umi_length = 12
    linker_sequence = 'GTAC'
    solution = []
    delta = lambda x, y, i, j: 1 if x[i] != y[j] else 0

    # Set up the repository path based on the script's location
    REPO_PATH = os.path.dirname(os.path.realpath(__file__))
    
    # Define the base command for running the Java ImmuneProfiler tool with a default memory allocation of 32GB
    workflow_base = 'java -Xmx32G -jar ' + os.path.join(REPO_PATH, 'src', 'ImmuneProfiler'+p_args.cogentIP_version+'.jar')
    
    # Adjust the memory allocation for the Java command if the user specifies a different memory size
    if p_args.memory_size != "32":
        workflow_base = 'java -Xmx' + p_args.memory_size + 'G -jar ' + os.path.join(REPO_PATH, 'src', 'ImmuneProfiler'+p_args.cogentIP_version+'.jar')
    
    # On Windows systems, simplify the Java command and limit the number of process threads if needed
    if platform.system() == 'Windows' :
        workflow_base = 'java -jar ' + os.path.join(REPO_PATH, 'src', 'ImmuneProfiler'+p_args.cogentIP_version+'.jar')
        if len(p_args.process_threads) > 0 and int(p_args.process_threads) > 4:
            p_args.process_threads = '4'
    
    # Modify the workflow command if using version 2 of the ImmuneProfiler
    if p_args.cogentIP_version.find("v2") >= 0 :
        workflow_base = workflow_base[:-6] + 'v2.jar'
    
    # Set up the links for various workflow steps
    workflow0_link = workflow_base + " " + p_args.out_dir
    workflow1_link = workflow_base + " " + p_args.out_dir
    workflow2_link = workflow_base + " " + p_args.out_dir
    
    # Enable AIRR (Adaptive Immune Receptor Repertoire) creation if using version 2 of the ImmuneProfiler
    if workflow2_link.find('ImmuneProfilerv2.jar') >= 0:
        p_args.airr_create = True
    
    # Initialize dictionaries for tracking alignment and clonotype statistics
    align_stats_dict = OrderedDict([
        ('Total sequencing reads', 'total reads'),
        ('Successfully aligned reads', 'aligned reads'), 
        ('Overlapped', 'overlapped'), 
        ('Overlapped and aligned', 'overlapped and aligned')
    ])
    clone_stats_dict = OrderedDict([('Final clonotype count', 'clonotype count')])
    
    # If running an internal test, use a more detailed set of clonotype statistics
    if p_args.internal_test:
        clone_stats_dict = OrderedDict([
            ('Average number of reads per clonotype', 'reads per clonotype'),
            ('Reads used in clonotypes, percent of total', 'Reads used in clonotypes'),
            ('Mapped low quality reads, percent of used', 'Mapped low quality reads'),
            ('Reads clustered in PCR error correction, percent of used', 'PCR error correction'),
            ('Reads dropped due to failed mapping, percent of total', 'failed mapping'),
            ('Final clonotype count', 'clonotype count')
        ])
    
    # ---------- Load barcodes and receptor information -------- #
    # Create a directory for storing reports
    create_dir(os.path.join(p_args.out_dir, 'report'))
    fd = OrderedDict()
    
    # ---------- Preprocessing Step -------- #
    # If using the Takara SmartSeq library kit, start the preprocessing workflow
    if (p_args.library_kit == 'takara_smartseq'):
        create_dir(os.path.join(p_args.out_dir, 'preprocess'))
        logger.info('Starting preprocessing')
        print(get_time() + ' [INFO] Starting preprocessing', flush=True)
        sd = OrderedDict()
    
        # Construct and run the preprocessing command
        prep_cmd = workflow0_link + ' workflow0 ' + p_args.receptor_kit + ' ' + p_args.fastq_dir \
                 + ' ' + p_args.meta_file + ' ' + p_args.out_dir + ' ' + p_args.species + ' ' + REPO_PATH
        if p_args.linker_correction : 
            prep_cmd += " -l"
        logger.info('Launching PREPROCESS using command:\n' + prep_cmd)
        log_name = 'run_prepro.log'
        out_log = os.path.join(p_args.out_dir, 'preprocess', log_name)
        run_prepro_err = os.path.join(p_args.out_dir, 'preprocess', 'run_prepro.error')
        run_system_cmd(prep_cmd, out_log, run_prepro_err, 'preprocess failed', log_name)
    
        # Load the sample data (sd) dictionary from the preprocessing results
        sd = load_sd(os.path.join(p_args.out_dir, 'preprocess', 'sd_dict.csv'))
    
        # Write the quality control (QC) reports for the samples
        write_sample_qc(p_args, sd)
        if p_args.airr_create and workflow2_link.find('ImmuneProfilerv1.6.jar') >= 0:
            create_dir(os.path.join(p_args.out_dir, 'airr_report'))
            write_airr_sample_qc(p_args, sd)
    
        # Create the final data (fd) dictionary from the sample data
        fd = create_fd(sd)
        logger.info('Completed preprocessing')
        print(get_time() + ' [INFO] Completed preprocessing', flush=True)
    
    # ---------- Run MIGEC (Molecular Identifier Guided Error Correction) ---------- #
    # If applicable, run the MIGEC analysis for UMI-guided correction
    if workflow2_link.find('ImmuneProfilerv1.6.jar') >= 0 and p_args.receptor_kit != 'TCRv1' and (p_args.receptor_kit != 'BCRv1' or p_args.species != 'mouse'):
        create_dir(os.path.join(p_args.out_dir, 'run_migec'))
        logger.info('Starting UMI guided analysis…')
        print(get_time() + ' [INFO] Starting UMI guided analysis…', flush=True)
        migec_start_time = datetime.datetime.now()
    
        # Prepare the input files for MIGEC
        print(get_time() + ' [INFO] Preparing MIGEC input files', flush=True)
        if p_args.receptor_kit != "TCRv1" and (p_args.receptor_kit != "BCRv1" or p_args.species != "mouse"):
            prep_migec_input(p_args, fd)
    
        # Run the MIGEC analysis
        print(get_time() + ' [INFO] Conducting MIGEC analysis')
        run_migec(p_args, 'run_migec', 'mig_run_migec.log', 'mig_run_migec.error')
    
        # Summarize the MIGEC results
        print(get_time() + ' [INFO] Summarizing MIGEC results', flush=True)
        summ_migec(p_args)
        logger.info('Finished collecting MIGEC statistics')
        logger.info('Finished MIGEC process')
        print(get_time() + ' [INFO] Finished MIGEC process', flush=True)
    
        # Clean up intermediate files if not needed
        if not p_args.keep_inter_file and not debug:
            pre_dir = os.path.join(p_args.out_dir, 'preprocess')
            pre_fqs = find_file('*.*', pre_dir)
            for fq in pre_fqs:
                os.remove(fq)
            os.rmdir(pre_dir)
    
    # ---------- Run MiXCR on MIGEC results or preprocess data ---------- #
    # If applicable, run the MiXCR analysis for read alignment and clonotype calling
    if os.path.isdir(os.path.join(p_args.out_dir, 'run_migec')) or p_args.receptor_kit == 'BCRv2hy':
        if os.path.exists(os.path.join(p_args.out_dir, 'run_mixcr')):
            change_run_mixcr_name(p_args, 'prePro')
        else: 
            create_dir(os.path.join(p_args.out_dir, 'run_mixcr'))
    
        # Set the read type and preprocessing directory for MiXCR
        read_type = 'mig'
        pre_dir = 'run_migec'
    
        # Run MiXCR for CDR3 and/or full-length regions
        print(get_time() + ' [INFO] Starting reads alignment, assembling, and clonotype calling', flush=True)
        if p_args.target_region == 'Both' or p_args.target_region == 'CDR3':
            mig_cdr3_fd = copy.deepcopy(fd)
            if p_args.receptor_kit == "BCRv2hy":
                run_mixcr_preIgGA(p_args, mig_cdr3_fd, read_type, 'cdr3', 'run_mixcr', pre_dir, 'run_mixcr_cdr3_sIg.log', 'run_mixcr_cdr3_sIg.error')
                change_run_mixcr_name(p_args, 'preSub')
            run_mixcr(p_args, mig_cdr3_fd, read_type, 'cdr3', 'run_mixcr', pre_dir, 'run_mixcr_cdr3.log', 'run_mixcr_cdr3.error')
    
        if p_args.target_region == 'Both' or p_args.target_region == 'Full_length':
            mig_fl_fd = copy.deepcopy(fd)
            if p_args.receptor_kit == "BCRv2hy":
                run_mixcr_preIgGA(p_args, mig_fl_fd, read_type, 'fl', 'run_mixcr', pre_dir, 'run_mixcr_cdr3_sIg.log', 'run_mixcr_cdr3_sIg.error')
                change_run_mixcr_name(p_args, 'preSub')
            run_mixcr(p_args, mig_fl_fd, 'mig', 'fl', 'run_mixcr', 'run_migec', 'run_mixcr_fl.log', 'run_mixcr_fl.error')
    
        logger.info('Finished MiXCR analysis')
    
    # ---------- Clean Up and Finish ---------- #
    # Remove unnecessary files and directories based on species and debugging settings
    if p_args.species == 'rhesus_monkey':
        os.remove(os.path.join(p_args.main_dir, 'imgt.201946-3.sv6.json'))
    else:
        pass
    logger.info('Analysis completed')
    print(get_time() + ' [INFO] Analysis completed', flush=True)
    sys.exit(0)
