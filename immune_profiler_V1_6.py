#!/usr/bin/env python3
# Author: Sean Sun

script_name = 'Cogent_NGS_immune_profiler'
version = 'CogentIP v1.6'
desc = """
       CogentIP v1.6
            Done:
            use new ImmuneProfilerv1.6.jar linked to MIGEC1.2.9 and MIXCR3
            add BCRv1,BCRv2,TCRv1,TCRv2 as choices in argument
            vdjviz input file generation
            add in customized UMI cutoff; not via metafile but by argument
            
            Todo:
            enable human, mouse
                tested: without -s argument(default human), -s human, -s mouse, 
                picocli.CommandLine$ExecutionException:
            remove analysis without MIGEC step
            rewrite run_mixcr and run_migec function
            double check error write for MIGEC and MiXCR
            catch error of MIGEC and MiXCR
            file finding error convert to user-understandable text

            In report_stats(file_dir, file_name, metaDict), the file_name should contain the out_name, thus in run_mixcr()
            and run_raw_mixcr() functions, changed following lines:
            report_stats(os.path.join(p_args.out_dir, mixcr_dir, process_type), p_args.out_name + '_' + process_type, cur_dict)
            report_stats(os.path.join(p_args.out_dir, mixcr_dir, process_type), p_args.out_name + '_' + process_type, cur_dict)
            basically, the stats filename now contains prefix with out_name: p_args.out_name + process_type

            check if out_name is string or directory string
"""

# ---------- import modules ---------- #

import argparse
#from curses import meta
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
#from Bio import pairwise2
from Bio import Align # Bio.Align.PairwiseAligner instead of pairwise2


def required_python_module_check():
    required_pkgs = ["argparse", "csv", "logging", "re", "collections",
               "datetime", "openpyxl", "subprocess", "contextlib", "os",
               "fnmatch", "sys", "copy", "gzip", "platform", "numpy",
               "time", "shutil", "matplotlib.pyplot", "Bio", "itertools"]
    # declare an empty list to track missing packages
    not_installed = []
    sys.stdout.write('\nChecking for python packages required by Immune ' +
                 'Profiler...\n')
    
    # import check. If an error is detected, it will add the name of
    # the package to the not_installed list
    for x in required_pkgs:
        try: 
            __import__(x)
        except ImportError as err:
            not_installed.append(x)
    
    # Will print out each module that is missing
    if len(not_installed) > 0:
        sys.stdout.write('\n')
        for package in not_installed:
            sys.stdout.write('\t-- ' + package + ' is not installed.\n')
        sys.stdout.write('\nSee https://docs.python.org/3/installing/ for ' +
                     'help with package installation.\n')
        sys.exit(0)
    # if no packages failed to import, writes a success message to the screen.
    else:
        sys.stdout.write('\nAll required packages are installed.\n\n')
    
    # check whether java version fit cogentIP version
    vjava = subprocess.check_output(['java', '-version'], stderr=subprocess.STDOUT)
    vjava = vjava.decode("utf-8").split("\n")[0]
    segs = vjava.split('"')
    vjava = segs[1].split('"')[0]
    vm = int(vjava.split('.')[0])
    match = re.match(r'([a-z]+)([0-9]+)', p_args.cogentIP_version, re.I)
    if match and int(match.groups()[1]) < 2 and vm > 15: #11
        print('your java version: %s, required v1.8 or v8-15' % vjava)
        sys.exit(0)
    
    # check installed full path with/without space
    if len( p_args.meta_file.split(' ') ) > 1:
        print('The address for CogentIP installation should be string without space')
        sys.exit(0)
    return vjava

# ---------- fxn | check directory, file, get time & set timer ------- #
def check_dir(d, exp):
    if os.path.isdir(d) is False:
        print(get_time() + ' [ERROR] ' + d + ' ' + exp, flush=True)
        sys.exit(1)
    else:
        return True


def check_file(f, exp):
    if os.path.isfile(f) is False:
        print(get_time() + ' [ERROR] ' + f + ' ' + exp, flush=True)
        sys.exit(1)
    else:
        return True


def get_time():
    ts = time.time()
    cur_time = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return cur_time


def timer(start_time):
    end_time = datetime.datetime.now()
    script_time = end_time - start_time
    secs = round(script_time.total_seconds() % 60, 2)
    mins = int(script_time.total_seconds() // 60) % 60
    hour = int(script_time.total_seconds() // 3600)
    diff = ':'.join([str(hour) + 'h', str(mins) + 'm', str(secs) + 's'])
    return diff


# ---------- fxn | load meta data ---------- #
def load_meta(meta_file,fastq_dir):
    meta_dict = OrderedDict()
    try:
        infile = open(meta_file, 'r')
        next(infile)
        for line in infile:
            words = [word.strip() for word in line.rstrip().split(',')]
            sname = ''
            r1_fq = ''
            r2_fq = ''
            if len(words) == 3:
                sname, r1_fq, r2_fq = words[0], words[1], words[2]
            elif len(words) < 3:
                print( get_time() + '[ERROR] sample ' + sname + ' does not have paired-fastqs')
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
            for i in range(len(meta_dict[sname])):
                fq = meta_dict[sname][i]
                if find_file(fq,fastq_dir) == []:
                    print(get_time() + '[ERROR] The input FASTQ: ' + fq + ' is not found. \n' +
                        'Please make sure the FASTQ directory is correct & ' + 
                        'FASTQ name in metadata file is correct (case sensitive) ', flush=True)
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
    bcs = []
    with open(bc_file) as f:
        for line in f:
            l = line.rstrip().split(',')[0]
            bcs.append(l)
    return bcs


def load_pseq(ps_file):
    with open(ps_file) as f:
        ps_dict = OrderedDict()
        for line in f:
            (cur_seq, ori_seq) = re.split(',', line.strip())
            ps_dict[cur_seq] = ori_seq
    return ps_dict


def update_meta(bcs, meta_dict, nchain):
    idx = 0
    for k in meta_dict.keys():
        meta_dict[k].append(bcs[idx:(idx+nchain)])
        idx += nchain
    return

def load_sd(sd_fname):
    sd = OrderedDict()
    with open(sd_fname) as f:
        for line in f:
            segs = line.rstrip().split('\t')
            key = segs[0][:-1]
            values = []
            for i in range(1, len(segs)):
                items = segs[i].split(',')
                for j in range(len(items)):
                    if j==0:
                        item0 = items[j]
                    if j==1 and items[j].isnumeric():
                        item1 = int(items[j])
                    elif j==2 and items[j].replace(".","").isnumeric():
                        item2 = float(items[j])
                mytuple = (item0, item1, item2)
                values.append(mytuple)
            sd[key] = values
    return sd


def change_dir(dirname, distname):
    if os.path.isdir(dirname):
        try:
            os.rename(dirname, distname)
        except OSError as err:
            print('ERROR: Unable to change directory from:' + dirname + " to " + distname, flush=True)
            sys.exit(1)
    else:
        print('ERROR: The folder ' + dirname + 'does not exists', flush=True)
        sys.exit(1)


def create_dir(dirname):
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
    if f.endswith('.gz'):
        f = gzip.open(f, mode = 'rt')
    else:
        f = open(f)
    return(f)


def get_sample_bcs(sample, meta_dict):
    bcs = meta_dict[sample][2]
    return bcs


def ham_check(hh, bc):
    try:
        return hh[bc]
    except KeyError:
        return ''

def get_errorbs_bytwoseqs(seq1, seq2, errbases):
#pre: seq and seq2 have the same length
    errbcount = 0
    j=0
    seqL = len(seq1)
    while errbcount <= errbases and j <seqL:
        if seq1[j]!= seq2[j]:
            errbcount += 1
        j += 1
    if errbcount <= errbases:
        return seq2, errbcount
    return '', errbcount

def get_alia_seqs(seq):
#pred: seq contains character 'R' or 'S' or 'I' or 'Y' or 'K' or 'M'
#post: resturn seqs list including all possible alias
    seqs = []
    seqs.append(seq)
    while len(seqs)>0 and (seqs[0].find('R')>=0 or seqs[0].find('S')>=0 or seqs[0].find('I')>=0 or seqs[0].find('Y')>=0 or seqs[0].find('K')>=0 or seqs[0].find('M')>=0):
        del_num = len(seqs)
        for i in range(len(seqs)):
            if seqs[i].find('R') >= 0:
                myseq = seqs[i].replace('R', 'A', 1)
                seqs.append(myseq)
                myseq = seqs[i].replace('R', 'G', 1)
                seqs.append(myseq)
            elif seqs[i].find('S') >= 0:
                myseq = seqs[i].replace('S', 'C', 1)
                seqs.append(myseq)
                myseq = seqs[i].replace('S', 'G', 1)
                seqs.append(myseq)
            elif seqs[i].find('I') >= 0:
                myseq = seqs[i].replace('I', 'A', 1)
                seqs.append(myseq)
                myseq = seqs[i].replace('I', 'C', 1)
                seqs.append(myseq)
                myseq = seqs[i].replace('I', 'T', 1)
                seqs.append(myseq)
            elif seqs[i].find('Y') >= 0:
                myseq = seqs[i].replace('Y', 'C', 1)
                seqs.append(myseq)
                myseq = seqs[i].replace('Y', 'T', 1)
                seqs.append(myseq)
            elif seqs[i].find('K') >= 0:
                myseq = seqs[i].replace('K', 'G', 1)
                seqs.append(myseq)
                myseq = seqs[i].replace('K', 'T', 1)
                seqs.append(myseq)
            elif seqs[i].find('M') >= 0:
                myseq = seqs[i].replace('M', 'A', 1)
                seqs.append(myseq)
                myseq = seqs[i].replace('M', 'C', 1)
                seqs.append(myseq)
        if len(seqs) > del_num:
            del seqs[0:del_num]
    return seqs

# pre: prs[0] and bc should be same length
def ham_check_flex(prs, bc):
    errbs = int(len(bc)//4)
    seqs = list(prs.keys())
    for i in range(len(seqs)):
        # no insert/delete
        seq = seqs[i]
        mappedseq, mappederr = get_errorbs_bytwoseqs(bc, seq, errbs)
        if len(mappedseq)>0:
            return mappedseq        
        if seq.find('R') >= 0 or seq.find('S') >= 0 or seq.find('I') >= 0 or seq.find('Y') >= 0 or seq.find('K') >= 0 or seq.find('M') >= 0:
            seqalias = get_alia_seqs(seq)
            for i in range(len(seqalias)):
                mappedseq, mappederr = get_errorbs_bytwoseqs(bc, seqalias[i], errbs)
                if len(mappedseq)>0:
                    return seq
            
    return ''

def get_dyn_align_score(bc, seq, match, mismatch, opengap, extgap):
#pre: the function to replace pairwise2.align.globalms(bc, seq, match, mismatch, opengap, extgap)[0][2]
#post: return alignment score
    #score = pairwise2.align.globalms(bc, seq, match, mismatch, opengap, extgap)[0][2]
    #alignments = aligner.align(bc, seq)
    maxscore = aligner.align(bc, seq)[0].score
    """
    maxscore = 0
    for alignment in sorted(alignments):
        print("Score = %.1f:" % alignment.score)
        print(alignment)
        if alignment.score > maxscore:
            maxscore = alignment.score
    if debug and maxscore > len(bc):
        print(" bc=%s\nseq=%s\noldscore=%f\nmaxscore=%f" % (bc, seq, score, maxscore) )
    """
    return maxscore-2, maxscore

def best_dynamic_score(bc, seq, match, mismatch, opengap, extgap):
    score = -1
    seqs = get_alia_seqs(seq)
    for count in range(len(seqs)):
        myseq=seqs[count]
        #myscore = pairwise2.align.globalms(bc, myseq, match, mismatch, opengap, extgap)[0][2]
        myscore, maxscore = get_dyn_align_score(bc, myseq, match, mismatch, opengap, extgap)
        #if debug and maxscore > len(bc):
        #    print(" bc=%s\nseq=%s\noldscore=%f\nmaxscore=%f" % (bc, seq, score, maxscore) )
        if myscore > len(bc):
            return myscore
        elif myscore > score:
            score = myscore
    return score

def save_to_pass_ambiguous_chain_map(bc, seq, score, p_args, sample, prs):
    if p_args.internal_test and score == len(bc)+1:
        fname = os.path.join(p_args.out_dir, 'pass_ambiguous_chain_map_' + sample + '.csv')
        f = open(fname, 'a')
        if not os.path.exists(fname) or os.stat(fname).st_size == 0:
            f.write("seq,chain\n")
        f.write(bc + ',' + seq + '(' + prs[seq] + ')\n')
        f.close()

def align_prs_bestscore(sample, prs, bc):
    match=2
    mismatch=-2
    opengap=-2
    extgap=-1
    if bc.find('GGGGGG') >= 0:
        return ''
    seqs = list(prs.keys())
    best_score = 0
    best_id = -1
    for i in range(len(seqs)):
        seq = seqs[i]
        score, maxscore = get_dyn_align_score(bc, seq, match, mismatch, opengap, extgap)
        if score > len(bc):
            save_to_pass_ambiguous_chain_map(bc, seq, score, p_args, sample, prs)
            return seq
        elif seq.find('R')>=0 or seq.find('S')>=0 or seq.find('I')>=0 or seq.find('Y')>=0 or seq.find('K')>=0 or seq.find('M')>=0:
            score = best_dynamic_score(bc, seq, match, mismatch, opengap, extgap)
            if score > len(bc):
                save_to_pass_ambiguous_chain_map(bc, seq, score, p_args, sample, prs)
                return seq
        if score > best_score:
            best_score = score
            best_id = i
    if p_args.internal_test and best_score == len(bc):
        fname = os.path.join(p_args.out_dir, 'fail_ambiguous_chain_map_' + sample + '.csv')
        f = open(fname, 'a')
        if not os.path.exists(fname) or os.stat(fname).st_size == 0:
            f.write("seq,chain\n")
        f.write(bc + ',' + seqs[best_id] + '(' + prs[seqs[best_id]] + ')\n')
        f.close()
    return ''

def get_partial_primer_start_end(p_args, prs):
    if p_args.receptor_kit == 'TCRv1':
        start, end = 18, 24
        if p_args.species == 'mouse':
            start, end = 7, 13
    elif p_args.receptor_kit == 'TCRv2':
        start, end = 5, 11
        if p_args.species == 'mouse':
            start, end = 7, 13
    elif p_args.receptor_kit == 'BCRv1':
        start, end = 16, 22
        if p_args.species == 'mouse':
            start, end = 16, 23
    elif p_args.receptor_kit == 'BCRv2' or p_args.receptor_kit == 'BCRv2sub' or p_args.receptor_kit == "BCRv2hy":
        start, end = 17, 24
        if p_args.species == 'human' and len(list(prs.keys())[0])>9:
            start, end = 17-3, 26
        if p_args.species == 'mouse':
            start, end = 16-3, 27
    
    if end-start != len(list(prs.keys())[0]): # mouse BCRv1
        #print("check and adjust the length of bc")
        if len(list(prs.keys())[0])>11 and p_args.receptor_kit == 'BCRv1' and p_args.species == 'mouse':
            start = start - 2
        end = start + len(list(prs.keys())[0])
    return start, end

def get_subseq(sample, read, repo_type, species, pseq_dict, prs, len_cutoff=30, mismatch=1):
    start, end = get_partial_primer_start_end(p_args, prs)
    if len(read[1]) >= len_cutoff:
        subseq = read[1][start: end]
    else:
        subseq = 'short'
    if mismatch == 0 or subseq in prs.keys() or subseq == 'short':
        pass
    else:
        #rawseq = ham_check(pseq_dict, subseq)
        rawseq = ham_check_flex(prs, subseq)
        if len(rawseq) == 0 and p_args.species == 'mouse' and p_args.receptor_kit.find('BCR')>=0:
            rawseq = align_prs_bestscore(sample, prs, subseq)
        if (len(rawseq) > 0):
            subseq = rawseq
    return subseq

def get_umi(read, umi_start=0, len_umi=12):
    if len(read[1]) < len_umi:
        umi = 'na'
    else:
        umi = read[1][umi_start : umi_start + len_umi] 
    return umi


def get_index(read):
    line = read[0].split()[0]
    segs = line.split(":")
    num = len(segs)
    return int(segs[num-3]), int(segs[num-2]), int(segs[num-1])


def process_read(r1, r2, bcs, bc_idx, umi, p_args, workflow2_link):
    if p_args.receptor_kit in ['TCRv2', 'BCRv2', 'BCRv2sub', 'BCRv2hy'] or (p_args.receptor_kit=="BCRv1" and p_args.species=="human"):
        if workflow2_link.find('ImmuneProfiler.jar') >= 0 and p_args.library_kit == 'takara_smartseq':
            bc = bcs[bc_idx].replace('"', '') 
            tmp = r1[0].split()  
            tmp[0] += '_bc_' + bc + '_umi_' + umi 
            r1[0] = ' '.join(tmp)
            
            tmp = r2[0].split()
            tmp[0] += '_bc_' + bc + '_umi_' + umi
            r2[0] = ' '.join(tmp)
            r2[1] = bc + r2[1]
            r2[3] = 'I' * len(bc) + r2[3]
        elif workflow2_link.find("ImmuneProfilerv2.jar") >= 0 and p_args.receptor_kit == 'BCRv2hy':
            tmp = r1[0].split()
            segs = tmp[0].split(':')
            if len(segs)>=7:
                segs[2] += '_umi_' + umi
            else:
                segs[2] += '_umi_' + umi
            tmp[0] = ':'.join(segs)
            r1[0] = ' '.join(tmp)
            
            tmp = r2[0].split()
            segs = tmp[0].split(':')
            if len(segs)>=7:
                segs[2] += '_umi_' + umi
            else:
                segs[2] += '_umi_' + umi
            tmp[0] = ':'.join(segs)
            r2[0] = ' '.join(tmp)
        elif workflow2_link.find("ImmuneProfilerv2.jar") >= 0 and p_args.library_kit == 'takara_smartseq':
            tmp = r1[0].split()  
            tmp[0] += '_umi_' + umi 
            r1[0] = ' '.join(tmp)
            
            tmp = r2[0].split()
            tmp[0] += '_umi_' + umi
            r2[0] = ' '.join(tmp)
    else:
        pass
    return


def correct_linker(read, start=12, linker='GTAC'):
    seq = read[1][start:start + len(linker)]
    if seq != linker:
        return True
    else:
        return False

def get_subigindex(lines, p_s, p_e, umi):
    pos = -1
    mid = int((p_s + p_e)/2)
    mid = (mid//4)*4
    
    if p_e > p_s :
        myumi = lines[mid+1][:12]
        if umi > myumi:
            pos = get_subigindex(lines, mid, p_e, umi)
        elif umi < myumi:
            pos = get_subigindex(lines, p_s, mid, umi)
        else: # umi == myumi
            pos = mid
    return pos

def get_umi_index(lines, pos):
    umi = lines[pos+1][:12]
    line = lines[pos].split()[0]
    segs = line.split(":")
    l = len(segs)
    return [umi, int(segs[l-3]), int(segs[l-2]), int(segs[l-1])]

def do_merge2(read2, r2):
    return read2

def do_merge1(read1, ref_subr1, pos):
    with open(ref_subr1, "r") as fr:
         lines = fr.readlines()
    r1 = lines[pos:pos+4]
    read1 = do_merge2(read1, r1)
    return read1

def merge_subIG_head(numl, fr1, fr2, sample, subseq, read1, read2, repo_type, species, pseq_dict, prs, len_cutoff=30, mismatch=1):
    infos = get_umi_index(read2, 0)
    lines2 = []
    lines1 = []
    for i in range(4):
         lines2.append(fr2.readline())
         lines1.append(fr1.readline())
    myinfos = get_umi_index(lines2, 0)
    #pos = get_subigindex(lines, 0, len(lines)-1, infos[0])
    while myinfos[0] < infos[0]:
        lines2.clear()
        lines1.clear()
        for i in range(4):
             lines2.append(fr2.readline())
             lines1.append(fr1.readline())
        myinfos = get_umi_index(lines2, 0)
    
    if myinfos[0] < infos[0]:
        print("exact same to do merging")
    elif myinfos[0] > infos:
        print("give to do merging")
    elif myinfos[0] == infos:
        print("exact same to do merging")
        do_merge2(read2, lines2)
        do_merge2(read1, lines1)
    else:
        print("get most similar r2")
    return read1, read2

def merge_subIG_mid(sample, subseq, read1, read2, repo_type, species, pseq_dict, prs, len_cutoff=30, mismatch=1 ):
    ref_subr1 = os.path.join(p_args.out_dir, "preprocess", sample+"_subIGG_R1_sorted.fastq")
    ref_subr2 = os.path.join(p_args.out_dir, "preprocess", sample+"_subIGG_R2_sorted.fastq")
    if subseq == "GGGAAGA":
        ref_subr1 = os.path.join(p_args.out_dir, "preprocess", sample+"_subIGA_R1_sorted.fastq")
        ref_subr2 = os.path.join(p_args.out_dir, "preprocess", sample+"_subIGA_R2_sorted.fastq")
    infos = get_umi_index(read2, 0)
    with open(ref_subr2, "r") as fr:
         lines = fr.readlines()
    pos = get_subigindex(lines, 0, len(lines)-1, infos[0])
    if pos >=0 and pos < len(lines):
        while get_umi_index(lines, pos) < infos:
            pos += 4
        while get_umi_index(lines, pos) > infos:
            pos -= 4
        
        if get_umi_index(lines, pos) < infos:
            pos += 4
        elif get_umi_index(lines, pos) > infos:
            pos -= 4
        elif get_umi_index(lines, pos) == infos:
            print("exact same to do merging")
        else:
            print("get most similar r2")
        do_merge2(read2, lines[pos:pos+4])
        do_merge1(read1, ref_subr1, pos)
    else:
        print("have no the corresponding the sub r2; give up")
    return read1, read2

def save_sortedfastq(umis, r1_fq, r2_fq, p_args):
    for key in range(2):
        fpath = r2_fq
        if key == 1:
            fpath = r1_fq
        fdname = os.path.dirname(fpath)
        fname = os.path.basename(fpath)
        fname_sorted = fname.replace(".fastq", "_sorted.fastq")
        with open(fpath, "r") as fr:
            lines = fr.readlines()
        #with open(os.path.join(p_args.fastq_dir, fname), 'w') as fw:
        with open(os.path.join(fdname, fname_sorted), 'w') as fw:
            ncols = len(umis[0])
            for i in range(len(umis)):
                index = umis[i][ncols-1]
                for j in range(4):
                    fw.write('%s\n' % lines[index*4+j].rstrip())  

def write_sample_qc(p_args, sd):
    file_dir = p_args.out_dir
    file_name = p_args.out_name
    with open(os.path.join(file_dir, 'report', file_name + '_sample_QC_stats.csv'), 'w', newline='') as qc_file:
        writer = csv.writer(qc_file, delimiter=',', lineterminator="\n")
        col_names = ['Sample_ID']
        key1 = list(sd)[0]
        for item in sd[key1]:
            col_names.append(item[0].split('_')[0])
            col_names.append(item[0].split('_')[0]+'%')
        writer.writerow(col_names)
        for key, val in sd.items():
            sample_id = key
            # stats = [ str(item[1]) + '(' + '{:.1%}'.format(item[2])+ ')' for item in val ]
            l1 = [str(item[1]) for item in val ]
            l2 = ['{:.1%}'.format(item[2]) for item in val]
            stats = list(chain.from_iterable(zip(l1, l2)))
            writer.writerow([sample_id] + stats)

def write_airr_sample_qc(p_args, sd):
    file_dir = p_args.out_dir
    file_name = p_args.out_name
    with open(os.path.join(file_dir, 'airr_report', file_name + '_sample_QC_stats.csv'), 'w', newline='') as qc_file:
        writer = csv.writer(qc_file, delimiter=',', lineterminator="\n")
        col_names = ['organism_id', 'sample_processing_id']
        key1 = list(sd)[0]
        for item in sd[key1]:
            col_names.append(item[0].split('_')[0])
            col_names.append(item[0].split('_')[0]+'%')
        writer.writerow(col_names)
        organism = p_args.species
        for key, val in sd.items():
            sample_id = key
            # stats = [ str(item[1]) + '(' + '{:.1%}'.format(item[2])+ ')' for item in val ]
            l1 = [str(item[1]) for item in val ]
            l2 = ['{:.1%}'.format(item[2]) for item in val]
            stats = list(chain.from_iterable(zip(l1, l2)))
            writer.writerow([organism] + [sample_id] + stats)

def create_fd(sd):
    sublist = ['short', 'undetermined', 'flc', 'total']
    if workflow2_link.find('ImmuneProfilerv2.jar')>=0 and p_args.internal_test:
        sublist = ['short', 'total']
    pd = OrderedDict((k, sd[k]) for k in sd.keys())
    for k, v in pd.items():
        l1 = []
        l2 = []
        for i in range(len(v)):
            if v[i][1] > 0:  # cutoff to discard chain category, used to be v[i][2] >= 0.05
                l1.append(i)
            if v[i][0] not in sublist:
            #if v[i][0] not in ['short', 'undetermined', 'flc', 'total']:
                l2.append(i)
        l = [value for value in l1 if value in l2]
        pd[k] = [v[idx] for idx in l]
    fd = OrderedDict()
    for k, v in pd.items():
        sample = k
        bc = ''
        for i in range(len(v)):
            chain = v[i][0].split('_')[0]
            if chain == 'undetermined' or chain == 'flc':
                bc = chain #bc = v[0][0].split('_')[1]
            else:
                bc = v[i][0].split('_')[1]
            name = sample + '_' + chain
            fd[name] = [ name + '_R1.fastq', name + '_R2.fastq', bc]
    return fd


def dos2unix(fq_file):
    """ Function: convert dos format file with \r\n to Unix with \n
    This function is writen for MiXCR, as MiXCR takes input fastqs
    with \n as EOL
    """
    try:
        text = open(fq_file).read()
        open(fq_file, 'w',newline='\n').write(text)
    except EnvironmentError as err:
        print('error', 'Unable to open file: ' + fq_file + '\n', flush=True)
    return


def prep_migec_input(p_args, fd):
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
    """Function: Run a system command.
    Stores stdout in outfile and stderr in errfile
    err_message is a message to be logged if the command fails
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
    with open(fq, 'r') as fp:
        num_reads = len(fp.readlines())//4
        return num_reads

def run_migec(p_args, migec_dir, log_name, errfile_name):
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
        print("copy file into folder assemble")
        for key, value in fd.items():
            src = os.path.join(p_args.out_dir, "preprocess", value[0])
            dst = os.path.join(p_args.out_dir, "run_migec", "assemble", value[0])
            if not os.path.isdir(os.path.join(p_args.out_dir, "run_migec", "assemble")):
                os.makedirs(os.path.join(p_args.out_dir, "run_migec", "assemble"))
            copyfile(src, dst)
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
    if os.path.isfile(os.path.join(p_args.out_dir, 'run_migec', 'histogram', 'estimates.txt')):
        infile = open(os.path.join(p_args.out_dir, 'run_migec', 'histogram', 'estimates.txt'), 'r')
        next(infile)  
        for line in infile:
            items = line.rstrip().split('\t')
            sample = items[0]
            total_reads = items[2]
            total_migs = items[3]
            overseq_threshold = items[4]
            collison_threshold = items[5]
            umi_qual_threshold = items[6]
            if sample in fd.keys():
                fd[sample].append(total_reads)
                fd[sample].append(total_migs)
                if p_args.umi_cutoff == '':
                    fd[sample].append(overseq_threshold)
                else:
                    fd[sample].append(str(p_args.umi_cutoff))
                fd[sample].append(collison_threshold)
                fd[sample].append(umi_qual_threshold)
        infile.close()


def find_file(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def getsubdirsandfiles(distfd):
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
    umi_cutoff = 1
    f = open(refinefile, 'r')
    lines = f.readlines()
    f.close()
    line = lines[len(lines)-2]
    if line.find("Effective threshold:") >= 0:
        umi_cutoff = str(int(float(line.split(": ")[1].split('\n')[0])))
    return umi_cutoff

def summ_results(analysis_dir, folder, inputDict, align_stats_dict, clone_stats_dict):
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
                if len(files)>0 :
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
        # need to rewrite this so that even no clonotype reported, 0 as stats could be reported
        tmp = find_file(key + '_' + folder + clonesall, full_path)
        if len(tmp) == 0:
            logger.info('No clonos_all file generated stop analysis \n')
            next
        else:
            for i in range(len(tmp)):
                clones_infile = tmp[i]
                #clones_outfile = clones_infile.split('clones')[0] + '_result.csv'
                clones_outfile = clones_infile[:len(clones_infile)-4] + '.csv' # output as csv format
                if p_args.cogentIP_version.find('v2')<0:
                    clones_outfile = clones_infile.replace('_clones_all.txt','_clones_result.csv')
                parse_clones_report(clones_infile, clones_outfile, key, inputDict)
    return inputDict


def report_stats(r_type, file_dir, file_name, metaDict):
    subheads = ['total MIG', 'UMI threshold', 'number of reads after MIG collapse']
    if list(metaDict.values())[0][5] == None or p_args.receptor_kit=='TCRv1' or (p_args.receptor_kit=='BCRv1' and p_args.species=='mouse'):
        subheads = ['aligned reads', 'without UMI', 'number of reads available in file']
    elif workflow2_link.find('ImmuneProfilerv2.jar')>=0:
        subheads = ['aligned reads', 'UMI threshold', 'number of reads available in file']
    with open(os.path.join(file_dir,file_name + '_mapping_stats.csv'), 'w', newline='') as stats_file:
        writer = csv.writer(stats_file, delimiter=',', lineterminator="\n")  # write as csv format
        heads = ['sample type','total read'] + subheads + ['aligned', 'pair-read overlap', 'overlapped and aligned']
        if r_type in ['BCRv1','BCRv2','BCRv2sub','BCRv2hy']:
            if p_args.internal_test:
                heads += ['reads per clonotype','reads in clonotypes','mapped low quality reads','pcr error correction','failed mapping',
                          'clonotype count','IgG','IgM','IgK','IgL','IgA','IgD','IgE','IgH(lack constant region)']
            else:
                heads += ['clonotype count','IgG','IgM','IgK','IgL','IgA','IgD','IgE','IgH(lack constant region)']
        #elif r_type == 'TCRv2':
        elif r_type in ['TCRv1', 'TCRv2']:
            if p_args.internal_test:
                heads += ['reads per clonotype','reads in clonotypes','mapped low quality reads','pcr error correction','failed mapping',
                          'clonotype count','TRA','TRB']
            else:
                heads += ['clonotype count','TRA','TRB']
        writer.writerow(heads)
        lheads = len(heads)
        for key, val in metaDict.items():
            sample_id = key
            stats2report = val[3:6] + [item.split(':')[1].split('(')[0].replace(' ', '') for item in val[8:]]
            #if len(stats2report) != 16:
            if len(stats2report) != lheads-1:
                for i in range(len(stats2report)):
                    #print( stats2report[i])
                    if stats2report[i] == None or not stats2report[i].isnumeric():
                        stats2report[i]='0'
                stats2report = stats2report + [0]*(lheads-1-len(stats2report))
            writer.writerow([sample_id] + stats2report)

def report_stats_old(r_type, file_dir, file_name, metaDict):
    subheads = ['total MIG', 'UMI threshold', 'number of reads after MIG collapse']
    if list(metaDict.values())[0][5] == None or p_args.receptor_kit=='TCRv1' or (p_args.receptor_kit=='BCRv1' and p_args.species=='mouse'):
        subheads = ['aligned reads', 'without UMI', 'number of reads available in file']
    elif workflow2_link.find('ImmuneProfilerv2.jar')>=0:
        subheads = ['aligned reads', 'UMI threshold', 'number of reads available in file']
    with open(os.path.join(file_dir,file_name + '_mapping_stats.csv'), 'w', newline='') as stats_file:
        writer = csv.writer(stats_file, delimiter=',')  # write as csv format
        heads = ['sample type','total read'] + subheads
        if r_type in ['BCRv1','BCRv2','BCRv2sub','BCRv2hy']:
            if p_args.internal_test:
                writer.writerow(['sample type','total read', 
                            subheads[0], subheads[1],subheads[2],
                            'aligned', 'pair-read overlap',
                            'overlapped and aligned',
                            'reads per clonotype','reads in clonotypes','mapped low quality reads','pcr error correction','failed mapping',
                            'clonotype count','IgG','IgM','IgK','IgL','IgA','IgD','IgE','IgH(lack constant region)'])
            else:
                writer.writerow(['sample type','total read', 
                             subheads[0], subheads[1],subheads[2],
                             'aligned', 'pair-read overlap',
                             'overlapped and aligned', 
                             'clonotype count','IgG','IgM','IgK','IgL','IgA','IgD','IgE','IgH(lack constant region)'])
        #elif r_type == 'TCRv2':
        elif r_type in ['TCRv1', 'TCRv2']:
            if p_args.internal_test:
                writer.writerow(['sample type','total read', 
                            subheads[0], subheads[1],subheads[2],
                            'aligned', 'pair-read overlap',
                            'overlapped and aligned',
                            'reads per clonotype','reads in clonotypes','mapped low quality reads','pcr error correction','failed mapping',
                            'clonotype count','TRA','TRB'])
            else:
                writer.writerow(['sample type','total read', 
                             subheads[0], subheads[1],subheads[2],
                             'aligned', 'pair-read overlap',
                             'overlapped and aligned',
                             'clonotype count','TRA','TRB'])
        for key, val in metaDict.items():
            sample_id = key
            stats2report = val[3:6] + [item.split(':')[1].split('(')[0].replace(' ', '') for item in val[8:]]
            if len(stats2report) != 16:
                stats2report = stats2report + [0]*(16-len(stats2report))
            writer.writerow([sample_id] + stats2report)

def stats_csv_to_airr( in_csv, airr_csv):
    rows = []
    with open(in_csv, 'rt', encoding='utf8') as f:
        reader = csv.reader(f, delimiter=',')
        count = 0
        for row in reader:
            if count == 0:
                row = ["organism_id", "sample_processing_id"] + airrhead(row)
            else:
                sid = row[0].split("_")[0]
                row = [p_args.species, sid] + row
            count += 1
            rows.append(row)
    save_csv(airr_csv, rows)

def save_csv(file_path, datarows):
    with open(file_path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',', lineterminator="\n")  # write as csv format
        for row in datarows:
            writer.writerow(row)

def airrhead(items):
    items = [sub.replace(' ', '_') for sub in items]
    items = [sub.lower() for sub in items]
    for i in range(len(items)):
        if (items[i] == "amino_acid_length"):
            items[i] = "junction_length_aa"
        elif (items[i] == "igh(lack_constant_region)"):
            items[i] = "IGH(lack_constant_region)"
        if (items[i] == "igh") or (items[i] == "igi") or (items[i] == "igk") or (items[i] == "igl"):
            items[i] = items[i].upper()
        elif items[i] == "iga" or items[i] == "igd" or items[i] == "ige" or items[i] == "igg" or items[i] == "igm":
            items[i] = items[i].upper()
        elif items[i] == "tra" or items[i] == "trb" or items[i] == "trd" or items[i] == "trg":
              items[i] = items[i].upper()
    return items


def merge_airr_csv(sample_id, process_type, report_name, file_dir, output_dir, airr_dir):
    mergedflist = []
    stats_file = find_file('*_mapping_stats.csv',file_dir)[0]
    with open(stats_file, 'rt', encoding='utf8') as f:
        reader = csv.reader(f, delimiter=',')
        file_prefix_list = []
        ct = 0
        for row in reader:
            if ct == 0:
                pass#print("keep header row")
            elif row[0].split('_')[0] == sample_id:
                file_prefix_list.append(row[0])
            else:
                pass
            ct += 1
    for key in file_prefix_list:
        sid = key.split('_')[0]
        if os.path.isdir(os.path.join(p_args.out_dir, 'airr_report', sid)) is False:
            os.makedirs(os.path.join(p_args.out_dir, 'airr_report', sid))
        else:
            pass
        #infile = find_file(key + '_' + process_type + '_clones_result.csv',output_dir)[0]
        infile = find_file(key + '_' + process_type + '_clones_*.csv',output_dir)[0]
        frowlist = []
        with open(infile, 'rt', encoding='utf8') as f:
            reader = csv.reader(f, delimiter=',')
            i = 0
            for row in reader:
                if i==0 :
                    row = ["organism_id", "sample_processing_id"] + airrhead(row)
                else:
                    row = [p_args.species, sid] + row
                frowlist.append(row)
                i += 1
        save_csv( os.path.join(p_args.out_dir, 'airr_report', sid, key + '_' + process_type + '_clones_result.csv'), frowlist)
        if len(mergedflist) == 0:
            mergedflist += frowlist
        else :
            mergedflist += frowlist[1:]
        #os.rename(infile, os.path.join(p_args.out_dir, 'report', sid, os.path.basename(infile)))
    #re-calculate the fraction based on sum all reads for each chain
    sum = 0
    for i in range(1,len(mergedflist)):
        sum += float(mergedflist[i][2])
    for i in range(1,len(mergedflist)):
        mergedflist[i][3] = float(mergedflist[i][2])/sum
    save_csv( os.path.join(p_args.out_dir, 'airr_report', p_args.out_name + '_' +sid + '_' + process_type + '_report.csv'), mergedflist)


def merge_csv(sample_id, process_type, report_name, file_dir, output_dir):
    wb = openpyxl.Workbook()
    ws = wb.active
    stats_file = find_file('*_mapping_stats.csv',file_dir)[0]
    ws.title = 'stats'
    with open(stats_file, 'rt', encoding='utf8') as f:
        reader = csv.reader(f, delimiter=',')
        file_prefix_list = []
        ct = 0
        for row in reader:
            if ct == 0:
                ws.append(row)  # keep header row
            elif row[0].split('_')[0] == sample_id:
                ws.append(row)
                file_prefix_list.append(row[0])
            else:
                pass
            ct += 1
    for key in file_prefix_list:
        sid = key.split('_')[0]
        if os.path.isdir(os.path.join(p_args.out_dir, 'report', sid)) is False:
            os.makedirs(os.path.join(p_args.out_dir, 'report', sid))
        else:
            pass
        pattern = '_clones_*.csv'
        if workflow2_link.find("ImmuneProfilerv2.jar") >= 0:
            pattern = '.clones_*.csv'
        infiles = find_file(key + '_' + process_type + pattern, file_dir)
        for i in range(len(infiles)):
            infile = infiles[i]  # add '_' for name seperation ****
            sheet_name = key + '_clone'
            if key.find("_TR") < 0 and key.find("_IG") < 0:
                sheet_name = key + infile[len(infile)-9:len(infile)-4] + '_clone'
            ws = wb.create_sheet(sheet_name)
            with open(infile, 'rt', encoding='utf8') as f:
                reader = csv.reader(f, delimiter=',')
                for row in reader:
                    ws.append(row)
            os.rename(infile, os.path.join(p_args.out_dir, 'report', sid, os.path.basename(infile)))
    wb.save(os.path.join(output_dir, report_name + '_report.xlsx'))
    return stats_file

def getsubIg_dict(input_dict, subIg):
    temp_dict = {}
    for key, val in input_dict.items():
        if key.find(subIg.upper()) >= 0:
            temp_dict[key] = val
    return temp_dict        

def findmatchedpos(seg, subs, fixpos):
    poslist=[]
    pos = seg.find(subs)
    while pos >=0:
        poslist.append(pos)
        pos = seg.find(subs, pos+1)
    if len(poslist)==0:
        l_subs = len(subs)
        temp = subs
        atcg = ["A", "T", "C", "G"]
        count = 0
        while count < l_subs:
            if count != fixpos:
                for i in range(4):
                    temp = subs[:count] + atcg[i] + subs[count+1:]
                    pos = seg.find(temp)
                    while pos>=0:
                        #break
                        if pos not in poslist:
                            poslist.append(pos)
                        pos = seg.find(temp, pos+1)
            count += 1
    poslist.sort()
    if len(poslist) == 0:
        return -1
    else: return poslist[len(poslist)-1]

def get_subIgA(seg, count): #, clonesall):
    print("This is subIgA")
    subIgA = 'IGA'
    return subIgA

def inbin(v0, v1, s, l):
    if (v0 >= 0 and v1 >= 0 and v1-v0 >= s and v1-v0 <= l):
        return True
    else: return False

def adjustPs(p0, p1, p2, p3, p4):
    if p0 >= 0 and p1 >= 0:
        if inbin(p1, p2, 44, 46) or inbin(p1, p3, 44, 46) or inbin(p1, p4, 44, 46):
            p0=-1
        elif inbin(p0, p2, 44, 46) or inbin(p0, p3, 44, 46) or inbin(p0, p4, 44, 46):
            p1 = -1
    elif p0 >= 0: #p1<0
        if not inbin(p0, p2, 44, 46) and not inbin(p0, p3, 44, 46) and not inbin(p0, p4, 44, 46):
            p0 = -1
    elif p1 >= 0: #p0<0
        if not inbin(p1, p2, 44, 46) and not inbin(p1, p3, 44, 46) and not inbin(p1, p4, 44, 46):
            p1 = -1
    #else: #p0<0 p1<0
    #    print("do nothing")       
    return p0, p1, p2, p3, p4

def get_subIgG(seg, count):
    subIgG = 'IGG'
    labs = ['GCTTCCAC','GCCTCCAC','AGGAGCACCTCT','AAGAGCACCTCT','AGGAGCACCTCC']
    p0 = findmatchedpos(seg, labs[0], 2) #seg.find(labs[0])
    p1 = findmatchedpos(seg, labs[1], 2) #seg.find(labs[1])
    p2 = findmatchedpos(seg, labs[2], 1) #seg.find(labs[2])
    p3 = findmatchedpos(seg, labs[3], 1) #seg.find(labs[3])
    p4 = findmatchedpos(seg, labs[4], 11) #seg.find(labs[4])
    #if p0==-1 and p1>=0 and p2==-1 and p3>=0 and p4==-1:
    #    print("p4>0")
    p0, p1, p2, p3, p4 = adjustPs(p0, p1, p2, p3, p4)
    if (p0>=0):
        subIgG = 'IGG3/4'
    elif (p1>=0):
        subIgG = 'IGG1/2'
    if (p0>=0 or p1>=0) and (p2>=0 or p3>=0 or p4>=0):
        if inbin(p0, p4, 44, 46):
            subIgG = 'IGG4'
        elif inbin(p1, p4, 44, 46):
            subIgG = 'IGG2'
        elif inbin(p0, p3, 44, 46) or inbin(p1, p3, 44, 46):
            subIgG = 'IGG1'
        elif inbin(p0, p2, 44, 46) or inbin(p1, p2, 44, 46):
            subIgG = 'IGG3'
        else:
            subIgG = 'IGG' #"undet"
    return subIgG

def analyzeallCAlignments(p_args, segs):
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
    for i in range(1,len(subIgAtype)):
        if subIgAvalue[i] == subIgAvalue[i-1]:
            subIgA += ','+subIgAtype[i]
    return subIgA

def get_subIg(p_args, aligns, key, read_type, target_region, mixcr_dir, clonesall):
    process_type = str(read_type) + '_' + str(target_region)
    subIgfile = os.path.join(p_args.out_dir, mixcr_dir, process_type,
                             key + '_' + process_type + '_subIg' + key[-1] + '.tsv')
    subIG = 'subIGG'
    if key.find('IGA') >= 0:
        subIG = 'subIGA'
    if not os.path.exists(aligns) and p_args.cogentIP_version == 'v2':
        aligns = aligns.replace('_alignments.tsv', '_readIds.tsv')
    fr = open(aligns, 'r')
    line = fr.readline()
    line = line.replace('targetSequences', subIG)
    line = line.replace('targetQualities', 'targetSequences')
    fw = open(subIgfile, 'w')
    fw.write(line)
    count = 0
    while True:
        line = fr.readline()
        if ("" == line):
            break;
        segs = line.split('\t')
        if key.find('IGG') >= 0:
            subIG = get_subIgG(segs[0], count)
        elif key.find('IGA') >= 0:
            #subIG = get_subIgA(segs[0], count)
            subIG = analyzeallCAlignments(p_args, segs)
        line = subIG + '\t' + segs[0]
        for i in range(2, len(segs)):
            line += '\t' + segs[i]
        fw.write(line)
        count += 1
    fr.close()
    fw.close()

def get_rightindex(miglines, mig_pos, asslines, ass_pos):
    print("False:")
    print("%s" % miglines[mig_pos+1][:30])
    print("%s" % asslines[ass_pos][:30])
    print("mig_pos=", mig_pos)
    print("ass_pos=", ass_pos)
    if asslines[ass_pos][:8] != "GTACGGGG":
        return mig_pos+4, ass_pos+1
    mig_p = mig_pos
    ass_p = ass_pos
    while miglines[mig_p+1][:30] != asslines[ass_p][:30] and ass_p < ass_pos+30 and ass_p < len(asslines):
        ass_p += 1
    if ass_p < ass_pos+30 and ass_p < len(asslines):
        ass_pos = ass_p
    else:
        while miglines[mig_p+1][:30] != asslines[ass_p][:30] and mig_p < mig_pos+30 and mig_p < len(miglines):
            mig_p += 4
        if mig_p < mig_pos+30 and mig_p < len(miglines):
            mig_pos = mig_p
        else:
            ass_pos += 1
            mig_pos += 4
    print("new mig_pos=", mig_pos)
    print("new ass_pos=", ass_pos)
    return mig_pos, ass_pos

def scanindex(mig1line, mig2line, asslines, p_thread):
    ass_pos =-1
    m2line = mig2line.split('\n')[0]
    m2line = m2line[len(m2line)-16:]
    m1line = mig1line.split('\n')[0]
    m1line = m1line[:16]
    for i in range(1, len(asslines)):
        assline = asslines[i].split('\t')[0]
        pe=findmatchedpos(m2line, assline[len(assline)-16:], -1)
        ps=findmatchedpos(m1line, assline[:16], -1)
        if pe >= 0 and ps >= 0:
            ass_pos = i
            #if ass_pos > p_thread:
            break
    return ass_pos

def get_umis_by_readIds_oldright(p_args, readIds_fpath, sample, subig):
#pre: only for BCRv2hy
    readIdslines = []
    umis = []
    if not os.path.exists(readIds_fpath) :
        return umis
    with open(readIds_fpath, 'r') as readIds_fr:
        for readIdsline in readIds_fr:
            segs = readIdsline.split('\t')
            readIdslines.append(segs[len(segs)-1])
    prefd = os.path.join(p_args.out_dir, 'run_migec', 'assemble')
    prelab = '_R2.t*'
    if p_args.cogentIP_version == 'v2':
        prefd = os.path.join(p_args.out_dir, 'preprocess')
        prelab = '_R1.fastq'
    ig_fpath = find_file(sample + '_' + subig + prelab, prefd)[0]
    with open(ig_fpath, 'r') as ig_fr:
        igfqlines = ig_fr.readlines()
    for i in range(1, len(readIdslines)):
        index_fq = int( readIdslines[i].rstrip() )
        #if index != i-1: print("stop here")
        if p_args.cogentIP_version == 'v2':
           umi_fq = igfqlines[4*index_fq].split(':')[2].split('_umi_')[1]
        else:
           umi_fq = igfqlines[4*index_fq].split(':')[1]
        umis.append(umi_fq)
    return umis

def get_umis_by_readIds(p_args, readIds_fpath, sample, subig):
#pre: only for BCRv2hy
    readIdslines = []
    umis = []
    if not os.path.exists(readIds_fpath) :
        return umis
    with open(readIds_fpath, 'r') as readIds_fr:
        for readIdsline in readIds_fr:
            segs = readIdsline.split('\t')
            readIdslines.append(segs[len(segs)-1])
    prefd = os.path.join(p_args.out_dir, 'run_migec', 'assemble')
    prelab = '_R2.t*'
    if p_args.cogentIP_version == 'v2':
        prefd = os.path.join(p_args.out_dir, 'preprocess')
        prelab = '_R1.fastq'
    ig_fpath = find_file(sample + '_' + subig + prelab, prefd)[0]
    count_fqline = -1
    igfqline = ''
    ig_fr = open(ig_fpath, 'r')
    for i in range(1, len(readIdslines)):
        index_fq = int( readIdslines[i].rstrip() )
        while count_fqline < 4*index_fq:
            igfqline = ig_fr.readline()
            count_fqline += 1
        if p_args.cogentIP_version == 'v2':
           umi_fq = igfqline.split(':')[2].split('_umi_')[1]
        else:
           umi_fq = igfqline.split(':')[1]
        umis.append(umi_fq)
    ig_fr.close()
    return umis

def get_umis_seqs_by_readIds(p_args, readIds_fpath, sample, subig):
#pre: only for BCRv2hy
    readIdslines = []
    assemseqs = []
    umis = []
    if not os.path.exists(readIds_fpath) :
        return umis, assemseqs 
    with open(readIds_fpath, 'r') as readIds_fr:
        for readIdsline in readIds_fr:
            segs = readIdsline.split('\t')
            if segs[0] != 'targetSequences':
                assemseqs.append(segs[0])
                readIdslines.append(segs[len(segs)-1])
    prefd = os.path.join(p_args.out_dir, 'run_migec', 'assemble')
    prelab = '_R2.t*'
    if p_args.cogentIP_version == 'v2':
        prefd = os.path.join(p_args.out_dir, 'preprocess')
        prelab = '_R1.fastq'
    ig_fpath = find_file(sample + '_' + subig + prelab, prefd)[0]
    
    count_fqline = -1
    igfqline = ''
    ig_fr = open(ig_fpath, 'r')
    # tiny file 1.readlines() whole lines. 2. based on the readIdslines each value to extract umi
    for i in range(len(readIdslines)):
        index_fq = int( readIdslines[i].rstrip() )
        while count_fqline < 4*index_fq:
            igfqline = ig_fr.readline()
            count_fqline += 1
        if p_args.cogentIP_version == 'v2':
           umi_fq = igfqline.split(':')[2].split('_umi_')[1]
        else:
           umi_fq = igfqline.split(':')[1]
        umis.append(umi_fq)
    
    ig_fr.close()
    return umis, assemseqs

def alignmentsPrettySeq(assemPretty_f, index, migIgAGr1_f, migIgAGr2_f, count):
    seq = ''
    return seq


def replace_mig_subIG_R1_by_assembledseq(p_args, sample, subig, process_type, mixcr_dir): #, sub_dir  #rep_IG_R2_tx_cf_by_subIG_R2( IG_R2, subIG_R2)
#pre: mixcr4+ analysis result
    mydir = os.path.join(p_args.out_dir, mixcr_dir, process_type)
    readIds_fpath = find_file(sample + '_' + subig + '_' + process_type + '_readIds.tsv', mydir)[0]
    assem_fpath = find_file(sample + '_' + subig + '_' + process_type + '_alignments.tsv', mydir)[0]
    assemPretty_fpath = find_file(sample + '_' + subig + '_' + process_type + '_alignmentsPretty.tsv', mydir)[0]
    umis_readids = get_umis_by_readIds(p_args, readIds_fpath, sample, subig)
    umis_unique = list(set(umis_readids))
    umis_unique.sort()
    umis_fpath = os.path.join(p_args.out_dir, mixcr_dir, process_type, sample + "_" + subig + "_"  + process_type + '_umis.txt')
    if not os.path.exists( umis_fpath):
        save_list_to_file( umis_fpath, umis_readids, True)
    with open(assem_fpath, 'r') as assem_fr:
        assemlines = assem_fr.readlines()
    mydir = os.path.join(p_args.out_dir, 'run_migec', 'assemble')
    migsubr1_fpath = find_file(sample + '_sub' + subig + '_R1*.fastq', mydir)[0]
    migsubr2_fpath = find_file(sample + '_sub' + subig + '_R2*.fastq', mydir)[0]
    #migIgAGr1_fpath = find_file(sample + '_' + subig + '_R1*.fastq', mydir)[0]
    #migIgAGr2_fpath = find_file(sample + '_' + subig + '_R2*.fastq', mydir)[0]
    with open(migsubr1_fpath, 'r') as migsub_fr1:
        migsub1lines = migsub_fr1.readlines()
    with open(migsubr2_fpath, 'r') as migsub_fr2:
        migsub2lines = migsub_fr2.readlines()
    migsub_fr1.close()
    migsub_fr2.close()
    subr1 = []
    subr2 = []
    for count in range(0, len(migsub1lines), 4): #= mig_pos = 0
        umi = migsub1lines[count].split(":")[1]
        if umi in umis_readids:
            index = umis_readids.index(umi)
            #readIds(umis) has the same lines with the mapped alignments
            seq = assemlines[index+1].split('\t')[0]
            #if seq.find(',') >= 0:
            #    print('check alignmentsPretty.txt with umis_readIds=%d' % index)
            #    alignmentsPrettySeq(assemPretty_fpath, index, migsub1lines, migsub2lines, count)
            if seq.find(',') < 0:
                subr1.append( migsub1lines[count])
                subr1.append( seq )
                subr1.append( migsub1lines[count+2])
                subr1.append( 'I' * len(seq))
                subr2.append( migsub2lines[count])
                subr2.append( migsub2lines[count+1])
                subr2.append( migsub2lines[count+2])
                subr2.append( migsub2lines[count+3])
    print('new migec sub R1 done is ready to save')
    save_list_to_file(migsubr1_fpath, subr1, True)
    save_list_to_file(migsubr2_fpath, subr2, True)

def replace_pre_subIG_R2_by_assembledseq(p_args, sample, subig, process_type, mixcr_dir): #, sub_dir  #rep_IG_R2_tx_cf_by_subIG_R2( IG_R2, subIG_R2)
    mydir = os.path.join(p_args.out_dir, mixcr_dir, process_type)
    readIds_fpath = find_file(sample + '_' + subig + '_' + process_type + '_readIds.tsv', mydir)[0]
    readIds_subfpath = readIds_fpath.replace(sample+'_'+subig, sample+'_sub'+subig)
    assem_fpath = '' + readIds_fpath
    assemPretty_fpath = find_file(sample + '_' + subig + '_' + process_type + '_alignmentsPretty.tsv', mydir)[0]
    
    umis_igreadids, assemlines = get_umis_seqs_by_readIds(p_args, readIds_fpath, sample, subig)
    
    # read index from subIGA_prePro_cdr3_readIds.tsv generated by mixcr to map back high quality reads in _subIGA/G_R1/2.fastq
    subreadids = []
    if os.path.exists(readIds_subfpath):
        subreadIds_fr = open(readIds_subfpath, 'r')
        subreadIdsline = subreadIds_fr.readline()
        while subreadIdsline != '':
            subreadIdsline = subreadIds_fr.readline().rstrip()
            if subreadIdsline != '':
                segs = subreadIdsline.split('\t')
                subreadids.append(int(segs[len(segs)-1]))
        subreadIds_fr.close()
    
    mydir = os.path.join(p_args.out_dir, 'preprocess')
    migsubr1_fpath = find_file(sample + '_sub' + subig + '_R1*.fastq', mydir)[0]
    migsubr2_fpath = find_file(sample + '_sub' + subig + '_R2*.fastq', mydir)[0]
    migsubr1_f =  migsubr1_fpath.replace('_R1.fastq', '_matched_R1.fastq')
    migsubr2_f =  migsubr2_fpath.replace('_R2.fastq', '_matched_R2.fastq')
    subr1 = []
    subr2 = []
    count = 0
    
    if not os.path.exists(migsubr1_fpath) or os.stat(migsubr1_fpath).st_size == 0:
        print("The file does not exist or size 0, the replacement can not do it")
        return
    with open(migsubr1_fpath, 'r') as migsub_fr1, open(migsubr2_fpath, 'r') as migsub_fr2, ExitStack() as cm:
        fhs1 = cm.enter_context(open(migsubr1_f, 'a'))
        fhs2 = cm.enter_context(open(migsubr2_f, 'a'))
        for x, y in zip(migsub_fr1, migsub_fr2):
            subr1.append(x.rstrip())
            subr2.append(y.rstrip())
            count += 1
            if len(subr1) == 4:
                umi_subig = subr1[0].split(':')[2].split('_umi_')[1] #umi = migsub1lines[count].split(":")[1]
                if count//4-1 in subreadids and umi_subig in umis_igreadids:
                    index = umis_igreadids.index(umi_subig)
                    #readIds(umis) has the same lines with the mapped alignments
                    print('umi:%s' % umi_subig)
                    print('index in mixcr umis_readIds:', index)
                    print('count in migec _subIgG/A fq:', count)
                    print('subIgG/A_R1.t.cf:%s' % subr1[1][:30])
                    print('assem_alignments:%s' % assemlines[index+1][:30])
                    seq = assemlines[index+1].split('\t')[0]
                    if seq.find(',') < 0:
                        subr1[1] = seq 
                        subr1[3] = ('I' * len(seq))
                        for item in subr1:
                            fhs1.write('%s\n' % item)
                        for item in subr2:
                            fhs2.write('%s\n' % item)
                del subr1[:]
                del subr2[:]
    print('new sub R1/2 done is ready to rename')
    migsubr1_ftemp = migsubr1_fpath.replace('_R1.fastq','_temp_R1.fastq')
    migsubr2_ftemp = migsubr1_fpath.replace('_R2.fastq','_temp_R2.fastq')
    os.rename(migsubr1_fpath, migsubr1_ftemp)
    os.rename(migsubr2_fpath, migsubr2_ftemp)
    os.rename(migsubr1_f, migsubr1_fpath)
    os.rename(migsubr2_f, migsubr2_fpath)

def repone_subIG_Rn_tx_cf_by_alignassem(Rn, p_args, sample, subig, process_type, mixcr_dir):   #rep_IG_R2_tx_cf_by_subIG_R2( IG_R2, subIG_R2)
    #Rn = "subR2"
    assem_f = "" #subig = 'IGG'
    readIds_f = ""
    migsubr1_f = ""
    migsubr2_f = ""
    for mixf in os.listdir(os.path.join(p_args.out_dir, mixcr_dir, process_type)):
        if mixf.endswith(sample + '_' + subig + '_' + process_type + '_alignments.txt'):
            assem_f = mixf
        elif mixf.endswith(sample + '_' + subig + '_' + process_type + '_readIds.txt'):
            readIds_f = mixf
    
    umis = get_umis_by_readIds(p_args, mixcr_dir, process_type, readIds_f, sample, subig)
    print("umi done to save")
    umi_f = sample + "_" + subig + "_"  + process_type + '_umis.txt' #sample
    if not os.path.exists(os.path.join(p_args.out_dir, mixcr_dir, process_type, umi_f)):
        save_list_to_file(os.path.join(p_args.out_dir, mixcr_dir, process_type, umi_f), umis, True)
            
    with open(os.path.join(p_args.out_dir, mixcr_dir, process_type, assem_f), 'r') as assem_fr:
        asslines = assem_fr.readlines()
    #new_migr1_f = sample + 'merged_' + subig + '_R1.t1.cf.fastq' #sample
    sample = assem_f.split('_')[0]
    nonsubig = subig.replace("sub", "")
    for migsubfs in os.listdir(os.path.join(p_args.out_dir, 'run_migec', 'assemble')):
        if migsubfs.find(sample + "_sub" + nonsubig + "_R1") >= 0:
            migsubr1_f = migsubfs
        elif migsubfs.find(sample + "_sub" + nonsubig + "_R2") >= 0:
            migsubr2_f = migsubfs
            #break
    with open(os.path.join(p_args.out_dir, 'run_migec', 'assemble', migsubr1_f), 'r') as migsub_fr1:
        migsub1lines = migsub_fr1.readlines()
    with open(os.path.join(p_args.out_dir, 'run_migec', 'assemble', migsubr2_f), 'r') as migsub_fr2:
        migsub2lines = migsub_fr2.readlines()
    migsub_fr1.close()
    migsub_fr2.close()
    subr1 = []
    subr2 = []
    for count in range(0, len(migsub1lines), 4): #= mig_pos = 0
        umi = migsub1lines[count].split(":")[1]
        if umi in umis:
            index = umis.index(umi)
            #if count == 23072:
            #    print("stop here")
            print("count:", count)
            print("%s" % migsub1lines[count+1][:30])
            print("%s" % asslines[index+1][:30])
            seq = asslines[index+1].split('\t')[0]
            if seq.find(",") < 0:
                #migsub1lines[count+1] = asslines[index+1].split('\t')[0]
                #migsub1lines[count+3] = 'I' * len(migsub1lines[count+1])
                if subig.find("sub") < 0:
                    subr1.append( migsub1lines[count])
                    subr1.append( seq )
                    subr1.append( migsub1lines[count+2])
                    subr1.append( 'I' * len(seq))
                    for i in range(4):
                        subr2.append( migsub2lines[count+i])
                else: #seq.find(",") < 0:
                    subr2.append( migsub2lines[count])
                    subr2.append( seq )
                    subr2.append( migsub2lines[count+2])
                    subr2.append( 'I' * len(seq))
                    for i in range(4):
                        subr1.append( migsub1lines[count+i])                    
    print("new migec sub R1 done is ready to save")
    new_migsubr1_f =  migsubr1_f
    new_migsubr2_f =  migsubr2_f
    #new_migsubr1_f =  migsubr1_f.replace(sample, sample+"merged")
    save_list_to_file(os.path.join(p_args.out_dir, "run_migec", "assemble", new_migsubr1_f), subr1, True)
    save_list_to_file(os.path.join(p_args.out_dir, "run_migec", "assemble", new_migsubr2_f), subr2, True)

def replace_mig_ig2(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    process_type = str(read_type) + '_' + str(target_region)
    subig ="IGG"
    for key, value in input_dict.items():
        if key.find("_subIGG") >= 0:
            subig ="subIGG"
        elif key.find("_subIGA") >= 0:
            subig ="subIGA"
        sample = key[:len(key)-7]
        repone_subIG_Rn_tx_cf_by_alignassem("subR2", p_args, sample, subig, process_type, mixcr_dir)

def replace_subig_by_assembledseq(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    logger.info('Launching function replace_subig_by_assembledseq to execute assembled seq to replace _subIg?_R1/2 \n')
    process_type = str(read_type) + '_' + str(target_region)
    subig ="IGG"
    for key, value in input_dict.items():
        if key.find("_IGG") >= 0 or key.find("_IGA") >= 0 :
            subig = key[-3:]
            sample = key[:len(key)-4]
            if p_args.cogentIP_version == 'cenegtv1.5':
                replace_mig_subIG_R1_by_assembledseq(p_args, sample, subig, process_type, mixcr_dir) #, sub_dir
            else:
                replace_pre_subIG_R2_by_assembledseq(p_args, sample, subig, process_type, mixcr_dir)

def save_list_to_file(filename, mylist, overwrite):
    if overwrite or not os.path.exists(filename):
        with open( filename, 'w') as fr:
            for i in range(len(mylist)):
                item = mylist[i].split('\n')[0]
                fr.write("%s\n" % item)
        if platform.system() == 'Windows':
            dos2unix( filename)

def getumidict(reads_fq):
    umidect = {}
    l_read = len(reads_fq)
    for i in range(0, l_read, 4):
        umi = reads_fq[i].split(":")[1]
        if umi in umidect.keys():
            print("error")
        umidect[umi] = i
    return umidect

def reverse_complement(seq):
    l_seq = len(seq)
    newseq = ''
    for i in range(l_seq):
        ch = seq[l_seq-1-i]
        if ch == 'A':
            ch = 'T'
        elif ch == 'T':
            ch = 'A'
        elif ch == 'C':
            ch = 'G'
        elif ch == 'G':
            ch = 'C'
        newseq += ch
    return newseq

def match_R2umi_non_sub(fnonr2, fsubr2, p_args, sub_dir, assemble_folder):
    uminonr2 = []
    umisubr2 = []
    with open(fnonr2, 'r') as fnonr2_fr:
        fnonr2lines = fnonr2_fr.readlines()
    for i in range(0, len(fnonr2lines), 4):
        umi = fnonr2lines[i].split(":")[1]
        uminonr2.append(umi)
    with open(fsubr2, 'r') as fsubr2_fr:
        fsubr2lines = fsubr2_fr.readlines()
    for i in range(0, len(fsubr2lines), 4):
        umi = fsubr2lines[i].split(":")[1]
        umisubr2.append(umi)
    non_r2lines = []
    sub_r2lines = []
    for i in range(len(umisubr2)):
        umi = umisubr2[i]
        if umi in uminonr2:
            index = uminonr2.index(umi)
            for j in range(4):
                if j == 0:
                    line = fsubr2lines[4*i+j]
                    non_r2lines.append(line.replace("R2","R1"))
                else:
                    non_r2lines.append(fnonr2lines[4*index+j])
                sub_r2lines.append(fsubr2lines[4*i+j])
        else: #umi not in uminonr2:
            for j in range(4):
                if j == 0:
                    line = fsubr2lines[4*i+j]
                    non_r2lines.append(line.replace("R2","R1"))
                else:
                    non_r2lines.append(fsubr2lines[4*i+j])
                sub_r2lines.append(fsubr2lines[4*i+j])     
    segs1 = fnonr2.split("_R2")
    #segs2 = fpath2.split("_R2")
    non_fq2 = segs1[0] + "_match_R1" + segs1[1] #fpath1
    sub_fq2 = fsubr2 #segs2[0] + "_match_R2" + segs2[1] #fpath2
    save_list_to_file(non_fq2, non_r2lines, True)
    save_list_to_file(sub_fq2, sub_r2lines, True)
    return non_fq2, sub_fq2

def run_mixcr_sub2(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    # first create the sample_IGG_R1.t1.cf.fastq
    out_log = os.path.join(p_args.out_dir, mixcr_dir, log_name)
    run_mixcr_err = os.path.join(p_args.out_dir, mixcr_dir, errfile_name)
    process_type = str(read_type) + '_' + str(target_region)
    create_dir(os.path.join(p_args.out_dir, mixcr_dir, process_type))
    sub_dict = OrderedDict()
    for key, value in input_dict.items():
        if key.find("_subIGG") >= 0 or key.find("_subIGA") >= 0:
            sub_dict[key] = value
    for key, value in input_dict.items():
        if key.find("_IGG") >= 0 or key.find("_IGA") >= 0:
            fq1 = find_file(key + '_R2*.fastq', os.path.join(p_args.out_dir, sub_dir, 'assemble'))[0]
            fq2 = find_file(key.replace("_","_sub") + '_R2*.fastq', os.path.join(p_args.out_dir, sub_dir, 'assemble'))[0]
            match_fq1, match_fq2 = match_R2umi_non_sub(fq1, fq2, p_args, sub_dir, 'assemble')
            subkey = key[:len(key)-3]+"sub"+key[len(key)-3:]
            sub_dict[subkey][0] = os.path.basename(match_fq1) #value[1]
            sub_dict[subkey][1] = os.path.basename(match_fq2) #value[1]    
    run_mixcr_base(2, p_args, sub_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)
    replace_mig_ig2(p_args, sub_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)

def matched_subIg_fsize(key, fd):
    size = 0
    segs = key.split('_')
    key = segs[0] + '_sub' + segs[1]
    fnames = find_file(key + '*.fastq', fd)
    if len(fnames) > 0:
        size = os.stat(fnames[0]).st_size
    return size

def sort_fastqs_in_list(sub_dict, p_args, sub_dir):
    for key, value in sub_dict.items():
        forward_R1fq = os.path.join(p_args.out_dir, sub_dir, value[0])
        reverse_R2fq = os.path.join(p_args.out_dir, sub_dir, value[1])
        forward_R1out = os.path.join(p_args.out_dir, sub_dir, value[0].replace('R1.fastq', 'sorted_R1.fastq'))
        reverse_R2out = os.path.join(p_args.out_dir, sub_dir, value[1].replace('R2.fastq', 'sorted_R2.fastq'))
        
        forward_file = open(forward_R1fq, "r")
        reverse_file = open(reverse_R2fq, "r")
        
        forward_out = open(forward_R1out, "w")
        reverse_out = open(reverse_R2out, "w")
        
        forward_records = sorted(SeqIO.parse(forward_file,"fastq"), key=lambda x: x.id)
        reverse_records = sorted(SeqIO.parse(reverse_file,"fastq"), key=lambda x: x.id)
        
        SeqIO.write(forward_records, forward_out, "fastq")
        SeqIO.write(reverse_records, reverse_out, "fastq")
        
        forward_file.close()
        reverse_file.close()
        forward_out.close()
        reverse_out.close()
        if platform.system() == 'Windows':
            dos2unix(forward_R1out)
            dos2unix(reverse_R2out)

def sort_fastqs_in_list_new(sub_dict, p_args, sub_dir):
    for key, value in sub_dict.items():
        readsR1f = os.path.join(p_args.out_dir, sub_dir, value[0])
        readsR2f = os.path.join(p_args.out_dir, sub_dir, value[1])
        
        forward_reads = SeqIO.parse(readsR1f,"fastq")
        reverse_reads = SeqIO.parse(readsR2f,"fastq")
                
        sorted_forward_reads = sorted(forward_reads, key=lambda read: read.id)
        sorted_reverse_reads = sorted(reverse_reads, key=lambda read: read.id)
        
        #outR1f = readsR1f.replace('R1.fastq', 'sorted_R1.fastq')
        #outR2f = readsR2f.replace('R2.fastq', 'sorted_R2.fastq')
        
        with open(readsR1f, "w") as forward_outfile:
            SeqIO.write(sorted_forward_reads, forward_outfile, "fastq")
        with open(readsR2f, "w") as reverse_outfile:
            SeqIO.write(sorted_reverse_reads, reverse_outfile, "fastq")
        
        if platform.system() == 'Windows':
            dos2unix(readsR1f)
            dos2unix(readsR2f)

def run_mixcr_preIgGA(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    # first create the sample_IGG_R1.t1.cf.fastq
    sub_dict = OrderedDict()
    for key, value in input_dict.items():
        if key.find('_IGG') >= 0 or key.find('_IGA') >= 0:
            sub_dict[key] = value.copy()
    if workflow2_link.find('ImmuneProfilerv2.jar') >= 0:
        for key, value in input_dict.items():
            if key.find('_subIGG') >= 0 or key.find('_subIGA') >= 0:
                sub_dict[key] = value.copy()
    sort_fastqs_in_list_new(sub_dict, p_args, sub_dir)
    run_mixcr(p_args, sub_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)
    replace_subig_by_assembledseq(p_args, sub_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)
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


#def run_mixcr_base(lab, p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
def run_mixcr_base(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    out_log = os.path.join(p_args.out_dir, mixcr_dir, log_name)
    run_mixcr_err = os.path.join(p_args.out_dir, mixcr_dir, errfile_name)
    process_type = str(read_type) + '_' + str(target_region)
    tmp_dict = copy.deepcopy(input_dict)
    is_initUnd = False
    preset = 'takara'
    if p_args.species == 'human':
        preset += '-human'
    elif p_args.species == 'mouse':
        preset += '-mouse'
    if p_args.receptor_kit == 'TCRv1':
        if  p_args.species == 'human': 
            preset += '-tcr-V1'
        elif p_args.species == 'mouse':
            preset += '-tcr'
    elif p_args.receptor_kit == 'TCRv2':
        if  p_args.species == 'human': 
            preset += '-tcr-V2'
        elif p_args.species == 'mouse':
            preset += '-tcr-V2'
    elif p_args.receptor_kit.find( 'BCR') >= 0:
        preset += '-bcr'
        if p_args.species == 'mouse' and p_args.receptor_kit == 'BCRv2':
            preset += '-V2'
    if target_region.find('cdr3') >= 0:
        preset += '-cdr3'
    else:
        preset += '-full-length'
    samps_clns = os.path.join(p_args.out_dir, mixcr_dir, process_type, '*.clns')
    alignQc_pdf = os.path.join(p_args.out_dir, 'report', 'alignQc.pdf')
    chainUsage_pdf = os.path.join(p_args.out_dir, 'report', 'chainUsage.pdf')
    createPDF = False
    for key, val in input_dict.items():
        if target_region == 'cdr3':
            reg_str = 'CDR3'
        elif target_region == 'fl':
            reg_str = 'Full_length'
        print(get_time() + ' [INFO] Processing ' + reg_str + ' region of sample ' + key + '...', flush=True)
        readslayout = 'Opposite'
        fq1 = ''
        fq2 = ''
        #if val[0].find("_TR") < 0 and val[0].find("_IG") < 0 and val[0].find("_undetermined") < 0 : # read_type == "rawUnd"
        if sub_dir == 'rawUnd' or read_type == 'rawUnd':
            fq1 = os.path.join(p_args.fastq_dir, val[0])
            fq2 = os.path.join(p_args.fastq_dir, val[1])
        #if p_args.receptor_kit != 'TCRv1' and (p_args.receptor_kit != "BCRv1" or p_args.species != 'mouse') and workflow2_link.find('ImmuneProfiler.jar') >= 0:
        elif sub_dir == 'run_migec':
            fq1 = find_file(key + '_R2*.fastq', os.path.join(p_args.out_dir, sub_dir, 'assemble'))[0]
            fq2 = find_file(key + '_R1*.fastq', os.path.join(p_args.out_dir, sub_dir, 'assemble'))[0]
            readslayout = 'Collinear'
        else:
            fq1 = find_file(key + '_R1*.fastq', os.path.join(p_args.out_dir, 'preprocess'))[0]
            fq2 = find_file(key + '_R2*.fastq', os.path.join(p_args.out_dir, 'preprocess'))[0]
        #if lab == 2: fq2 = find_file(key.replace("_sub","_") + '_match_R1*.fastq', os.path.join(p_args.out_dir, sub_dir, 'assemble'))[0]
        vdjca = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '.vdjca')
        align_report = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_align_report.txt')
        clns = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '.clns')
        clone_report = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_clones_report.txt')
        clone_all = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_clones_all.txt')
        vdjviz_input = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_vdjviz.txt')
        airr_tsv = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_airr.tsv')
        output_id = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type)
        na_R1 = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_notAligned_R1.fastq')
        na_R2 = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_notAligned_R2.fastq')
        exportaligns = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_alignments.tsv')
        exportPretty = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_alignmentsPretty.tsv')
        exportreadIds = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_readIds.tsv')
        exportclonesall = os.path.join(p_args.out_dir, mixcr_dir, process_type, key + '_' + process_type + '_clones_all.tsv')
        if p_args.receptor_kit in ['TCRv1', 'TCRv2']:
            chains = 'TRA,TRB'
        elif p_args.receptor_kit in ['BCRv1', 'BCRv2', 'BCRv2sub', 'BCRv2hy']:
            chains = 'IG'
        if p_args.species == 'rhesus_monkey':
            # sp_cmd = ' --library imgt -s 9544 '
            mjson=os.path.join(REPO_PATH, 'src', 'config', 'imgt.201946-3.sv6.json')
            # copyfile(mjson, os.path.join(p_args.out_dir, 'run_migec', 'assemble', 'imgt.201946-3.sv6.json'))
            copyfile(mjson,os.path.join(p_args.main_dir, 'imgt.201946-3.sv6.json'))
            sp_cmd = ' --library imgt.201946-3.sv6 -s rhesus_monkey '
        elif p_args.species == 'mouse':
            sp_cmd = ' -s mouse '
        else:
            sp_cmd = ' -s human '
        threads = ''
        if len(p_args.process_threads) > 0:
            threads = ' --threads ' + p_args.process_threads
        if workflow2_link.find('ImmuneProfilerv1.6.jar') >= 0:
            align_cmd = workflow2_link + ' workflow2 align ' + threads + ' -p kaligner2 -OreadsLayout=' + readslayout + \
                    ' -OvParameters.geneFeatureToAlign=Vtranscript ' + \
                    sp_cmd + fq1 + ' ' + fq2 + ' ' + vdjca + ' --report ' + align_report
            if p_args.internal_test :
                align_cmd = workflow2_link + ' workflow2 align -f ' + threads + ' -p kaligner2 -OreadsLayout=' + readslayout + \
                    ' -OvParameters.geneFeatureToAlign=Vtranscript -OsaveOriginalReads=True ' + \
                    ' --not-aligned-R1 ' + na_R1 + ' --not-aligned-R2 ' + na_R2 + \
                    sp_cmd  + fq1 + ' ' + fq2 + ' ' + vdjca + ' --report ' + align_report
                if p_args.allow_partial:
                    align_cmd = workflow2_link + ' workflow2 align -f ' + threads + ' -p kaligner2 -OreadsLayout=' + readslayout + \
                        ' -OvParameters.geneFeatureToAlign=Vtranscript -OsaveOriginalReads=True ' + \
                        ' -OallowPartialAlignments=True -OallowNoCDR3PartAlignments=True ' + \
                        ' --not-aligned-R1 ' + na_R1 + ' --not-aligned-R2 ' + na_R2 + \
                        sp_cmd  + fq1 + ' ' + fq2 + ' ' + vdjca + ' --report ' + align_report    
            logger.info('Launching MiXCR alignment using command:\n' + align_cmd)
            run_system_cmd(align_cmd, out_log, run_mixcr_err, 'alignment failed', log_name)
            if target_region == 'cdr3':
                assem_cmd = workflow2_link +  ' workflow2 assemble ' + vdjca + ' ' + clns + ' --report ' + clone_report
            elif target_region == 'fl':
                assem_cmd = workflow2_link +  ' workflow2 assemble ' + ' -OassemblingFeatures=VDJRegion ' + \
                        vdjca + ' ' + clns + ' --report ' + clone_report
            logger.info('Launching MiXCR assembling using command:\n' + assem_cmd)
            run_system_cmd(assem_cmd, out_log, run_mixcr_err, 'assembling failed', log_name)
            
            expt_aligns_cmd = workflow2_link +  ' workflow2 exportAlignments ' + vdjca + ' ' + exportaligns
            logger.info('Launching MiXCR exportAlignments using command:\n' + expt_aligns_cmd)
            run_system_cmd(expt_aligns_cmd, out_log, run_mixcr_err, 'exportAlignments failed', log_name)
            #expt_cmd = workflow2_link +  ' workflow2 exportClones --preset full -uniqueTagCount UMI ' + \
            expt_cmd = workflow2_link +  ' workflow2 exportClones --preset full ' + ' -c ' + chains + ' ' + clns + ' ' + clone_all
            logger.info('Launching MiXCR clone exporting using command:\n' + expt_cmd)
            run_system_cmd(expt_cmd, out_log, run_mixcr_err, 'clone exporting failed', log_name)
        else:#workflow2_link.find('ImmuneProfilerv2.jar')>=0
            cmd = workflow2_link + ' workflow2 analyze ' + preset + ' '
            if p_args.umi_cutoff != '' :#and p_args.species == 'human':
                cmd += ' -Massemble.consensusAssemblerParameters.assembler.minRecordsPerConsensus=' 
                cmd += p_args.umi_cutoff + ' -MrefineTagsAndSort.parameters.postFilter=null '
            cmd += fq1 + ' ' + fq2 + ' ' + output_id + threads
            logger.info('Launching MiXCR alignment using command:\n' + cmd)
            run_system_cmd(cmd, out_log, run_mixcr_err, 'alignment failed', log_name)
            airr_cmd = workflow2_link + ' workflow2 exportAirr ' + clns + ' ' + airr_tsv
            logger.info('Launching MiXCR exportAirr using command:\n' + airr_cmd)
            run_system_cmd(airr_cmd, out_log, run_mixcr_err, 'exportAirr failed', log_name) 
        if (p_args.receptor_kit == 'BCRv2hy' and (key.find('IGG') >= 0 or key.find('IGA') >= 0) ) or p_args.internal_test :
            expt_alignsPretty_cmd = workflow2_link +  ' workflow2 exportAlignmentsPretty ' + vdjca + ' ' + exportPretty
            logger.info('Launching MiXCR exportAlignmentsPretty using command:\n' + expt_alignsPretty_cmd)
            run_system_cmd(expt_alignsPretty_cmd, out_log, run_mixcr_err, 'exportAlignmentsPretty failed', log_name)
            expt_readIds_cmd = workflow2_link +  ' workflow2 exportAlignments -readIds ' + vdjca + ' ' + exportreadIds
            logger.info('Launching MiXCR exportAlignmentsreadIds using command:\n' + expt_readIds_cmd)
            run_system_cmd(expt_readIds_cmd, out_log, run_mixcr_err, 'exportAlignmentsreadIds failed', log_name)
            if key.find("_subIGG") >= 0 or key.find("_subIGA") >= 0: 
                get_subIg(p_args, exportaligns, key, read_type, target_region, mixcr_dir, exportclonesall) 
        if p_args.vdjviz_inputs:
            vdjviz_convert(clone_all,vdjviz_input)
    if workflow2_link.find('ImmuneProfilerv2.jar')>=0 and createPDF:
        alignQc_cmd = workflow2_link + ' workflow2 exportQc align '  + samps_clns + ' ' + alignQc_pdf
        logger.info('Launching MiXCR exportQc using command:\n' + alignQc_cmd)
        run_system_cmd(alignQc_cmd, out_log, run_mixcr_err, 'exportQc align failed', log_name)
        chainUsage_cmd = workflow2_link + ' workflow2 exportQc chainUsage '  + samps_clns + ' ' + chainUsage_pdf
        logger.info('Launching MiXCR exportQc using command:\n' + chainUsage_cmd)
        run_system_cmd(chainUsage_cmd, out_log, run_mixcr_err, 'exportQc chainUsage failed', log_name)
    if not p_args.keep_inter_file:
        vdjca_list = find_file('*.vdjca', os.path.join(p_args.out_dir, mixcr_dir, process_type))
        clns_list = find_file('*.clns', os.path.join(p_args.out_dir, mixcr_dir, process_type))
        del_list = vdjca_list + clns_list
        if len(del_list) > 0:
            for item in del_list:
                os.remove(item)

def run_mixcr(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    #out_log = os.path.join(p_args.out_dir, mixcr_dir, log_name)
    #run_mixcr_err = os.path.join(p_args.out_dir, mixcr_dir, errfile_name)
    process_type = str(read_type) + '_' + str(target_region)
    tmp_dict = copy.deepcopy(input_dict)
    create_dir(os.path.join(p_args.out_dir, mixcr_dir, process_type))
    run_mixcr_base(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)
    
    cur_dict = summ_results(p_args.out_dir, process_type, tmp_dict, align_stats_dict, clone_stats_dict)
    report_stats(p_args.receptor_kit, os.path.join(p_args.out_dir, mixcr_dir, process_type), p_args.out_name + '_' + process_type, cur_dict)
    sample_list = list(set([k.split("_")[0] for k in fd.keys()]))
    for sample_id in sample_list:
        stats_file = merge_csv(sample_id, process_type,
                  p_args.out_name + '_' + sample_id + '_' + process_type,
                  os.path.join(p_args.out_dir, mixcr_dir, process_type),
                  os.path.join(p_args.out_dir, 'report'))
        if p_args.receptor_kit == 'BCRv2hy':
            subIgfs = find_file(sample_id+'_subIG*_subIg?.tsv', os.path.join(p_args.out_dir,mixcr_dir,process_type))
            for subIgf in subIgfs:
                os.rename(subIgf, os.path.join(p_args.out_dir, 'report', os.path.basename(subIgf)))
        if p_args.airr_create:
            merge_airr_csv(sample_id, process_type, 
                  p_args.out_name + '_' + sample_id + '_' + process_type,
                  os.path.join(p_args.out_dir, mixcr_dir, process_type),
                  os.path.join(p_args.out_dir, 'report'),
                  os.path.join(p_args.out_dir, 'airr_report'))     
    os.rename(stats_file, os.path.join(p_args.out_dir, 'report', os.path.basename(stats_file)))
    reg_str = 'CDR3'
    if target_region == 'fl':
        reg_str = 'Full_length'
    print(get_time() + ' [INFO] Finished ' + reg_str + ' region processing for all samples', flush=True)
    if p_args.airr_create:
        stats_csv_to_airr( os.path.join(p_args.out_dir, 'report', os.path.basename(stats_file)), os.path.join(p_args.out_dir, 'airr_report', os.path.basename(stats_file)))

def run_raw_mixcr(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name):
    #out_log = os.path.join(p_args.out_dir, mixcr_dir, log_name)
    #run_mixcr_err = os.path.join(p_args.out_dir, mixcr_dir, errfile_name)
    process_type = str(read_type) + '_' + str(target_region)
    tmp_dict = copy.deepcopy(input_dict)
    create_dir(os.path.join(p_args.out_dir, mixcr_dir, process_type))
    run_mixcr_base(p_args, input_dict, read_type, target_region, mixcr_dir, sub_dir, log_name, errfile_name)

    cur_dict = summ_results(p_args.out_dir, process_type, tmp_dict, align_stats_dict, clone_stats_dict)
    report_stats(p_args.receptor_kit, os.path.join(p_args.out_dir, mixcr_dir, process_type), p_args.out_name + '_' + process_type, cur_dict)
    sample_list = list(set([k.split("_")[0] for k in fd.keys()]))
    for sample_id in sample_list:
        stats_file = merge_csv(sample_id,
                  process_type,
                  p_args.out_name + '_' + sample_id + '_' + process_type,
                  os.path.join(p_args.out_dir, mixcr_dir, process_type),
                  os.path.join(p_args.out_dir, 'report'))
        if workflow2_link.find('ImmuneProfilerv2.jar')>=0 and p_args.airr_create:
            airrfiles = find_file(sample_id+'_*_airr.tsv', os.path.join(p_args.out_dir, mixcr_dir, process_type))
            for airr_orig in airrfiles:
                airr_dist = os.path.join(p_args.out_dir, "report", sample_id, os.path.basename(airr_orig))
                os.rename(airr_orig, airr_dist)
    os.rename(stats_file, os.path.join(p_args.out_dir, 'report', os.path.basename(stats_file)))
    reg_str = 'CDR3'
    if target_region == 'fl':
        reg_str = 'Full_length'
    print(get_time() + ' [INFO] Finished raw ' + reg_str + ' region processing for all samples', flush=True)


def change_run_mixcr_name(p_args, mytype):
    #if p_args.receptor_kit != 'TCRv1' and (p_args.receptor_kit != 'BCRv1' or p_args.species != 'mouse'):
    #if p_args.library_kit == 'takara_smartseq' and (workflow2_link.find('ImmuneProfiler.jar')>=0 or p_args.internal_test):
    change_dir(os.path.join(p_args.out_dir, 'run_mixcr'), os.path.join(p_args.out_dir, mytype+'_mixcr'))
    create_dir(os.path.join(p_args.out_dir, 'run_mixcr'))
    change_dir(os.path.join(p_args.out_dir, 'report'), os.path.join(p_args.out_dir, mytype+'_report'))
    create_dir(os.path.join(p_args.out_dir, 'report'))
    #if not read_type.endswith('rawUnd'):
    src = os.path.join(p_args.out_dir, mytype+'_report', p_args.out_name+"_sample_QC_stats.csv")
    dst = os.path.join(p_args.out_dir, 'report', p_args.out_name+"_sample_QC_stats.csv")
    shutil.copy(src, dst)
    if p_args.airr_create and workflow2_link.find('ImmuneProfiler.jar')>=0:
        change_dir(os.path.join(p_args.out_dir, 'airr_report'), os.path.join(p_args.out_dir, mytype+'_airr_report'))
        create_dir(os.path.join(p_args.out_dir, 'airr_report'))
        if mytype.find('raw')>=0:
            src = os.path.join(p_args.out_dir, mytype+'_airr_report', p_args.out_name+"_sample_QC_stats.csv")
            dst = os.path.join(p_args.out_dir, 'airr_report', p_args.out_name+"_sample_QC_stats.csv")
            shutil.copy(src, dst)
    if p_args.receptor_kit == 'BCRv2hy' and not p_args.internal_test and not debug:
        pre_dir =os.path.join(p_args.out_dir, mytype+'_mixcr')
        shutil.rmtree(pre_dir)
        pre_dir =os.path.join(p_args.out_dir, mytype+'_report')
        shutil.rmtree(pre_dir)
        if p_args.airr_create:
            pre_dir =os.path.join(p_args.out_dir, mytype+'_airr_report')
            shutil.rmtree(pre_dir)

def parse_clones_report(in_file, out_file, key, inputDict):
    with open(in_file, 'r') as infile, open(out_file, 'a') as outfile:
        next(infile)
        header = ['Read Count', 'Fraction',
                  'Clonal Sequence', 'Clonal Sequence Quality',
                  'CDR3 Min Quality',
                  'CDR3 Sequence', 'CDR3 Amino Acid Sequence',
                  'Clonal Type',
                  'Frame Shift', 'Stop Codon',
                  'Amino Acid Length',
                  'V segment', 'all V hits',
                  'D segment', 'all D hits',
                  'J segment', 'all J hits',
                  'C segment', 'all C hits']
        outfile.write(','.join(header) + '\n')
        chains = set()
        clones = set()
        chain_types = set()
        countNoC = 0
        countAll = 0
        countHnoC = 0
        countG = 0
        countM = 0
        countK = 0
        countL = 0
        countA = 0
        countD = 0
        countE = 0
        countTRA = 0
        countTRB = 0
        countClnNoC = 0
        countClnAll = 0
        countClnHnoC = 0
        countClnG = 0
        countClnM = 0
        countClnK = 0
        countClnL = 0
        countClnA = 0
        countClnD = 0
        countClnE = 0
        countClnTRA = 0
        countClnTRB = 0
        shift = 0
        if workflow2_link.find('ImmuneProfilerv2.jar') >= 0:
            shift =2
            if p_args.species == 'mouse' and p_args.receptor_kit in ['BCRv1', 'TCRv1']:
                shift = 0
        for line in infile:
            line = line.rstrip().split('\t')
            cloneId = line[0]
            cloneCount = line[1]
            cloneFraction = line[2]
            clonalSequence = line[3+shift]
            clonalSequenceQuality = line[4+shift]
            
            nSeqCDR3 = line[13+shift]
            minQualCDR3 = line[14+shift]
            aaSeqCDR3 = line[15+shift]
            if workflow2_link.find('ImmuneProfiler.jar') >= 0:
                nSeqCDR3 = line[23]
                minQualCDR3 = line[24]
                aaSeqCDR3 = line[32]
            if '_' in aaSeqCDR3:
                aaframeshift = 'Frame Shift'
            else:
                aaframeshift = ''
            if '*' in aaSeqCDR3:
                aaStopCodon = 'Stop Codon'
            else:
                aaStopCodon = ''
            aaLength = str(len(aaSeqCDR3))
            vlist = line[5+shift].split(',')
            Vsegments = [x.split('*')[0] for x in vlist][0]
            vmulti = ';'.join([x.split('*')[0] for x in vlist])
            dlist = line[6+shift].split(',')
            Dsegments = [x.split('*')[0] for x in dlist][0]
            dmulti = ';'.join([x.split('*')[0] for x in dlist])
            jlist = line[7+shift].split(',')
            Jsegments = [x.split('*')[0] for x in jlist][0]
            jmulti = ';'.join([x.split('*')[0] for x in jlist])
            c_list = line[8+shift].split(',')
            Csegments = [x.split('*')[0] for x in c_list][0]
            cmulti = ';'.join([x.split('*')[0] for x in c_list])
            
            # use V region to find if is IG or TR
            chain_type = Vsegments.split(',')[0][:2]
            if chain_type == 'TR':
                chain = Vsegments.split(',')[0][2]
                clonalType = Vsegments.split(',')[0][:3]
                if Csegments.split(',')[0][:2] == '':
                    countNoC += int(float(cloneCount))
                    countClnNoC += 1
            elif chain_type == 'IG':
                if Csegments.split(',')[0][:2] != '':
                    # if any constant region, use for heavy chain seperation
                    if Csegments.split(',')[0][2:3] == 'H':
                        chain = Csegments.split(',')[0][3]  # G/M
                        clonalType = Csegments.split(',')[0][:2] + \
                                     Csegments.split(',')[0][3]  # IGG, IGM
                    else:
                        chain = Csegments.split(',')[0][2]  # K/L
                        clonalType = Csegments.split(',')[0][:3]
                else:  # no constant region, use V region
                    chain = Vsegments.split(',')[0][2]  # K/L/H
                    clonalType = Vsegments.split(',')[0][:3]
                    countNoC += int(float(cloneCount))
                    countClnNoC += 1
            chain_types.add(chain_type)  # either IG or TR
            chains.add(chain)  # H,G,M,K,L,D,E,A,B
            countAll += int(float(cloneCount))
            countClnAll += 1
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
            if 'TR' in chain_types:  # align to TRA,TRB only
                if chain == 'A':
                    countTRA += int(float(cloneCount))
                    countClnTRA += 1
                elif chain == 'B':
                    countTRB += int(float(cloneCount))
                    countClnTRB += 1
            result = [cloneCount, cloneFraction, clonalSequence,
                      clonalSequenceQuality, minQualCDR3,
                      nSeqCDR3, aaSeqCDR3, clonalType,
                      aaframeshift, aaStopCodon, aaLength,
                      Vsegments, vmulti, Dsegments, dmulti,
                      Jsegments, jmulti, Csegments, cmulti
                      ]
            outfile.write(','.join(result))  # write as csv format
            outfile.write('\n')
        # write stats as csv file
        # merge this portion into summary stats table
        # outfile.write(',No.Reads, No.Clone\n')
        length = len(inputDict[key])
        if 'IG' in chain_types:
            if inputDict[key][length-8]==None or inputDict[key][length-8].find('IgG:')<0:
                inputDict[key].append('IgG:' + str(countClnG))
                inputDict[key].append('IgM:' + str(countClnM))
                inputDict[key].append('IgK:' + str(countClnK))
                inputDict[key].append('IgL:' + str(countClnL))
                inputDict[key].append('IgA:' + str(countClnA))
                inputDict[key].append('IgD:' + str(countClnD))
                inputDict[key].append('IgE:' + str(countClnE))
                inputDict[key].append('IgH(lack constant region):' + str(countClnHnoC))
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
            if inputDict[key][length-2].find('TRA:')<0:
                inputDict[key].append('TRA:' + str(countClnTRA))
                inputDict[key].append('TRB:' + str(countClnTRB))
            else:
                merge_chain_value(inputDict, key, 'TRA:', countClnTRA, length-2)
                merge_chain_value(inputDict, key, 'TRB:', countClnTRB, length-1)


def merge_chain_value(inputDict, key, chain, value, index):
    old_value = int(inputDict[key][index].split(chain)[1])
    inputDict[key][index] = chain + str(old_value + value)
    return

def vdjviz_convert(file_clones_all, chord_file):
    """ Function: vdjviz_convert, convert mixcr clones_all.txt
        into vdjviz input file
    """
    if os.path.getsize(file_clones_all) > 0:
        with open(chord_file, 'w') as out_file:
            #headers = ['#count','freq','cdr3nt','cdr3aa','v','d','j']
            headers = ['count','frequency','CDR3nt','CDR3aa','V','D','J']
            out_file.write('\t'.join(headers)+'\n')
        with open(file_clones_all) as all_clones:
                all_clones_line_num = 0
                for line in all_clones:
                    all_clones_line_num += 1
                    if all_clones_line_num != 1:
                        all_clones_line_content = line.rstrip()
                        all_clones_line_content = \
                            re.split(r"\t", 
                                all_clones_line_content.rstrip("\t"))
                        all_clones_line_content[5] = \
                            all_clones_line_content[5].split("*")[0]
                        all_clones_line_content[6] = \
                            all_clones_line_content[6].split("*")[0]
                        all_clones_line_content[7] = \
                            all_clones_line_content[7].split("*")[0]
                        if len(all_clones_line_content[5]) == 0:
                            all_clones_line_content[5]=str(r".")
                        if len(all_clones_line_content[6]) == 0:
                            all_clones_line_content[6]=str(r".")
                        if len(all_clones_line_content[7]) == 0:
                            all_clones_line_content[7]=str(r".") 
                        newline = [str(all_clones_line_content[1]),
                                   str(all_clones_line_content[2]),
                                   str(all_clones_line_content[23]),
                                   str(all_clones_line_content[32]),
                                   str(all_clones_line_content[5]),
                                   str(all_clones_line_content[6]), 
                                   str(all_clones_line_content[7])]
                        with open(chord_file, "a") as out_file:
                            out_file.write("\t".join(newline[0:])+"\n")
                if all_clones_line_num < 3:
                    os.remove(chord_file)
                else:
                    pass
    else:
        pass

def get_mergedfqs(p_args, log_name, errfile_name):
    out_log = os.path.join(p_args.out_dir, log_name)
    merge_err = os.path.join(p_args.out_dir, errfile_name)
    print("merged fastq for hybird kit done")

def split_align_unalign_by_readIds(fq_fn, id_fn):
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
    fq_fn_aligned = fq_fn[:len(fq_fn)-6] + "_aligned.fastq"
    fq_fn_unaligned = fq_fn[:len(fq_fn)-6] + "_unaligned.fastq"
    save_list_to_file(fq_fn_aligned, readsaligned, True)
    save_list_to_file(fq_fn_unaligned, readsunaligned, True)

def umi_polyGNs(umi):
    ispolyGs = False
    i = 0
    while i <len(umi) and not ispolyGs:
        if umi[i]=='G' and i+4<len(umi) and umi[i+1] =='G' and umi[i+2] =='G' and umi[i+3] =='G' and umi[i+4] =='G':
            ispolyGs = True
        if umi[i]=='N' and i+4<len(umi) and umi[i+1] =='N' and umi[i+2] =='N' and umi[i+3] =='N' and umi[i+4] =='N':
            ispolyGs = True
        i += 1
    return ispolyGs

def get_adapter_args(p_args):
    adapter_start = 0
    adapter_len = 0
    if p_args.library_kit == 'takara_smartseq':
        if p_args.species == 'human':
            if p_args.receptor_kit == 'TCRv1':
                adapter_start = 8
                adapter_len = 18 #+1
            elif p_args.receptor_kit == 'TCRv2':
                adapter_start = 5
                adapter_len = 21
            elif p_args.receptor_kit == 'BCRv1':
                adapter_start = 0
                adapter_len = 0
            elif p_args.receptor_kit == 'BCRv2':
                adapter_start = 7
                adapter_len = 19 #+5
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
                adapter_len = 20 #+2
    return adapter_start, adapter_len 

def remove_umi_rawUnd(p_args, valfqs):
    raw_noumi_fq = os.path.join(p_args.fastq_dir, valfqs[0])
    linker_range = 0
    umi_dict = OrderedDict()
    umi_start = 0
    len_umi = 12
    len_umi_link = len_umi +7
    raw_umi_fq = os.path.join(p_args.fastq_dir, valfqs[1]) #r2 for 'takara_smartseq'
    if p_args.library_kit == 'cellecta_air':
        len_umi = 14
        len_umi_link = 14
    r1fq = []
    r2fq = []
    umiread_count = 0
    total = 0
    # remove adapter from r1
    adapter_start = 0
    adapter_len = 0
    if p_args.rm_adapter:
        adapter_start, adapter_len = get_adapter_args(p_args)
    with load_fastq(raw_noumi_fq) as f1, load_fastq(raw_umi_fq) as f2 :
        prefix1 = os.path.basename(raw_noumi_fq).split('.fastq')[0]
        out_noumi_fq = os.path.join(p_args.fastqRMumi_dir,  prefix1 + '.fastq')
        prefix2 = os.path.basename(raw_umi_fq).split('.fastq')[0]
        out_rmumi_fq = os.path.join(p_args.fastqRMumi_dir,  prefix2 + '.fastq')
        fw1 = open(out_noumi_fq, 'w')
        fw2 = open(out_rmumi_fq, 'w')
        for x, y in zip(f1,f2):
            r1fq.append(x.rstrip())
            r2fq.append(y.rstrip())
            if len(r1fq) == 4:
                total += 1
                umi = get_umi(r2fq, umi_start, len_umi)
                if not umi_polyGNs(umi):
                    if (p_args.library_kit == 'takara_smartseq' and r2fq[1][len_umi:len_umi_link] == 'GTACGGG') or p_args.library_kit != 'takara_smartseq':
                        umiread_count += 1
                        if umi in umi_dict.keys():
                            umi_dict[umi] += 1
                        else:
                            umi_dict[umi] = 1
                r2fq[1] = r2fq[1][umi_start+len_umi_link:]
                r2fq[3] = r2fq[3][umi_start+len_umi_link:]
                
                r1fq[1] = r1fq[1][adapter_start + adapter_len:]
                r1fq[3] = r1fq[3][adapter_start + adapter_len:]
                for item in r1fq:
                    fw1.write('%s\n' % item)
                for item in r2fq:
                    fw2.write('%s\n' % item)
                del r1fq[:]
                del r2fq[:]
        fw1.close()
        fw2.close()
    if platform.system() == 'Windows':
        dos2unix(out_noumi_fq)
        dos2unix(out_rmumi_fq)
    #if raw_noumi_fq.endswith('.gz'):
    #    #cmd = 'tar -vczf ' + out_rmumi_fq + '.gz ' + out_rmumi_fq
    #    os.system('gunzip ' + out_noumi_fq)
    #    run_system_cmd(cmd, 'run_gzip.log', 'run_gzip.error', 'gzip failed', 'NA')
    #    os.remove(out_rmumi_fq)
    print("total reads=", total, "\numi reads=", umiread_count)
    return umi_dict, total, umiread_count

def insertmax_returnmin(list_key5max, list_val5max, key_seq, val):
    len_5top = 5
    index = 0
    while index < len_5top and val > list_val5max[index]:
        index += 1
    if index > 0:
        index -= 1
        for i in range(index):
            if i<len_5top:
                list_key5max[i] = list_key5max[i+1]
                list_val5max[i] = list_val5max[i+1]
        list_key5max[index] = key_seq
        list_val5max[index] = val
    return list_key5max[0], list_val5max[0]

def visualize_umi_distribution(p_args, target_range, umi_dict, total, umireads_count):
    numkeys = 5
    list_key5max = ['NNNNNNNNNNNN']*numkeys
    list_val5max = [0]*numkeys 
    keymax = list(umi_dict.keys())[0]
    valmax = umi_dict[keymax]
    keymax, valmax = insertmax_returnmin(list_key5max, list_val5max, keymax, valmax)
    for key_umi, val_umicount in umi_dict.items():
        if val_umicount > valmax:
            valmax = val_umicount
            keymax = key_umi
            keymax, valmax = insertmax_returnmin(list_key5max, list_val5max, keymax, valmax)
    umi_summary_f = os.path.join(p_args.fastqRMumi_dir,target_range + '_umi_summary.csv')
    with open(umi_summary_f, 'w') as fw:
        print('total reads: %s' % str(total))
        fw.write('total reads:' + str(total) + '\n')
        print('umi reads: %s' % str(umireads_count))
        fw.write('umi reads:' + str(umireads_count) + '\n')
        print('number of umis: %s' % len(umi_dict.keys()))
        fw.write('number of umis:' + str(len(umi_dict.keys()))+'\n')
        print('top 5 umis=', list_key5max)
        fw.write('top 5 umis=')
        fw.write(','.join(str(item) for item in list_key5max))
        fw.write('\n')
        print('top 5 umis with count=', list_val5max)
        fw.write('top 5 umis with count=')
        fw.write(','.join(str(item) for item in list_val5max))
        fw.write('\n')
    valmax = list_val5max[-1]
    umi_counts = [0] * (valmax+1)
    for key_umi, val_umicount in umi_dict.items():
        umi_counts[ val_umicount ] += 1
    x = [*range(0,valmax+1,1)]
    y = umi_counts
    umi_range_f = os.path.join(p_args.fastqRMumi_dir,target_range + '_umi_distribution.csv')
    with open(umi_range_f, 'w') as fp:
        fp.write('umi_type_numbers,')
        fp.write(','.join(str(item) for item in x))
        fp.write('\n')
        fp.write('frequency,')
        fp.write(','.join(str(item) for item in y))    
    drawplot(x, y, umi_range_f[:-4]+'.png')
    return

def drawplot(x, y, fname):
    len_exp = len(y)
    if y[len(y)-1]==1:
        len_exp -=1
    while len_exp>50 and (y[len_exp-1]==0 or y[len_exp-2]==0 or y[len_exp-3]==0 or y[len_exp-4]==0):
        len_exp -=1
    x = x[:len_exp]
    y = y[:len_exp]
    basename = os.path.basename(fname)
    sampleID = basename.split('_')[0]
    setname = 'takara TCR '
    if fname.find('outCellecta') > 0:
        setname = 'cellecta '
    else:
        setname = 'takara '
    if fname.find('TCR') > 0:
        setname += 'TCR '
    else:
        setname += 'BCR '
    for i in range(4):
        plt.plot(x, y)
        plt.xlabel('umi_type_numbers')
        plt.ylabel('frequency of umi_types')
        plt.title('x=umi_types and y=frequency of umi_types with ' + setname + sampleID)
        #plt.xticks(np.arange(0, len_exp, step=20))
        #plt.yticks(np.arange(0, 8000, step=1000))
        if i == 1:
            plt.xlim([0,50])
        if i == 2:
            plt.xlim([0,50])
            plt.ylim([0,20000])
        if i == 3:
            plt.ylim([0,2000])
        fwname = fname[:-4] + str(i) +'.png'
        if i == 0:
            fwname = fname[:-4] + '.png'
        f = plt.savefig(fwname)
        plt.close(f)


def drawUMIdistribution():
    for i in range(9,11):
        frname = 'C:\\Users\\suns\\Downloads\\John_peak\\outmTCRv2ilku1\\John_peak_rmumi\\S'+str(i)+'_umi_distribution.csv'
        fr = open(frname, 'r')
        a = fr.readline().replace('umi_count,', '').split(',')
        b = fr.readline().replace('frequency,', '').split(',')
        fr.close()
        x = [eval(i) for i in a]
        y = [eval(i) for i in b]
        drawplot(x, y, frname[:-4] + '.png')

def get_meta_file(p_args):
    basename = os.path.basename(p_args.fastq_dir)
    p_args.fastq_dir = os.path.dirname(p_args.fastq_dir)
    fastqs = find_file(basename, p_args.fastq_dir)
    fastqs.sort()
    #if p_args.meta_file == '':
    if True:
        p_args.meta_file = os.path.join(p_args.fastq_dir, p_args.receptor_kit + '_meta.csv')
        f = open(p_args.meta_file,'w')
        f.write('sampleID,read1_file_name,read2_file_name\n')
        sampleID = ''
        for i in range(len(fastqs)):
            fqname = os.path.basename(fastqs[i])
            if fqname.split('_R')[0] != sampleID:
                sampleID = fqname.split('_R')[0]
                f.write(sampleID + ',' + fqname)
            else:
                f.write(',' + fqname +'\n')
        f.close()
        if platform.system() == 'Windows':
            dos2unix(p_args.meta_file)
    return p_args.fastq_dir, p_args.meta_file

def comp_mixcr_report_v3_v4_matrix(v3reportf, v4reportf, receptor):
    v3matrix = []
    v4matrix = []
    results = []
    chains = ['TRA','TRB']
    id_chain = 7
    id_seq = 5
    if receptor == 'BCRv2':
        chains = ['IGG','IGM','IGK','IGL','IGA','IGD','IGE','IGH']
    for i in range(len(chains)) :
        v3matrix.append([])
        v4matrix.append([])
    for i in range(2):
        vf =v3reportf
        if i == 1:
            vf =v4reportf
        with open(vf, mode ='r')as file:
            csvFile = csv.reader(file)
            count = 0
            for line in csvFile:
                if count > 0:
                    chain_type = line[id_chain]
                    seq = line[id_seq]
                    for id in range(len(chains)):
                        if chain_type == chains[id]:
                            if i == 0:
                                v3matrix[id].append(seq)
                            else:
                                v4matrix[id].append(seq)
                count += 1
    for i in range(len(chains)) :
        v3set = set(v3matrix[i])
        v4set = set(v4matrix[i])
        interset = v3set.intersection(v4set)
        overlap = len(interset)
        results.append([len(v3matrix[i]), len(v3set), overlap, len(v4matrix[i]), len(v4set)])
        print(len(v3matrix[i]), len(v3set), overlap, len(v4matrix[i]), len(v4set))
    return results

def plot_multi_bar(results, fb):
    N = 3 * len(results) - 1
    bases = (results[0][2], results[0][2])
    diffs = (results[0][1]-results[0][2], results[0][4]-results[0][2])
    for i in range(1, len(results)):
        bases = bases + (0, results[i][2], results[i][2])
        diffs = diffs + (0, results[i][1]-results[i][2], results[i][4]-results[i][2])
    ind = np.arange(N)  
    width = 0.35 
    
    fig = plt.subplots(figsize =(1*N, 7))
    p1 = plt.bar(ind, bases, width)#, yerr = boyStd)
    p2 = plt.bar(ind, diffs, width, bottom = bases)#, yerr = girlStd)
    
    plt.ylabel('clonotype count')
    tname = os.path.basename(fb).split('_bar.png')[0]
    segs = tname.split('_')
    tname = segs[len(segs)-1]
    if len(results) == 2:
        tname += ' TRA,TRB clonotype count MIXCR v3 vs v4'
        xlab = ('v3TRA','v4TRA','','v3TRB','v4TRB')
    else:
        tname += ' IgG,IgM,IgK,IgL,IgA,IgD,IgE,IgH clonotype count MIXCR v3 vs v4'
        xlab = ('v3IgG','v4IgG','','v3IgM','v4IgM','','v3IgK','v4IgK','','v3IgL','v4IgL','','v3IgA','v4IgA','','v3IgD','v4IgD','','v3IgE','v4IgE','','v3IgH','v4IgH')
    plt.title(tname)
    plt.xticks(ind, xlab)
    #plt.yticks(np.arange(0, 81, 10))
    plt.legend((p1[0], p2[0]), ('bases', 'diffs'))
    
    f = plt.savefig(fb)
    plt.close(f)

def get_modify_fastq_head(ifd, ifr, barcode, ofd):
    fwname = os.path.joun(ofd, ifr)
    fw = open(fwname, 'w')
    fr = os.path.join(ifd, ifr)
    with open(fr, 'r') as f:
        count = 0
        for line in f:
            count += 1
            if count % 4 == 1:
                segs = line.split(' ')
                newline = segs[0] + '_' + barcode + ' ' + segs[1]
                fw.write(newline)
            else:
                fw.write(line)
    fw.close()

def get_demux_counts_all_csv_fastq(ifpath, ofd):
    ifd = os.path.dirname(ifpath)
    with open(ifpath) as f:
        for line in f:
            segs = line.split(',')
            if len(segs[0]) == 16:
                key = segs[1].replace('_','')
                ifr = 'demux_' + segs[0]+'_'+key+'_R1.fastq.gz'
                get_modify_fastq_head(ifd, ifr, segs[0], ofd)
                ifr = 'demux_' + segs[0]+'_'+key+'_R2.fastq.gz'
                get_modify_fastq_head(ifd, ifr, segs[0], ofd)

def splitfastq_by_adapter(sample, fastq_dir, adapterseq):
    r1 = []
    r2 = []
    raw_r1_fq = os.path.join(fastq_dir, sample + "_L001_R1_001.fastq.gz")
    raw_r2_fq = os.path.join(fastq_dir, sample + "_L001_R2_001.fastq.gz")
    out_r1_fqs = []
    out_r2_fqs = []
    prs = ["with-adapter","without-adapter"]
    
    for directory in prs:
        if not os.path.isdir( os.path.join(fastq_dir, directory) ):
            create_dir(os.path.join(fastq_dir, directory))
        out_r1_fq = os.path.join(fastq_dir, directory, sample + '_L001_R1_001.fastq')
        out_r2_fq = os.path.join(fastq_dir, directory, sample + '_L001_R2_001.fastq')
        out_r1_fqs.append(out_r1_fq)
        out_r2_fqs.append(out_r2_fq)
    for item in out_r1_fqs:
        open(item, 'w+').close()
    for item in out_r2_fqs:
        open(item, 'w+').close()
    
    fhs1 = []
    fhs2 = []
    with load_fastq(raw_r1_fq) as f1, load_fastq(raw_r2_fq) as f2, ExitStack() as cm:
        for name in out_r1_fqs:
            fhs1.append(cm.enter_context(open(name, 'a')))
        for name in out_r2_fqs:
            fhs2.append(cm.enter_context(open(name, 'a')))
        
        for x, y in zip(f1, f2):
            r1.append(x.rstrip())
            r2.append(y.rstrip())
            if len(r1) == 4:
                if r1[1].find(adapterseq)>=0 or r2[1].find(adapterseq)>=0:
                    for item in r1:
                        fhs1[0].write('%s\n' % item)
                    for item in r2:
                        fhs2[0].write('%s\n' % item)
                else:
                    for item in r1:
                        fhs1[1].write('%s\n' % item)
                    for item in r2:
                        fhs2[1].write('%s\n' % item)
                del r1[:]
                del r2[:]
    if platform.system() == 'Windows':
        for name in out_r1_fqs:
            dos2unix(name)
        for name in out_r2_fqs:
            dos2unix(name)
    #tar in windows could not correctly compress fastq into fastq.gz
    for directory in prs:
        for f in os.listdir(os.path.join(fastq_dir, directory)):
            if f.endswith('.fastq'):
                fq = os.path.join(fastq_dir, directory, f)
                cmd = 'tar -vczf ' + fq + '.gz ' + fq
                run_system_cmd(cmd, 'run_gzip.log', 'run_gzip.error', 'gzip failed', 'NA')
                os.remove(fq)
    return

def clonetypes_comparison(fp1, fp2):
    f1 = open(fp1, 'r')
    f2 = open(fp2, 'r')
    flist1 = f1.readlines()
    flist2 = f2.readlines()
    fwlist = flist1.copy()
    count_same = 0 
    for j in range (1,len(flist2)):
        segs2 = flist2[j].split(",")
        i = 1
        is_same = False
        while not is_same and i < len(flist1):
            segs1 = flist1[i].split(",")
            if segs2[2]==segs1[2] and segs2[11]==segs1[11] and segs2[13]==segs1[13] and segs2[15]==segs1[15] and segs2[17]==segs1[17]:
                count_same += 1
                is_same = True
            i += 1
        if not is_same:
            fwlist.append(flist2[j])
    fw = fp1.replace('_IGG_', '_IGGMerged_')
    with open(fw, 'w') as f:
        for item in fwlist:
            f.write("%s" % item)
    return count_same, len(flist2)-1-count_same

def dosformat(fpath):
    is_dosformat = False
    f = open(fpath, 'rb')
    line = f.readline()
    if line.find(b'\r\n') >= 0 or line.find(b'\n\r') >= 0 or line.find(b'\r') >= 0:
        is_dosformat = True
    f.close()
    return is_dosformat

# plot pie for all QC chains(flc 6%)
def plot_QC_chains(fname):
    f = open(fname, 'r')
    dirname = os.path.dirname(fname)
    flist = f.readlines()
    labels = []
    for i in range (0,len(flist)):
        segs = flist[i].split(",")
        y = []
        mylabels = []
        for j in range(1, len(segs)-2):
            if i == 0:
                if j%2 == 1:
                    labels.append(segs[j])
            else:
                if j%2 == 1:
                    value = int(segs[j])
                    if value > 1:
                        y.append(value)
                if j%2 == 0 and value > 1:
                    mylabels.append(labels[j//2-1]+'_'+segs[j])
        if i > 0:
            ofname = os.path.join(os.path.dirname(dirname), 'QC_pie_plot_' + segs[0] + ".png")
            y[1], y[2] = y[2], y[1]
            mylabels[1], mylabels[2] = mylabels[2], mylabels[1]
            if len(y)>7:
                y[5], y[6] = y[6], y[5]
                mylabels[5], mylabels[6] = mylabels[6], mylabels[5]
            myexplode = [0]*len(y)
            myexplode[len(myexplode)-1]=0.2
            plt.pie(y, labels=mylabels, explode = myexplode)
            plt.title("QC plot sample_" + segs[0])
            #plt.show()
            of = plt.savefig(ofname)
            #of = plt.savefig(ofname, dpi=100)
            plt.close(of)
    return

def flc_analysis(fname):
    dirname = os.path.dirname(fname)
    count = 0
    heads = []
    data = []
    with open(fname, 'r') as fp:
        for line in fp:
            line = line.strip()
            segs = line.split(',')
            if count == 0:
                heads = segs.copy()
            else:
                data.append(segs[1:len(segs)-2])
            count += 1
    df = pd.DataFrame(data=data, columns=heads[1:len(heads)-2])
    total = len(df)
    df['pos_GTACGGG'] = df['pos_GTACGGG'].astype('int32')
    df['p0_40'] = df['p0_40'].astype('int32')
    df['sim(GTACGGG)'] = df['sim(GTACGGG)'].astype('float32')
    df['p0_40_sim_max'] = df['p0_40_sim_max'].astype('float32')
    
    for key in ['pos', 'sim_max']:
        y=[]
        mylabels = []
        if key == 'pos':
            df = df.sort_values(by=['pos_GTACGGG', 'p0_40_sim_max'], ascending=False)
            labels = ['p-1(0_40sim<3):','p-1(0_40sim=3):','p-1(0_40sim=4):','p-1(0_40sim=5):','p-1(0_40sim=6):','p0-10(sim=7):','p11(sim=7):','p13(sim=7):','p14-40(sim=7):','p41-117(sim=7):' ]
            for i in range(len(labels)):
                if i == 0:
                    sub_df = df.loc[(df['pos_GTACGGG'] ==-1) & (df['p0_40_sim_max'] < 3)]
                elif i == 1:
                    sub_df = df.loc[(df['pos_GTACGGG'] ==-1) & (df['p0_40_sim_max'] == 3)]
                elif i == 2:
                    sub_df = df.loc[(df['pos_GTACGGG'] ==-1) & (df['p0_40_sim_max'] == 4)]
                elif i == 3:
                    sub_df = df.loc[(df['pos_GTACGGG'] ==-1) & (df['p0_40_sim_max'] == 5)]
                elif i == 4:
                    sub_df = df.loc[(df['pos_GTACGGG'] ==-1) & (df['p0_40_sim_max'] == 6)]
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
                #print(sub_df.shape[0])
                y.append(sub_df.shape[0])
                mylabels.append(labels[i]+str(sub_df.shape[0])+"("+str(int(100*len(sub_df)/total+0.5))+"%)")
        else:
            df = df.sort_values(by=['p0_40_sim_max', 'pos_GTACGGG'], ascending=False)
            labels = ['0_40sim_max<3:', '0_40sim_max=3:', '0_40sim_max=4:','0_40sim_max=5:','0_40sim_max=6:','0_40sim_max=7:' ]
            for i in range(len(labels)):
                if i == 0:
                    sub_df = df.loc[ df['p0_40_sim_max'] < 3 ]
                elif i == 1:
                    sub_df = df.loc[ df['p0_40_sim_max'] == 3 ]
                elif i == 2:    
                    sub_df = df.loc[ df['p0_40_sim_max'] == 4 ]
                elif i == 3:
                    sub_df = df.loc[ df['p0_40_sim_max'] == 5 ]
                elif i == 4:
                    sub_df = df.loc[ df['p0_40_sim_max'] == 6 ]
                elif i == 5:
                    sub_df = df.loc[ df['p0_40_sim_max'] == 7 ]
                y.append(sub_df.shape[0])
                mylabels.append(labels[i]+str(sub_df.shape[0]))            
        ofname = os.path.join(dirname, 'flc_pie_plot_' + key + '_' + segs[0] + ".png")
        #y[1], y[2] = y[2], y[1]
        #mylabels[1], mylabels[2] = mylabels[2], mylabels[1]
        plt.pie(y, labels=mylabels)
        plt.title("Pie plot sample_" +  key + '_' + segs[0])
        #plt.show()
        of = plt.savefig(ofname, dpi=250)
        #of = plt.savefig(ofname)
        plt.close(of)
    return

# ---------- main ---------- #
if __name__ == '__main__':

    # ---------- setup | parser ---------- #

    parser = argparse.ArgumentParser(description=desc)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)

    # ---------- setup | user options ---------- #

    # add required arguments
    required.add_argument('-r', '--receptor_kit', dest='receptor_kit',
                          help='specify receptor kit: TCRv1, TCRv2, BCRv1 or BCRv2', choices=['TCRv1', 'TCRv2', 'BCRv1','BCRv2'],
                          required=True)
    required.add_argument('-f', '--fastq_dir', dest='fastq_dir',
                          help='a folder stores  all input FASTQs', required=True)
    required.add_argument('-m', '--meta_file', dest='meta_file',
                          help='a file contains sample ID and corresponding FASTQ pair',
                          required=True)
    required.add_argument('-o', '--output_name', dest='out_name',
                          help='Name an output directory to be created to store results ' + \
                               'and use as file prefix; ' + \
                               'directory name should be less than 20 characters',
                          required=True)
    required.add_argument('-t', '--target_region', dest='target_region',
                          help='specify target regions reads should map to',
                          choices=['CDR3', 'Full_length','Both'],
                          required=True)

    # add optional arguments
    optional.add_argument('-k', '--keep_inter_file', dest='keep_inter_file',
                          help='decide if keep intermediate files, including MiXCR files & ' +
                          'preprocessed FASTQs [Default: False]',
                          action='store_true', default=False)
    optional.add_argument('-l', '--linker_correction', dest='linker_correction',
                          help='decide if remove reads based on sequence match of linker' +
                               ' [Default: False]',
                          action='store_true', default=False)
    optional.add_argument('-s', '--species',  dest = 'species',
                          help = 'specify the genome species: human, mouse' + ' [Default: human]',
                          default='human')
    optional.add_argument('-u', '--umi_cutoff',  dest = 'umi_cutoff',
                          help = 'specify an integer to use as the UMI cutoff' + ' [Default: \'\']',
                          default='')
    optional.add_argument('-e', '--memory_size',  dest = 'memory_size',
                          help = 'decide how many memory cab be assigned for current project',
                          default='32')

    # ---------- setup | parse/check user options ---------- #
    p_args = parser.parse_args()
    if p_args.receptor_kit == 'TCRv1' or (p_args.receptor_kit == 'BCRv1' and p_args.species == 'mouse'):
        p_args.linker_correction = False
        p_args.umi_cutoff = ''
    if p_args.fastq_dir[len(p_args.fastq_dir)-1] == '/' or p_args.fastq_dir[len(p_args.fastq_dir)-1] == '\\':
        p_args.fastq_dir = p_args.fastq_dir[:-1]
    p_args.internal_test = False #True
    p_args.cogentIP_version = 'v1.6'
    p_args.process_threads = ''
    p_args.library_kit = 'takara_smartseq'
    p_args.rm_adapter = False
    p_args.vdjviz_inputs = False
    p_args.allow_partial = False
    
    if p_args.cogentIP_version != 'v2':
        p_args.airr_create = True
    else:
        p_args.airr_create = False
    if p_args.umi_cutoff == '':
        pass
    elif not p_args.umi_cutoff.isdigit():
        print('The umi_cutoff ' + str(p_args.umi_cutoff) + ' specified by user is not integer, please modify and relaunch')
        sys.exit(1)
    if p_args.species == 'mouse' and p_args.receptor_kit != 'TCRv2' :
        print('In the case of mouse species, support will only be provided for the receptor_kit TCRv2')
        sys.exit(1)
    if p_args.receptor_kit == 'BCRv2hy' and p_args.species == 'mouse':
        print('The receptor_kit ' + str(p_args.receptor_kit) + ' specified by human, please modify and relaunch')
        sys.exit(1)
    
    # p_args = p_args_indebug(p_args)
    # ---------- setup | start ---------- #
    vjava = required_python_module_check()
    start_time = datetime.datetime.now()
    # ---------- setup | end ---------- #
    if '/' in p_args.out_name or '\\' in p_args.out_name:
        print('out_name specified by -o should be a string, not directory path; ' + 
            'the string is served as user-defined name, a new folder with this ' + 
            'name would be created at the folder store metadata file')
        sys.exit(1)
    
    # --------- setup aligner | end ---------- #
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    
    # extract fastq fils by type of fastq file format
    basename = os.path.basename(p_args.fastq_dir)
    if basename.find('*')>=0 or basename.find('?')>=0:
        p_args.fastq_dir, p_args.meta_file = get_meta_file(p_args)

    p_args.main_dir = os.path.split(os.path.abspath(p_args.meta_file))[0]
    p_args.out_dir = os.path.join(p_args.main_dir, p_args.out_name)

    check_dir(p_args.fastq_dir, 'does not exist!\n')
    if len(p_args.meta_file) > 0:
        check_file(p_args.meta_file, 'does not exist!\n')

    # ---------- setup | load meta, create outdir ---------- #
    meta_dict = load_meta(p_args.meta_file, p_args.fastq_dir)
    print(get_time() + ' [INFO] Loading meta data file specified by user', flush=True)

    # ---------- delete the folder only for internal_test or debug ----- #
    debug = False #True False
    logger = logging.getLogger('temp_recoder')
    #debugfunction(p_args)
    if os.path.isdir(p_args.out_dir) and (p_args.internal_test or debug):
        logging.shutdown()
        shutil.rmtree(p_args.out_dir)

    if os.path.isdir(p_args.out_dir):
        print(get_time() + ' [INFO] Analysis dir already exists: ' + p_args.out_dir, flush=True)
        sys.exit(1)
    elif p_args.out_name == os.path.basename(os.path.normpath(p_args.main_dir)):
        print('The output folder name you defined is identical to its upper folder, please rename', flush=True)
        sys.exit(1)
    else:
        try:
            os.makedirs(p_args.out_dir)
        except OSError as err:
            print(get_time() + ' [ERROR] Unable to create directory: ' + p_args.out_dir, flush=True)
            sys.exit(1)

    # ---------- setup | logger ---------- #
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
    umi_length = 12
    linker_sequence = 'GTAC'
    solution = []
    delta = lambda x,y,i,j: 1 if x[i] != y[j] else 0
    
    #cwd = os.getcwd()
    #__file__ = os.path.join(cwd, 'immune_profiler.py')
    REPO_PATH = os.path.dirname(os.path.realpath(__file__))
    workflow_base = 'java -Xmx32G -jar ' + os.path.join(REPO_PATH, 'src', 'ImmuneProfiler'+p_args.cogentIP_version+'.jar')
    if p_args.memory_size != "32":
        workflow_base = 'java -Xmx' + p_args.memory_size + 'G -jar ' + os.path.join(REPO_PATH, 'src', 'ImmuneProfiler'+p_args.cogentIP_version+'.jar')
    if platform.system() == 'Windows' :
        workflow_base = 'java -jar ' + os.path.join(REPO_PATH, 'src', 'ImmuneProfiler'+p_args.cogentIP_version+'.jar')
        if len(p_args.process_threads) > 0 and int(p_args.process_threads) > 4:
            p_args.process_threads = '4'
    if p_args.cogentIP_version.find("v2") >= 0 :
        workflow_base = workflow_base[:-6] + 'v2.jar'
    workflow0_link = workflow_base + " " + p_args.out_dir
    workflow1_link = workflow_base + " " + p_args.out_dir
    workflow2_link = workflow_base + " " + p_args.out_dir
    
    if workflow2_link.find('ImmuneProfilerv2.jar') >= 0:
        p_args.airr_create = True
    
    align_stats_dict = OrderedDict([('Total sequencing reads', 'total reads'),
        ('Successfully aligned reads', 'aligned reads'), 
        ('Overlapped', 'overlapped'), 
        ('Overlapped and aligned', 'overlapped and aligned')])
    clone_stats_dict = OrderedDict([('Final clonotype count', 'clonotype count')])
    if p_args.internal_test:# and p_args.receptor_kit == 'TCRv2':
        clone_stats_dict = OrderedDict([('Average number of reads per clonotype', 'reads per clonotype'),
            ('Reads used in clonotypes, percent of total', 'Reads used in clonotypes'),
            ('Mapped low quality reads, percent of used', 'Mapped low quality reads'),
            ('Reads clustered in PCR error correction, percent of used', 'PCR error correction'),
            ('Reads dropped due to failed mapping, percent of total', 'failed mapping'),
            ('Final clonotype count', 'clonotype count')])

    # ---------- load bcs & receptor info -------- #
    create_dir(os.path.join(p_args.out_dir, 'report'))
    fd = OrderedDict()
    
    # ---------- preprocess -------- #
    #if (p_args.library_kit == 'takara_smartseq') and (workflow2_link.find('ImmuneProfiler.jar')>=0 or p_args.internal_test):
    if (p_args.library_kit == 'takara_smartseq'):
        create_dir(os.path.join(p_args.out_dir, 'preprocess'))
        logger.info('Starting preprocessing')
        print(get_time() + ' [INFO] Starting preprocessing', flush=True)
        sd = OrderedDict()
        # preprocess start
        prep_cmd = workflow0_link + ' workflow0 ' + p_args.receptor_kit + ' ' + p_args.fastq_dir \
                 + ' ' + p_args.meta_file + ' ' + p_args.out_dir + ' ' + p_args.species + ' ' + REPO_PATH
        if p_args.linker_correction : 
            prep_cmd += " -l"
        logger.info('Launching PREPROCESS using command:\n' + prep_cmd)
        log_name = 'run_prepro.log'
        out_log = os.path.join(p_args.out_dir, 'preprocess', log_name)
        run_prepro_err = os.path.join(p_args.out_dir, 'preprocess', 'run_prepro.error')
        run_system_cmd(prep_cmd, out_log, run_prepro_err, 'preprocess failed', log_name)
        sd = load_sd( os.path.join(p_args.out_dir, 'preprocess', 'sd_dict.csv') )
        # preprocess end
        write_sample_qc(p_args, sd)
        if p_args.airr_create and workflow2_link.find('ImmuneProfilerv1.6.jar')>=0:
            create_dir(os.path.join(p_args.out_dir, 'airr_report'))
            write_airr_sample_qc(p_args, sd)
        fd = create_fd(sd)
        logger.info('Completed preprocessing')
        print(get_time() + ' [INFO] Completed preprocessing', flush=True)
    
    # ---------- Run MIGEC ---------- #
    if workflow2_link.find('ImmuneProfilerv1.6.jar')>=0 and p_args.receptor_kit != 'TCRv1' and (p_args.receptor_kit != 'BCRv1' or p_args.species !='mouse'):
        create_dir(os.path.join(p_args.out_dir, 'run_migec'))
        logger.info('Starting UMI guided analysis…')
        print(get_time() + ' [INFO] Starting UMI guided analysis…', flush=True)
        migec_start_time = datetime.datetime.now()
        print(get_time() + ' [INFO] Preparing MIGEC input files', flush=True)
        if p_args.receptor_kit != "TCRv1" and (p_args.receptor_kit != "BCRv1" or p_args.species !="mouse"):
            prep_migec_input(p_args, fd)
        print(get_time() + ' [INFO] Conducting MIGEC analysis')
        run_migec(p_args, 'run_migec', 'mig_run_migec.log', 'mig_run_migec.error')
        print(get_time() + ' [INFO] Summarizing MIGEC results', flush=True)
        summ_migec(p_args)                     # *****
        logger.info('Finished collecting MIGEC statistics')
        logger.info('Finished MIGEC process')
        print(get_time() + ' [INFO] Finished MIGEC process', flush=True)
        
        # delete preprocessing folder as finished raw read MiXCR analysis & MIGEC analysis 
        if not p_args.keep_inter_file and not debug:
            pre_dir =os.path.join(p_args.out_dir, 'preprocess')
            pre_fqs = find_file('*.*', pre_dir) #'*.fastq'
            for fq in pre_fqs:
                os.remove(fq)
            os.rmdir(pre_dir)

    # ---------- Run MiXCR on MIGEC results or preprocess data---------- #
    if os.path.isdir(os.path.join(p_args.out_dir, 'run_migec')) or p_args.receptor_kit == 'BCRv2hy':
        if os.path.exists(os.path.join(p_args.out_dir, 'run_mixcr')) :
            change_run_mixcr_name(p_args, 'prePro')
        else: 
            create_dir(os.path.join(p_args.out_dir, 'run_mixcr'))
        read_type = 'mig'
        pre_dir = 'run_migec'                
        print(get_time() + ' [INFO] Starting reads alignment, assembling, and clonotype calling', flush=True)
        if p_args.target_region == 'Both' or p_args.target_region == 'CDR3':
            mig_cdr3_fd = copy.deepcopy(fd)
            if p_args.receptor_kit == "BCRv2hy" :
                run_mixcr_preIgGA(p_args, mig_cdr3_fd, read_type, 'cdr3', 'run_mixcr', pre_dir, 'run_mixcr_cdr3_sIg.log', 'run_mixcr_cdr3_sIg.error')
                change_run_mixcr_name(p_args, 'preSub')
            run_mixcr(p_args, mig_cdr3_fd, read_type, 'cdr3', 'run_mixcr', pre_dir, 'run_mixcr_cdr3.log', 'run_mixcr_cdr3.error')

        if p_args.target_region == 'Both' or p_args.target_region == 'Full_length':
            mig_fl_fd = copy.deepcopy(fd)
            if p_args.receptor_kit == "BCRv2hy" :
                run_mixcr_preIgGA(p_args, mig_fl_fd, read_type, 'fl', 'run_mixcr', pre_dir, 'run_mixcr_cdr3_sIg.log', 'run_mixcr_cdr3_sIg.error')
                change_run_mixcr_name(p_args, 'preSub')
            run_mixcr(p_args, mig_fl_fd, 'mig', 'fl', 'run_mixcr',
                  'run_migec', 'run_mixcr_fl.log', 'run_mixcr_fl.error')
        logger.info('Finish MiXCR analysis')
    
    # ---------- clean up & finish ---------- #
    if p_args.species == 'rhesus_monkey':
        os.remove(os.path.join(p_args.main_dir, 'imgt.201946-3.sv6.json'))
    else:
        pass
    logger.info('Analysis completed')
    print(get_time() + ' [INFO] Analysis completed', flush=True)
    sys.exit(0)
