#!/usr/bin/env python3

#-------------------------------------------------------------------------------
#   align.py: alignment functions for genotyping.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   This file is part of arcasHLA.
#
#   arcasHLA is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   arcasHLA is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with arcasHLA.  If not, see <https://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------

import os
import sys
import re
import json
import pickle
import argparse
import logging as log

import pandas as pd

from datetime import date
from argparse import RawTextHelpFormatter
from textwrap import wrap
from collections import Counter, defaultdict
from itertools import combinations

from reference import check_ref, get_exon_combinations
from arcas_utilities import *

#-------------------------------------------------------------------------------

__version__     = '0.4.0'
__date__        = '2022-01-27'

#-------------------------------------------------------------------------------
#   Paths and filenames
#-------------------------------------------------------------------------------

rootDir = os.path.dirname(os.path.realpath(__file__)) + '/../'

#-----------------------------------------------------------------------------
# Process and align FASTQ input
#-----------------------------------------------------------------------------

def pseudoalign(fqs, sample, paired, reference, outdir, temp, threads, avg, std):
    '''Calls Kallisto to pseudoalign reads.'''

    log.info('[alignment] Analyzing read length')

    awk = "| awk '{if(NR%4==2) print length($1)}'"

    if fqs[0].endswith('.gz'):
        cat = 'zcat'
    else:
        cat = 'cat'

    # Get read length stats
    reads_file = ''.join([temp, sample, '.reads.txt'])

    for fq in fqs:
        command = [cat, '<', fq, awk, '>>', reads_file]
        run_command(command)

    num = 0
    with open(reads_file, 'r') as f:
        for _ in f:
            num += 1
    if num == 0:
        sys.exit('[genotype] Error: FASTQ files are empty; check arcasHLA extract for issues.')

    command = ['kallisto pseudo -i', reference, '-t', threads, '-o', temp]

    if not paired:
        command.extend(['--single -l', str(avg), '-s', str(std)])

    command.extend(fqs)
    run_command(command, '[alignment] Pseudoaligning with Kallisto: ')

    return num

#-----------------------------------------------------------------------------
# Process transcript assembly output
#-----------------------------------------------------------------------------

def process_counts(count_file, eq_file, gene_list, allele_idx, allele_lengths): 
    '''Processes pseudoalignment output, returning compatibility classes.'''
    log.info('[alignment] Processing pseudoalignment')
    # Process count information
    counts = dict()
    with open(count_file, 'r', encoding='UTF-8') as file:
        for line in file:
            line = line.rstrip('\n')
            eq, count = line.split('\t')
            counts[eq] = float(count)


    # Process compatibility classes
    eqs = dict()
    with open(eq_file, 'r', encoding='UTF-8') as file:
        for line in file:
            line = line.rstrip('\n')
            eq, indices = line.split('\t')
            eqs[eq] = indices.split(',')

    # Set up compatibility class index
    eq_idx = defaultdict(list)

    count_unique = 0
    count_multi = 0

    # Pre-compute gene for each index to avoid repeated string splits
    idx_gene = {}
    for idx, alleles in allele_idx.items():
        if alleles:
            idx_gene[idx] = alleles[0].split('*')[0]

    for eq, indices in eqs.items():
        if any(idx not in idx_gene for idx in indices):
            continue

        # Determine gene set with early exit on multi-gene
        gene = None
        multi_gene = False
        for idx in indices:
            g = idx_gene[idx]
            if gene is None:
                gene = g
            elif g != gene:
                multi_gene = True
                break

        count = counts[eq]

        if not multi_gene and count > 0:
            eq_idx[gene].append((indices, count))
            count_unique += count
        else:
            count_multi += count

    # Alleles mapping to their respective compatibility classes
    allele_eq = defaultdict(set)
    for eqs in eq_idx.values():
        for eq,(indices,_) in enumerate(eqs):
            for idx in indices:
                allele_eq[idx].add(eq)
    
    return eq_idx, allele_eq, [count_unique, count_multi]

def process_partial_counts(count_file, eq_file, allele_idx, allele_lengths, 
                           exon_idx, exon_combos):
    '''Processes pseudoalignment output, returning compatibility classes.'''
    
    log.info('[alignment] Processing pseudoalignment')
    counts_index = dict()
    with open(count_file,'r', encoding='UTF-8') as file:
        for line in file:
            line = line.rstrip('\n')
            eq, count = line.split('\t')
            counts_index[eq] = float(count)

    eqs = dict()
    with open(eq_file,'r', encoding='UTF-8') as file:
        for line in file:
            line = line.rstrip('\n')
            eq, indices = line.split('\t')
            eqs[eq] = indices.split(',')

    eq_idx = {str(i):defaultdict(list) for i in exon_combos}
    count_unique = 0
    count_multi = 0

    # Pre-compute gene for each index to avoid repeated string splits
    idx_gene = {}
    for idx, alleles in allele_idx.items():
        if alleles:
            idx_gene[idx] = alleles[0].split('*')[0]

    for eq, indices in eqs.items():
        if any(index not in idx_gene for index in indices):
            continue

        # Determine gene set with early exit on multi-gene
        gene = None
        multi_gene = False
        for index in indices:
            g = idx_gene[index]
            if gene is None:
                gene = g
            elif g != gene:
                multi_gene = True
                break

        count = counts_index[eq]

        if not multi_gene and count > 0:
            # Single-pass exon grouping with defaultdict
            exon_groups = defaultdict(list)
            for index in indices:
                exon_groups[exon_idx[index]].append(index)

            for exon, exon_indices in exon_groups.items():
                eq_idx[exon][gene].append((exon_indices, count))
            count_unique += count
        else:
            count_multi += count
                
    return eq_idx, [count_unique, count_multi]
           

def get_count_stats(eq_idx, gene_length):
    '''Returns counts and relative abundance of genes.'''
    stats = {gene:[0,0,0.] for gene in eq_idx}

    abundances = defaultdict(float)
    for gene, eqs in eq_idx.items():
        count = sum(count for eq,count in eqs)
        abundances[gene] = count / gene_length[gene]
        stats[gene][0] = count
        stats[gene][1] = len(eqs)
        
    total_abundance = sum(abundances.values())

    for gene, abundance in abundances.items():
        stats[gene][2] = abundance / total_abundance

    return stats

def alignment_summary(align_stats, partial = False):
    '''Prints alignment summary to log.'''

    count_unique, count_multi, total, _, _ = align_stats
    log.info('[alignment] Processed {:.0f} reads, {:.0f} pseudoaligned '
             .format(total, count_unique + count_multi)+
             'to HLA reference')
              
    log.info('[alignment] {:.0f} reads mapped to a single HLA gene'
             .format(count_unique))

def gene_summary(gene_stats):
    '''Prints gene read count and relative abundance to log.'''

    log.info('[alignment] Observed HLA genes:')

    log.info('\t\t{: <10}    {}    {}    {}'
             .format('gene','abundance','read count','classes'))

    for g,(c,e,a) in sorted(gene_stats.items()):
        log.info('\t\tHLA-{: <6}    {: >8.2f}%    {: >10.0f}    {: >7.0f}'
                 .format(g, a*100, c, e))

def get_alignment(fqs, sample, reference, reference_info, outdir,
        temp, threads, single, partial = False, avg = 200, std = 20):
    '''Runs pseudoalignment and processes output.'''
    paired = not single
        
    count_file = ''.join([temp, 'pseudoalignments.tsv'])
    eq_file = ''.join([temp, 'pseudoalignments.ec'])

    total = pseudoalign(fqs,
                sample,
                paired,
                reference,
                outdir,
                temp,
                threads,
                avg,
                std)
    
    # Process partial genotyping pseudoalignment
    if partial:
        (commithash, (gene_set, allele_idx, exon_idx,
            lengths, partial_exons, partial_alleles)) = reference_info

        gene_set = set(gene_set)
        lengths = {int(a): int(x) for a, x in lengths.items()}
        partial_alleles = set(partial_alleles)
        
        exon_combos = get_exon_combinations() 
        
        eq_idx, align_stats = process_partial_counts(count_file,
                                               eq_file,
                                               allele_idx, 
                                               lengths,
                                               exon_idx,
                                               exon_combos)
        align_stats.extend([total, avg, std])
        
        alignment_summary(align_stats, True)
        
        with open(''.join([outdir,sample,'.partial_alignment.p']),'wb') as file:
            alignment_info = [commithash, eq_idx, [], paired, 
                              align_stats, []]
            pickle.dump(alignment_info, file)
            
    # Process regular pseudoalignment
    else:
        (commithash, (gene_set, allele_idx,
             lengths, gene_length)) = reference_info

        gene_set = set(gene_set)
        gene_length = {a: int(x) for a, x in gene_length.items()}
        lengths = {int(a): int(x) for a, x in lengths.items()}

        eq_idx, allele_eq, align_stats = process_counts(count_file,
                                                        eq_file, 
                                                        gene_set, 
                                                        allele_idx, 
                                                        lengths)
        
        align_stats.extend([total, avg, std])
        
        alignment_summary(align_stats)
        
        gene_stats = get_count_stats(eq_idx, gene_length)
        gene_summary(gene_stats)

        with open(''.join([outdir, sample, '.alignment.p']), 'wb') as file:
            alignment_info = [commithash, eq_idx, allele_eq, paired,
                    align_stats, gene_stats]
            pickle.dump(alignment_info, file)
            
        with open(''.join([outdir,sample,'.genes.json']), 'w') as file:
            json.dump(gene_stats, file)
            
    return alignment_info

def load_alignment(file, commithash, partial = False):
    '''Loads previous pseudoalignment.'''
    
    log.info('[alignment] Loading previous alignment %s', file)
    
    with open(file, 'rb') as file:
        alignment_info = pickle.load(file)
    
    # Compatibility with arcasHLA 1.0
    if len(alignment_info) != 5:
        if partial:
            (commithash_alignment, eq_idx, paired,
                 _, _, _, _) = alignment_info
            alignment_info = [commithash, eq_idx, None, paired, 
                              None, None]
        else:
            (commithash_alignment, eq_idx, allele_eq, paired, 
                 align_stats, gene_stats, num, avg, std) = alignment_info
            align_stats.extend([num, avg, std])
            alignment_info = [commithash, eq_idx, allele_eq, paired, 
                              align_stats, gene_stats]
    
    commithash_alignment, _,_,_, align_stats, gene_stats = alignment_info
        
    if commithash != commithash_alignment:
        sys.exit('[alignment] Error: reference used for alignment ' +
                 'different than the one in the database')
        
    if align_stats: alignment_summary(align_stats)
    if not partial: gene_summary(gene_stats)
        
    return alignment_info

#-----------------------------------------------------------------------------
