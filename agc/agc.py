#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw
import numpy as np
np.int = int

import sys

__author__ = "FAUCHOIS Antoine"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["FAUCHOIS Antoine"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "FAUCHOIS Antoine"
__email__ = "antoine.fauchois92@gmail.com"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """    
    with gzip.open(amplicon_file, "rt") as fasta_input:
        #Initialize sequence
        sequence = ""
        #Looking length for each sequence
        for ligne in fasta_input:
            if ligne.startswith(">"):
                if len(sequence) >= minseqlen:
                    yield sequence
                sequence = ""
            else:
                sequence += ligne.strip()
        #Return the last sequence        
        if len(sequence) >= minseqlen:
            yield sequence

def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    #Call generator from read_fasta
    reads = read_fasta(amplicon_file, minseqlen)
    #Get occurence
    occurence = Counter(reads)
    #Parse occurrence
    for sequence, count in occurence.most_common():
        if count >= mincount:
            yield [sequence, count] 
    

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    identical_base = sum([1 for i in range(len(alignment_list[0])) if alignment_list[0][i] == alignment_list[1][i]])
    identical_rate = identical_base / len(alignment_list[0]) * 100
    return identical_rate

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int = 50, kmer_size: int = 22) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    #Call dereplicated reads
    derepli_reads = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    #Build OTU
    otu_list = [] 
    for derepli_read in list(derepli_reads):
        is_new_otu = True
        for j, otu in enumerate(otu_list):
            #Extract query and target sequence
            query_seq = derepli_read[0]
            target_seq = otu[0]
            #Make a global alignment
            alignment = nw.global_align(query_seq, target_seq, 
            gap_open = -3, gap_extend = -1, matrix = str(Path(__file__).parent / "MATCH"))
            #Compute identical rate
            identical_rate = get_identity(alignment)
            #Add occurence if there is a identity
            if identical_rate >= 97:
                otu_list[j][1] += derepli_read[1]
                is_new_otu = False
                break
        if is_new_otu:
            otu_list.append(derepli_read)
    return otu_list   
            


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w") as fasta_output:
        for i, otu in enumerate(OTU_list):
            fasta_output.write(f">OTU_{i + 1} occurrence:{otu[1]}\n")
            fasta_output.write(f"{textwrap.fill(otu[0], width=80)}\n")


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    #get reads and their occurence from fasta file
    OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount)

    write_OTU(OTU_list, args.output_file)

if __name__ == '__main__':
    main()
