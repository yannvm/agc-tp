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
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Vander Meersche Yann"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Yann Vander Meersche"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Vander Meersche Yann"
__email__ = "yann-vm@hotmail.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

#==============================================================
# Main program
#==============================================================
def read_fasta(amplicon_file, minseqlen):
    """Reads a compressed Fasta file.

    Arguments
    ---------
    amplicon_file: Path to the Fasta file

    Yield
    -----
    Generator yielding Fasta reads
    """
    if amplicon_file.endswith("gz"):
        with gzip.open(amplicon_file, "rb") as f_in:
            seq = b""
            for line in f_in:
                if line.startswith(b">"):
                    if len(seq) >= minseqlen:
                        yield seq.decode('ascii')
                    seq = b""
                else:
                    seq += line.strip()
            yield seq

    else:
        with open(amplicon_file, "r") as f_in:
            seq = ""
            for line in f_in:
                if line.startswith(">"):
                    if len(seq) >= minseqlen:
                        yield seq
                    seq = ""
                else:
                    seq += line.strip()
            yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """ Good Dockstring """
    gen_fasta = read_fasta(amplicon_file, minseqlen)

    seq_count = Counter()
    for seq in gen_fasta:
        seq_count[seq] += 1

    for count in seq_count.most_common():
        if count[1] >= mincount:
            yield count


def get_chunks(sequence, chunk_size):
    """ Good Dockstring """
    list_chunk = []
    for i in range(0, len(sequence), chunk_size):
        chunk = sequence[i:i+chunk_size]
        if len(chunk) == chunk_size:
            list_chunk.append(chunk)
    if len(list_chunk) >= 4:
        return list_chunk


def get_unique(ids):
    """ Good Dockstring """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    """ Good Dockstring """
    return list(set(lst1) & set(lst2))


def cut_kmer(sequence, kmer_size):
    """Cut a sequence in k-mers.

    Arguments
    ---------
    sequence: Fasta sequence
    kmer_size: K-mer size

    Yield
    -----
    Generator yielding Fasta sequence k-mers
    """
    for i in range(len(sequence) - kmer_size+1):
        yield sequence[i:i+kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """ Good Dockstring """
    gen_kmer = cut_kmer(sequence, kmer_size)

    for kmer in gen_kmer:
        if kmer not in kmer_dict:
            kmer_dict[kmer] = [id_seq]
        else:
            kmer_dict[kmer].append(id_seq)

    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """ Good Dockstring """
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size)
        if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list):
    """ Good Dockstring """
    nb_same = 0

    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            nb_same += 1

    return nb_same / len(alignment_list[0]) * 100


def detect_chimera(perc_identity_matrix):
    """ Good Dockstring """
    list_std = []
    bool_seq0 = False
    bool_seq1 = False

    for liste in perc_identity_matrix:
        list_std.append(statistics.stdev(liste))
        if liste[0] > liste[1]:
            bool_seq0 = True
        if liste[0] < liste[1]:
            bool_seq1 = True

    return bool(statistics.mean(list_std) > 5 and bool_seq0 and bool_seq1)



def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ Good Dockstring """

    list_nonchim = []
    kmer_dict = {}
    chimera = False
    id_kmer = 0

    for i, sequence in enumerate(dereplication_fulllength(amplicon_file, minseqlen, mincount)):

        if i%100 == 0:
            print(i)


        chunk_list = get_chunks(sequence[0], chunk_size)

        mates_list = []
        for chunk in chunk_list:
            mates_list.append(search_mates(kmer_dict, chunk, kmer_size))

        common_list = []
        for j in range(len(mates_list)):
            common_list = common(common_list, mates_list[j])


        if len(common_list) > 1:
            for c in common_list[0:2]:
                chunk_ref = get_chunks(list_nonchim[c], chunk_size)

                identity_matrix = [[]*4]
                for k in range(len(chunk_list)):
                    align = nw.global_align(chunk_list[k], chunk_ref[k],
                        gap_open=-1, gap_extend=-1, matrix="agc/MATCH")
                    identity_matrix[k].append(get_identity(align))

            chimera = detect_chimera(identity_matrix)

        if not chimera:
            kmer_dict = get_unique_kmer(kmer_dict, sequence[0], id_kmer, kmer_size)
            list_nonchim.append(sequence[0])
            id_kmer += 1
            yield sequence

"""
def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    OTU = []
    for sequence in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
        OTU.append(sequence)

    return OTU
"""

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ Good Dockstring """
    gen_chimerar_rem = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)

    OTU = [next(gen_chimerar_rem)]

    for sequence in gen_chimerar_rem:
        bool_id = True
        for i in range(len(OTU)):
            align = nw.global_align(OTU[i][0], sequence[0],
                gap_open=-1, gap_extend=-1, matrix="agc/MATCH")
            if get_identity(align) < 97:
                bool_id = False
                break

        if bool_id:
            OTU.append(sequence)

    return OTU




def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    """ Good Dockstring """
    with open(output_file, "wt") as f_out:
        for i, OTU in enumerate(OTU_list):
            f_out.write(f">OTU_{i+1} occurrence:{OTU[1]}\n")
            f_out.write(f"{fill(OTU[0])}\n")


def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen,
        args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(OTU_list, args.output_file)

if __name__ == '__main__':
    main()
