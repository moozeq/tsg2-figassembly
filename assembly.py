#!/usr/bin/env python3

import argparse
import itertools
from pprint import pprint
from typing import List, Dict


def read_records(filename: str):
    """Get records from file as Seqs objects"""
    from Bio import SeqIO
    print(f'[*] Reading sequences from {filename}')
    seqs = [record for record in SeqIO.parse(filename, 'fasta')]
    print(f'[+] Read {len(seqs)} sequences')
    return seqs


def read_sequences(filename: str) -> List[str]:
    """Get sequences from file as list of strings"""
    seqs = [str(record.seq) for record in read_records(filename)]
    return seqs


class Seq:
    seq_id = 0

    def __init__(self, seq: str):
        self.seq = seq
        self.id = Seq.seq_id
        Seq.seq_id += 1

    def _single_longest_prefix(self, sec_seq: 'Seq'):
        s1 = self.seq
        s2 = sec_seq.seq

        i = s1.find(s2[0])
        while i >= 0:
            if s2.startswith(s1[i:]):
                return s1[i:]

            i = s1.find(s2[0], i + 1)

        return ''

    def _single_longest_suffix(self, sec_seq: 'Seq'):
        s1 = sec_seq.seq
        s2 = self.seq

        i = s1.find(s2[0])
        while i >= 0:
            if s2.startswith(s1[i:]):
                return s1[i:]

            i = s1.find(s2[0], i + 1)

        return ''

    def longest_fix(self, seqs: Dict[int, 'Seq']) -> (int, str, str):
        lp = ''
        lp_seq_id = -1
        fix = ''

        for seq_id, seq in seqs.items():
            slp = self._single_longest_prefix(seq)
            if len(slp) > len(lp):
                lp = slp
                lp_seq_id = seq_id
                fix = 'prefix'

        for seq_id, seq in seqs.items():
            slp = self._single_longest_suffix(seq)
            if len(slp) > len(lp):
                lp = slp
                lp_seq_id = seq_id
                fix = 'suffix'

        return lp_seq_id, lp, fix

    def merged_best_seq(self, seqs: Dict[int, 'Seq']) -> 'Seq':
        best_seq_id, longest_prefix, fix = self.longest_fix(seqs)
        if best_seq_id != -1:
            best_seq = seqs[best_seq_id]
            self.id = best_seq_id

        if fix == 'prefix':
            self.seq = f'{self.seq}{best_seq.seq[len(longest_prefix):]}'
        elif fix == 'suffix':
            self.seq = f'{best_seq.seq[:len(longest_prefix)]}{self.seq}'
        else:
            print(f'[*] Skip seq with length = {len(self.seq)}')

        return self


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assembly sequences')
    parser.add_argument('file', type=str, help='fasta file with reads')
    parser.add_argument('-n', '--num', type=int, default=10, help='k for greedy scs')
    parser.add_argument('-o', '--output', type=str, help='output fasta')
    args = parser.parse_args()

    seqs = [Seq(s) for s in read_sequences(args.file)]
    seqs = {s.id: s for s in seqs}

    run = 1
    while len(seqs) > 1:
        p_seqs_len = len(seqs)
        seqs_keys = list(seqs.keys())
        print(f'[*] Run = {run}; keys = {len(seqs_keys)}')
        for i, s_id in enumerate(seqs_keys):
            if i & 0xf == 0:
                print(f'[*] Iteration = {i}')
            if s_id in seqs:
                start_seq_obj: Seq = seqs.pop(s_id)
                seq_obj: Seq = start_seq_obj.merged_best_seq(seqs)
                seqs[seq_obj.id] = seq_obj
        run += 1
        if p_seqs_len == len(seqs):  # nothing changed, merge all left sequences
            break

    super_seq = ''.join(s.seq for s in seqs.values())
    print(f'[+] Done, super sequence length = {len(super_seq)}')

    filename = args.output if args.output else f'training/super.fasta'

    with open(filename, 'w') as rf:
        rf.write('>a\n')
        rf.write(super_seq)
