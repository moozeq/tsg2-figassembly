#!/usr/bin/env python3

import argparse
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


def sim_score(s1: str, s2: str) -> float:
    score = 0.0
    for l1, l2 in zip(s1, s2):
        if l1 == l2:
            score += 1.0
    return score / len(s1)


def sim_exact(s1: str, s2: str) -> float:
    return 1.0 if s1 == s2 else 0.0


SIM_SCORE = 1.0
SIM_FUNC = sim_exact
START_SEQ_LEN = 0
FIX_LONG = 4


class Seq:
    seq_id = 0

    def __init__(self, seq: str):
        self.seq = seq
        self.id = Seq.seq_id
        Seq.seq_id += 1

    @staticmethod
    def get_longest_fix(s1: str, s2: str):
        i = s1.find(s2[0])
        while i >= 0:
            s1_post = s1[i:]
            s2_pre = s2[:len(s1_post)]
            if SIM_FUNC(s2_pre, s1_post) >= SIM_SCORE:
                return s1_post

            i = s1.find(s2[0], i + 1)

        return ''

    @staticmethod
    def _longest_fix(sel_seq: 'Seq', seqs: Dict[int, 'Seq']) -> (int, str, str):
        lp_seq_id = -1
        lp = ''
        fix = ''

        for seq_id, seq in seqs.items():
            slp = Seq.get_longest_fix(sel_seq.seq, seq.seq)
            if len(slp) > len(lp):
                lp = slp
                lp_seq_id = seq_id
                fix = 'prefix'

        for seq_id, seq in seqs.items():
            slp = Seq.get_longest_fix(seq.seq, sel_seq.seq)
            if len(slp) > len(lp):
                lp = slp
                lp_seq_id = seq_id
                fix = 'suffix'

        return lp_seq_id, lp, fix

    def merged_best_seq(self, seqs: Dict[int, 'Seq']) -> 'Seq':
        best_seq_id, longest_fix, fix = Seq._longest_fix(self, seqs)
        if len(longest_fix) < START_SEQ_LEN // FIX_LONG or best_seq_id < 0:
            return self

        best_seq = seqs[best_seq_id]
        self.id = best_seq_id

        print(f'Longest fix: {fix}, {len(longest_fix)}')
        if fix == 'prefix':
            self.seq = f'{self.seq}{best_seq.seq[len(longest_fix):]}'
        elif fix == 'suffix':
            self.seq = f'{best_seq.seq[:-len(longest_fix)]}{self.seq}'
        else:
            raise Exception('Wrong path')

        return self


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assembly sequences')
    parser.add_argument('file', type=str, help='fasta file with reads')
    parser.add_argument('-s', '--score', type=float, default=1.0, help='similarity score')
    parser.add_argument('-o', '--output', type=str, help='output fasta')
    args = parser.parse_args()

    seqs = [Seq(s) for s in read_sequences(args.file)]
    START_SEQ_LEN = len(seqs[0].seq)
    seqs = {s.id: s for s in seqs}

    run = 1
    while len(seqs) > 1:
        seqs_keys = list(seqs.keys())
        print(f'[*] Run = {run}; keys = {len(seqs_keys)}')
        for k, s_id in enumerate(seqs_keys):
            if k & 0xff == 0:
                print(f'[*] Done keys = {k}, seqs left = {len(seqs)}')
                pass
            if s_id in seqs:
                start_seq_obj: Seq = seqs.pop(s_id)
                seq_obj = start_seq_obj.merged_best_seq(seqs)
                seqs[seq_obj.id] = seq_obj
        run += 1
        SIM_SCORE *= 0.8
        SIM_FUNC = sim_score

    super_seq = ''.join(s.seq for s in seqs.values())
    print(f'[+] Done, super sequence length = {len(super_seq)}')

    filename = args.output if args.output else f'training/super.fasta'

    with open(filename, 'w') as rf:
        rf.write('>super_sequence\n')
        rf.write(super_seq)
