import argparse
import itertools
from typing import List


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1  # move just past previous match


def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss) - 1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i + 1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i + 1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest


def pick_maximal_overlap(reads, k):
    """ Return a pair of reads from the list with a
        maximal suffix/prefix overlap >= k.  Returns
        overlap length 0 if there are no such overlaps."""
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen


def greedy_scs(reads, k):
    """ Greedy shortest-common-superstring merge.
        Repeat until no edges (overlaps of length >= k)
        remain. """
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[-(len(read_b) - olen):])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assembly sequences')
    parser.add_argument('file', type=str, help='fasta file with reads')
    parser.add_argument('-n', '--num', type=int, default=20, help='k for greedy scs')
    args = parser.parse_args()

    seqs = read_sequences(args.file)
    print(len(gscs := greedy_scs(seqs, args.num)))

    with open(f'training/r{args.num}.fasta', 'w') as rf:
        rf.write('>aaa\n')
        rf.write(gscs)


"""
Results:

=========================================
r5.fasta:

1 reads; of these:
  1 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    1 (100.00%) aligned >1 times
100.00% overall alignment rate
Pokrycie referencji: 0.371875
Pokrycie odczytów: 0.0514999913445
Ocena identyczności: 0.5
Ocena rozdrobnienia: 0.732486760359
Łączna ocena: 0.00701413180687
=========================================

=========================================
r10.fasta:

1 reads; of these:
  1 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    1 (100.00%) aligned >1 times
100.00% overall alignment rate
Pokrycie referencji: 0.445375
Pokrycie odczytów: 0.0601522799791
Ocena identyczności: 0.5
Ocena rozdrobnienia: 0.732486760359
Łączna ocena: 0.00981177797392
=========================================

=========================================
r15.fasta:

1 reads; of these:
  1 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    1 (100.00%) aligned >1 times
100.00% overall alignment rate
Pokrycie referencji: 0.093375
Pokrycie odczytów: 0.0124928922634
Ocena identyczności: 0.5
Ocena rozdrobnienia: 1.0
Łączna ocena: 0.000583261907549
=========================================

=========================================
r20.fasta:

1 reads; of these:
  1 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    1 (100.00%) aligned >1 times
100.00% overall alignment rate
Pokrycie referencji: 0.1585
Pokrycie odczytów: 0.0209375670811
Ocena identyczności: 0.5
Ocena rozdrobnienia: 0.898244401704
Łączna ocena: 0.00149045890396
=========================================
"""
