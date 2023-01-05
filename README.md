# TSG2 Assembly project

## Description

**F**ast-**I**terative-**G**reedy **A**ssembly algorithm for turning multiple FASTA reads into
one super sequence.

## Usage

Script requires only `Python >= 3.8` and does not use any
external libraries.

```bash
# assembly
./assembly.py training/reads/reads1.fasta -o training/super1.fasta
./assembly.py training/reads/reads2.fasta -o training/super2.fasta
./assembly.py training/reads/reads3.fasta -o training/super3.fasta
```

## Algorithm

Algorithm in each iteration is going through all loaded
sequences, finding the longest prefix or suffix, then merging them
in place.

First run is performed finding exact -fix of specified minimum length
without any miss-matches.

In every next iteration, similarity function is used to finding the
longest -fix, loosening similarity threshold at the same time in each run.

## Results

On training data following results were achieved:

|data|assembly length [bp]|assembly time [s]|total score|
|:---:|:---:|:---:|:---:|
|reads1.fasta|4122|116|0.112|
|reads2.fasta|6247|136|0.116|
|reads3.fasta|8205|170|0.067|

Used commands:

```bash
# evaluation
./evaluate.sh super1.fasta
./evaluate.sh super2.fasta
./evaluate.sh super3.fasta
```

```bash
# measure assembly time
python3 -m cProfile -s time ./assembly.py training/reads/reads1.fasta -o training/super1.fasta
python3 -m cProfile -s time ./assembly.py training/reads/reads2.fasta -o training/super2.fasta
python3 -m cProfile -s time ./assembly.py training/reads/reads3.fasta -o training/super3.fasta
```
