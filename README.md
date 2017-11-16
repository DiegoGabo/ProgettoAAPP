# ProgettoAAPP

## Description
We developed a software that counts the number of occurrences of every k-mer in a long string. With k-mer we mean a substring of 
length k.

## Documentation
You can read the file presentation in which there are all details about our implementation and the results of our tests.
You can also read the paper "A fast, lock-free approach for efficient parallel counting of occurrences of k-mers" written by
Guillaume Marçais1, and Carl Kingsford  in which there is the report of a similar algorithm.

## Datasets used
We used several datasets to test our implementation. They consist in real DNA sequences taken from the National Center for Biotechnology Information. Here there are our datasets that can be found in the folder dna_sequences:
* DNA randomly generated;
* gbbct155.seq: complete genome of Escherichia Coli;
* gbbct367.seq: complete genome of Lactococcus lactis;
* gbgss116.seq: Chlorocebus aethiops genomic clone CH252-491O17, genomic surveey sequence;
* hs_alt_CHM1_1.1_chr22.gbk: Homo sapiens chromosome 22 genomic scaffold.


## Versions
We developed three versions of the program:
* Serial version (folder 1_serial_version);
* Lock-parallel version (folder 2_parallel_lock_version);
* Optimized-parallel version (folder 3_parallel_version_with_optimisation).

### Serial Version
Implementation of an hash table with entries key-values that represent the following features:
* Key = the substring corresponding to the k-mer;
* Value = the number of times the k-mer appeasr in the input string.

### Lock-parallel version
We obtained parallelism using OpenMp. The main cycle is executed in parallel by different threads (SIMD parallelism)
and conflicts are solved using the critical openMP pragma.

### Optimized-parallel version
We optimized the Lock-parallel version using CAS instruction to improve performances and applying other optimisations.

## Tests
Tests can be seen in the folder test and in the presentation. We tested the following things:
* execution time of all version with differents datasets;
* execution time of all version with increasing length of the dataset;
* execution time of all version with increasing K;
* execution time of all version with increasing number of threads;
* memory occupied with respect to dataset length;
* memory occupied with respect to K.
