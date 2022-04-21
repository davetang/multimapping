Minimap2 usage.

    minimap2 || true

    ## Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]
    ## Options:
    ##   Indexing:
    ##     -H           use homopolymer-compressed k-mer (preferrable for PacBio)
    ##     -k INT       k-mer size (no larger than 28) [15]
    ##     -w INT       minimizer window size [10]
    ##     -I NUM       split index for every ~NUM input bases [4G]
    ##     -d FILE      dump index to FILE []
    ##   Mapping:
    ##     -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [0.0002]
    ##     -g NUM       stop chain enlongation if there are no minimizers in INT-bp [5000]
    ##     -G NUM       max intron length (effective with -xsplice; changing -r) [200k]
    ##     -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]
    ##     -r NUM[,NUM] chaining/alignment bandwidth and long-join bandwidth [500,20000]
    ##     -n INT       minimal number of minimizers on a chain [3]
    ##     -m INT       minimal chaining score (matching bases minus log gap penalty) [40]
    ##     -X           skip self and dual mappings (for the all-vs-all mode)
    ##     -p FLOAT     min secondary-to-primary score ratio [0.8]
    ##     -N INT       retain at most INT secondary alignments [5]
    ##   Alignment:
    ##     -A INT       matching score [2]
    ##     -B INT       mismatch penalty (larger value for lower divergence) [4]
    ##     -O INT[,INT] gap open penalty [4,24]
    ##     -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]
    ##     -z INT[,INT] Z-drop score and inversion Z-drop score [400,200]
    ##     -s INT       minimal peak DP alignment score [80]
    ##     -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]
    ##   Input/Output:
    ##     -a           output in the SAM format (PAF by default)
    ##     -o FILE      output alignments to FILE [stdout]
    ##     -L           write CIGAR with >65535 ops at the CG tag
    ##     -R STR       SAM read group line in a format like '@RG\tID:foo\tSM:bar' []
    ##     -c           output CIGAR in PAF
    ##     --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]
    ##     --MD         output the MD tag
    ##     --eqx        write =/X CIGAR operators
    ##     -Y           use soft clipping for supplementary alignments
    ##     -t INT       number of threads [3]
    ##     -K NUM       minibatch size for mapping [500M]
    ##     --version    show version number
    ##   Preset:
    ##     -x STR       preset (always applied before other options; see minimap2.1 for details) []
    ##                  - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
    ##                  - map-hifi - PacBio HiFi reads vs reference mapping
    ##                  - ava-pb/ava-ont - PacBio/Nanopore read overlap
    ##                  - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    ##                  - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
    ##                  - sr - genomic short-read mapping
    ## 
    ## See `man ./minimap2.1' for detailed description of these and other advanced command-line options.

Generate 100 unique reads with read lengths of 100 bp and repeat each
unique read incrementally, i.e. first read is unique, second read is
repeated twice, third read is repeated thrice, and so on. Next use these
reads to generate a reference sequence where each read is separated by a
spacer sequence that is two times the length of a read, which is 200 bp
in this case.

> Without any options, minimap2 takes a reference database and a query
> sequence file as input and produce approximate mapping, without
> base-level alignment (i.e. coordinates are only approximate and no
> CIGAR in output).

    l=100
    r=100
    ../script/create_ref_se.pl -l ${l} -r ${r}
    minimap2 -a l${l}_r${r}_ref.fa l${l}_r${r}_reads.fa 2> /dev/null > default.sam

    ## FASTA reads written to l100_r100_reads.fa
    ## FASTA reference written to l100_r100_ref.fa
    ## Done

Multimapping reads are reported with one read set as the primary
alignment and the others set as secondary alignment/s.

    head default.sam

    ## @SQ  SN:ref  LN:1515200
    ## @PG  ID:minimap2 PN:minimap2 VN:2.24-r1122   CL:minimap2 -a l100_r100_ref.fa l100_r100_reads.fa
    ## 1    0   ref 201 60  100M    *   0   0   AGAAAGGCGATTGGGATAAAAGACACAACCTTGTCGTCACGATGCGGTCCAGTTTCTGGCACAATCCCTAGTTCGCGCTAGGTCGTCGTAATGAGCTCAT    *   NM:i:0  ms:i:200    AS:i:200    nn:i:0  tp:A:P  cm:i:13 s1:i:87 s2:i:0  de:f:0  rl:i:0
    ## 2    0   ref 501 0   100M    *   0   0   GGATTGGGTCGGAATTTTCCTGCGGAGGTCCCTCCGCACAGGGCCGGCAACTGGCGGTTCTGAAGCACCCCTCCAGATTGTTTAAGCAGCCTGAAACGCA    *   NM:i:0  ms:i:200    AS:i:200    nn:i:0  tp:A:P  cm:i:14 s1:i:85 s2:i:85 de:f:0  rl:i:0
    ## 2    256 ref 801 0   100M    *   0   0   *   *   NM:i:0  ms:i:200    AS:i:200    nn:i:0  tp:A:S  cm:i:14 s1:i:85 de:f:0  rl:i:0
    ## 3    0   ref 1701    0   100M    *   0   0   CCAGCCACATAAACCATTGCTGTGCAGCAGGGGCGTGAACCGGCGTGCGTACTGAGACCCGAAGTATACTTTTATATTTGAACGCTCTGTGTGCGTCTTG    *   NM:i:0  ms:i:200    AS:i:200    nn:i:0  tp:A:P  cm:i:15 s1:i:91 s2:i:91 de:f:0  rl:i:0
    ## 3    256 ref 1401    0   100M    *   0   0   *   *   NM:i:0  ms:i:200    AS:i:200    nn:i:0  tp:A:S  cm:i:15 s1:i:91 de:f:0  rl:i:0
    ## 3    256 ref 1101    0   100M    *   0   0   *   *   NM:i:0  ms:i:200    AS:i:200    nn:i:0  tp:A:S  cm:i:15 s1:i:91 de:f:0  rl:i:0
    ## 4    0   ref 2001    0   100M    *   0   0   CCGAGGCCCACATGACCTCGTGAACTCTAAGAGGGCTTCTTAAGTAGGCTGCTCGGGGACTTAGACCAATACCCGTGGGAGATCGATTTGCGCGAATTAC    *   NM:i:0  ms:i:200    AS:i:200    nn:i:0  tp:A:P  cm:i:19 s1:i:97 s2:i:97 de:f:0  rl:i:0
    ## 4    256 ref 2901    0   100M    *   0   0   *   *   NM:i:0  ms:i:200    AS:i:200    nn:i:0  tp:A:S  cm:i:19 s1:i:97 de:f:0  rl:i:0

The default setting is to keep only 5 secondary alignments.

> -N INT retain at most INT secondary alignments \[5\]

    cat default.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2n | head

    ##       1 1
    ##       2 2
    ##       3 3
    ##       4 4
    ##       5 5
    ##       6 6
    ##       6 7
    ##       6 8
    ##       6 9
    ##       6 10

Report up to 100 secondary alignments.

    l=100
    r=100
    ../script/create_ref_se.pl -l ${l} -r ${r}
    minimap2 -N ${r} -a l${l}_r${r}_ref.fa l${l}_r${r}_reads.fa 2> /dev/null \
       | grep -v "^@" \
       | cut -f1 \
       | sort \
       | uniq -c \
       | sort -k2rn \
       | head

    ## FASTA reads written to l100_r100_reads.fa
    ## FASTA reference written to l100_r100_ref.fa
    ## Done
    ##     100 100
    ##      99 99
    ##      98 98
    ##      97 97
    ##      96 96
    ##      95 95
    ##      94 94
    ##      93 93
    ##      92 92
    ##      91 91

Report up to 1,000 secondary alignments.

    l=100
    r=1000
    ../script/create_ref_se.pl -l ${l} -r ${r}
    minimap2 -N ${r} -a l${l}_r${r}_ref.fa l${l}_r${r}_reads.fa 2> /dev/null \
       | grep -v "^@" \
       | cut -f1 \
       | sort \
       | uniq -c \
       | sort -k2rn \
       | head

    ## FASTA reads written to l100_r1000_reads.fa
    ## FASTA reference written to l100_r1000_ref.fa
    ## Done
    ##    1000 1000
    ##     999 999
    ##     998 998
    ##     997 997
    ##     996 996
    ##     995 995
    ##     994 994
    ##     993 993
    ##     992 992
    ##     991 991

Report up to 10,000 secondary alignments fails (and returns with an
error but `rmarkdown::render` continues).

    l=100
    r=10000
    ../script/create_ref_se.pl -l ${l} -r ${r}
    minimap2 -t 8 -N ${r} -a l${l}_r${r}_ref.fa l${l}_r${r}_reads.fa

    echo $?

    ## FASTA reads written to l100_r10000_reads.fa
    ## FASTA reference written to l100_r10000_ref.fa
    ## Done
    ## minimap2: bseq.c:95: mm_bseq_read3: Assertion `ks->seq.l <= INT32_MAX' failed.
    ## bash: line 3:   179 Aborted                 (core dumped) minimap2 -t 8 -N ${r} -a l${l}_r${r}_ref.fa l${l}_r${r}_reads.fa
    ## 134

Report up to 5,000 secondary alignments also fails (but with with a
return code of 0) and all reads are unmapped.

    l=100
    r=5000
    ../script/create_ref_se.pl -l ${l} -r ${r}
    minimap2 -t 8 -N ${r} -a l${l}_r${r}_ref.fa l${l}_r${r}_reads.fa > r${r}.sam

    echo $?

    head r${r}.sam

    ## FASTA reads written to l100_r5000_reads.fa
    ## FASTA reference written to l100_r5000_ref.fa
    ## Done
    ## [WARNING][1;31m failed to parse the first FASTA/FASTQ record. Continue anyway.[0m
    ## [M::mm_idx_gen::5.240*1.00] collected minimizers
    ## [M::mm_idx_gen::5.242*1.00] sorted minimizers
    ## [M::main::5.242*1.00] loaded/built the index for 0 target sequence(s)
    ## [M::mm_mapopt_update::5.242*1.00] mid_occ = 10
    ## [M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 0
    ## [M::mm_idx_stat::5.242*1.00] distinct minimizers: 0 (-nan% are singletons); average occurrences: -nan; average spacing: -nan; total length: 0
    ## [M::worker_pipeline::5.259*1.01] mapped 5000 sequences
    ## [M::main] Version: 2.24-r1122
    ## [M::main] CMD: minimap2 -t 8 -N 5000 -a l100_r5000_ref.fa l100_r5000_reads.fa
    ## [M::main] Real time: 5.840 sec; CPU: 5.872 sec; Peak RSS: 3.496 GB
    ## 0
    ## @PG  ID:minimap2 PN:minimap2 VN:2.24-r1122   CL:minimap2 -t 8 -N 5000 -a l100_r5000_ref.fa l100_r5000_reads.fa
    ## 1    4   *   0   0   *   *   0   0   AGTGTACAACCCTCTTCGTACACCGATGCCGAAGGCGCCTTAGTGGAGGCGTCCGATTCAGGAGCGTGCTCTTGTATGAGTACATGGTGCGGTTTAATGT    *   rl:i:0
    ## 2    4   *   0   0   *   *   0   0   ATCGATGTTCCTCGAACGCTGACTCTCAGTCTCGGCACGCAATAACACACTGTGACACCGTACTAGGACCGTCAATTAGATTGGGTGATTGGGACGAGAT    *   rl:i:0
    ## 3    4   *   0   0   *   *   0   0   GTAGTTACGTGAAACGCCTCTCGGGGAGTTCACATTTCTTTCAGCCTCGACCGTACAACCAATGTAAGGTTATGGAGCCATGTTGTTTCAAGAATGTGTT    *   rl:i:0
    ## 4    4   *   0   0   *   *   0   0   CGTAAATTGTCAGGATAGCCCAACGCACGTCATTGCTAACTGACGACTAGGTGACCACAGGATAACGCTATGATCAGCATTGCAAATCCCATTGATCGCC    *   rl:i:0
    ## 5    4   *   0   0   *   *   0   0   CTCGATTTGCACAGCCTTGCTGGCTGAAAACCGGGGACTTTGAGTTGAGACCGAGCCGCGGCCTGTGTGCGCGGTGAAAATGGTCTGGTCTATAACCGCT    *   rl:i:0
    ## 6    4   *   0   0   *   *   0   0   ATGTCGGATCTTTGCGCGCAAAAGGTCAACTATCTGGGTAGAATCTCCCGACATCCGAACACCTATTCTCTCCGTAACGCAGGCATAAATCACTGGTGCG    *   rl:i:0
    ## 7    4   *   0   0   *   *   0   0   AAAGGCAAAACATAATAAGCGGATTACCGAGAAGCAGTAACATTCCCCATGGTCAGGAAACTGCCTGCATAACCCTTTCGTTACAGTATAACAGTAATGT    *   rl:i:0
    ## 8    4   *   0   0   *   *   0   0   TTGGAAACCACGGTAGGCACTTGCACAAATTTACACTCGGTCTACCCGATTTTAATCCCTGTTTCATACTCAAATGTTAATATGAAGGTCGCCTACATTG    *   rl:i:0
    ## 9    4   *   0   0   *   *   0   0   CTTTTGGCTGGATTAAGCCCTCAAACCTGTCTGTAGTGGATTCACTGTCAATTAAATAAACACATCTTGATCGAGCGCACTTAACCGTGACGTTAGCGGG    *   rl:i:0

Report up to 3,500 secondary alignments works!

    l=100
    r=3500
    ../script/create_ref_se.pl -l ${l} -r ${r}
    minimap2 -N ${r} -a l${l}_r${r}_ref.fa l${l}_r${r}_reads.fa 2> /dev/null \
       | grep -v "^@" \
       | cut -f1 \
       | sort \
       | uniq -c \
       | sort -k2rn \
       | head

    ## FASTA reads written to l100_r3500_reads.fa
    ## FASTA reference written to l100_r3500_ref.fa
    ## Done
    ##    3500 3500
    ##    3499 3499
    ##    3498 3498
    ##    3497 3497
    ##    3496 3496
    ##    3495 3495
    ##    3494 3494
    ##    3493 3493
    ##    3492 3492
    ##    3491 3491

We will repeat the mapping step 100 times to test whether reported loci
are randomly selected. We have also set `-N` to 100 to report up to 100
additional loci, to provide more loci to (potentially) randomly select
from.

However, the reported mapping of read `20` is always the same, so even
though the reported loci seems to be randomly selected, the same loci is
chosen.

    l=100
    r=100
    for i in {1..100}; do
       minimap2 -N ${r} -a l${l}_r${r}_ref.fa l${l}_r${r}_reads.fa 2> /dev/null
    done | grep $'^20\t0' | sort | uniq -c

    ##     100 20   0   ref 61401   0   100M    *   0   0   CTCAGCGGCGATGCCGAGTAAGGGGATTCGAGAAGCGATTTTGCATACTAGGCGTTCATCCTTCCGTCACCTGCTGTCCGGGTCTATAAGGGTAGCGTGC    *   NM:i:0  ms:i:200    AS:i:200    nn:i:0  tp:A:P  cm:i:16 s1:i:94 s2:i:94 de:f:0  rl:i:0

Mapping repetitive reads to chromosome X, retaining up to 500 additional
loci. However, some reads are uniquely mapped such as read `5`.

    minimap2 -t 8 -N 500 -a ../data/chrX.fa.gz ../data/chrx_kmer.fa.gz 2> /dev/null \
       | grep -v "^@" \
       | cut -f1 \
       | sort -n \
       | uniq -c \
       | head

    ##       1 1
    ##       2 2
    ##       4 3
    ##       4 4
    ##       1 5
    ##       6 6
    ##       7 7
    ##       8 8
    ##       1 9
    ##       1 10

Re-mapping with an index generates the same results for the first couple
of reads.

    minimap2 -d chrx.mmi ../data/chrX.fa.gz 2> /dev/null
    minimap2 -t 8 -N 500 -a chrx.mmi ../data/chrx_kmer.fa.gz 2> /dev/null \
       | grep -v "^@" \
       | cut -f1 \
       | sort -n \
       | uniq -c \
       | head

    ##       1 1
    ##       2 2
    ##       4 3
    ##       4 4
    ##       1 5
    ##       6 6
    ##       7 7
    ##       8 8
    ##       1 9
    ##       1 10

Note that using the short read mode eliminates secondary alignments!

    minimap2 -ax sr -t 8 -N 500 chrx.mmi ../data/chrx_kmer.fa.gz 2> /dev/null \
       | grep -v "^@" \
       | cut -f1 \
       | sort -n \
       | uniq -c \
       | head

    ##       1 1
    ##       1 2
    ##       1 3
    ##       1 4
    ##       1 5
    ##       1 6
    ##       1 7
    ##       1 8
    ##       1 9
    ##       1 10

Clean up.

    rm *.sam *.fa core* *.mmi

Minimap2 version used for this README.

    minimap2 --version

    ## 2.24-r1122
