HISAT2 usage.

    hisat2 || true

    ## No index, query, or output file specified!
    ## HISAT2 version 2.2.1 by Daehwan Kim (infphilo@gmail.com, www.ccb.jhu.edu/people/infphilo)
    ## Usage: 
    ##   hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
    ## 
    ##   <ht2-idx>  Index filename prefix (minus trailing .X.ht2).
    ##   <m1>       Files with #1 mates, paired with files in <m2>.
    ##              Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
    ##   <m2>       Files with #2 mates, paired with files in <m1>.
    ##              Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
    ##   <r>        Files with unpaired reads.
    ##              Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
    ##   <sam>      File for SAM output (default: stdout)
    ## 
    ##   <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
    ##   specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.
    ## 
    ## Options (defaults in parentheses):
    ## 
    ##  Input:
    ##   -q                 query input files are FASTQ .fq/.fastq (default)
    ##   --qseq             query input files are in Illumina's qseq format
    ##   -f                 query input files are (multi-)FASTA .fa/.mfa
    ##   -r                 query input files are raw one-sequence-per-line
    ##   -c                 <m1>, <m2>, <r> are sequences themselves, not files
    ##   -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
    ##   -u/--upto <int>    stop after first <int> reads/pairs (no limit)
    ##   -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
    ##   -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
    ##   --phred33          qualities are Phred+33 (default)
    ##   --phred64          qualities are Phred+64
    ##   --int-quals        qualities encoded as space-delimited integers
    ## 
    ##  Presets:                 Same as:
    ##    --fast                 --no-repeat-index
    ##    --sensitive            --bowtie2-dp 1 -k 30 --score-min L,0,-0.5
    ##    --very-sensitive       --bowtie2-dp 2 -k 50 --score-min L,0,-1
    ## 
    ##  Alignment:
    ##   --bowtie2-dp <int> use Bowtie2's dynamic programming alignment algorithm (0) - 0: no dynamic programming, 1: conditional dynamic programming, and 2: unconditional dynamic programming (slowest)
    ##   --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
    ##   --ignore-quals     treat all quality values as 30 on Phred scale (off)
    ##   --nofw             do not align forward (original) version of read (off)
    ##   --norc             do not align reverse-complement version of read (off)
    ##   --no-repeat-index  do not use repeat index
    ## 
    ##  Spliced Alignment:
    ##   --pen-cansplice <int>              penalty for a canonical splice site (0)
    ##   --pen-noncansplice <int>           penalty for a non-canonical splice site (12)
    ##   --pen-canintronlen <func>          penalty for long introns (G,-8,1) with canonical splice sites
    ##   --pen-noncanintronlen <func>       penalty for long introns (G,-8,1) with noncanonical splice sites
    ##   --min-intronlen <int>              minimum intron length (20)
    ##   --max-intronlen <int>              maximum intron length (500000)
    ##   --known-splicesite-infile <path>   provide a list of known splice sites
    ##   --novel-splicesite-outfile <path>  report a list of splice sites
    ##   --novel-splicesite-infile <path>   provide a list of novel splice sites
    ##   --no-temp-splicesite               disable the use of splice sites found
    ##   --no-spliced-alignment             disable spliced alignment
    ##   --rna-strandness <string>          specify strand-specific information (unstranded)
    ##   --tmo                              reports only those alignments within known transcriptome
    ##   --dta                              reports alignments tailored for transcript assemblers
    ##   --dta-cufflinks                    reports alignments tailored specifically for cufflinks
    ##   --avoid-pseudogene                 tries to avoid aligning reads to pseudogenes (experimental option)
    ##   --no-templatelen-adjustment        disables template length adjustment for RNA-seq reads
    ## 
    ##  Scoring:
    ##   --mp <int>,<int>   max and min penalties for mismatch; lower qual = lower penalty <6,2>
    ##   --sp <int>,<int>   max and min penalties for soft-clipping; lower qual = lower penalty <2,1>
    ##   --no-softclip      no soft-clipping
    ##   --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
    ##   --rdg <int>,<int>  read gap open, extend penalties (5,3)
    ##   --rfg <int>,<int>  reference gap open, extend penalties (5,3)
    ##   --score-min <func> min acceptable alignment score w/r/t read length
    ##                      (L,0.0,-0.2)
    ## 
    ##  Reporting:
    ##   -k <int>           It searches for at most <int> distinct, primary alignments for each read. Primary alignments mean 
    ##                      alignments whose alignment score is equal to or higher than any other alignments. The search terminates 
    ##                      when it cannot find more distinct valid alignments, or when it finds <int>, whichever happens first. 
    ##                      The alignment score for a paired-end alignment equals the sum of the alignment scores of 
    ##                      the individual mates. Each reported read or pair alignment beyond the first has the SAM ‘secondary’ bit 
    ##                      (which equals 256) set in its FLAGS field. For reads that have more than <int> distinct, 
    ##                      valid alignments, hisat2 does not guarantee that the <int> alignments reported are the best possible 
    ##                      in terms of alignment score. Default: 5 (linear index) or 10 (graph index).
    ##                      Note: HISAT2 is not designed with large values for -k in mind, and when aligning reads to long, 
    ##                      repetitive genomes, large -k could make alignment much slower.
    ##   --max-seeds <int>  HISAT2, like other aligners, uses seed-and-extend approaches. HISAT2 tries to extend seeds to 
    ##                      full-length alignments. In HISAT2, --max-seeds is used to control the maximum number of seeds that 
    ##                      will be extended. For DNA-read alignment (--no-spliced-alignment), HISAT2 extends up to these many seeds
    ##                      and skips the rest of the seeds. For RNA-read alignment, HISAT2 skips extending seeds and reports 
    ##                      no alignments if the number of seeds is larger than the number specified with the option, 
    ##                      to be compatible with previous versions of HISAT2. Large values for --max-seeds may improve alignment 
    ##                      sensitivity, but HISAT2 is not designed with large values for --max-seeds in mind, and when aligning 
    ##                      reads to long, repetitive genomes, large --max-seeds could make alignment much slower. 
    ##                      The default value is the maximum of 5 and the value that comes with -k times 2.
    ##   -a/--all           HISAT2 reports all alignments it can find. Using the option is equivalent to using both --max-seeds 
    ##                      and -k with the maximum value that a 64-bit signed integer can represent (9,223,372,036,854,775,807).
    ##   --repeat           report alignments to repeat sequences directly
    ## 
    ##  Paired-end:
    ##   -I/--minins <int>  minimum fragment length (0), only valid with --no-spliced-alignment
    ##   -X/--maxins <int>  maximum fragment length (500), only valid with --no-spliced-alignment
    ##   --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
    ##   --no-mixed         suppress unpaired alignments for paired reads
    ##   --no-discordant    suppress discordant alignments for paired reads
    ## 
    ##  Output:
    ##   -t/--time          print wall-clock time taken by search phases
    ##   --un <path>           write unpaired reads that didn't align to <path>
    ##   --al <path>           write unpaired reads that aligned at least once to <path>
    ##   --un-conc <path>      write pairs that didn't align concordantly to <path>
    ##   --al-conc <path>      write pairs that aligned concordantly at least once to <path>
    ##   (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g.
    ##   --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)
    ##   --summary-file <path> print alignment summary to this file.
    ##   --new-summary         print alignment summary in a new style, which is more machine-friendly.
    ##   --quiet               print nothing to stderr except serious errors
    ##   --met-file <path>     send metrics to file at <path> (off)
    ##   --met-stderr          send metrics to stderr (off)
    ##   --met <int>           report internal counters & metrics every <int> secs (1)
    ##   --no-head             suppress header lines, i.e. lines starting with @
    ##   --no-sq               suppress @SQ header lines
    ##   --rg-id <text>        set read group id, reflected in @RG line and RG:Z: opt field
    ##   --rg <text>           add <text> ("lab:value") to @RG line of SAM header.
    ##                         Note: @RG line only printed when --rg-id is set.
    ##   --omit-sec-seq        put '*' in SEQ and QUAL fields for secondary alignments.
    ## 
    ##  Performance:
    ##   -o/--offrate <int> override offrate of index; must be >= index's offrate
    ##   -p/--threads <int> number of alignment threads to launch (1)
    ##   --reorder          force SAM output order to match order of input reads
    ##   --mm               use memory-mapped I/O for index; many 'hisat2's can share
    ## 
    ##  Other:
    ##   --qc-filter        filter out reads that are bad according to QSEQ filter
    ##   --seed <int>       seed for random number generator (0)
    ##   --non-deterministic seed rand. gen. arbitrarily instead of using read attributes
    ##   --remove-chrname   remove 'chr' from reference names in alignment
    ##   --add-chrname      add 'chr' to reference names in alignment 
    ##   --version          print version information and quit
    ##   -h/--help          print this usage message
    ## (ERR): hisat2-align exited with value 1

Generate 100 unique reads with read lengths of 100 bp and repeat each
unique read incrementally, i.e. first read is unique, second read is
repeated twice, third read is repeated thrice, and so on. Next use these
reads to generate a reference sequence where each read is separated by a
spacer sequence that is two times the length of a read, which is 200 bp
in this case.

The reference sequence is first indexed and then the reads are mapped to
the reference using default settings. The `-f` parameter specifies that
input reads are FASTA files.

    l=100
    r=100
    ../script/create_ref_se.pl -l ${l} -r ${r}
    hisat2-build l${l}_r${r}_ref.fa l${l}_r${r}_ref 2> /dev/null > /dev/null
    hisat2 -f -x l${l}_r${r}_ref l${l}_r${r}_reads.fa > default.sam 2> /dev/null

    ## FASTA reads written to l100_r100_reads.fa
    ## FASTA reference written to l100_r100_ref.fa
    ## Done

The `-k` and `--max-seeds` parameters control the reporting of
multimappers.

> -k <int> It searches for at most <int> distinct, primary alignments
> for each read. Primary alignments mean alignments whose alignment
> score is equal to or higher than any other alignments. The search
> terminates when it cannot find more distinct valid alignments, or when
> it finds <int>, whichever happens first. The alignment score for a
> paired-end alignment equals the sum of the alignment scores of the
> individual mates. Each reported read or pair alignment beyond the
> first has the SAM ‘secondary’ bit (which equals 256) set in its FLAGS
> field. For reads that have more than <int> distinct, valid alignments,
> hisat2 does not guarantee that the <int> alignments reported are the
> best possible in terms of alignment score. Default: 5 (linear index)
> or 10 (graph index). Note: HISAT2 is not designed with large values
> for -k in mind, and when aligning reads to long, repetitive genomes,
> large -k could make alignment much slower.
>
> –max-seeds <int> HISAT2, like other aligners, uses seed-and-extend
> approaches. HISAT2 tries to extend seeds to full-length alignments. In
> HISAT2, –max-seeds is used to control the maximum number of seeds that
> will be extended. For DNA-read alignment (–no-spliced-alignment),
> HISAT2 extends up to these many seeds and skips the rest of the seeds.
> For RNA-read alignment, HISAT2 skips extending seeds and reports no
> alignments if the number of seeds is larger than the number specified
> with the option, to be compatible with previous versions of HISAT2.
> Large values for –max-seeds may improve alignment sensitivity, but
> HISAT2 is not designed with large values for –max-seeds in mind, and
> when aligning reads to long, repetitive genomes, large –max-seeds
> could make alignment much slower. The default value is the maximum of
> 5 and the value that comes with -k times 2.

Reads are mapped up to 5 times as stated above.

    cat default.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2n | head

    ##       1 1
    ##       2 2
    ##       3 3
    ##       4 4
    ##       5 5
    ##       5 6
    ##       5 7
    ##       5 8
    ##       5 9
    ##       5 10

The `-a` parameter will report all alignments that can be found!

> -a/–all HISAT2 reports all alignments it can find. Using the option is
> equivalent to using both –max-seeds and -k with the maximum value that
> a 64-bit signed integer can represent (9,223,372,036,854,775,807).

    l=100
    r=100
    hisat2 -a -f -x l${l}_r${r}_ref l${l}_r${r}_reads.fa 2> /dev/null \
       | grep -v "^@" \
       | cut -f1 \
       | sort \
       | uniq -c \
       | sort -k2rn \
       | head

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

We will repeat the mapping step 100 times to test whether reported loci
are randomly selected. The reported mapping of read `100` is always the
same.

    l=100
    r=100
    for i in {1..100}; do
       hisat2 -a -f -x l${l}_r${r}_ref l${l}_r${r}_reads.fa 2> /dev/null
    done | grep $'^100\t0' | sort | uniq -c

    ##     100 100  0   ref 1507401 1   100M    *   0   0   CTAGCGACAATTGTGGGATGGCCCCTGGGGATCACTGGAGACGAAAAGATGCTATGATGCACCAGTGCGTGCTTGACACTAGTGTATATAAAACTGTACA    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100    YT:Z:UU NH:i:100

HISAT2 can map up to 3,500 loci (and more!).

    l=100
    r=3500
    ../script/create_ref_se.pl -l ${l} -r ${r}
    hisat2-build -p 8 l${l}_r${r}_ref.fa l${l}_r${r}_ref 2> /dev/null > /dev/null
    hisat2 -a -p 8 -f -x l${l}_r${r}_ref l${l}_r${r}_reads.fa > r${r}.sam 2> /dev/null

    cat r${r}.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2rn | head

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

Mapping repetitive reads to chromosome X.

    gunzip -c ../data/chrX.fa.gz > chrX.fa
    hisat2-build chrX.fa chrx 2> /dev/null > /dev/null
    gunzip -c ../data/chrx_kmer.fa.gz > chrx_kmer.fa
    hisat2 -f -a -x chrx chrx_kmer.fa 2> /dev/null > chrx.sam

Reads repeating up to 10 times are mapped as expected.

    cat chrx.sam \
       | grep -v "^@" \
       | cut -f1 \
       | sort -n \
       | uniq -c \
       | head

    ##       1 1
    ##       2 2
    ##       3 3
    ##       4 4
    ##       5 5
    ##       6 6
    ##       7 7
    ##       8 8
    ##       9 9
    ##      10 10

Reads repeating up to *n* times are also mapped as expected!

    cat chrx.sam \
       | grep -v "^@" \
       | cut -f1 \
       | sort -n \
       | uniq -c \
       | sort -k2rn \
       | head

    ##    1470 1470
    ##    1469 1469
    ##    1468 1468
    ##    1467 1467
    ##    1465 1465
    ##    1463 1463
    ##    1452 1452
    ##    1449 1449
    ##    1448 1448
    ##    1447 1447

Clean up.

    rm *.fa *.ht2 *.sam

HISAT2 version used for this README.

    hisat2 --version

    ## /opt/hisat2-2.2.1/hisat2-align-s version 2.2.1
    ## 64-bit
    ## Built on 2fc4bb9bdc13
    ## Tue 19 Apr 2022 02:47:26 AM UTC
    ## Compiler: gcc version 9.4.0 (Ubuntu 9.4.0-1ubuntu1~20.04.1) 
    ## Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY -std=c++11
    ## Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
