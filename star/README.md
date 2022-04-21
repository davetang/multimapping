STAR usage can be accessed using `--help`. (Since the help page is quite
long, only the first 10 lines are shown below.)

    STAR --help | head

    ## Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq
    ## Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2020
    ## 
    ## STAR version=2.7.10a
    ## STAR compilation time,server,dir=2022-04-19T02:52:25+00:00 c2f985c2d07e:/usr/src/STAR-2.7.10a/source
    ## For more details see:
    ## <https://github.com/alexdobin/STAR>
    ## <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>
    ## ### versions
    ## versionGenome           2.7.4a

Generate 100 unique reads with read lengths of 100 bp and repeat each
unique read incrementally, i.e. first read is unique, second read is
repeated twice, third read is repeated thrice, and so on. Next use these
reads to generate a reference sequence where each read is separated by a
spacer sequence that is two times the length of a read, which is 200 bp
in this case.

The reference sequence is first indexed and then the reads are mapped to
the reference using default settings. We will set
`--genomeSAindexNbases` based on the recommendation:

> For small genomes, the parameter –genomeSAindexNbases must be scaled
> down to min(14, log2(GenomeLength)/2 - 1).

    l=100
    r=100
    ../script/create_ref_se.pl -l ${l} -r ${r}

    genome=$(cat l${l}_r${r}_ref.fa | grep -v "^>" | wc -m)
    param=$(printf "%.f" $(bc -l<<<"(l(${genome}) / l(2)) / 2 - 1"))

    STAR \
       --runMode genomeGenerate \
       --genomeDir star_idx \
       --genomeFastaFiles l${l}_r${r}_ref.fa \
       --genomeSAindexNbases ${param} \
       --runThreadN 8 > /dev/null

    STAR --runMode alignReads \
       --genomeDir star_idx \
       --readFilesIn l${l}_r${r}_reads.fa \
       --outFileNamePrefix default. > /dev/null

    ## FASTA reads written to l100_r100_reads.fa
    ## FASTA reference written to l100_r100_ref.fa
    ## Done

STAR will map to up to 10 loci by default.

> outFilterMultimapNmax 10 int: maximum number of loci the read is
> allowed to map to. Alignments (all of them) will be output only if the
> read maps to no more loci than this value. Otherwise no alignments
> will be output, and the read will be counted as “mapped to too many
> loci” in the Log.final.out .

Any read mapping to more than 10 loci will become unmapped and will be
excluded from the output SAM by default as per the `outSAMunmapped`
parameter. To include unmapped reads use `outSAMunmapped Within` (not
used in this README).

> outSAMunmapped None string(s): output of unmapped reads in the SAM
> format 1st word: None … no output Within … output unmapped reads
> within the main SAM file (i.e. Aligned.out.sam) 2nd word: KeepPairs …
> record unmapped mate for each alignment, and, in case of unsorted
> output, keep it adjacent to its mapped mate. Only affects
> multi-mapping reads.

    cat default.Aligned.out.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2n

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

Even with `outFilterMultimapNmax` set to 100, the maximum number of loci
is 50.

    l=100
    r=100

    STAR --runMode alignReads \
       --genomeDir star_idx \
       --readFilesIn l${l}_r${r}_reads.fa \
       --outFilterMultimapNmax 100 \
       --outFileNamePrefix m100. > /dev/null

    cat m100.Aligned.out.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2rn | head

    ##      50 50
    ##      49 49
    ##      48 48
    ##      47 47
    ##      46 46
    ##      45 45
    ##      44 44
    ##      43 43
    ##      42 42
    ##      41 41

Setting `winAnchorMultimapNmax` to 100 does not increase the number of
multimappers.

> winAnchorMultimapNmax 50 int&gt;0: max number of loci anchors are
> allowed to map to

    l=100
    r=100

    STAR --runMode alignReads \
       --genomeDir star_idx \
       --readFilesIn l${l}_r${r}_reads.fa \
       --outFilterMultimapNmax 100 \
       --winAnchorMultimapNmax 100 \
       --outFileNamePrefix ma100. > /dev/null

    cat ma100.Aligned.out.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2rn | head

    ##      50 50
    ##      49 49
    ##      48 48
    ##      47 47
    ##      46 46
    ##      45 45
    ##      44 44
    ##      43 43
    ##      42 42
    ##      41 41

Increasing `limitOutSAMoneReadBytes` also does not increase the number
of multimappers.

> limitOutSAMoneReadBytes 100000 int&gt;0: max size of the SAM record
> (bytes) for one read. Recommended value:
> &gt;(2*(LengthMate1+LengthMate2+100)*outFilterMultimapNmax

    l=100
    r=100

    STAR --runMode alignReads \
       --genomeDir star_idx \
       --readFilesIn l${l}_r${r}_reads.fa \
       --outFilterMultimapNmax 100 \
       --winAnchorMultimapNmax 100 \
       --limitOutSAMoneReadBytes 1000000 \
       --outFileNamePrefix ml100. > /dev/null

    cat ml100.Aligned.out.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2rn | head

    ##      50 50
    ##      49 49
    ##      48 48
    ##      47 47
    ##      46 46
    ##      45 45
    ##      44 44
    ##      43 43
    ##      42 42
    ##      41 41

Increasing `--winAnchorMultimapNmax` and `seedPerWindowNmax` to 100
increases the number of loci a read can map to.

> seedPerWindowNmax 50 int&gt;0: max number of seeds per window

    l=100
    r=100

    STAR --runMode alignReads \
       --genomeDir star_idx \
       --readFilesIn l${l}_r${r}_reads.fa \
       --outFilterMultimapNmax 100 \
       --winAnchorMultimapNmax 100 \
       --seedPerWindowNmax 100 \
       --outFileNamePrefix ms100. > /dev/null

    cat ms100.Aligned.out.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2rn | head

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

Regarding the order of multimappers, there is the `outMultimapperOrder`
parameter:

> outMultimapperOrder Old\_2.4 string: order of multimapping alignments
> in the output files Old\_2.4 … quasi-random order used before 2.5.0
> Random … random order of alignments for each multi-mapper. Read mates
> (pairs) are always adjacent, all alignment for each read stay
> together. This option will become default in the future releases.

However it seems like the same loci is chosen as the primary alignment.

    l=100
    r=100

    STAR --runMode alignReads \
       --genomeDir star_idx \
       --readFilesIn l${l}_r${r}_reads.fa \
       --outFileNamePrefix default2. > /dev/null

    STAR --runMode alignReads \
       --genomeDir star_idx \
       --readFilesIn l${l}_r${r}_reads.fa \
       --outFileNamePrefix default3. > /dev/null

    STAR --runMode alignReads \
       --genomeDir star_idx \
       --readFilesIn l${l}_r${r}_reads.fa \
       --outFileNamePrefix default4. > /dev/null

    cat default*.sam | grep "^5"

    ## 5    0   ref 3201    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:1  AS:i:98 nM:i:0
    ## 5    256 ref 3501    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:2  AS:i:98 nM:i:0
    ## 5    256 ref 3801    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:3  AS:i:98 nM:i:0
    ## 5    256 ref 4101    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:4  AS:i:98 nM:i:0
    ## 5    256 ref 4401    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:5  AS:i:98 nM:i:0
    ## 5    0   ref 3201    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:1  AS:i:98 nM:i:0
    ## 5    256 ref 3501    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:2  AS:i:98 nM:i:0
    ## 5    256 ref 3801    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:3  AS:i:98 nM:i:0
    ## 5    256 ref 4101    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:4  AS:i:98 nM:i:0
    ## 5    256 ref 4401    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:5  AS:i:98 nM:i:0
    ## 5    0   ref 3201    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:1  AS:i:98 nM:i:0
    ## 5    256 ref 3501    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:2  AS:i:98 nM:i:0
    ## 5    256 ref 3801    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:3  AS:i:98 nM:i:0
    ## 5    256 ref 4101    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:4  AS:i:98 nM:i:0
    ## 5    256 ref 4401    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:5  AS:i:98 nM:i:0
    ## 5    0   ref 3201    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:1  AS:i:98 nM:i:0
    ## 5    256 ref 3501    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:2  AS:i:98 nM:i:0
    ## 5    256 ref 3801    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:3  AS:i:98 nM:i:0
    ## 5    256 ref 4101    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:4  AS:i:98 nM:i:0
    ## 5    256 ref 4401    0   100M    *   0   0   GGGGCGCCGCGATCTTCCGACGCAGGGTTTTATTACTCGGGCTATATAGGTTAGCGACGATGTGTAATTACTAGATTACGACATATATAAGACCAGGCCG    *   NH:i:5  HI:i:5  AS:i:98 nM:i:0

Map repetitive reads to chrX.

    gunzip -c ../data/chrX.fa.gz > chrX.fa

    STAR \
       --runMode genomeGenerate \
       --genomeDir chrx \
       --genomeFastaFiles chrX.fa \
       --genomeSAindexNbases 12 \
       --runThreadN 8 > /dev/null

    STAR --runMode alignReads \
       --genomeDir chrx \
       --outFilterMultimapNmax 1000 \
       --winAnchorMultimapNmax 1000 \
       --seedPerWindowNmax 1000 \
       --readFilesCommand "gunzip -c" \
       --readFilesIn ../data/chrx_kmer.fa.gz \
       --outFileNamePrefix chrx. \
       --runThreadN 8 > /dev/null

Reads repeating up to 10 times are mapped as expected.

    cat chrx.Aligned.out.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2n | head

    ##       1 1
    ##       2 2
    ##       3 3
    ##       4 4
    ##       6 5
    ##       6 6
    ##       7 7
    ##       8 8
    ##       9 9
    ##      10 10

However reads repeating up to 1,000 times are not mapped to all their
loci. In addition, note that reads mapping over 1,000 times are
discarded (therefore STAR does map up to 1,000 times but it just does
not report all the loci with the provided parameters).

    cat chrx.Aligned.out.sam | grep -v "^@" | cut -f1 | sort | uniq -c | sort -k2rn | head

    ##     120 1000
    ##     127 999
    ##     119 998
    ##     120 997
    ##     101 995
    ##     103 993
    ##     120 992
    ##     103 991
    ##     127 990
    ##     101 989

Clean up.

    rm -rf *.sam *.fa star_idx chrx *.out *.out.tab

STAR version used for this README.

    STAR --version

    ## 2.7.10a
