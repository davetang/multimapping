Download chromosome X (if necessary).

    if [[ ! -e chrX.fa.gz ]]; then
       wget --quiet https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrX.fa.gz
    fi

    if [[ ! -e chrX.fa ]]; then
       gunzip chrX.fa.gz
    fi

Run Jellyfish with the following parameters:

-   `-m` - Length of mer
-   `-s` - Initial hash size
-   `-t` - Number of threads
-   `-C` - Count both strand, canonical representation (false)

and output as FASTA.

    jellyfish count -m 100 -s 100M -t 8 -C chrX.fa
    jellyfish dump mer_counts.jf > mer_counts_dumps.fa

Create FASTA with one representative for *k*-mer number.

    ../script/jellyfish_rep.pl mer_counts_dumps.fa | gzip > chrx_kmer.fa.gz

Clean up.

    rm mer_counts_dumps.fa mer_counts.jf
    gzip chrX.fa

Jellyfish version used for this README.

    jellyfish --version

    ## jellyfish 2.3.0
