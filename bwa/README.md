Prepare reference and map.

    l=100
    r=100
    ../script/create_ref_se.pl -l ${l} -r ${r}
    bwa index l${l}_r${r}_ref.fa
    bwa mem l${l}_r${r}_ref.fa l${l}_r${r}_reads.fa 2> /dev/null > default.sam

    ## FASTA reads written to l100_r100_reads.fa
    ## FASTA reference written to l100_r100_ref.fa
    ## Done
    ## [bwa_index] Pack FASTA... 0.05 sec
    ## [bwa_index] Construct BWT for the packed sequence...
    ## [bwa_index] 0.25 seconds elapse.
    ## [bwa_index] Update BWT... 0.05 sec
    ## [bwa_index] Pack forward-only FASTA... 0.02 sec
    ## [bwa_index] Construct SA from BWT and Occ... 0.19 sec
    ## [main] Version: 0.7.17-r1188
    ## [main] CMD: bwa index l100_r100_ref.fa
    ## [main] Real time: 0.661 sec; CPU: 0.558 sec

Examine results.

    head default.sam

    ## @SQ  SN:ref  LN:3030500
    ## @PG  ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem l100_r100_ref.fa l100_r100_reads.fa
    ## 1    0   ref 501 60  100M    *   0   0   CCTCGGCTTTTAATTGGCGCCAATTAGAACTCACCTATAATTTACATACGCGGTCGTGCGAATACAGCAGCATTATGGTACCGTAAACGTAGAATTTGGG    *   NM:i:0  MD:Z:100    AS:i:100    XS:i:0
    ## 2    0   ref 1101    0   100M    *   0   0   TAAGATGGTCTAGTGCTTCTCGGCCATTAGTCTGGGCCTTCACGCGTTACTGAAGGATCCTTCAAACTATTTTATCTTTCAACGCAGGTGGAGCCCGATT    *   NM:i:0  MD:Z:100    AS:i:100    XS:i:100    XA:Z:ref,+1701,100M,0;
    ## 3    0   ref 3501    0   100M    *   0   0   AGACATCTACATGAAGGGAGGGTGAGCTCTCGTCGAACTTGTGCGGCCAGCCAGCGCTGGTTTTTTACAAACAGTCCTTGTAGACGGGATACCTCAAGCG    *   NM:i:0  MD:Z:100    AS:i:100    XS:i:100    XA:Z:ref,+2301,100M,0;ref,+2901,100M,0;
    ## 4    0   ref 5901    0   100M    *   0   0   TCTACCTCACATCGCAAGACCAGCCCCAGAATGGTTAAGTTATAACGGAATTATAAAAGCGGTTGGACGGGTAGACGTTAAAGGTTTTGACCTGCTAAAC    *   NM:i:0  MD:Z:100    AS:i:100    XS:i:100    XA:Z:ref,+4701,100M,0;ref,+5301,100M,0;ref,+4101,100M,0;
    ## 5    0   ref 7701    0   100M    *   0   0   AGATCCTTGATTCTCGTTATACTTTACCCAGTTCGTCAGAACTCTTCTTTGCATAAACCTAATCCGCCGGGAACAGGCTGGCAGTCAACCTGTAAAATGC    *   NM:i:0  MD:Z:100    AS:i:100    XS:i:100    XA:Z:ref,+6501,100M,0;ref,+8901,100M,0;ref,+7101,100M,0;ref,+8301,100M,0;
    ## 6    0   ref 11901   0   100M    *   0   0   GAGATTCGCAGATAAAGCGCACAATCGTGTAGAGTTGGCCCTCCTCTCAAGGCGCAACATTCGCATTAACCGCCTGTAAGCCTGTATAATAGGGTTCCAT    *   NM:i:0  MD:Z:100    AS:i:100    XS:i:100    XA:Z:ref,+10101,100M,0;ref,+11301,100M,0;ref,+9501,100M,0;ref,+12501,100M,0;ref,+10701,100M,0;
    ## 7    0   ref 14901   0   100M    *   0   0   CCCCTTCCCGTAGCATAACCCTCCATTCCGCCTTCTGGTTGAGTATGCTAGGGCAGCCGATCTATACCGACATGAGCCACTGTAAGGTGGAACATACAAT    *   NM:i:0  MD:Z:100    AS:i:100    XS:i:100
    ## 8    0   ref 18501   0   100M    *   0   0   GGATTGGGTCGGAATTTTCCTGCGGAGGTCCCTCCGCACAGGGCCGGCAACTGGCGGTTCTGAAGCACCCCTCCAGATTGTTTAAGCAGCCTGAAACGCA    *   NM:i:0  MD:Z:100    AS:i:100    XS:i:100

Clean up.

    rm *.sam *.fa *.fa.amb *.fa.ann *.fa.bwt *.fa.pac *.fa.sa
