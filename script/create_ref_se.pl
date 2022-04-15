#!/usr/bin/env perl
#
# Create -r number of reads of length -l where reads are repeated incrementally from 1 to -r.
# For example, the first read is present only once, the second read is repeated twice, and so on.
# The script will overwrite files, which is fine because they can be regenerated with the same parameters.
#

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('h:l:r:', \%opts);

srand(1984);

if ($opts{'h'} ||
    !exists $opts{'l'} ||
    !exists $opts{'r'}
){
   usage();
}

my $infile = $opts{'f'};
my $rl = $opts{'l'};
my $repeat = $opts{'r'};

sub make_read {
   my ($len) = @_;
   my @nuc = qw/A C G T/;
   my $read = '';
   foreach(1..$len){
      $read .= $nuc[int(rand($#nuc+1))];
   }
   return($read)
}

my %reads = ();
while(scalar keys %reads < $repeat){
   my $read = make_read($rl);
   if (! exists $reads{$read}){
      $reads{$read} = 1;
   }
}

my $spacer = 'A' x int($rl * 2);
my $ref = $spacer;
my $i = 1;

my $read_fasta = "l${rl}_r${repeat}_reads.fa";
my $ref_fasta = "l${rl}_r${repeat}_ref.fa";

open(my $read_fh, '>', $read_fasta) || die "Could not open $read_fasta for writing: $!\n";
open(my $ref_fh, '>', $ref_fasta) || die "Could not open $ref_fasta for writing: $!\n";

foreach my $read (keys %reads){
   print $read_fh ">$i\n$read\n";
   my $ins = ($read . $spacer) x $i;
   $ref .= $ins;
   ++$i;
}

print $ref_fh ">ref\n$ref\n";

close($read_fh);
close($ref_fh);

warn("FASTA reads written to $read_fasta\nFASTA reference written to $ref_fasta\nDone\n");

sub usage {
print STDERR <<EOF;
Usage: $0 < -l read length > < -r max repeat >

Where:   -l         generate reads of this length
         -r         maximum number of times a read repeats
         -h         this helpful usage message

EOF
exit();
}

