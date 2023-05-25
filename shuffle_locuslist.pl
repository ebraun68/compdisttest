#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use List::Util qw/shuffle/;

############################################################################
# Set the global variables
############################################################################
my($progname) = $0;

my($hter);
my($iter);
my($jter);
my($kter);
my($lter);
my($mter);
my($nter);
my($qter);
my($zter);

my($tempch1);
my($tempch2);
my($tempvar);
my($tempstr);
my @temparray;

my($cmd);
my($verbose) = 0;

my($expneg1) = 0.367879441171442;

if ( @ARGV != 4 ) {
	print STDERR "Usage:\n  \$ $progname <locuslist> <outlist> <samplesize> <randseed>\n";
	print STDERR "  locuslist  = text file listing translated loci to analyze\n";
	print STDERR "  outlist    = rearranged text file listing proteins seqeunces\n";
	print STDERR "  samplesize = number of proteins seqeunces to sample\n";
	print STDERR "       -1 -- jackknife with deletion probability of exp(-1)\n";
	print STDERR "        0 -- output the full number of proteins\n";
	print STDERR "        n -- output a list of n files (if n > the total number\n";
	print STDERR "             of file then n = total number of sequences)\n";
	print STDERR "  randseed   = random integer to seed the random number generator\n";
	print STDERR "        use 0 for implicit srand() call \n";
	print STDERR "exiting...\n";
	exit;
}

my($locuslist) = $ARGV[0];
my($outlist)   = $ARGV[1];
my($ss)        = $ARGV[2];
my($seed)      = $ARGV[3];
if ( $seed > 0 ) { srand($seed); }

############################################################################
# Read the locus list
############################################################################
print "Reading the locus list $locuslist... ";
open (my $INF, $locuslist) or die "Could not open file $locuslist for input.\n";
my @locusdata = shuffle <$INF>; # Read and shuffle the input file
close($INF) or die "Could not close file $locuslist.\n";
my($locuslines) = $#locusdata + 1;
print "$locuslines file names read...\n";

my($samplesize) = $locuslines;
my($jackknife) = 0;
if ( $ss < $locuslines ) { $samplesize = $ss; }
if ( $ss == 0 ) { $samplesize = $locuslines; }
if ( $ss == -1 ) {
	$samplesize = $locuslines;
	$jackknife = 1;
}

open (my $OUTF, ">$outlist") or die "Could not open file $outlist for output.\n";

for ($iter=0; $iter<$samplesize; $iter++) {
	chomp($locusdata[$iter]);
	if ( $jackknife == 0 ) { $tempvar = 1; }
	else { $tempvar = rand(); }
	if ( $tempvar > $expneg1 ) { print $OUTF "$locusdata[$iter]\n"; }
}

close($OUTF) or die "Could not close file $outlist.\n";

