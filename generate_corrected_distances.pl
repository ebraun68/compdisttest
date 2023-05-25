#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# generate_corrected_distances.pl
#
# This program uses the output from generate_dstance_matrices.pl and applies
# corrections to those distances, outputting a modified distance matrix in
# nexus format.
############################################################################

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

if ( @ARGV != 4 ) {
	print STDERR "Usage:\n  \$ $progname <distlist> <alpha> <outfile> <prefix>\n";
	print STDERR "  distlist = text file listing evolutionary distances\n";
	print STDERR "  alpha    = alpha parameter for distance correction\n";
	print STDERR "    -- use negative number for no correction\n";
	print STDERR "    -- use 0 for Poisson-type correction\n";
	print STDERR "  outfile  = nexus format distance matrix\n";
	print STDERR "  prefix   = prefix for the treefiles\n";
	print STDERR "exiting...\n";
	exit;
}

my($distlist) = $ARGV[0];
my($alpha)    = $ARGV[1];
my($outfile)  = $ARGV[2];
my($prefix)   = $ARGV[3];

############################################################################
# Read the taxon list
############################################################################
print STDERR "Reading the distances from file $distlist... ";
open (my $INF, $distlist) or die "Could not open file $distlist for input.\n";
my @distdata = <$INF>; # Read the input file
close($INF) or die "Could not close file $distlist.\n";
my($distlines) = $#distdata + 1;
for ($iter=0; $iter<$distlines; $iter++) {
	chomp($distdata[$iter]);
}

my($ntax) = $distdata[0];
print STDERR "pairwise distances for $ntax taxa read...\n";
my($arraysize) = ( $ntax * $ntax ) + 1;
my($distdataline) = $ntax + 1;

my @taxdata;
for ($iter=0; $iter<$ntax; $iter++) {
	$jter = $iter + 1;
	$taxdata[$iter] = $distdata[$jter];
	print "$taxdata[$iter]\n";
}

############################################################################
# Generate the corrected distances
############################################################################
my @aa_distarray;
my @day_distarray;
my @han_distarray;
my @hp_distarray;
my @size_distarray;
my @gc_distarray;

# iterate over the distance matrix elements for the raw amino acid coding
(@temparray) = split(/\t/, $distdata[$distdataline]);
for ($iter=0; $iter<$arraysize; $iter++) {
	if ( $iter == 0 ) { $aa_distarray[$iter] = $temparray[$iter]; }
	else { $aa_distarray[$iter] = corrD($alpha,$temparray[$iter]); }
}

# iterate over the distance matrix elements for Dayhoff coding
$distdataline++;
(@temparray) = split(/\t/, $distdata[$distdataline]);
for ($iter=0; $iter<$arraysize; $iter++) {
	if ( $iter == 0 ) { $day_distarray[$iter] = $temparray[$iter]; }
	else { $day_distarray[$iter] = corrD($alpha,$temparray[$iter]); }
}

# iterate over the distance matrix elements for Hanada coding
$distdataline++;
(@temparray) = split(/\t/, $distdata[$distdataline]);
for ($iter=0; $iter<$arraysize; $iter++) {
	if ( $iter == 0 ) { $han_distarray[$iter] = $temparray[$iter]; }
	else { $han_distarray[$iter] = corrD($alpha,$temparray[$iter]); }
}

# iterate over the distance matrix elements for HP coding
$distdataline++;
(@temparray) = split(/\t/, $distdata[$distdataline]);
for ($iter=0; $iter<$arraysize; $iter++) {
	if ( $iter == 0 ) { $hp_distarray[$iter] = $temparray[$iter]; }
	else { $hp_distarray[$iter] = corrD($alpha,$temparray[$iter]); }
}

# iterate over the distance matrix elements for size coding
$distdataline++;
(@temparray) = split(/\t/, $distdata[$distdataline]);
for ($iter=0; $iter<$arraysize; $iter++) {
	if ( $iter == 0 ) { $size_distarray[$iter] = $temparray[$iter]; }
	else { $size_distarray[$iter] = corrD($alpha,$temparray[$iter]); }
}

# iterate over the distance matrix elements for GC coding
$distdataline++;
(@temparray) = split(/\t/, $distdata[$distdataline]);
for ($iter=0; $iter<$arraysize; $iter++) {
	if ( $iter == 0 ) { $gc_distarray[$iter] = $temparray[$iter]; }
	else { $gc_distarray[$iter] = corrD($alpha,$temparray[$iter]); }
}

############################################################################
# Output the distances
############################################################################
open (my $NEXF, ">$outfile") or die "Could not open file $outfile for output.\n";

print $NEXF "#NEXUS\n\n";

#
# Amino acid distances
#
print $NEXF "BEGIN DISTANCES; \[ Amino acid tiled compression distances with alpha = $alpha \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$ntax";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$ntax; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$ntax; $jter++) {
		$zter = $jter + ($iter*$ntax);
		print $NEXF "\t$aa_distarray[$zter+1]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "log file=$prefix.LS_log.txt replace;\n";
print $NEXF "\[! Amino acid tiled compression distances with alpha = $alpha \]\n";
print $NEXF "hsearch;\n";
print $NEXF "describetrees 1/ plot=phylogram;\n";
print $NEXF "savetrees file=$prefix.LS_aa_compdist.tre format=altnexus brlens=yes replace;\n";
print $NEXF "END;\n\n";

#
# Dayhoff distances
#
print $NEXF "BEGIN DISTANCES; \[ Dayhoff tiled compression distances with alpha = $alpha \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$ntax";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$ntax; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$ntax; $jter++) {
		$zter = $jter + ($iter*$ntax);
		print $NEXF "\t$day_distarray[$zter+1]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! Dayhoff tiled compression distances with alpha = $alpha \]\n";
print $NEXF "hsearch;\n";
print $NEXF "describetrees 1/ plot=phylogram;\n";
print $NEXF "savetrees file=$prefix.LS_Dayhoff_compdist.tre format=altnexus brlens=yes replace;\n";
print $NEXF "END;\n\n";

#
# Hanada distances
#
print $NEXF "BEGIN DISTANCES; \[ Hanada tiled compression distances with alpha = $alpha \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$ntax";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$ntax; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$ntax; $jter++) {
		$zter = $jter + ($iter*$ntax);
		print $NEXF "\t$han_distarray[$zter+1]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! Hanada tiled compression distances with alpha = $alpha \]\n";
print $NEXF "hsearch;\n";
print $NEXF "describetrees 1/ plot=phylogram;\n";
print $NEXF "savetrees file=$prefix.LS_Hanada_compdist.tre format=altnexus brlens=yes replace;\n";
print $NEXF "END;\n\n";

#
# HP distances
#
print $NEXF "BEGIN DISTANCES; \[ HP tiled compression distances with alpha = $alpha \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$ntax";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$ntax; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$ntax; $jter++) {
		$zter = $jter + ($iter*$ntax);
		print $NEXF "\t$hp_distarray[$zter+1]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! HP tiled compression distances with alpha = $alpha \]\n";
print $NEXF "hsearch;\n";
print $NEXF "describetrees 1/ plot=phylogram;\n";
print $NEXF "savetrees file=$prefix.LS_HP_compdist.tre format=altnexus brlens=yes replace;\n";
print $NEXF "END;\n\n";

#
# Amino acid size distances
#
print $NEXF "BEGIN DISTANCES; \[ Amino acid size tiled compression distances with alpha = $alpha \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$ntax";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$ntax; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$ntax; $jter++) {
		$zter = $jter + ($iter*$ntax);
		print $NEXF "\t$size_distarray[$zter+1]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! Amino acid size tiled compression distances with alpha = $alpha \]\n";
print $NEXF "hsearch;\n";
print $NEXF "describetrees 1/ plot=phylogram;\n";
print $NEXF "savetrees file=$prefix.LS_aasize_compdist.tre format=altnexus brlens=yes replace;\n";
print $NEXF "END;\n\n";

#
# GC distances
#
print $NEXF "BEGIN DISTANCES; \[ GC tiled compression distances with alpha = $alpha \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$ntax";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$ntax; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$ntax; $jter++) {
		$zter = $jter + ($iter*$ntax);
		print $NEXF "\t$gc_distarray[$zter+1]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! GC tiled compression distances with alpha = $alpha \]\n";
print $NEXF "hsearch;\n";
print $NEXF "describetrees 1/ plot=phylogram;\n";
print $NEXF "savetrees file=$prefix.LS_GC_compdist.tre format=altnexus brlens=yes replace;\n";
print $NEXF "log stop;\n";
print $NEXF "END;\n\n";


close($NEXF) or die "Could not close file $outfile\n";

exit;


#########################
## END OF MAIN PROGRAM ##
## SUBROUTINES FOLLOW ###
#########################


############################################################################
# Subroutine to correct distances
############################################################################
sub corrD {

	my($a,$P) = @_;
	my($localdist);
#	print "alpha = $a -- raw dist = $P\n";
	
	if ( $P == 0 ) { $localdist = 0; }
	elsif ( $a < 0 ) { $localdist = $P; }
	elsif ( $a == 0 ) { $localdist = -1 * log( 1 - $P ); }
	else { $localdist = $a * (((1-$P)**(-1/$a))-1); }
	
	return $localdist;
}
