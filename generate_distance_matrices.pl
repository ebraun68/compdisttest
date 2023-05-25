#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# generate_distance_matrices.pl
#
# This program takes a list of taxa and a list of file identifiers as input
# (the file identifier list should point to protein fasta files that include
# the taxa of interest).
# 
# It requires a path to the pairwise_compression_dist.pl program and it will
# generate a nexus distance matrix and list of distances.
############################################################################

############################################################################
# Set the global variables
############################################################################
my($progname) = $0;
my($pairwise_exec) = "./pairwise_compression_dist.pl";

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

if ( @ARGV != 3 ) {
	print STDERR "Usage:\n  \$ $progname <taxalist> <prefix> <outfile>\n";
	print STDERR "  taxalist = text file listing taxa to analyze\n";
	print STDERR "  prefix   = prefix for file list\n";
	print STDERR "  outfile  = prefix for output files\n";
	print STDERR "Expected format for list of files:\n";
	print STDERR "  <prefix>.sequences.txt = list of sequence files\n";
	print STDERR "Format of the output files:\n";
	print STDERR "  <outfile>.nex      = nexus distance matrices\n";
	print STDERR "  <outfile>.dist.txt = list of evolutionary distances\n";
	print STDERR "exiting...\n";
	exit;
}

my($taxlist) = $ARGV[0];
my($prefix)  = $ARGV[1];
my($outfile) = $ARGV[2];

############################################################################
# Read the taxon list
############################################################################
print "Reading the taxon list $taxlist... ";
open (my $INF, $taxlist) or die "Could not open file $taxlist for input.\n";
my @taxdata = <$INF>; # Read the input file
close($INF) or die "Could not close file $taxlist.\n";
my($taxlines) = $#taxdata + 1;
print "$taxlines taxon names read...\n";

############################################################################
# Calculate the distances
############################################################################
my($pairwise_file);
my @distfiledata;

# Initialize and zero out the distance arrays; also, strip endlines from taxon names
my @aa_distarray;
my @day_distarray;
my @han_distarray;
my @hp_distarray;
my @size_distarray;
my @gc_distarray;

for ($iter=0; $iter<$taxlines; $iter++) {
	chomp($taxdata[$iter]);
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		$aa_distarray[$zter] = 0;
		$day_distarray[$zter] = 0;
		$han_distarray[$zter] = 0;
		$hp_distarray[$zter] = 0;
		$size_distarray[$zter] = 0;
		$gc_distarray[$zter] = 0;		
	}
}

# Iterate over the taxa, calculate distances
for ($iter=0; $iter<$taxlines; $iter++) {
	for ($jter=($iter+1); $jter<$taxlines; $jter++) {
		print "$taxdata[$iter] -- $taxdata[$jter]\n";
		
		# Generate the distance estimates
		$pairwise_file = "compression_distances." . "$taxdata[$iter]" . "_" . "$taxdata[$jter]" . ".txt";
		$cmd = "$pairwise_exec $prefix $taxdata[$iter] $taxdata[$jter] > " . "$pairwise_file";
		print "$cmd\n";
		system("$cmd");
		
		# Read the distance file
		open (my $DISTF, $pairwise_file) or die "Could not open file $pairwise_file for input.\n";
		my @distfiledata = <$DISTF>; # Read the distances file
		close($DISTF) or die "Could not close file $pairwise_file.\n";
		my($tempvar) = $#distfiledata + 1;

		# Extract the distance information for amino acid distances
		for ($kter=0; $kter<$tempvar; $kter++) {
			chomp($distfiledata[$kter]);
			(@temparray) = split(/\s+/, $distfiledata[$kter]);
			$temparray[0] //= "foo"; # placeholder for undefined $temparray[0]
			if ( $temparray[0] eq "corrdist" && $temparray[2] eq "aa_tile" ) {
				$zter = $jter + ($iter*$taxlines);
				$aa_distarray[$zter] = $temparray[3];
				$zter = $iter + ($jter*$taxlines);
				$aa_distarray[$zter] = $temparray[3];
				$kter=$tempvar;
			}
		} # end amino acid distances - for ($kter=0; $kter<$tempvar...

		# Extract the distance information for dayhoff distances
		for ($kter=0; $kter<$tempvar; $kter++) {
			chomp($distfiledata[$kter]);
			(@temparray) = split(/\s+/, $distfiledata[$kter]);
			$temparray[0] //= "foo"; # placeholder for undefined $temparray[0]
			if ( $temparray[0] eq "corrdist" && $temparray[2] eq "day_tile" ) {
				$zter = $jter + ($iter*$taxlines);
				$day_distarray[$zter] = $temparray[3];
				$zter = $iter + ($jter*$taxlines);
				$day_distarray[$zter] = $temparray[3];
				$kter=$tempvar;
			}
		} # end dayhoff distances - for ($kter=0; $kter<$tempvar...
		
		# Extract the distance information for hanada distances
		for ($kter=0; $kter<$tempvar; $kter++) {
			chomp($distfiledata[$kter]);
			(@temparray) = split(/\s+/, $distfiledata[$kter]);
			$temparray[0] //= "foo"; # placeholder for undefined $temparray[0]
			if ( $temparray[0] eq "corrdist" && $temparray[2] eq "han_tile" ) {
				$zter = $jter + ($iter*$taxlines);
				$han_distarray[$zter] = $temparray[3];
				$zter = $iter + ($jter*$taxlines);
				$han_distarray[$zter] = $temparray[3];
				$kter=$tempvar;
			}
		} # end hanada distances - for ($kter=0; $kter<$tempvar...
		
		# Extract the distance information for HP distances
		for ($kter=0; $kter<$tempvar; $kter++) {
			chomp($distfiledata[$kter]);
			(@temparray) = split(/\s+/, $distfiledata[$kter]);
			$temparray[0] //= "foo"; # placeholder for undefined $temparray[0]
			if ( $temparray[0] eq "corrdist" && $temparray[2] eq "hp_tile" ) {
				$zter = $jter + ($iter*$taxlines);
				$hp_distarray[$zter] = $temparray[3];
				$zter = $iter + ($jter*$taxlines);
				$hp_distarray[$zter] = $temparray[3];
				$kter=$tempvar;
			}
		} # end HP distances - for ($kter=0; $kter<$tempvar...
		
		# Extract the distance information for aa size distances
		for ($kter=0; $kter<$tempvar; $kter++) {
			chomp($distfiledata[$kter]);
			(@temparray) = split(/\s+/, $distfiledata[$kter]);
			$temparray[0] //= "foo"; # placeholder for undefined $temparray[0]
			if ( $temparray[0] eq "corrdist" && $temparray[2] eq "size_tile" ) {
				$zter = $jter + ($iter*$taxlines);
				$size_distarray[$zter] = $temparray[3];
				$zter = $iter + ($jter*$taxlines);
				$size_distarray[$zter] = $temparray[3];
				$kter=$tempvar;
			}
		} # end aa size distances - for ($kter=0; $kter<$tempvar...
		
		# Extract the distance information for GC distances
		for ($kter=0; $kter<$tempvar; $kter++) {
			chomp($distfiledata[$kter]);
			(@temparray) = split(/\s+/, $distfiledata[$kter]);
			$temparray[0] //= "foo"; # placeholder for undefined $temparray[0]
			if ( $temparray[0] eq "corrdist" && $temparray[2] eq "gc_tile" ) {
				$zter = $jter + ($iter*$taxlines);
				$gc_distarray[$zter] = $temparray[3];
				$zter = $iter + ($jter*$taxlines);
				$gc_distarray[$zter] = $temparray[3];
				$kter=$tempvar;
			}
		} # end GC distances - for ($kter=0; $kter<$tempvar...

	}
}

############################################################################
# Output the distances
############################################################################
open (my $NEXF, ">$outfile.nex") or die "Could not open file $outfile.nex for output.\n";

print $NEXF "#NEXUS\n\n";

#
# Amino acid distances
#
print $NEXF "BEGIN DISTANCES; \[ Amino acid tiled compression distances \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$taxlines";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$taxlines; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $NEXF "\t$aa_distarray[$zter]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! Amino acid tiled compression distances \]\n";
print $NEXF "nj treefile=$prefix.NJ_aa_compdist.tre replace;\n";
print $NEXF "END;\n\n";

#
# Dayhoff distances
#
print $NEXF "BEGIN DISTANCES; \[ Dayhoff tiled compression distances \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$taxlines";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$taxlines; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $NEXF "\t$day_distarray[$zter]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! Dayhoff tiled compression distances \]\n";
print $NEXF "nj treefile=$prefix.NJ_Dayhoff_compdist.tre replace;\n";
print $NEXF "END;\n\n";

#
# Hanada distances
#
print $NEXF "BEGIN DISTANCES; \[ Hanada tiled compression distances \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$taxlines";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$taxlines; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $NEXF "\t$han_distarray[$zter]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! Hanada tiled compression distances \]\n";
print $NEXF "nj treefile=$prefix.NJ_Hanada_compdist.tre replace;\n";
print $NEXF "END;\n\n";

#
# HP distances
#
print $NEXF "BEGIN DISTANCES; \[ HP tiled compression distances \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$taxlines";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$taxlines; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $NEXF "\t$hp_distarray[$zter]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! HP tiled compression distances \]\n";
print $NEXF "nj treefile=$prefix.NJ_HP_compdist.tre replace;\n";
print $NEXF "END;\n\n";

#
# Amino acid size distances
#
print $NEXF "BEGIN DISTANCES; \[ Amino acid size tiled compression distances \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$taxlines";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$taxlines; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $NEXF "\t$size_distarray[$zter]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! Amino acid size tiled compression distances \]\n";
print $NEXF "nj treefile=$prefix.NJ_aasize_compdist.tre replace;\n";
print $NEXF "END;\n\n";

#
# GC distances
#
print $NEXF "BEGIN DISTANCES; \[ GC tiled compression distances \]\n";
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$taxlines";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$taxlines; $iter++) {
	print $NEXF "$taxdata[$iter]";
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $NEXF "\t$gc_distarray[$zter]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "\[! GC tiled compression distances \]\n";
print $NEXF "nj treefile=$prefix.NJ_GC_compdist.tre replace;\n";
print $NEXF "END;\n\n";


close($NEXF) or die "Could not close file $outfile.nex.\n";

open (my $OUTF, ">$outfile.dist.txt") or die "Could not open file $outfile.dist.txt for output.\n";

print $OUTF "$taxlines\n";
for ($iter=0; $iter<$taxlines; $iter++) { print $OUTF "$taxdata[$iter]\n"; }
print $OUTF "aa_comp_dist";
for ($iter=0; $iter<$taxlines; $iter++) {
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $OUTF "\t$aa_distarray[$zter]";
	}
}
print $OUTF "\nday_comp_dist";
for ($iter=0; $iter<$taxlines; $iter++) {
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $OUTF "\t$day_distarray[$zter]";
	}
}
print $OUTF "\nhan_comp_dist";
for ($iter=0; $iter<$taxlines; $iter++) {
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $OUTF "\t$han_distarray[$zter]";
	}
}
print $OUTF "\nHP_comp_dist";
for ($iter=0; $iter<$taxlines; $iter++) {
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $OUTF "\t$hp_distarray[$zter]";
	}
}
print $OUTF "\nsize_comp_dist";
for ($iter=0; $iter<$taxlines; $iter++) {
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $OUTF "\t$size_distarray[$zter]";
	}
}
print $OUTF "\nGC_comp_dist";
for ($iter=0; $iter<$taxlines; $iter++) {
	for ($jter=0; $jter<$taxlines; $jter++) {
		$zter = $jter + ($iter*$taxlines);
		print $OUTF "\t$gc_distarray[$zter]";
	}
}

close($OUTF) or die "Could not close file $outfile.dist.txt.\n";

print "Data written to $outfile.nex and $outfile.dist.txt\n";
print "Exiting...\n\n";
