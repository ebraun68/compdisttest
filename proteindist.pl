#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# Initialize variables
############################################################################

my($progname) = $0;
my($iter);
my($jter);
my($kter);
my($lter);
my($mter);
my($zter);
my($tvar1);
my($tvar2);
my($tch);

my @tmparr1;
my @tmparr2;

if ( @ARGV != 4 ) {
	print "Usage:\n  \$ $progname <prefix> <alphabet> <mode> <distance>\n";
	print "  prefix   = prefix for input and output files\n";
	print "      use '.phy' as the extension for the input file\n";
	print "  alphabet = amino acid coding methods\n";
	print "      stdAA    -- standard amino acid alphabet\n";
	print "      univRY   -- universal code, binary nucleotides\n";
#	print "      univNT   -- universal code (A, C, G, T, and ambiguities)\n";
	print "      Dayhoff  -- six-state Dayhoff recoding\n";
	print "      SR6      -- six-state S&R recoding\n";
	print "      KGB6pam  -- six-state KGB recoding from PAM matrix\n";
	print "      KGB6wag  -- six-state KGB recoding from WAG matrix\n";
	print "      Hanada   -- four-state Hanada recoding\n";
	print "      HP       -- two-state hydrophobic-polar recoding\n";
	print "      size     -- two-state side chain size (large vs small) recoding\n";
	print "      GC       -- three-state GARP, FYMINK, other recoding\n";
	print "  mode     = run mode\n";
	print "      dist     -- generate nexus distance matrix and recoded alignment\n";
	print "      distonly -- generate nexus distance matrix only\n";
	print "      recode   -- generate recoded alignment (phylip format)\n";
	print "               -- NOTE: to use a separate authority file (taxon list) instead\n";
	print "                        of the input file simply append \"auth\" to the mode\n";
	print "                        (i.e, \"distauth\", \"distonlyauth\", and \"recodeauth\")\n";
	print "  distance = authority file name (optional, default to input file)\n";
	print "      pdist   -- p-distance\n";
	print "      Poisson -- Poisson correction\n";
	print "      gamma   -- Gamma correction, must provide alpha parameter\n";
	print "         use gamma=alpha where alpha is >0\n";
	print "exiting...\n";
	exit;
}

my($prefix)   = $ARGV[0];
my($suffix)   = $ARGV[1]; # the alphabet is stored as $suffix
my($mode)     = $ARGV[2];
my($disttype) = $ARGV[3];

# Set up the distance type
(@tmparr1) = split(/=/, $disttype);
my($distance) = -1; # -1 is used for p-distance; p-distance is the default
if ( lc($tmparr1[0]) eq "poisson" ) { $distance = 0; }
elsif ( lc($tmparr1[0]) eq "gamma" ) { $distance = $tmparr1[1]; }

# deal with the authority file
my($seqfname)  = "$prefix" . ".phy";
my($authfname) = "$prefix" . ".phy";
if ( lc($mode) eq "distauth" ) { $authfname = "$prefix" . ".authority.txt"; $mode = "dist"; }
if ( lc($mode) eq "distonlyauth" ) { $authfname = "$prefix" . ".authority.txt"; $mode = "distonly"; }
if ( lc($mode) eq "recodeauth" ) { $authfname = "$prefix" . ".authority.txt"; $mode = "recode"; }

# Set up the run mode flags
my($outputnewalingment) = 1;
my($outputdistmatrix)   = 1;
if ( lc($mode) eq "recode" ) {
	$outputdistmatrix = 0;
}
if ( lc($mode) eq "distonly" ) {
	$outputnewalingment = 0;
}
if ( lc($suffix) eq "stdaa" && lc($mode) ne "distonly" ) {
	print "Alphabet \"stdAA\" should only be used with \"distonly\" mode\n";
	print "exiting...\n\n";
	exit;
}

############################################################################
# Initialize arrays for codons and alternative alphabets
############################################################################

# Flag for whether the output will be codons or recoded amino acids
my($tripletflag) = 3; 
# amino acid alphabet = 1
# codons              = 3

# Amino acid numbers and universal code
my @aminoacid = ("F", "L", "I", "M", "V", "S", "P", "T", "A", "Y", "H", "Q", "N", "K", "D", "E", "C", "W", "R", "G");

# Universal code, with  R = 0; Y = 1
#  AA order:   F      L      I      M      V      S      P      T      A      Y      H      Q      N      K      D      E      C      W      R      G
my @univRY = ("111", "11?", "01?", "010", "01?", "???", "11?", "01?", "01?", "101", "101", "100", "001", "000", "001", "000", "101", "100", "?0?", "00?");
#  RY codon:   YYY    YYN    RYN    RYR    RYN    NNN    YYN    RYN    RYN    YRY    YRY    YRR    RRY    RRR    RRY    RRR    YRY    YRR    MRN    RRN
# 0.  F -- YYY (UUU, UUC)
# 1.  L -- YYN (UUA, UUG, CUN)
# 2.  I -- RYN (AUU, AUC, AUA)
# 3.  M -- RYR (ATG)
# 4.  V -- RYN (GUN)
		
# 5.  S -- NNN (UCN, AGU, AGC)
# 6.  P -- YYN (CCN)
# 7.  T -- RYN (ACN)
# 8.  A -- RYN (GCN)
		
# 9.  Y -- YRY (UAU, UAC)
# 10. H -- YRY (CAU, CAC)
# 11. Q -- YRR (CAA, CAG)
# 12. N -- RRY (AAU, AAC)
# 13. K -- RRR (AAA, AAG)
# 14. D -- RRY (GAU, GAC)
# 15. E -- RRR (GAA, GAG)
		
# 16. C -- YRY (UGU, UGC)
# 17. W -- YRR (UGG)
# 18. R -- NRN (CGN, AGA, AGG)
# 19. G -- RRN (GGN)

# Universal code (nucleotides with ambiguity)
#  AA order:   F      L      I      M      V      S      P      T      A      Y      H      Q      N      K      D      E      C      W      R      G
my @univNT = ("TTY", "YT?", "AT?", "ATG", "GT?", "???", "CC?", "AC?", "GC?", "TAY", "CAY", "CAR", "AAY", "AAR", "GAY", "GAR", "TGY", "TGG", "?G?", "GG?");

# stdAA: no recoding
#  AA order:    F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @stdaa   = ("F", "L", "I", "M", "V", "S", "P", "T", "A", "Y", "H", "Q", "N", "K", "D", "E", "C", "W", "R", "G");

# Dayhoff codes (six state)
#  AA order:    F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @dayhoff = ("5", "4", "4", "4", "4", "1", "1", "1", "1", "5", "3", "2", "2", "3", "2", "2", "0", "5", "3", "1");
#   0 = C
#   1 = AGPST
#   2 = NDEQ
#   3 = RHK
#   4 = ILMV
#   5 = FWY

# S&R6 saturation codes (six state) 
# from Susko, E., and Roger, A.J. (2007). On reduced amino acid alphabets for phylogenetic 
# inference. Mol. Biol. Evol. 24, 2139–2150.
#  AA order:    F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @sr6sat  = ("5", "3", "3", "3", "3", "0", "0", "0", "0", "5", "0", "2", "1", "2", "1", "1", "4", "4", "2", "1");
#   0 = APST
#   1 = DENG
#   2 = QKR
#   3 = ILMV
#   4 = WC
#   5 = FYH

# KGB PAM codes (six state) 
# from Kosiol, C., Goldman, N., and Buttimore, N.H. (2004). A new criterion and method for 
# amino acid classification. J. Theor. Biol. 228, 97–106.
#  AA order:     F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @kgb6pam  = ("4", "2", "2", "2", "5", "0", "0", "1", "0", "4", "1", "1", "1", "1", "1", "1", "5", "3", "1", "0");
#   0 = AGPS
#   1 = DENQHKRT
#   2 = MIL
#   3 = W
#   4 = FY
#   5 = CV

# KGB WAG codes (six state) 
# from Kosiol, C., Goldman, N., and Buttimore, N.H. (2004). A new criterion and method for 
# amino acid classification. J. Theor. Biol. 228, 97–106.
#  AA order:     F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @kgb6wag  = ("4", "3", "3", "3", "2", "0", "0", "0", "0", "4", "0", "0", "0", "0", "0", "0", "2", "5", "0", "1");
#   0 = HRKQNEDSTPA
#   1 = G
#   2 = CV
#   3 = IML
#   4 = FY
#   5 = W

# Hanada codes (four state, recoded as nucleotides)
#  AA order:   F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @hanada = ("G", "C", "C", "C", "C", "A", "A", "A", "A", "G", "G", "G", "A", "G", "T", "T", "A", "G", "G", "A");
#   A = ANCGPST
#   C = ILMV
#   G = RQHKFWY
#   T = DE

# Two-state HP (hydrophobic-polar) alphabet
#   Braun, E. L. (2018). An evolutionary model motivated by physicochemical properties of
#   amino acids reveals variation among proteins. Bioinformatics, 34(13), i350-i356. (Fig. 1)
#  AA order:   F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @hp     = ("1", "1", "1", "1", "1", "0", "1", "0", "1", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "1");
#   0 = RNDCQEHKSTWY
#   1 = AGILMFPV

# Two-state size (large-small) alphabet
#   Braun, E. L. (2018). An evolutionary model motivated by physicochemical properties of
#   amino acids reveals variation among proteins. Bioinformatics, 34(13), i350-i356. (Fig. 1)
#
#  AA order:       F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @largesmall = ("0", "0", "0", "0", "1", "1", "1", "1", "1", "0", "0", "0", "1", "0", "1", "0", "1", "0", "0", "1");
#   0 = RQEHILKMFWY
#   1 = ANDCGPSTV

# Three state GC-AT alphabet
#   Singer, G. A., & Hickey, D. A. (2000). Nucleotide bias causes a genomewide bias in the 
#   amino acid composition of proteins. Molecular biology and evolution, 17(11), 1581-1588.
#
#  AA order:       F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @gcat       = ("0", "1", "0", "0", "1", "1", "2", "1", "2", "0", "1", "1", "0", "0", "1", "1", "1", "1", "2", "2");
#   0 = FYMINK (AT-rich)
#	1 = LVSTHQDECW
#	2 = GARP (GC-rich)

my @alphabet;	my($nstates) = 20;
# NOTE: default is no recoding (stdAA)
if ( lc($suffix) eq "univry" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $univRY[$iter]; } $nstates = 2; }
elsif ( lc($suffix) eq "dayhoff" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $dayhoff[$iter]; } $tripletflag = 1; $nstates = 6; }
elsif ( lc($suffix) eq "sr6" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $sr6sat[$iter]; } $tripletflag = 1; $nstates = 6; }
elsif ( lc($suffix) eq "kgb6pam" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $kgb6pam[$iter]; } $tripletflag = 1; $nstates = 6; }
elsif ( lc($suffix) eq "kgb6wag" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $kgb6wag[$iter]; } $tripletflag = 1; $nstates = 6; }
elsif ( lc($suffix) eq "hanada" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $hanada[$iter]; } $tripletflag = 1; $nstates = 4; }
elsif ( lc($suffix) eq "hp" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $hp[$iter]; } $tripletflag = 1; $nstates = 2; }
elsif ( lc($suffix) eq "size" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $largesmall[$iter]; } $tripletflag = 1; $nstates = 2; }
elsif ( lc($suffix) eq "gc" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $gcat[$iter]; } $tripletflag = 1;$nstates = 3; }
else { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $stdaa[$iter]; } $tripletflag = 1; $suffix = "stdAA"; } # default


############################################################################
# Read the authority file
############################################################################
my @authlist;
open (my $AUTHF, "$authfname") or die "Could not open file < $authfname > for input.\n";
@authlist = <$AUTHF>; # Read the authority file
close($AUTHF) or die "Could not close file < $authfname >\n";
my($authlistnum) = $#authlist + 1;

my($ntax);
my @taxname;
my($maxnamelen) = 10;

chomp($authlist[0]);
($ntax) = split(/\s+/, $authlist[0]);
# print "$ntax\n";
for ($iter=0; $iter<$ntax; $iter++) {
	$jter = $iter + 1;
	chomp($authlist[$jter]);
	($taxname[$iter]) = split(/\s+/, $authlist[$jter]);
	if ( length($taxname[$iter]) > $maxnamelen ) { $maxnamelen = length($taxname[$iter]); }
#	print "  $taxname[$iter] -- $maxnamelen\n";
}

############################################################################
# Read the sequence file
############################################################################
my @seqlist;
open (my $SEQF, "$seqfname") or die "Could not open file < $seqfname > for input.\n";
@seqlist = <$SEQF>; # Read the relaxed phylip format input file
close($SEQF) or die "Could not close file < $seqfname >\n";
my($seqlistnum) = $#seqlist + 1;
for ($iter=0; $iter<$seqlistnum; $iter++) { chomp($seqlist[$iter]); }

my($nchar);

($tvar1,$nchar) = split(/\s+/, $seqlist[0]);
my($totalchar) = $nchar * $tripletflag;

############################################################################
# Generate the recoded data, output the data if $outputnewalingment == 1
############################################################################
my($taxonfound);	my($padding);
my($seqname);		my($sequence);
my @recodedseq;

# output the information about the alphabet
if ( lc($suffix) eq "stdaa" ) { print "Distance analysis conducted without recoding amino acid data\n"; }
else { print "Distance analysis conducted using $suffix recoding of amino acid data\n"; }

# set up counters for the state frequencies
my($definedcharacters) = 0;
my @statecounts; 	my @statefreqs; 
for ($iter=0; $iter<$nstates; $iter++) { $statecounts[$iter] = 0; }

my($outfname) = "$prefix" . "." . "$suffix" . ".phy";
#open (my $OUTF, ">$outfname") or die "Could not open file < $outfname > for output.\n";
#print $OUTF "$ntax $totalchar\n";

# iterate through taxon names from authority file
for ($iter=0; $iter<$ntax; $iter++) {
	
	# output the name and pad with spaces
#	print $OUTF "$taxname[$iter]";
#	$padding = $maxnamelen - length($taxname[$iter]) + 2;
#	for ($jter=0; $jter<$padding; $jter++) { print $OUTF " "; }
	
	# search for the taxon name
	$taxonfound = 0;
	for ($kter=1; $kter<$seqlistnum; $kter++) {
		($seqname,$sequence) = split(/\s+/, $seqlist[$kter]);
		if ( $seqname eq $taxname[$iter] ) {
			$recodedseq[$iter] = recodeprotein($sequence);
			
			$taxonfound = 1;
			# count the occurrences of each state
			$tvar1 = length($recodedseq[$iter]);
			for ($zter=0; $zter<$tvar1; $zter++) {
				$tch = substr($recodedseq[$iter],$zter,1);
				$tvar2 = state2value($tch);
				if ( $tvar2 > -1 ) {
					$statecounts[$tvar2]++;
					$definedcharacters++;
				}
			}
			$kter = $seqlistnum;
		}
	}
	
	# generate data for taxa that were not found
	if ( $taxonfound == 0 ) {
		$recodedseq[$iter] = "?";
		for ($lter=1; $lter<$totalchar; $lter++) { 
			$recodedseq[$iter] = "$recodedseq[$iter]" . "?";
#			print $OUTF "?"; 
		}
#		print $OUTF "\n";
	}
}

# Output recoded matrix in relaxed phylip format if $outputnewalingment == 1
if ( $outputnewalingment == 1 ) {
#	my($outfname) = "$prefix" . "." . "$suffix" . ".phy";
	open (my $OUTF, ">$outfname") or die "Could not open file < $outfname > for output.\n";
	print $OUTF "$ntax $totalchar\n";

	
	for ($iter=0; $iter<$ntax; $iter++) {
	
		# output the name and pad with spaces
		print $OUTF "$taxname[$iter]";
		$padding = $maxnamelen - length($taxname[$iter]) + 2;
		for ($jter=0; $jter<$padding; $jter++) { print $OUTF " "; }
		
		# output the sequence
		print $OUTF "$recodedseq[$iter]\n";
	}

	close($OUTF) or die "Could not close file < $outfname >\n";
	print "Recoded data file exported as $outfname\n";
}

my($b) = 1;
print "Data matrix includes $definedcharacters defined characters\n";
print "State\tState_count\tState_freq\n";
for ($iter=0; $iter<$nstates; $iter++) { 
	$statefreqs[$iter] = $statecounts[$iter] / $definedcharacters;
	print "$iter\t$statecounts[$iter]\t$statefreqs[$iter]\n";
	$b = $b - ( $statefreqs[$iter] * $statefreqs[$iter] );
}
print "Equal input saturation parameter b = $b\n";


############################################################################
# Calculate and output the distance estimates ( if $outputdistmatrix == 0 )
############################################################################
# exit if $outputdistmatrix == 0
if ( $outputdistmatrix == 0 ) { exit; }

print "\nEstimating distances: ";
if ( $distance == -1 ) { print "p-distances\n"; }
elsif ( $distance == 0 ) { print "Poisson distances\n"; }
else { print "gamma distances with alpha = $distance\n"; }


my @pd; # array for p-distance [0] and proportion of defined sites [1]
my @pdistmatrix;
my @defmatrix;
my @corrdistmat;	my($p_over_b);
my($probdist) = 0; # counter for undefined distances

# Generate the distance matrix
for ($iter=0; $iter<$ntax; $iter++) {
	for ($jter=0; $jter<=$iter; $jter++) {
		@pd = pdist($recodedseq[$iter],$recodedseq[$jter]);
		$pdistmatrix[matcoord($iter,$jter)] = $pd[0];
		$pdistmatrix[matcoord($jter,$iter)] = $pd[0];
		$defmatrix[matcoord($iter,$jter)]   = $pd[1];
		$defmatrix[matcoord($jter,$iter)]   = $pd[1];
		if ( $pd[0] != -1 ) {
			if ( $distance == -1 ) { 
				$corrdistmat[matcoord($iter,$jter)] = $pd[0];
				$corrdistmat[matcoord($jter,$iter)] = $pd[0];
			}
			elsif ( $distance == 0 ) { # Dist: d = -b ln( 1 - p/b )
				$p_over_b = $pd[0] / $b;
				if ( $p_over_b >= 1.0 ) { $p_over_b = 0.9999999999999999; } # set maximum for p/b			
				$corrdistmat[matcoord($iter,$jter)] = -1.0 * $b * log(1.0-$p_over_b);
				$corrdistmat[matcoord($jter,$iter)] = $corrdistmat[matcoord($iter,$jter)]
			}
			else { # Gamma distance: d = ba [ (1-p/b)^-1/a - 1 ]
				$p_over_b = $pd[0] / $b;
				if ( $p_over_b >= 1.0 ) { $p_over_b = 0.9999999999999999; } # set maximum for p/b			
				$corrdistmat[matcoord($iter,$jter)] = $b * $distance * ( ((1.0-$p_over_b)**(-1.0/$distance)) - 1.0 );
				$corrdistmat[matcoord($jter,$iter)] = $corrdistmat[matcoord($iter,$jter)]
			}
		}
		else { 
			$corrdistmat[matcoord($iter,$jter)] = -1;
			$corrdistmat[matcoord($jter,$iter)] = -1;
			$probdist++; # increment counter for undefined distances
		}
	}
}
if ( $probdist > 0 ) {
	print "\nWARNING: $probdist distance estimates are undefined\n";
	print "         distance matrix cannot be used as is for tree estimation unless\n";
	print "         some taxa are excluded\n\n";
}

# Output a nexus distance matrix
my($distoutfile);
if ( $distance == -1 ) { $distoutfile = "$prefix" . "." . "$suffix" . ".pdist.nex"; }
elsif ( $distance == 0 ) { $distoutfile = "$prefix" . "." . "$suffix" . ".Poisson.nex"; }
else { $distoutfile = "$prefix" . "." . "$suffix" . ".gamma. " . "$distance" . ".nex"; $distoutfile =~ tr/ //ds; }

open (my $NEXF, ">$distoutfile") or die "Could not open file $distoutfile for output.\n";

print $NEXF "#NEXUS\n\n";
if ( $distance == -1 ) { print $NEXF "\[ Amino acid p-distances with $suffix coding \]\n"; }
elsif ( $distance == 0 ) { print $NEXF "\[ Amino acid Poisson distances with $suffix coding \]\n"; }
else { print $NEXF "\[ Amino acid gamma distances with $suffix coding, alpha = $distance \]\n"; }
print $NEXF "BEGIN DISTANCES;\n"; 
print $NEXF "\tDIMENSIONS NEWTAXA NTAX=$ntax";
print $NEXF ";\n\tFORMAT TRIANGLE=Both;\n";
print $NEXF "MATRIX\n";

for ($iter=0; $iter<$ntax; $iter++) {
	print $NEXF "$taxname[$iter]";
	for ($jter=0; $jter<$ntax; $jter++) {
		$tvar2=matcoord($iter,$jter);
		print $NEXF "\t$corrdistmat[$tvar2]";
	}
	print $NEXF "\n";
}

print $NEXF "\t;\nEND;\n\n";

my($disttreefile);
if ( $distance == -1 ) { $disttreefile = "$prefix" . "." . "$suffix" . ".pdist.tre"; }
elsif ( $distance == 0 ) { $disttreefile = "$prefix" . "." . "$suffix" . ".Poisson.tre"; }
else { $disttreefile = "$prefix" . "." . "$suffix" . ".gamma. " . "$distance" . ".tre"; $disttreefile =~ tr/ //ds; }

print $NEXF "[ This block can be commented out by adding a square bracket before  ]\n";
print $NEXF "[ 'BEGIN PAUP' and eliminating the square bracket before the 'end of ]\n";
print $NEXF "[ PAUP block'                                                        ]\n";
print $NEXF "BEGIN PAUP;\n";
print $NEXF "hsearch;\n";
print $NEXF "describetrees 1/ plot=phylogram;\n";
print $NEXF "savetrees file=$disttreefile format=altnexus brlens=yes replace;\n";
print $NEXF "log stop;\n";
print $NEXF "END; [end of PAUP block]\n\n";

close($NEXF) or die "Could not close file $distoutfile\n";

# Output a phylip distance matrix
if ( $distance == -1 ) { $distoutfile = "$prefix" . "." . "$suffix" . ".pdist.phydist.txt"; }
elsif ( $distance == 0 ) { $distoutfile = "$prefix" . "." . "$suffix" . ".Poisson.phydist.txt"; }
else { $distoutfile = "$prefix" . "." . "$suffix" . ".gamma. " . "$distance" . ".phydist.txt"; $distoutfile =~ tr/ //ds; }

open (my $PHYF, ">$distoutfile") or die "Could not open file $distoutfile for output.\n";

print $PHYF "$ntax\n";

for ($iter=0; $iter<$ntax; $iter++) {
	print $PHYF "$taxname[$iter]";
	for ($jter=0; $jter<$ntax; $jter++) {
		$tvar2=matcoord($iter,$jter);
		print $PHYF "\t$corrdistmat[$tvar2]";
	}
	print $PHYF "\n";
}

close($PHYF) or die "Could not close file $distoutfile\n";

# Output a matrix with information on the number of defined sites
$distoutfile = "$prefix" . "." . "$suffix" . ".defined_sites.txt";

open (my $DEFF, ">$distoutfile") or die "Could not open file $distoutfile for output.\n";

$tvar1 = length($recodedseq[0]);
print $DEFF "$ntax\t$tvar1\n"; # output number of taxa and number of sites

for ($iter=0; $iter<$ntax; $iter++) {
	print $DEFF "$taxname[$iter]";
	for ($jter=0; $jter<$ntax; $jter++) {
		$tvar2=matcoord($iter,$jter);
		print $DEFF "\t$defmatrix[$tvar2]";
	}
	print $DEFF "\n";
}

close($DEFF) or die "Could not close file $distoutfile\n";

print "Data written to output files. Run complete...\n\n";

exit;

############################################################################
# Subroutine to recode the protein sequence, either by reverse translating
# or using a reduced alphabet
############################################################################
sub recodeprotein {
	my($localinseq) = @_;
	my($aalen) = length($localinseq);

	my $localoutseq;
	my($aa);	my($codon);
	my($locfirstaa) = 1;
	
	my($allmissing) = "???";
	my($allgap) = "---";
	if ( $tripletflag == 1 ) { $allmissing = "?"; $allgap = "-"; }

	for ($zter=0; $zter<$aalen; $zter++) {

		$aa = substr($localinseq,$zter,1);
		
		# then back translate the amino acid
		$codon = "$allmissing";
		if ( $aa eq "-" ) { $codon = "$allgap"; }
		elsif ( uc($aa) eq "F" ) { $codon = "$alphabet[0]";  } # amino acid 0
		elsif ( uc($aa) eq "L" ) { $codon = "$alphabet[1]";  } # amino acid 1
		elsif ( uc($aa) eq "I" ) { $codon = "$alphabet[2]";  } # amino acid 2
		elsif ( uc($aa) eq "M" ) { $codon = "$alphabet[3]";  } # amino acid 3
		elsif ( uc($aa) eq "V" ) { $codon = "$alphabet[4]";  } # amino acid 4
		
		elsif ( uc($aa) eq "S" ) { $codon = "$alphabet[5]";  } # amino acid 5
		elsif ( uc($aa) eq "P" ) { $codon = "$alphabet[6]";  } # amino acid 6
		elsif ( uc($aa) eq "T" ) { $codon = "$alphabet[7]";  } # amino acid 7
		elsif ( uc($aa) eq "A" ) { $codon = "$alphabet[8]";  } # amino acid 8
		
		elsif ( uc($aa) eq "Y" ) { $codon = "$alphabet[9]";  } # amino acid 9
		elsif ( uc($aa) eq "H" ) { $codon = "$alphabet[10]"; } # amino acid 10
		elsif ( uc($aa) eq "Q" ) { $codon = "$alphabet[11]"; } # amino acid 11
		elsif ( uc($aa) eq "N" ) { $codon = "$alphabet[12]"; } # amino acid 12
		elsif ( uc($aa) eq "K" ) { $codon = "$alphabet[13]"; } # amino acid 13
		elsif ( uc($aa) eq "D" ) { $codon = "$alphabet[14]"; } # amino acid 14
		elsif ( uc($aa) eq "E" ) { $codon = "$alphabet[15]"; } # amino acid 15
		
		elsif ( uc($aa) eq "C" ) { $codon = "$alphabet[16]"; } # amino acid 16
		elsif ( uc($aa) eq "W" ) { $codon = "$alphabet[17]"; } # amino acid 17
		elsif ( uc($aa) eq "R" ) { $codon = "$alphabet[18]"; } # amino acid 18
		elsif ( uc($aa) eq "G" ) { $codon = "$alphabet[19]"; } # amino acid 19
		
		if ( $locfirstaa == 1 ) { 
			$localoutseq = "$codon";
			$locfirstaa = 0;
		}
		else { 
			$localoutseq = "$localoutseq" . "$codon";
		}
	}
	
	return $localoutseq;
}

############################################################################
# Subroutine to recode the protein sequence, either by reverse translating
# or using a reduced alphabet
############################################################################
sub state2value {
	my($instate) = $_[0];
	
	my($stateval) = -1;
	
	if ( lc($suffix) eq "stdaa" ) {
		if    ( uc($instate) eq "F" ) { $stateval = 0;  } # amino acid 0
		elsif ( uc($instate) eq "L" ) { $stateval = 1;  } # amino acid 1
		elsif ( uc($instate) eq "I" ) { $stateval = 2;  } # amino acid 2
		elsif ( uc($instate) eq "M" ) { $stateval = 3;  } # amino acid 3
		elsif ( uc($instate) eq "V" ) { $stateval = 4;  } # amino acid 4
		
		elsif ( uc($instate) eq "S" ) { $stateval = 5;  } # amino acid 5
		elsif ( uc($instate) eq "P" ) { $stateval = 6;  } # amino acid 6
		elsif ( uc($instate) eq "T" ) { $stateval = 7;  } # amino acid 7
		elsif ( uc($instate) eq "A" ) { $stateval = 8;  } # amino acid 8
		
		elsif ( uc($instate) eq "Y" ) { $stateval = 9;  } # amino acid 9
		elsif ( uc($instate) eq "H" ) { $stateval = 10; } # amino acid 10
		elsif ( uc($instate) eq "Q" ) { $stateval = 11; } # amino acid 11
		elsif ( uc($instate) eq "N" ) { $stateval = 12; } # amino acid 12
		elsif ( uc($instate) eq "K" ) { $stateval = 13; } # amino acid 13
		elsif ( uc($instate) eq "D" ) { $stateval = 14; } # amino acid 14
		elsif ( uc($instate) eq "E" ) { $stateval = 15; } # amino acid 15
		
		elsif ( uc($instate) eq "C" ) { $stateval = 16; } # amino acid 16
		elsif ( uc($instate) eq "W" ) { $stateval = 17; } # amino acid 17
		elsif ( uc($instate) eq "R" ) { $stateval = 18; } # amino acid 18
		elsif ( uc($instate) eq "G" ) { $stateval = 19; } # amino acid 19
	}
	elsif ( lc($suffix) eq "hanada" ) {
		if    ( uc($instate) eq "A" ) { $stateval = 0;  } # Hanada state 0
		elsif ( uc($instate) eq "C" ) { $stateval = 1;  } # Hanada state 1
		elsif ( uc($instate) eq "G" ) { $stateval = 2;  } # Hanada state 2
		elsif ( uc($instate) eq "T" ) { $stateval = 3;  } # Hanada state 3
	}
	else { 
		if ( $instate ne "-" && $instate ne "?" ) { $stateval = int($instate); }
	}
	
	return $stateval;
}

############################################################################
# Convert coordinates for 2D array to use with 1D array
############################################################################
sub matcoord {
	my($coord1,$coord2) = @_;
	return (($coord1*$ntax)+$coord2);
}

############################################################################
# Subroutine to calculate the p-distance and number of defined amino acids
# for two protein sequences
############################################################################
sub pdist {
	my($sequence0,$sequence1) = @_;
	my($seqlen) = length($sequence1);
	my @locaa;
	
	my($locdef);
	my($locdefnum)  = 0;
	my($locdiff)    = 0;
	
	for ($mter=0; $mter<$seqlen; $mter++) {
		$locaa[0] = substr($sequence0,$mter,1);
		$locaa[1] = substr($sequence1,$mter,1);
		$locdef = 1;
	#	print "$mter-$locaa[0]-$locaa[1] ";
		if ( $locaa[0] eq "?" || $locaa[0] eq "-" ) { $locdef = 0; }
		if ( $locaa[1] eq "?" || $locaa[1] eq "-" ) { $locdef = 0; }
		if ( $locdef == 1 ) {
			if ( $locaa[0] ne $locaa[1] ) { 
				$locdiff++;
			}
			$locdefnum++;
		}
	}
	
	my @retp;
	if ( $locdefnum == 0 ) {
		$retp[0] = -1; # flag for undefined p-distance
		$retp[1] = 0;  # proportion of defined sites
	}
	else {
		$retp[0] = $locdiff / $locdefnum; # p-distance
		$retp[1] = $locdefnum / $seqlen;  # proportion of defined sites
	}
	
	return @retp;
}

