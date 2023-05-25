# compdisttest
Programs to test the behavior of compression distances in phylogenomic analyses of proteins

```
################################################################################
#
# Programs for compression distance calculations
#	- all programs written in perl
#	- PAUP required for distance analyses
#	- the goal of these programs is proof of concept and method exploration so
#     they have not been optimized for speed
#	- these programs do limited checking of the command line arguments and this
#	  can lead to errors; check input arguments carefully when run failures occur
#
################################################################################
```

Workflow --	the programs are run sequentially with the output from one programs used as
			input for the next program:
	1. shuffle_locuslist.pl (optional)
	2. generate_distance_matrices.pl (requires pairwise_compression_dist.pl)
		- pairwise_compression_dist.pl calls the compressor
	3. generate_corrected_distances.pl
	4. Estimate trees in PAUP

Input --	requires this information:
	1. a text file with paths to a set of protein fasta files for analysis
		- single line protein fasta files
	2. a list of taxa to be included in analyses

```
################################################################################
# Detailed workflow for each program
################################################################################
```
**********
1. shuffle_locuslist.pl
	- purpose: to shuffle, subsample, and jackknife a list of fasta files
	- the fasta file list should be named <prefix>.sequences.txt
```	
Usage:
  $ ./shuffle_locuslist.pl <locuslist> <outlist> <samplesize> <randseed>
  locuslist  = text file listing translated loci to analyze
  outlist    = rearranged text file listing proteins seqeunces
  samplesize = number of proteins seqeunces to sample
       -1 -- jackknife with deletion probability of exp(-1)
        0 -- output the full number of proteins
        n -- output a list of n files (if n > the total number
             of file then n = total number of sequences)
  randseed   = random integer to seed the random number generator
        use 0 for implicit srand() call
```
	
The format of the locus list is simple - a text file listing fasta files for analysis. The
full path to each fasta file is listed on each line, like this example:

/Users/edwardbraun/OrthoMam_align_v10c/omm_filtered_AA_CDS/ENSG00000000457_SCYL3_AA.fasta
/Users/edwardbraun/OrthoMam_align_v10c/omm_filtered_AA_CDS/ENSG00000001084_GCLC_AA.fasta
/Users/edwardbraun/OrthoMam_align_v10c/omm_filtered_AA_CDS/ENSG00000001460_STPG1_AA.fasta
/Users/edwardbraun/OrthoMam_align_v10c/omm_filtered_AA_CDS/ENSG00000001561_ENPP4_AA.fasta
etc.

The output file is the same type of list, but it is shuffled and truncated based on the
the <samplesize> variable. Note that jackknifing (samplesize == -1) produces files with
different sizes when the program is run with different random number seeds.

The output file has the same format as the input file; the output file is intended to be
used as the input for generate_distance_matrices.pl (step 2 of this workflow) so it should
end with ".sequences.txt"

**********
2. generate_distance_matrices.pl
	- purpose: to generate compression distances for all pairs of taxa
	- requires the path to an executable version of pairwise_compression_dist.pl
	- the path to pairwise_compression_dist.pl is a variable ($pairwise_exec) near
	  the beginning of this program
	- the compressor is set by a variable ($compexec) near the beginning of the
	  pairwise_compression_dist.pl program

```
Usage:
  $ ./generate_distance_matrices.pl <taxalist> <prefix> <outfile>
  taxalist = text file listing taxa to analyze
  prefix   = prefix for file list
  outfile  = prefix for output files
Expected format for list of files:
  <prefix>.sequences.txt = list of sequence files
Format of the output files:
  <outfile>.nex      = nexus distance matrices
  <outfile>.dist.txt = list of evolutionary distances
```
	  
The nexus outfile has embedded PAUP blocks that generate neighbor joining trees from
each distance matrix when the file is executed in PAUP.

The  ".dist.txt" outfile contains the distance matrices for all amino acid codings in
a semi-tabular format:

number_of_taxa
taxname_1
taxname_2
taxname_3
...
taxname_N
aa_comp_dist	0	dist1	dist2 ...
etc. (total of six compression distances)

The distances are tab-delimited lines with all elements of a square matrix listed after
the distance type (e.g., if your analyses produce a four-taxon distance matrix, the first
numerical value - element 1 if the line is split into an array - corresponds to element 
0,0 of a square 4x4 matrix, this is followed by element 2 = 0,1; element 3 = 0,2; element
4 = 0,3; element 5 = 1,0; element 6 = 1,1; and so forth up to element 16 = 3,3).

The following is an example of a ".dist.txt" outfile for a four-taxon analysis:

```
4
Monodelphis_domestica
Canis_familiaris
Mus_musculus
Homo_sapiens
aa_comp_dist		0	0.474326888779173	0.505570201770356	0.479603580165895	0.474326888779173	0	0.38080570297089	0.294460853905096	0.505570201770356	0.38080570297089	0	0.360436073449617	0.479603580165895	0.294460853905096	0.360436073449617	0
day_comp_dist		0	0.412633799990905	0.440005109480969	0.420277853701423	0.412633799990905	0	0.323514774315049	0.252166972225877	0.440005109480969	0.323514774315049	0	0.306376331184798	0.420277853701423	0.252166972225877	0.306376331184798	0
han_comp_dist		0	0.396327922954512	0.421498045872989	0.40630872911347	0.396327922954512	0	0.304957146406577	0.238909628050055	0.421498045872989	0.304957146406577	0	0.289205671836408	0.40630872911347	0.238909628050055	0.289205671836408	0
HP_comp_dist		0	0.531125985722307	0.553484836682509	0.532334374037801	0.531125985722307	0	0.425349559710654	0.33869615060523	0.553484836682509	0.425349559710654	0	0.406302891660676	0.532334374037801	0.33869615060523	0.406302891660676	0
size_comp_dist		0	0.541964872449069	0.573561299086416	0.542881158546323	0.541964872449069	0	0.447652590105229	0.347168265158443	0.573561299086416	0.447652590105229	0	0.421254181278427	0.542881158546323	0.347168265158443	0.421254181278427	0
GC_comp_dist		0	0.597126465763397	0.624693355443062	0.599907492568409	0.597126465763397	0	0.503542380154037	0.410043650499122	0.624693355443062	0.503542380154037	0	0.484965459693378	0.599907492568409	0.410043650499122	0.484965459693378	0
```
	  
The ".dist.txt" file is input for the distance correction program.

**********
3. generate_corrected_distances.pl
	- purpose: apply a gamma-correction to distances in a ".dist.txt" file

```
Usage:
  $ generate_corrected_distances.pl <distlist> <alpha> <outfile> <prefix>
  distlist = text file listing evolutionary distances
  alpha    = alpha parameter for distance correction
    -- use negative number for no correction
    -- use 0 for Poisson-type correction
  outfile  = nexus format distance matrix
  prefix   = prefix for the treefiles
```

<distlist> is the ".dist.txt" file written by generate_distance_matrices.pl. The alpha
parameter for a gamma correction is passed on the command line. This program generates a
nexus format outfile identical in format to the generate_distance_matrices.pl nexus 
outfile (i.e., square distance matrices with embedded PAUP blocks).

If we call the distances output by generate_distance_matrices.pl "$P" (note that $P is 
expected to range from 0 to 1) then passing 0 will result in a Poisson-like correction:

	$dist = -1 * log( 1 - $P )

otherwise a gamma correction will be applied:

	$dist = $alpha * (((1-$P)**(-1/$alpha))-1)
	
PAUP blocks that perform a heuristic search for the optimal tree are embedded in the file.


-------------
# Additional program

proteindist.pl 
	- purpose: generate nexus distance matrix for protein alignments
	- used for comparisons with the compression distances
	- can calculate distances for many alternative amino acid alphabets
	- can be used to generate phylip files with the alternative alphabets 
	
```
Usage:
  $ proteindist.pl <prefix> <alphabet> <mode> <distance>
  prefix   = prefix for input and output files
      use '.phy' as the extension for the input file
  alphabet = amino acid coding methods
      stdAA    -- standard amino acid alphabet
      univRY   -- universal code, binary nucleotides
      Dayhoff  -- six-state Dayhoff recoding
      SR6      -- six-state S&R recoding
      KGB6pam  -- six-state KGB recoding from PAM matrix
      KGB6wag  -- six-state KGB recoding from WAG matrix
      Hanada   -- four-state Hanada recoding
      HP       -- two-state hydrophobic-polar recoding
      size     -- two-state side chain size (large vs small) recoding
      GC       -- three-state GARP, FYMINK, other recoding
  mode     = run mode
      dist     -- generate nexus distance matrix and recoded alignment
      distonly -- generate nexus distance matrix only
      recode   -- generate recoded alignment (phylip format)
               -- NOTE: to use a separate authority file (taxon list) instead
                        of the input file simply append "auth" to the mode
                        (i.e, "distauth", "distonlyauth", and "recodeauth")
  distance = authority file name (optional, default to input file)
      pdist   -- p-distance
      Poisson -- Poisson correction
      gamma   -- Gamma correction, must provide alpha parameter
         use gamma=alpha where alpha is >0
```

