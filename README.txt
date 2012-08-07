This is a pipeline for processing amplicon sequence data using the DADA algorithm. 

I. UPSTREAM of DADA. the dada algorithm takes in data in a format we call .uniques. this a simple and useful way to succinctly store deep sequence data sets. it is:
READS_1 \t SEQUENCE_1 \n
READS_2 \t SEQUENCE_2 \n
...
where READS_i is the is the number of identical copies of each unique sequence SEQUENCE_i, given as an integer (\t and \n are tabs and newlines). We provide a few PERL scripts to get your data into this format:
	i. fasta2seqs.pl. strips away all header information and newlines that are found in .fasta files and creates files of the format: 
SEQUENCE1\n
SEQUENCE2\n
...
which we call .seqs file. fasta2seqs can be called with the name of a folder of .fasta files or a single .fasta file. for example:
"perl fasta2seqs.pl ~/myfastas" or "perl fasta2seqs.pl ~.myfasta.fasta"
A new file or folder will be created with a _seqs appended to its name
	ii. lengthTrim.pl: this takes in a folder of .seqs files, trims all sequences down to a particular length, and creates a folder of new .seqs files. sequences shorter than the given length are removed entirely. it is called in the following way:
"perl lengthTrim.pl ~/myfastas_seqs 100"
which would create the folder "~/myfastas_seqs100" containing in each .seqs file only 100bp sequences
	iii. clusterUniques.pl: this takes in a folder of .seqs files (or single such file, either length trimmed or not), and outputs a folder (or file) of unique sequences and the number of reads of each in the format described above for use in dada. this folder will have '_uniques' appended to its name.

II. DADA. 
Before running DADA it is necessary to compile the two .c files, align2dom_homo.c and align2dom.c in MATLAB. this is trivial -- open matlab and type 'mex align2dom.c' and 'mex align2dom_homo.c' into the command window, which will create executable files with a suffix dependent on your local architecture. these files must be placed in the same directory as dada.m, the MATLAB implementation of the DADA algorithm, which will be able to call these routines during sequence alignments. it's recommended to put both of these files in a folder somewhere on your MATLAB path (or in a folder that you subsequently add to the path).
DADA can be called in its simplest form by typing: "dada('inFolder')" where inFolder is a folder of .uniques files as discussed above. It can also be called with a number of options to alter its clustering, with the structure: 
"dada('inFolder','opt1',op1val,'opt2',op2val...)"
These options are:
'omegaA', the abundance p-value significance threshold, set to .01 by default
'omegaR', the read p-vale significance threshold, set to .01 by default
'G', the gap penalty (open and extension)
'GH', a homopolymer gap penalty
'context', whether or not to use context-dependent error probabilities during clustering

there are several other options for expert users described at the top of dada.m.

The output of this clustering will end up in the current folder and will begin with the date and time of calling and then list all options that DADA was called under. within this folder will be the results of clustering after each iteration of updated error probabilities, each residering a folder with the number of the iteration. within each iteration folder, there will be a .mat file for each clustered .uniques file as well as a file called 'ERR.mat'. the non 'ERR' .mat files contain the bin struct array, which contains information about each cluster. it also contains a 'reads' and 'reals' array containing the abundances and genotypes inferred as well as D, a pairwise distance matrix between these types. The ERR.mat file contains two matrices, ERR and ERRa, which are the context-dependent and context-independent error probabilities inferred by DADA. In the language of the DADA paper, ERR(L,i,j,R) = T^{L,R}_{ij} and ERR(i,j) = T_{ij}.

III. DOWNSTREAM of DADA.
If you wish to work with a single .fasta file of all inferred genotypes, we have included bin2fasta. This takes as its input the location of a folder of .mat files produced by DADA and creates a file all.fasta within the same folder, placing the abundance of each sequence in its header. 

We also include 'omegaA_rejoin.m', which takes a folder of .mat files produced by dada and returns a vector of the approximate omegaA values that would be required for each cluster to be reabsorbed into some other cluster. Making a histogram of these is a good way to see whether the results are likely to be sensitive to the choice of the omegaA parameter, and if more conservative values ought to be used.

Finally, we include 'discrimPlot.m', which takes various cluster parameters (see file for details) and significance levels omegaA and omegaR and returns the resulting discrimination plot for the two pvalue statistics. See the dada manuscipt for examples of these plots.