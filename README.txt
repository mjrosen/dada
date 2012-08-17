This is a pipeline for processing amplicon sequence data using the DADA algorithm. 

I. UPSTREAM of DADA. 

DADA takes in data in a format we call .uniques. This a succinct way to store deep sequence data sets. For example:

435	GGCGGCGATATCCTGC
384	GGCTGCGATATCCTGC
286	GGCGGCGATATCATGC
199	GGCGGCGATAGCCTGC
...

Each line consists of a number of identical copies and a sequence, which are separated by a tab. We provide the following PERL scripts to easily get your data into .uniques format from .fasta.

	i. "fasta2seqs.pl" This takes a FASTA file and strips away all header information and newlines and creates a file that we call a '.seqs' file, with one read per line. For example:

GGCGGCGATATCCTGC
GGCTGCGATATCCTGC
GGCGGCGATATCCTGC
GGCGGCGATAGCCTGC

fasta2seqs.pl can be called along with a folder of .fasta files or a single .fasta file. For example:

"perl fasta2seqs.pl ~/myfastas" 

or 

"perl fasta2seqs.pl ./myfasta.fasta"

are both acceptable.

A new file or folder will be created with "_seqs" appended to its name

	ii. "lengthTrim.pl" This takes in a folder of .seqs files, trims all sequences down to a particular length, and creates a folder of new .seqs files. Sequences shorter than the given length are removed. It can be called as:

"perl lengthTrim.pl ~/myfastas_seqs 100"

This creates the folder "~/myfastas_seqs100", which will contain .seqs files of 100bp sequences.

	iii. "clusterUniques.pl" This takes in a folder of .seqs files (or single .seqs file) and outputs a folder (or file) in the .unique format described above. The folder will have "_uniques" appended to its name.

II. DADA. 

Before running DADA it is necessary to compile the MEX (Matlab Executable) C programs "align2dom_homo.c" and "align2dom.c". This is done within MATLAB -- just open MATLAB, go to the directory where the files are located, and type into the command window

"mex align2dom.c" 

and 

"mex align2dom_homo.c" 

This will create executable files with a suffix dependent on your local architecture. These files should then be placed in the same directory as "dada.m" because DADA will call these routines during sequence alignments. It's recommended that you place these files and "dada.m" somewhere on your MATLAB path.

DADA can be called most simply with 

"dada('inFolder')" 

where 'inFolder' is a folder of .uniques files (explained avove) It can also be called with a number of options via

"dada('inFolder','opt1',op1val,'opt2',op2val...)"

The main options are:
a. 'omegaA', the abundance p-value significance threshold, set to .01 by default
b. 'omegaR', the read p-vale significance threshold, set to .01 by default
c. 'G', the gap penalty (open and extension)
d. 'GH', a homopolymer gap penalty
e. 'context', whether or not to use context-dependent error probabilities during clustering, set to "false" by default. Use values "true" and "false" to turn context-dependence on or off
f. 'tol', a tolerance on the norm of the difference between T matrices on consecutive rounds for DADA to quit (see manuscript)
g. 'err', gives the user the option of providing an initial set of error probabilities for the first round of clustering. This may save time when clustering the same or similar data sets many times.
h. 'noUpdate', whether to not update error probabilities, set to "false" by default. Can be useful to set to true in combination with 'err' if you want to see the effect of clustering the data with a particular set of error probabilities.

The output folder is placed in the current folder when DADA is called, and contains the date and time of the call and all options that DADA was called under. 

Within this folder will be the results of clustering after each iteration of updated error probabilities, each residering a folder with the number of the iteration. within each iteration folder, there will be a .mat file for each clustered .uniques file as well as a file called 'ERR.mat'. the non 'ERR' .mat files contain the bin struct array, which contains information about each cluster. it also contains a 'reads' and 'reals' array containing the abundances and genotypes inferred as well as D, a pairwise distance matrix between these types. The ERR.mat file contains two matrices, ERR and ERRa, which are the context-dependent and context-independent error probabilities inferred by DADA. In the language of the DADA paper, ERR(L,i,j,R) = T^{L,R}_{ij} and ERR(i,j) = T_{ij}.

III. DOWNSTREAM of DADA.
If you wish to work with a single .fasta file of all inferred genotypes, we have included bin2fasta. This takes as its input the location of a folder of .mat files produced by DADA and creates a file all.fasta within the same folder, placing the abundance of each sequence in its header. 

We also include 'omegaA_rejoin.m', which takes a folder of .mat files produced by dada and returns a vector of the approximate omegaA values that would be required for each cluster to be reabsorbed into some other cluster. Making a histogram of these is a good way to see whether the results are likely to be sensitive to the choice of the omegaA parameter, and if more conservative values ought to be used.

Finally, we include 'discrimPlot.m', which takes various cluster parameters (see file for details) and significance levels omegaA and omegaR and returns the resulting discrimination plot for the two pvalue statistics. See the dada manuscipt for examples of these plots.