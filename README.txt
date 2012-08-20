This is a pipeline for processing amplicon sequence data using the DADA algorithm. 

I. UPSTREAM of DADA. 

DADA takes in data in a format we call .uniques, a succinct way that we store deep sequence data sets. Here is an example:

435	GGCGGCGATATCCTGC
384	GGCTGCGATATCCTGC
286	GGCGGCGATATCATGC
199	GGCGGCGATAGCCTGC
...

Each line consists of an integer and a string separated by a tab. The string is a nucleotide sequence and the integer is the number of times that this sequence was observed in the data set. We provide the following PERL scripts to easily get your data into .uniques format from .fasta.

i. "fasta2seqs.pl" takes a .fasta file, strips away the headers and newlines, and creates a file a .seqs file, which is just one complete sequence on each line. For example:

GGCGGCGATATCCTGC
GGCTGCGATATCCTGC
GGCGGCGATATCCTGC
GGCGGCGATAGCCTGC

"fasta2seqs.pl" can be called with a folder of .fasta files or a single .fasta file. For example:

"perl fasta2seqs.pl ~/myfastas" 

or 

"perl fasta2seqs.pl ./myfasta.fasta"

both work.

A new folder with "_seqs" appended to its name will be created with each .fasta turned into a .seqs

ii. "lengthTrim.pl" takes in a folder of .seqs files, trims sequences down to a given length, and creates a folder of new .seqs files. Sequences shorter than this length are removed. You can call it with

"perl lengthTrim.pl ~/myfastas_seqs 100"

This creates the folder "~/myfastas_seqs100", which will contain .seqs files of 100bp sequences.

iii. "clusterUniques.pl" takes a folder of .seqs files (or single .seqs file) and outputs a folder (or file) in the .unique format described above. This folder will have "_uniques" appended to its name.

II. DADA. 

Before running DADA it is necessary to compile the MEX (Matlab Executable) C programs "align2dom_homo.c" and "align2dom.c". This is done within MATLAB -- just go to the directory where the files are located, and type into the command window:

"mex align2dom.c" 

and 

"mex align2dom_homo.c" 

This creates executable files with a suffix dependent on your local architecture (you can find out what this extension will be by entering "mexext" into the command window). These files should then be placed in the same directory as "dada.m" because DADA calls these routines during sequence alignment. It is recommended that you place these files and "dada.m" in a folder together somewhere on your MATLAB path.

In its simplest form, DADA may be called by

"dada('inFolder')" 

where inFolder is the path for a folder of .uniques files. DADA can also be called with a number of options via

"dada('inFolder','opt1',op1val,'opt2',op2val...)"

The options are:
a. 'omegaA', the abundance p-value significance threshold, set to .01 by default
b. 'omegaR', the read p-value significance threshold, set to .01 by default
c. 'G', the gap penalty (open and extension), set to -4 by default
d. 'GH', a homopolymer gap penalty, set to -1 by default
e. 'context', whether or not to use context-dependent error probabilities during clustering, set to "false" by default. Use values "true" and "false" to turn context-dependence on or off
f. 'tol', a tolerance on the norm of the difference between T matrices on consecutive rounds for DADA to quit (see manuscript), set to 1e-6 by default.
g. 'err', optional initial set of error probabilities for the first round of clustering (instead of estimating these by assuming each file contains only a single genotype) -- may save time if clustering the same or similar data sets many times and can be useful if trying to determine what clustering a particular set of error probabilities would imply.
h. 'noUpdate', whether to not update error probabilities, set to "false" by default. It can be useful to set this to true in combination with 'err' if you want to see the effect of clustering the data with a particular set of error probabilities.

The output is placed in a new folder within the current folder when DADA is called that contains the date and time of the call and all options that DADA was called with. Within this folder will be folders named 0,1,2,... each contains the results of the clustering after each round of updating the error probabilities. The .mat files reflecting each input file contain the struct array "bin", which contains detailed information about the clustering, the "reads" and "reals" arrays containing the inferred abundances and genotypes, and D, a pairwise distance matrix between the inferred types. The folder may also contain (as long as 'noUpdate' was not set to "true") the file "ERR.mat", which contains the context-independent probabilities, ERRa, and context-dependent probilities, ERR. In the language of the methods of the DADA paper, ERR(L,i,j,R) = T^{L,R}_{ij} and ERRa(i,j) = T_{ij}.

III. DOWNSTREAM of DADA.

i. In case you wish to work with .fasta files of the inferred genotypes, we have included the MATLAB script "bin2fasta.m". This takes as its input the location of a folder of .mat files produced by DADA and creates fasta files of the inferred genotypes, placing the abundance of each sequence in its header. If you want all genotypes output to a single file called 'all.fasta', you call it as 

bin2fasta('inFolder')

you can also output a single .uniques file for each file, which you would achieve by calling

bin2fasta('inFolder','fasta',false)

ii. "omegaA_rejoin.m" takes a folder of .mat files produced by dada and returns a vector of approximate omegaA values that would be required for each cluster to be reabsorbed into some other cluster. Making a histogram of these is a good way to see whether the results are likely to be sensitive to the choice of the omegaA parameter, and if more conservative values ought to be used. These p-values are determined for the scenario that the reads of the reads of the error-free family of a given cluster were all errors away from some other cluster (for which all associated reads are included).

iii. "discrimPlot.m" shows a discrimination plot of the type shown in the DADA paper. Its arguments are:
a. omegaA
b. omegaR
c. T, a 4x4 matrix of context-independent error probabilities
d. rho, the number of reads in a cluster
e. bc, a 1x4 vector of the base composition of the sample genotype for a cluster
f. Nfam, the total number of families found for this cluster

For example, one could call

discrimPlot(.01,.01,ERRa,5000,[50 50 50 50],1000)

where ERRa came from the "ERR.mat" file of the final round of clustering for your data set. Because of the way the "effective Hamming distance" is defined (again see DADA paper), everything above the curved line or to the right of the vertical line would be rejected as an error by DADA at the given significance levels if the x-axis is interpreted as a true Hamming distance, so this plot will show you the minimum resolution you are getting by DADA. If you wanted finer resolution, you can run the algorithm again with one or the other significance level changed, but you may be getting more false positives in the process. 

