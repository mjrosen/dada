#!/usr/bin/perl -w

#args:
#either a single fasta or a directory of fastas
#if it is a folder, '_seqs' is appended to the folder name
#and .seqs if it is a file

my $in = shift @ARGV;

if (-d $in) { #user has input a directory
    #strip off final slash if present
    $in =~ s/\/$//;
    $out = $in.'_seqs'; #the output folder name
    system "rm -r -f $out";
    system "mkdir $out";
    my @files = <$in/*>;
    foreach $fasta (@files){
	($seqs = $fasta) =~ s/$in/$out/; #change directory
	$seqs =~ s/[^.]+$/seqs/ if $seqs =~ /[.]/;	
	&convert($fasta, $seqs);
    }
} elsif (-e $in) { #user has input a file
    ($fasta, $seqs) = ($in, $in);
    $seqs =~ s/[^.]+$/seqs/ if $seqs =~ /[.]/;
    &convert($fasta, $seqs);    
} else {
    die 'Input is neither a valid file or directory';
}

sub convert {    
    open IN, "<", shift @_;
    open OUT, ">", shift @_;

    my $read = '';
    while(<IN>) {
	chomp;
	if (/>/) { #read a new header
	    print OUT $read."\n" if $read;
	    $read = '';
	} else {#part of a read
	    $read .= $_;
	}
    }
    print OUT $read."\n" if $read; #dump the final line
}
