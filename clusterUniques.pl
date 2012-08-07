#!/usr/bin/perl -w

#read in a directory of files of amplicons, either length sorted or not, in .seqs format and output the unique sequences in .seqs format (that is, each line has "read num \t sequence \n").

#args:
#1. either a file of seqs or a directory of files of seqs. if it is a folder, '_uniques' is appended to the folder name, and likelise if a file

my $in = shift @ARGV;
$in =~ s/\/$//; #remove forward slash if one is present
my $out = $in.'_uniques';
if (-d $in) { #user has input a directory
    system "rm -r -f $out";
    system "mkdir $out"; #refresh the directory
    my @files = <$in/*>; #the input files
    foreach $seqs (@files){
	($uniques = $seqs) =~ s/$in/$out/; #change the directory
	$uniques =~ s/[^.]+$/uniques/ if $uniques =~ /[.]/;
	&sort($seqs,$uniques);
    }
} elsif (-e $in) { #user has input a file
    ($seqs, $uniques) = ($in, $in);
    $uniques =~ s/[^.]+$/uniques/ if $uniques =~ /[.]/;
    &sort($seqs,$uniques);
} else {
    die 'Input is neither a valid file or directory';
}

sub sort {
    open SEQS, "<", shift @_;    
    open UNIQUES, ">", shift @_;       
    my %seqHash = ();
    my @seqs = ();
    
    while(<SEQS>){
	chomp;
	push @seqs, $_;
    }
    
    foreach (@seqs){
	if (exists $seqHash{$_}) {
	    $seqHash{$_}++;
	} else {
	    $seqHash{$_} = 1;
	}
    }
    
    @sorted = sort {$seqHash{$b} <=> $seqHash{$a}} keys %seqHash;
    
    foreach (@sorted) {
	print UNIQUES $seqHash{$_}."\t".$_."\n";
    }   
}
