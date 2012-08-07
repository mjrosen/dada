#!/usr/bin/perl -w

#args:
#1. the name of a folder with files of amplicons in .seqs format
#2. a length to trim all reads (shorter reads are excluded)

($dir, $length) = @ARGV;
my @files = <$dir/*>;

system "rm -r -f $dir$length";
system "mkdir $dir$length";

foreach my $file (@files){
    open IN, $file;
    $file =~ s/$dir/$dir$length/g;
    open OUT, ">", $file or die "Could not open $file: $!";
    while(<IN>){
	chomp;
	print OUT substr($_,0,$length)."\n" if length >= $length;
    }
    close IN;
    close OUT;
}
