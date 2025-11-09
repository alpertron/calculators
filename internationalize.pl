#!/usr/bin/perl
use strict;
use warnings;

# Check arguments
if (@ARGV != 2) {
    die "Usage: $0 input_file output_file\n";
}

my ($infile, $outfile) = @ARGV;

open(my $in, "<", $infile) or die "Cannot open input file '$infile': $!";
open(my $out, ">", $outfile) or die "Cannot open output file '$outfile': $!";

while (my $line = <$in>) {
    chomp $line;
    # Match lines starting with LITERAL_, followed by identifier and a quoted string
    if ($line =~ /^(LITERAL_[A-Za-z0-9_]+)\s+(".*")/) {
        print $out "#define $1 $2\n";
    }
}

close($in);
close($out);
