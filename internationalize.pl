#!/usr/bin/perl
#
# This file is part of Alpertron Calculators.
#
# Copyright 2015-2021 Dario Alejandro Alpern
#
# Alpertron Calculators is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Alpertron Calculators is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
#
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
