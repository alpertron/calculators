#! /usr/bin/perl
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
my $fileNbr = $ARGV[0];
my $templateFile = $ARGV[1];
my $serviceWorkerFile = $ARGV[2];
my $oldNbr = "0000";
my $step = 1;
unlink $serviceWorkerFile;
open(templateFile, '<', $templateFile) or die "couldn't open HTML file";
open(serviceWorkerFile, '>', $serviceWorkerFile) or die "couldn't open JS file";

while (<templateFile>)
{
  s/$oldNbr/$fileNbr/g;
  print serviceWorkerFile;
}

close (templateFile); 
close (serviceWorkerFile);
