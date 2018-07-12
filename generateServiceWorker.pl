#! /usr/bin/perl
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
