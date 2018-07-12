#! /usr/bin/perl
use strict; 
use warnings;
my $commandLine = $ARGV[0];
my $htmlFile = $ARGV[1];
my $jsFile = $ARGV[2];
my $oldJS = "0000\.js";
my $oldWASM = "0000\.wasm";
my $newJS = $commandLine."\.js";
my $newWASM = $commandLine."\.wasm";
my $tempFile = "temp.tmp";
my $step = 1;
open(htmlFile, '<', $htmlFile) or die "couldn't open HTML file";
open(jsFile, '<', $jsFile) or die "couldn't open JS file";
open(tempFile, '>', $tempFile) or die "couldn't open temp file";

while (<htmlFile>)
{
  if ($step == 2)
  {
    if (/^\/\/-->/)
    {
      print tempFile;
      $step = 3;
    }
  }
  else
  {
    print tempFile;
    if (/^\<\!--/)
    {
      $step = 2;
      while (<jsFile>)
      {
        s/$oldJS/$newJS/g;
        s/$oldWASM/$newWASM/g;
        print tempFile;
      }
    }
  }
}

close (htmlFile); 
close (jsFile);
close (tempFile);
unlink $htmlFile;
unlink $jsFile;
rename $tempFile, $htmlFile;
