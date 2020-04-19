#! /usr/bin/perl
use strict; 
use warnings;
use MIME::Base64;

my $commandLine = $ARGV[0];
my $htmlFile = $ARGV[1];
my $jsFile = $ARGV[2];
my $wasmContents = "";
my $wasmb64 = "";
if ($#ARGV == 3)
{
  open my $wasmFile, '<', $ARGV[3] or die;
  binmode $wasmFile;
  while (1)
  {
    my $success = read $wasmFile, $wasmContents, 512, length($wasmContents);
    die $! if not defined $success;
    last if not $success;
  }
  close $wasmFile;
  $wasmb64 = encode_base64($wasmContents);
  $wasmb64 =~ s/\//!/g;
  $wasmb64 =~ s/[\n\r]//g;
  $wasmb64 .= "\n";
}
my $oldJS = "0000\.js";
my $newJS = $commandLine."\.js";
my $tempFile = "temp.tmp";
my $step = 1;
open(htmlFile, '<', $htmlFile) or die "couldn't open HTML file";
open(jsFile, '<', $jsFile) or die "couldn't open JS file";
open(tempFile, '>', $tempFile) or die "couldn't open temp file";

while (<htmlFile>)
{
  if (/^\<script type=\"text\/wasmb64\"/)
  {
    print tempFile;
    $step = 4;
    next;
  }
  elsif (/^\<script\>/)
  {
    $step = 2;
  }
  else 
  {
    if (/^\<\/script\>/)
    {
      $step = 1;
    }
  }
  if ($step == 3)
  {
    if (/^\/\/-->/)
    {
      print tempFile;
      $step = 1;
    }
  }
  elsif ($step == 4)
  {
    print tempFile $wasmb64;
    $step = 1;
  }
  else
  {
    print tempFile;
    if ($step == 2)
    {
      if (/^\<\!--/)
      {
        $step = 3;
        while (<jsFile>)
        {
          s/$oldJS/$newJS/g;
          print tempFile;
        }
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
