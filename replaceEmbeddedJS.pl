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
  open my $wasmFile, '<', $ARGV[3] or die "couldn't open wasm file ".$ARGV[3];
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
open(htmlFile, '<', $htmlFile) or die "couldn't open HTML file ".$htmlFile;
open(jsFile, '<', $jsFile) or die "couldn't open JS file ".$jsFile;
open(tempFile, '>', $tempFile) or die "couldn't open temp file ".$tempFile;

while (<htmlFile>)
{
  if (/id=\"wasmb64\"/)
  {
    print tempFile;
    print tempFile $wasmb64;
    $step = 4;
  }
  elsif (/^\<script\>/)
  {
    $step = 2;
  }
  elsif (/^\<\/script\>/)
  {
    $step = 1;
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
