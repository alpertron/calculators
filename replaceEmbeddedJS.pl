#! /usr/bin/perl
use strict; 
use warnings;
use MIME::Base64;

my $commandLine = $ARGV[0];
my $htmlFile = $ARGV[1];
my $jsFile = $ARGV[2];
my $workerContents = "";
my $wasmContents = "";
my $wasmb64 = "";
# Argument 3 or 4 can be optionally, WebAssembly file (file name ends with .wasm)
# and worker file (file name ends with .js).
for (my $cmdLineArgNbr = 3; $cmdLineArgNbr <= $#ARGV; $cmdLineArgNbr++)
{
  if ($ARGV[$cmdLineArgNbr] =~ /wasm$/)
  {          # WebAssembly file
    open my $wasmFile, '<', $ARGV[$cmdLineArgNbr] or die "couldn't open wasm file ".$ARGV[$cmdLineArgNbr];
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
  else
  {             # Worker file
    open my $workerFile, '<', $ARGV[$cmdLineArgNbr] or die "couldn't open worker file ".$ARGV[$cmdLineArgNbr];
    binmode $workerFile;
    while (1)
    {
      my $success = read $workerFile, $workerContents, 512, length($workerContents);
      die $! if not defined $success;
      last if not $success;
    }
    close $workerFile;
  }
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
  elsif (/id=\"worker\"/)
  {
    print tempFile;
    print tempFile $workerContents;
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
  elsif (/noand/)
  {
    $step = 5;
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
  elsif ($step == 5)
  {
    if (/\<\/div/)
    {
      $step = 1;
    }
    elsif (/noand/)
    {
    }
    else
    {
      print tempFile;
    }
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
