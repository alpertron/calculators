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
use bytes;

my $jsContents = "";
my $asmJSContents = "";
my $idxJS;
my $idxJSasm;
my $idxJSasmEnd;
my $idxStartAsmjsLoading;

open my $jsFile, '<', $ARGV[0] or die "couldn't open js file ".$ARGV[0]." for reading";
while (1)
{
  my $success = read $jsFile, $jsContents, 512, length($jsContents);
  die $! if not defined $success;
  last if not $success;
}
close $jsFile;

open my $asmJSFile, '<', $ARGV[1] or die "couldn't open asmJS file ".$ARGV[1]." for reading";
while (1)
{
  my $success = read $asmJSFile, $asmJSContents, 512, length($asmJSContents);
  die $! if not defined $success;
  last if not $success;
}
close $asmJSFile;

open(tempFile, '>', $ARGV[2]) or die "couldn't open temp file ".$ARGV[2]." for writing";

$idxStartAsmjsLoading = index($jsContents, "startLowLevelCodeCallback");
# Change fields of asm object according to asmJS file.
for (my $i=3; $i<=$#ARGV; $i++)
{
#For each function defined in asmJS:
  my $funcName = $ARGV[$i];
#Find function name (case insensitive).
  $idxJSasm = index(lc($asmJSContents), lc($funcName));
  if ($idxJSasm == -1)
  {
    die "couldn't find export function ".$funcName. " in asmJS file.";
  }
#Find next dot after the function name.
  $idxJSasm = index($asmJSContents, ".", $idxJSasm);
  if ($idxJSasm == -1)
  {
    die "couldn't find export function ".$funcName. " in asmJS file.";
  }
#Find function name (case insensitive).
  $idxJS = index(lc($jsContents), lc($funcName), $idxStartAsmjsLoading);
  if ($idxJS == -1)
  {
    die "couldn't find export function ".$funcName. " in js file.";
  }
#Find next dot after the function name.
  $idxJS = index($jsContents, ".", $idxJS);
  if ($idxJS == -1)
  {
    die "couldn't find export function ".$funcName. " in js file.";
  }
#Replace field name from the one found in the asmJS file.
  substr($jsContents, $idxJS+1, 1) = substr($asmJSContents, $idxJSasm+1, 1);
}

close($jsFile);
close($asmJSFile);

#Find and write first part of JS file.
$idxJS = index($jsContents, "//##");
print tempFile substr($jsContents, 0, $idxJS);

#Find and write asmJS from asmJS file.
$idxJSasm = index($asmJSContents, "// EMSCRIPTEN_START_ASM");
$idxJSasmEnd = index($asmJSContents, "// EMSCRIPTEN_END_ASM");
print tempFile substr($asmJSContents, $idxJSasm, $idxJSasmEnd - $idxJSasm);

#Write remainder of JS file.
print tempFile substr($jsContents, $idxJS);

close(tempFile);
unlink($ARGV[1]);
