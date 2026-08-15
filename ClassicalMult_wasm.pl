#
# This file is part of Alpertron Calculators.
#
# Copyright 2025 Dario Alejandro Alpern
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
#!/bin/perl
my $maxNbrLimbs       = 16;
my $FORMAT_STR        = "  %-25s# %s\n";
my $funcName          = "ClassicalMult";

my $indexMultiplicand = 0;
my $indexMultiplier   = 1;
my $numberLimbs       = 2;
my $counter           = 3;
my $Right0            = 4;
my $Right1            = 5;
my $Right2            = 6;
my $Right3            = 7;
my $ptrProduct        = 8;
my $productAccumul    = 9;
my $Right128Bits      = 10;
my $LeftWindow        = 11;
my $SumProd           = 12;
my $Left              = 13;
my $Product           = $Left + $maxNbrLimbs;

print <<"EOF";
  .extern arr
	.globl  $funcName
	.type	$funcName,\@function
$funcName:
  .functype $funcName (i32, i32, i32) -> ()
EOF
print "  .local i32, i64, i64, i64, i64, i32, i64, v128, v128, v128";

# Reserve space for local variables for multiplier.
for (my $idx=0; $idx < $maxNbrLimbs; $idx++) {
  print ", i64";
}
# Reserve space for local variables for product.
for (my $idx=0; $idx < $maxNbrLimbs; $idx++) {
  print ", i64";
}
print "\n";

#
# Multiply two groups of nbrLen limbs. The first one starts at idxFactor1
# and the second one at idxFactor2. The 2*nbrLen limb result is stored
# starting at idxFactor1.
# Accumulate products by result limb.
#
# On input:
# indexMultiplicand: Index to arr of multiplicand.
# indexMultiplier: Index to arr of multiplier.
# numberLimbs: Number of limbs.
#
#                                 A  B | C  D | E  F
#                                 H  I | J  K | L  M
#                                 ======+======+=====
#                                AM BM |CM DM |EM FM
#                            AL |BL CL |DL EL |FL
#                         ------+------+-----
#                         AK BK |CK DK |EK FK
#                     AJ |BJ CJ |DJ EJ |FJ
#                  ------+------+------
#                  AI BI |CI DI |EI FI
#              AH |BH CH |DH EH |FH
#              
#  Variables in the accumulators are named according to the group.
#  In this way, LeftH can be A, C or E, LeftL can be B, D or F,
#  Right1 can be H, J or L and Right0 can be I, K or M.
#
#  Use local variables Left0, Left1, Left2, ... for complete multiplicand,
#  Right1, Right0 for two limbs of multiplier and
#  Prod0, Prod1, Prod2, ... for running limbs of product. Lower-level
#  limbs are stored in memory.
#
#  Before the outer loop, the first two lines are computed in groups of two limbs.
#  The outer loop processes two lines at a time, in groups of two limbs.
#

sub store_mem_and_get_carry {
  my ($offset) = @_;
  printf($FORMAT_STR, "local.tee $productAccumul", "Store accumulator in local variable");
  printf($FORMAT_STR, "local.get $ptrProduct", "Get address of limb of product");
  printf($FORMAT_STR, "local.get $productAccumul", "Get accumulator again");
  printf($FORMAT_STR, "i64.const 2147483647", "");
  printf($FORMAT_STR, "i64.and", "");
  my $limbNbr = $offset/4;
  printf($FORMAT_STR, "i64.store32 $offset", "Store limb #$limbNbr of product to memory");
  printf($FORMAT_STR, "i64.const 31", "");
  printf($FORMAT_STR, "i64.shr_u", "Get carry from multiplication");
}

sub store_localvar_and_get_carry {
  my ($localVarIndex) = @_;
  printf($FORMAT_STR, "local.tee $productAccumul", "Store accumulator in local variable");
  printf($FORMAT_STR, "i64.const 2147483647", "");
  printf($FORMAT_STR, "i64.and", "");
  printf($FORMAT_STR, "local.set ${\ ($Product + $localVarIndex)}", "Store Prod$localVarIndex");
  printf($FORMAT_STR, "local.get $productAccumul", "Get accumulator");
  printf($FORMAT_STR, "i64.const 31", "");
  printf($FORMAT_STR, "i64.shr_u", "Get carry from multiplication");
}

sub add_product {
  my ($leftVarIndex, $rightVarIndex) = @_;
  printf($FORMAT_STR, "local.get ${\ ($Left + $leftVarIndex)}", "Get Left$leftVarIndex from local variable");
  printf($FORMAT_STR, "local.get ${\ ($Right0 + $rightVarIndex)}", "Get Right$rightVarIndex from local variable");
  printf($FORMAT_STR, "i64.mul", "Get Left$leftVarIndex * Right$rightVarIndex");
  printf($FORMAT_STR, "i64.add", "Add product");  
}

# nbrRows must be between 1 and 4.
sub multiply_multiple_rows {
  my ($nbrRows, $nbrLimbs, $initialRow) = @_;
  my $row, $col;
  # Get RightN.
  for ($row = $nbrRows - 1; $row >= 0; $row--) {
    printf($FORMAT_STR, "local.get $ptrProduct", "");
    my $currentOffset = ($nbrLimbs + $initialRow + $row) * 4;
    printf($FORMAT_STR, "i64.load32_u $currentOffset", "Get Right$row from memory");
    my $localVar = $Right0 + $row;
    if ($row == 0) {
      printf($FORMAT_STR, "local.tee $localVar", "Save Right$row to local variable");
    } else {
      printf($FORMAT_STR, "local.set $localVar", "Save Right$row to local variable");
    }
  }
  if ($nbrRows == 4) {
    printf($FORMAT_STR, "local.get $ptrProduct", "");
    my $currentOffset = ($nbrLimbs + $initialRow) * 4;
    printf($FORMAT_STR, "v128.load $currentOffset", "Get Right0-3 from memory");
    printf($FORMAT_STR, "local.set $Right128Bits", "Store Right0-3 in local variable");
  }
  # Add product of rightmost (incomplete) columns.
  printf($FORMAT_STR, "local.get $Left", "Get value of Left0");
  printf($FORMAT_STR, "i64.mul", "Compute product Left0 * Right0");
  if ($initialRow > 0) {
    printf($FORMAT_STR, "local.get $Product", "Get value of Prod0");
    printf($FORMAT_STR, "i64.add", "Add limb from product");
  }
  store_mem_and_get_carry($initialRow * 4);
  for ($col = 1; $col < $nbrRows - 1; $col++) {
    for ($row=0; $row <= $col; $row++) {
      add_product($col - $row, $row);
    }
    if ($initialRow > 0) {
      my $localVar = $Product + $col;
      printf($FORMAT_STR, "local.get $localVar", "Get value of Prod$col");
      printf($FORMAT_STR, "i64.add", "Add limb from product");
    }
    store_mem_and_get_carry(($initialRow + $col) * 4);
  }
  # Add products of complete columns.
  for (; $col < $nbrLimbs; $col++) {
    if ($nbrRows == 4) {
      if ($col == $nbrRows - 1) {
        printf($FORMAT_STR, "local.get ${\ ($Left + $col)}", "Get Left$col");
        printf($FORMAT_STR, "i32.wrap_i64", "Obtain least significant limb");
        printf($FORMAT_STR, "i32x4.splat",
               "Get Left$col: Left$col: Left$col: Left$col");
        printf($FORMAT_STR, "local.get ${\ ($Left + $col - 1)}", "Get Left${\ ($col - 1)}");
        printf($FORMAT_STR, "i32.wrap_i64", "Obtain least significant limb");
        printf($FORMAT_STR, "i32x4.replace_lane 1",
               "Get Left$col: Left${\ ($col - 1)}: Left$col: Left$col");
        printf($FORMAT_STR, "local.get ${\ ($Left + $col - 2)}", "Get Left${\ ($col - 2)}");
        printf($FORMAT_STR, "i32.wrap_i64", "Obtain least significant limb");
        printf($FORMAT_STR, "i32x4.replace_lane 2",
               "Get Left$col: Left${\ ($col - 1)}: Left${\ ($col - 2)}: Left$col");
        printf($FORMAT_STR, "local.get ${\ ($Left + $col - 3)}", "Get Left${\ ($col - 3)}");
        printf($FORMAT_STR, "i32.wrap_i64", "Obtain least significant limb");
        printf($FORMAT_STR, "i32x4.replace_lane 3",
               "Get Left$col: Left${\ ($col - 1)}: Left${\ ($col - 2)}: Left${\ ($col - 3)}");
      } else {
        printf($FORMAT_STR, "local.get $LeftWindow",
               "Get Left${\ ($col - 1)}: Left${\ ($col - 2)}: Left${\ ($col - 3)}: Left${\ ($col - 4)}");
        printf($FORMAT_STR, "local.get $LeftWindow", "shuffle requires the argument twice in stack");
        printf($FORMAT_STR, "i8x16.shuffle 12, 13, 14, 15,  0, 1, 2, 3,  4, 5, 6, 7,  8, 9, 10, 11", "");
        printf($FORMAT_STR, "",
               "Get Left${\ ($col - 4)}: Left${\ ($col - 1)}: Left${\ ($col - 2)}: Left${\ ($col - 3)}");
        printf($FORMAT_STR, "local.get ${\ ($Left + $col)}", "Get Left$col");
        printf($FORMAT_STR, "i32.wrap_i64", "Obtain least significant limb");
        printf($FORMAT_STR, "i32x4.replace_lane 0",
               "Get Left$col: Left${\ ($col - 1)}: Left${\ ($col - 2)}: Left${\ ($col - 3)}");
      }
      printf($FORMAT_STR, "local.tee $LeftWindow",
             "Store Left$col: Left${\ ($col - 1)}: Left${\ ($col - 2)}: Left${\ ($col - 3)}");
      printf($FORMAT_STR, "local.get $Right128Bits", "Get Right0: Right1: Right2: Right3");
      printf($FORMAT_STR, "i64x2.extmul_low_i32x4_u", "Get Left$col*Right0: Left${\ ($col - 1)}*Right1");
      printf($FORMAT_STR, "local.get $LeftWindow",
             "Get Left$col: Left${\ ($col - 1)}: Left${\ ($col - 2)}: Left${\ ($col - 3)}");
      printf($FORMAT_STR, "local.get $Right128Bits", "Get Right0: Right1: Right2: Right3");
      printf($FORMAT_STR, "i64x2.extmul_high_i32x4_u", "Get Left${\ ($col - 2)}*Right2: Left${\ ($col - 3)}*Right3");
      printf($FORMAT_STR, "i64x2.add", "Get Left$col*Right0 + Left${\ ($col - 2)}*Right2: Left${\ ($col - 1)}*Right1 + Left${\ ($col - 3)}*Right3");
      printf($FORMAT_STR, "local.tee $SumProd", "Store this vector");
      printf($FORMAT_STR, "i64x2.extract_lane 0", "Get Left$col*Right0 + Left${\ ($col - 2)}*Right2");
      printf($FORMAT_STR, "local.get $SumProd", "Get the vector again");
      printf($FORMAT_STR, "i64x2.extract_lane 1", "Get Left${\ ($col - 1)}*Right1 + Left${\ ($col - 3)}*Right3");
      printf($FORMAT_STR, "i64.add", "Get Left$col*Right0 + Left${\ ($col - 2)}*Right2 + Left${\ ($col - 1)}*Right1 + Left${\ ($col - 3)}*Right3");
      printf($FORMAT_STR, "i64.add", "Add carry");
    } else {
      for ($row=0; $row < $nbrRows; $row++) {
        add_product($col - $row, $row);
      }
    }
    my $localVar = $Product + $col;
    printf($FORMAT_STR, "local.get $localVar", "Get value of Prod$col");
    printf($FORMAT_STR, "i64.add", "Add limb from product");
    if ($col == $nbrRows - 1) {
      store_mem_and_get_carry(($initialRow + $col) * 4);
    } else {
      store_localvar_and_get_carry($col - $nbrRows);
    }
  }
  # Add products of leftmost (incomplete) columnsz.
  for ($col=$nbrLimbs; $col < $nbrRows + $nbrLimbs - 1; $col++) {
    for ($row=0; $row < $nbrRows; $row++) {
      if ($col - $row < $nbrLimbs) {  # Do not add out of range products.
        add_product($col - $row, $row);
      }
    }
    if (($initialRow > 0) && ($col < $nbrLimbs)) {
      my $localVar = $Product + $col;
      printf($FORMAT_STR, "local.get $localVar", "Get value of Prod$col");
      printf($FORMAT_STR, "i64.add", "Add limb from product");
    }
    store_localvar_and_get_carry($col - $nbrRows);    
  }
  my $localVarIndex = $nbrLimbs - 1;
  $localVar = $Product + $localVarIndex;
  printf($FORMAT_STR, "local.set $localVar", "Store Prod$localVarIndex in local variable");  
}

#
# Main routine
#
my $loopNbr;
for ($loopNbr=0; $loopNbr <= $maxNbrLimbs; $loopNbr++) {
  print "  block\n";
}
printf($FORMAT_STR, "i32.const arr", ""); 
printf($FORMAT_STR, "local.get $indexMultiplicand", "Get index of arguments and products");
printf($FORMAT_STR, "i32.const 2", "Convert index to offset to array arr"); 
printf($FORMAT_STR, "i32.shl", ""); 
printf($FORMAT_STR, "i32.add", "Get pointer to arguments and product"); 
printf($FORMAT_STR, "local.set $ptrProduct", "Save it in ptrProduct"); 
printf($FORMAT_STR, "i32.const 0", "Initialize counter");
printf($FORMAT_STR, "local.set $counter", "");
printf($FORMAT_STR, "local.get $numberLimbs", "Perform long switch for each value of nbrLimbs");
my $maxNbrLimbsMinus1 = $maxNbrLimbs - 1;
printf($FORMAT_STR, "i32.const 1", "Convert range [1, $maxNbrLimbs] to [0, $maxNbrLimbsMinus1]");
printf($FORMAT_STR, "i32.sub", "");
print "  br_table {";
for ($loopNbr=0; $loopNbr <= $maxNbrLimbs; $loopNbr++) {
  print "," if ($loopNbr ne 0);
  print " $loopNbr";
}
print "}\n";
for ($nbrLimbs=1; $nbrLimbs <= $maxNbrLimbs; $nbrLimbs++) {
  print "  end_block\n";
  if ($nbrLimbs > 1) {
    print "     /* Write $nbrLimbs-limb multiplicand to local variables */\n";
    for ($loopNbr=0; $loopNbr < $nbrLimbs; $loopNbr++) {
      printf($FORMAT_STR, "local.get $ptrProduct", "Get pointer to multiplicand");
      my $offset = 4 * $loopNbr;
      printf($FORMAT_STR, "i64.load32_u $offset", "Get Left$loopNbr from memory");
      $offset = $Left + $loopNbr;
      printf($FORMAT_STR, "local.set $offset", "Save Left$loopNbr to local variable");
    }
    print "     /* Classical multiplication for $nbrLimbs limbs */\n";
  }
  if ($nbrLimbs == 1) {
    print "     /* Classical multiplication for 1 limb */\n";
    printf($FORMAT_STR, "local.get $ptrProduct", "Get pointer to product");
    printf($FORMAT_STR, "local.get $ptrProduct", "Get pointer to multiplicand");
    printf($FORMAT_STR, "i64.load32_u 0", "Get multiplicand from memory");
    printf($FORMAT_STR, "local.get $ptrProduct", "Get pointer to multiplicand");
    printf($FORMAT_STR, "i64.load32_u 4", "Get multiplier from memory");
    printf($FORMAT_STR, "i64.mul", "Compute product");
    printf($FORMAT_STR, "local.tee $productAccumul", "Store product");
    printf($FORMAT_STR, "i32.wrap_i64", "Obtain least significant limb");
    printf($FORMAT_STR, "i32.const 2147483647", "");
    printf($FORMAT_STR, "i32.and", "");
    printf($FORMAT_STR, "i32.store 0", "Store limb #0 of product to memory");
    printf($FORMAT_STR, "local.get $ptrProduct", "Get pointer to product");
    printf($FORMAT_STR, "local.get $productAccumul", "Get accumulator from local variable");
    printf($FORMAT_STR, "i64.const 31", "");
    printf($FORMAT_STR, "i64.shr_u", "Get carry from multiplication");
    printf($FORMAT_STR, "i32.wrap_i64", "Obtain least significant limb");
    printf($FORMAT_STR, "i32.store 4", "Store limb #1 of product to memory");
  } elsif ($nbrLimbs == 2) {
    multiply_multiple_rows(2, 2, 0);
  } elsif ($nbrLimbs == 3) {
    multiply_multiple_rows(3, 3, 0);
  } else {
    multiply_multiple_rows(4, $nbrLimbs, 0);
    for ($row = 4; $row <= $nbrLimbs - 4; $row += 4) {
      multiply_multiple_rows(4, $nbrLimbs, $row);
    }
    if ($row < $nbrLimbs) {
      multiply_multiple_rows($nbrLimbs - $row, $nbrLimbs, $row);
    }
  }
  if ($nbrLimbs > 1) {
    my $doubleLimbs = 2 * $nbrLimbs;
    print "     /* Write upper half of $doubleLimbs-limb product to memory from local variables */\n";
    for ($loopNbr=0; $loopNbr < $nbrLimbs; $loopNbr++) {
      $offset = $Product + $loopNbr;
      printf($FORMAT_STR, "local.get $ptrProduct", "Get pointer to product");
      printf($FORMAT_STR, "local.get $offset", "Get Prod$loopNbr from local variable");
      my $currentProdLimb = $loopNbr + $nbrLimbs;
      my $offset = 4 * $currentProdLimb;
      printf($FORMAT_STR, "i64.store32 $offset", "Store it to limb #$currentProdLimb of product on memory");
    }
  }
  my $jumpTarget = $maxNbrLimbs - $nbrLimbs;
  printf($FORMAT_STR, "br $jumpTarget", "Exit all blocks");
}

printf($FORMAT_STR, "end_block", "");
printf($FORMAT_STR, "end_function", "");
