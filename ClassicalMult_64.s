/*
 This file is part of Alpertron Calculators.

 Copyright 2025 Dario Alejandro Alpern

 Alpertron Calculators is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Alpertron Calculators is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
*/

    .arch armv8-a
    .global ClassicalMult
    .type ClassicalMult, %function
    .text
    .align 4
/*
 Multiply two groups of nbrLen limbs. The first one starts at idxFactor1
 and the second one at idxFactor2. The 2*nbrLen limb result is stored
 starting at idxFactor1. Use arrayAux as temporary storage.
 Accumulate products by result limb.

 On input:
 w0: Index to arr of multiplicand.
 w1: Index to arr of multiplier.
 w2: Number of limbs.
 Other registers:
 x0: Pointer to multiplicand, then multiplier. Also pointer to product.
 x7: x1 = carry between columns.
 x4: x3 = Sum of products of double limbs in a column during long multiplication.
 x6: x5 = Product of double limbs in a column during long multiplication.
 x8, x9, ... = Double limbs of multiplicand and multiplier.
 */

    A0 .req x8
    A1 .req x10
    A2 .req x12
    A3 .req x14
    A4 .req x16
    A5 .req x19
    A6 .req x21
    A7 .req x23
    B0 .req x9
    B1 .req x11
    B2 .req x13
    B3 .req x15
    B4 .req x17
    B5 .req x20
    B6 .req x22
    B7 .req x24

    A0L .req w8
    A1L .req w10
    A2L .req w12
    A3L .req w14
    A4L .req w16
    A5L .req w19
    A6L .req w21
    A7L .req w23
    B0L .req w9
    B1L .req w11
    B2L .req w13
    B3L .req w15
    B4L .req w17
    B5L .req w20
    B6L .req w22
    B7L .req w24

    /* x3:x4 accumulator, x5:x6 temporary regs,
       x1:x7 carry (low:high), x0 output */

    .macro FIRST a, b
    mul   x3, \a, \b
    umulh x4, \a, \b
    .endm

    .macro MAC a, b
    mul   x5, \a, \b
    umulh x6, \a, \b
    adds  x3, x3, x5
    adc   x4, x4, x6
    .endm

    /* Add input carry, saves two 31-bit limbs and sets x1 to low carry */
    .macro CLOSE off, hi=xzr
    adds  x3, x3, x1
    adc   x4, x4, \hi
    and   x5, x3, #0x7fffffff
    ubfx  x6, x3, #31, #31
    stp   w5, w6, [x0, #\off]
    extr  x1, x4, x3, #62
    .endm
    
    .macro FINAL_CARRY off
    and   x5, x1, #0x7fffffff
    lsr   x6, x1, #31
    stp   w5, w6, [x0, #\off]
    .endm
    
ClassicalMult:
    /* Get external values and pointers to buffers. */
    adrp x3, :got:arr
    ldr x3, [x3, #:got_lo12:arr] /* x3 = Pointer to arr. */
    add x0, x3, x0, lsl #2       /* x0 = Pointer to multiplier. */
    adr x13, mm_jump_table       /* Jump to specific multiplication routine. */    
    add x13, x13, x12, lsl #2
    br x13
 
ClassicalMult1Limb:
    ldp w5, w6, [x0]             /* w5 = Multiplicand, w6 = Multiplier */
    umull x9, w5, w6             /* Compute product */
    lsr x10, x9, #31             /* Get most significant limb */
    bic w9, w9, #0x80000000      /* Get least significant limb */
    stp w9, w10, [x0]            /* Save product */
    ret

ClassicalMult2Limbs:
    ldp w5, w6, [x0]
    orr x5, x5, x6, lsl #31      /* x5 = Multiplicand */
    ldp w6, w7, [x0, #8]
    orr x6, x6, x7, lsl #31      /* x6 = Multiplier */
    mul x3, x5, x6
    umulh x4, x5, x6             /* x4: x3 = Product */
    and x9, x3, #0x7fffffff      /* Get limb #0 of product */
    ubfx x10, x3, #31, #31       /* Get limb #1 of product */
    stp w9, w10, [x0]            /* Store product in memory */
    extr x9, x5, x4, #62
    and x9, x9, #0x7fffffff      /* Get limb #2 of product */
    lsr x10, x4, #29             /* Get limb #3 of product */
    stp w9, w10, [x0, #8]        /* Store product in memory */
    ret

ClassicalMult3Limbs:
    ldp w5, w6, [x0]
    ldp w7, w8, [x0, #8]        /* w7: w6: w5 = Multiplicand */
    ldp w9, w10, [x0, #16]      /* w10: w9: w8 = Multiplier */
    umull x11, w5, w8
    lsr x12, x11, #31           /* Compute carry */
    bic w11, w11, #0x80000000   /* Get limb #0 of product */
    umaddl x12, w5, w9, x12
    umaddl x12, w6, w8, x12
    lsr x3, x12, #31            /* Compute carry */
    bic w12, w12, #0x80000000   /* Get limb #1 of product */
    stp w11, w12, [x0]          /* Store limbs #0 and #1 of product in memory */
    umaddl x3, w5, w10, x3
    umaddl x3, w6, w9, x3
    umaddl x3, w7, w8, x3
    lsr x4, x3, #31             /* Compute carry */
    bic w3, w3, #0x80000000     /* Get limb #2 of product */
    umaddl x4, w6, w10, x4
    umaddl x4, w7, w9, x4
    lsr x11, x4, #31            /* Compute carry */
    bic w4, w4, #0x80000000     /* Get limb #3 of product */
    stp w3, w4, [x0, #8]        /* Store limbs #2 and #3 of product in memory */
    umaddl x11, w7, w10, x11
    lsr x12, x11, #31           /* Get limb #5 of product */
    bic w11, w11, #0x80000000   /* Get limb #4 of product */
    stp w11, w12, [x0, #16]     /* Store limbs #4 and #5 of product in memory */
    ret
ClassicalMult4Limbs:
    /* Load multiplicand and multiplier into registers */
    ldp A0L, w3, [x0, #0]
    orr A0, A0, x3, lsl #31
    ldp A1L, w3, [x0, #8]
    orr A1, A1, x3, lsl #31
    ldp B0L, w3, [x0, #16]
    orr B0, x10, x3, lsl #31
    ldp B1L, w3, [x0, #24]
    orr B1, B1, x3, lsl #31
Multiply2DoubleLimbs:
    /* ---- Column #0: no input carry ---- */
    FIRST A0, B0
    and   x5, x3, #0x7fffffff
    ubfx  x6, x3, #31, #31
    stp   w5, w6, [x0, #0]
    extr  x1, x4, x3, #62

    /* ---- Column #1 (2 products) ---- */
    FIRST A1, B0
    MAC   A0, B1
    CLOSE 8

    /* ---- Column #2 (1 product) ---- */
    FIRST A1, B1
    CLOSE 16
    
    /* ---- Final carry (< 2^62): last two limbs ---- */
    FINAL_CARRY 24
EndMultiply2DoubleLimbs:
    ret
ClassicalMult5Limbs:
    /* Load multiplicand and multiplier into registers */
    ldr A2L, [x0, #16]
    ldp B0L, w3, [x0, #20]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #28]
    orr B1, B1, x3, lsl #31
    ldr B2L, [x0, #36]
    b 1f
ClassicalMult6Limbs:
    ldp A2L, w3, [x0, #16]
    orr A2, A2, x3, lsl #31
    ldp B0L, w3, [x0, #24]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #32]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #40]
    orr B2, B2, x3, lsl #31
1:    
    ldp A0L, w3, [x0, #0]
    orr A0, A0, x3, lsl #31
    ldp A1L, w3, [x0, #8]
    orr A1, A1, x3, lsl #31
Multiply3DoubleLimbs:
    /* ---- Column #0: no input carry ---- */
    FIRST A0, B0
    and   x5, x3, #0x7fffffff
    ubfx  x6, x3, #31, #31
    stp   w5, w6, [x0, #0]
    extr  x1, x4, x3, #62

    /* ---- Column #1 (2 products) ---- */
    FIRST A1, B0
    MAC   A0, B1
    CLOSE 8

    /* ---- Column #2 (3 products) ---- */
    FIRST A2, B0
    MAC   A1, B1
    MAC   A0, B2
    CLOSE 16

    /* ---- Column #3 (2 products) ---- */
    FIRST A2, B1
    MAC   A1, B2
    CLOSE 24

    /* ---- Column #4 (1 product) ---- */
    FIRST A2, B2
    CLOSE 32 
    tbnz w2, #0, EndMultiply3DoubleLimbs
  
    /* ---- Final carry (< 2^62): last two limbs ---- */
    FINAL_CARRY 40
EndMultiply3DoubleLimbs:
    ret
ClassicalMult7Limbs:
    /* Load multiplicand and multiplier into registers */
    ldr A3L, [x0, #24]
    ldp B0L, w3, [x0, #28]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #36]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #44]
    orr B2, B2, x3, lsl #31
    ldr B3L, [x0, #52]
    b 1f
ClassicalMult8Limbs:
    ldp A3L, w3, [x0, #24]
    orr A3, A3, x3, lsl #31
    ldp B0L, w3, [x0, #32]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #40]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #48]
    orr B2, B2, x3, lsl #31
    ldp B3L, w3, [x0, #56]
    orr B3, B3, x3, lsl #31
1:
    ldp A0L, w3, [x0, #0]
    orr A0, A0, x3, lsl #31
    ldp A1L, w3, [x0, #8]
    orr A1, A1, x3, lsl #31
    ldp A2L, w3, [x0, #16]
    orr A2, A2, x3, lsl #31
Multiply4DoubleLimbs:
    /* ---- Column #0: no input carry ---- */
    FIRST A0, B0
    and   x5, x3, #0x7fffffff
    ubfx  x6, x3, #31, #31
    stp   w5, w6, [x0, #0]
    extr  x1, x4, x3, #62

    /* ---- Column #1 (2 products) ---- */
    FIRST A1, B0
    MAC   A0, B1
    CLOSE 8

    /* ---- Column #2 (3 products) ---- */
    FIRST A2, B0
    MAC   A1, B1
    MAC   A0, B2
    CLOSE 16

    /* ---- Column #3 (4 products: the carry still fits in a word) ---- */
    FIRST A3, B0
    MAC   A2, B1
    MAC   A1, B2
    MAC   A0, B3
    CLOSE 24

    /* ---- Column #4 (3 products) ---- */
    FIRST A3, B1
    MAC   A2, B2
    MAC   A1, B3
    CLOSE 32
  
    /* ---- Column #5 (2 products) ---- */
    FIRST A3, B2
    MAC   A2, B3
    CLOSE 40
  
    /* ---- Column #6 (1 product) ---- */
    FIRST A3, B3
    CLOSE 48
    tbnz w2, #0, EndMultiply4DoubleLimbs
  
    /* ---- Final carry (< 2^62): last two limbs ---- */
    FINAL_CARRY 56
EndMultiply4DoubleLimbs:
    ret
ClassicalMult9Limbs:
    /* Load multiplicand and multiplier into registers */
    ldr A4L, [x0, #32]
    ldp B0L, w3, [x0, #36]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #44]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #52]
    orr B2, B2, x3, lsl #31
    ldp B3L, w3, [x0, #60]
    orr B3, B3, x3, lsl #31
    ldr B4L, [x0, #68]
    b 1f
ClassicalMult10Limbs:
    ldp A4L, w3, [x0, #32]
    orr A4, A4, x3, lsl #31
    ldp B0L, w3, [x0, #40]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #48]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #56]
    orr B2, B2, x3, lsl #31
    ldp B3L, w3, [x0, #64]
    orr B3, B3, x3, lsl #31
    ldp B4L, w3, [x0, #72]
    orr B4, B4, x3, lsl #31
1:
    ldp A0L, w3, [x0, #0]
    orr A0, A0, x3, lsl #31
    ldp A1L, w3, [x0, #8]
    orr A1, A1, x3, lsl #31
    ldp A2L, w3, [x0, #16]
    orr A2, A2, x3, lsl #31
    ldp A3L, w3, [x0, #24]
    orr A3, A3, x3, lsl #31
Multiply5DoubleLimbs:
    /* ---- Column #0: no input carry ---- */
    FIRST A0, B0
    and   x5, x3, #0x7fffffff
    ubfx  x6, x3, #31, #31
    stp   w5, w6, [x0, #0]
    extr  x1, x4, x3, #62

    /* ---- Column #1 (2 products) ---- */
    FIRST A1, B0
    MAC   A0, B1
    CLOSE 8

    /* ---- Column #2 (3 products) ---- */
    FIRST A2, B0
    MAC   A1, B1
    MAC   A0, B2
    CLOSE 16

    /* ---- Column #3 (4 products: the carry still fits in a word) ---- */
    FIRST A3, B0
    MAC   A2, B1
    MAC   A1, B2
    MAC   A0, B3
    CLOSE 24

    /* ---- Column #4: input carry in 1 word, output carry in 2 words ---- */
    FIRST A4, B0
    MAC   A3, B1
    MAC   A2, B2
    MAC   A1, B3
    MAC   A0, B4
    CLOSE 32
    lsr   x7, x4, #62

    /* ---- Column #5: input carry in 2 words, output carry in 1 word ---- */
    FIRST A4, B1
    MAC   A3, B2
    MAC   A2, B3
    MAC   A1, B4
    CLOSE 40, x7
  
    /* ---- Column #6 (3 products) ---- */
    FIRST A4, B2
    MAC   A3, B3
    MAC   A2, B4
    CLOSE 48
  
    /* ---- Column #7 (2 products) ---- */
    FIRST A4, B3
    MAC   A3, B4
    CLOSE 56
  
    /* ---- Column #8 (1 product) ---- */
    FIRST A4, B4
    CLOSE 64
    tbnz w2, #0, EndMultiply5DoubleLimbs
  
    /* ---- Final carry (< 2^62): last two limbs ---- */
    FINAL_CARRY 64
EndMultiply5DoubleLimbs:
    ret
ClassicalMult11Limbs:
    stp x19, x20, [sp, #-16]!
    /* Load multiplicand and multiplier into registers */
    ldr A5L, [x0, #40]
    ldp B0L, w3, [x0, #44]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #52]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #60]
    orr B2, B2, x3, lsl #31
    ldp B3L, w3, [x0, #68]
    orr B3, B3, x3, lsl #31
    ldp B4L, w3, [x0, #76]
    orr B4, B4, x3, lsl #31
    ldr B5L, [x0, #84]
    b 1f
ClassicalMult12Limbs:
    stp x19, x20, [sp, #-16]!
    /* Load multiplicand and multiplier into registers */
    ldp A5L, w3, [x0, #40]
    orr A5, A5, x3, lsl #31
    ldp B0L, w3, [x0, #48]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #56]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #64]
    orr B2, B2, x3, lsl #31
    ldp B3L, w3, [x0, #72]
    orr B3, B3, x3, lsl #31
    ldp B4L, w3, [x0, #80]
    orr B4, B4, x3, lsl #31
    ldp B5L, w3, [x0, #88]
    orr B5, B5, x3, lsl #31
1:
    ldp A0L, w3, [x0, #0]
    orr A0, A0, x3, lsl #31
    ldp A1L, w3, [x0, #8]
    orr A1, A1, x3, lsl #31
    ldp A2L, w3, [x0, #16]
    orr A2, A2, x3, lsl #31
    ldp A3L, w3, [x0, #24]
    orr A3, A3, x3, lsl #31
    ldp A4L, w3, [x0, #32]
    orr A4, A4, x3, lsl #31
Multiply6DoubleLimbs:
    /* ---- Column #0: no input carry ---- */
    FIRST A0, B0
    and   x5, x3, #0x7fffffff
    ubfx  x6, x3, #31, #31
    stp   w5, w6, [x0, #0]
    extr  x1, x4, x3, #62

    /* ---- Column #1 (2 products) ---- */
    FIRST A1, B0
    MAC   A0, B1
    CLOSE 8

    /* ---- Column #2 (3 products) ---- */
    FIRST A2, B0
    MAC   A1, B1
    MAC   A0, B2
    CLOSE 16

    /* ---- Column #3 (4 products: the carry still fits in a word) ---- */
    FIRST A3, B0
    MAC   A2, B1
    MAC   A1, B2
    MAC   A0, B3
    CLOSE 24

    /* ---- Column #4: input carry in 1 word, output carry in 2 words ---- */
    FIRST A4, B0
    MAC   A3, B1
    MAC   A2, B2
    MAC   A1, B3
    MAC   A0, B4
    CLOSE 32
    lsr   x7, x4, #62

    /* ---- Column #5 (6 products) ---- */
    FIRST A5, B0
    MAC   A4, B1
    MAC   A3, B2
    MAC   A2, B3
    MAC   A1, B4
    MAC   A0, B5
    CLOSE 40, x7
    lsr   x7, x4, #62

    /* ---- Column #6 (5 products) ---- */
    FIRST A5, B1
    MAC   A4, B2
    MAC   A3, B3
    MAC   A2, B4
    MAC   A1, B5
    CLOSE 48, x7
    lsr   x7, x4, #62

    /* ---- Column #7: input carry in 2 words, output carry in 1 word ---- */
    FIRST A5, B2
    MAC   A4, B3
    MAC   A3, B4
    MAC   A2, B5
    CLOSE 56, x7
  
    /* ---- Column #8 (3 products) ---- */
    FIRST A5, B3
    MAC   A4, B4
    MAC   A3, B5
    CLOSE 64
 
    /* ---- Column #9 (2 products) ---- */
    FIRST A5, B4
    MAC   A4, B5
    CLOSE 72
  
    /* ---- Column #10 (1 product) ---- */
    FIRST A5, B5
    CLOSE 80
    tbnz w2, #0, EndMultiply6DoubleLimbs
  
    /* ---- Final carry (< 2^62): last two limbs ---- */
    FINAL_CARRY 88
    /* ---- Restore registers and go out ---- */
EndMultiply6DoubleLimbs:
    ldp x19, x20, [sp], #16
    ret
ClassicalMult13Limbs:
    stp x19, x20, [sp, #-32]!
    stp x21, x22, [sp, #16]
    /* Load multiplicand and multiplier into registers */
    ldr A6L, [x0, #48]
    ldp B0L, w3, [x0, #52]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #60]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #68]
    orr B2, B2, x3, lsl #31
    ldp B3L, w3, [x0, #76]
    orr B3, B3, x3, lsl #31
    ldp B4L, w3, [x0, #84]
    orr B4, B4, x3, lsl #31
    ldp B5L, w3, [x0, #92]
    orr B5, B5, x3, lsl #31
    ldr B6L, [x0, #100]
    b 1f
ClassicalMult14Limbs:
    stp x19, x20, [sp, #-32]!
    stp x21, x22, [sp, #16]
    /* Load multiplicand and multiplier into registers */
    ldp A6L, w3, [x0, #48]
    orr A6, A6, x3, lsl #31
    ldp B0L, w3, [x0, #56]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #64]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #72]
    orr B2, B2, x3, lsl #31
    ldp B3L, w3, [x0, #80]
    orr B3, B3, x3, lsl #31
    ldp B4L, w3, [x0, #88]
    orr B4, B4, x3, lsl #31
    ldp B5L, w3, [x0, #96]
    orr B5, B5, x3, lsl #31
    ldp B6L, w3, [x0, #104]
    orr B6, B6, x3, lsl #31
1:
    ldp A0L, w3, [x0, #0]
    orr A0, A0, x3, lsl #31
    ldp A1L, w3, [x0, #8]
    orr A1, A1, x3, lsl #31
    ldp A2L, w3, [x0, #16]
    orr A2, A2, x3, lsl #31
    ldp A3L, w3, [x0, #24]
    orr A3, A3, x3, lsl #31
    ldp A4L, w3, [x0, #32]
    orr A4, A4, x3, lsl #31
    ldp A5L, w3, [x0, #40]
    orr A5, A5, x3, lsl #31
Multiply7DoubleLimbs:
    /* ---- Column #0: no input carry ---- */
    FIRST A0, B0
    and   x5, x3, #0x7fffffff
    ubfx  x6, x3, #31, #31
    stp   w5, w6, [x0, #0]
    extr  x1, x4, x3, #62

    /* ---- Column #1 (2 products) ---- */
    FIRST A1, B0
    MAC   A0, B1
    CLOSE 8

    /* ---- Column #2 (3 products) ---- */
    FIRST A2, B0
    MAC   A1, B1
    MAC   A0, B2
    CLOSE 16

    /* ---- Column #3 (4 products: the carry still fits in a word) ---- */
    FIRST A3, B0
    MAC   A2, B1
    MAC   A1, B2
    MAC   A0, B3
    CLOSE 24

    /* ---- Column #4: input carry in 1 word, output carry in 2 words ---- */
    FIRST A4, B0
    MAC   A3, B1
    MAC   A2, B2
    MAC   A1, B3
    MAC   A0, B4
    CLOSE 32
    lsr   x7, x4, #62

    /* ---- Column #5 (6 products) ---- */
    FIRST A5, B0
    MAC   A4, B1
    MAC   A3, B2
    MAC   A2, B3
    MAC   A1, B4
    MAC   A0, B5
    CLOSE 40, x7
    lsr   x7, x4, #62

    /* ---- Column #6 (7 products)---- */
    FIRST A6, B0
    MAC   A5, B1
    MAC   A4, B2
    MAC   A3, B3
    MAC   A2, B4
    MAC   A1, B5
    MAC   A0, B6
    CLOSE 48, x7
    lsr   x7, x4, #62

    /* ---- Column #7 (6 products) ---- */
    FIRST A6, B1
    MAC   A5, B2
    MAC   A4, B3
    MAC   A3, B4
    MAC   A2, B5
    MAC   A1, B6
    CLOSE 56, x7
    lsr   x7, x4, #62

    /* ---- Column #8 (5 products) ---- */
    FIRST A6, B2
    MAC   A5, B3
    MAC   A4, B4
    MAC   A3, B5
    MAC   A2, B6
    CLOSE 64, x7
    lsr   x7, x4, #62

    /* ---- Column #9: input carry in 2 words, output carry in 1 word ---- */
    FIRST A6, B3
    MAC   A5, B4
    MAC   A4, B5
    MAC   A3, B6
    CLOSE 72, x7
  
    /* ---- Column #10 (3 products) ---- */
    FIRST A6, B4
    MAC   A5, B5
    MAC   A4, B6
    CLOSE 80
  
    /* ---- Column #11 (2 products) ---- */
    FIRST A6, B5
    MAC   A5, B6
    CLOSE 88

    /* ---- Column #12 (1 product) ---- */
    FIRST A6, B6
    CLOSE 96
    tbnz w2, #0, EndMultiply7DoubleLimbs
   
    /* ---- Final carry (< 2^62): last two limbs ---- */
    FINAL_CARRY 104
    /* ---- Restore registers and go out ---- */
EndMultiply7DoubleLimbs:
    ldp   x21, x22, [sp, #16]
    ldp   x19, x20, [sp], #32
    ret
ClassicalMult15Limbs:
    stp x19, x20, [sp, #-48]!
    stp x21, x22, [sp, #16]
    stp x23, x24, [sp, #32]
    /* Load multiplicand and multiplier into registers */
    ldr A7L, [x0, #56]
    ldp B0L, w3, [x0, #60]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #68]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #76]
    orr B2, B2, x3, lsl #31
    ldp B3L, w3, [x0, #84]
    orr B3, B3, x3, lsl #31
    ldp B4L, w3, [x0, #92]
    orr B4, B4, x3, lsl #31
    ldp B5L, w3, [x0, #100]
    orr B5, B5, x3, lsl #31
    ldp B6L, w3, [x0, #108]
    orr B6, B6, x3, lsl #31
    ldr B7L, [x0, #116]
    b 1f
ClassicalMult16Limbs:
    stp x19, x20, [sp, #-48]!
    stp x21, x22, [sp, #16]
    stp x23, x24, [sp, #32]
    /* Load multiplicand and multiplier into registers */
    ldp A7L, w3, [x0, #56]
    orr A7, A7, x3, lsl #31
    ldp B0L, w3, [x0, #64]
    orr B0, B0, x3, lsl #31
    ldp B1L, w3, [x0, #72]
    orr B1, B1, x3, lsl #31
    ldp B2L, w3, [x0, #80]
    orr B2, B2, x3, lsl #31
    ldp B3L, w3, [x0, #88]
    orr B3, B3, x3, lsl #31
    ldp B4L, w3, [x0, #96]
    orr B4, B4, x3, lsl #31
    ldp B5L, w3, [x0, #104]
    orr B5, B5, x3, lsl #31
    ldp B6L, w3, [x0, #112]
    orr B6, B6, x3, lsl #31
    ldp B7L, w3, [x0, #120]
    orr B7, B7, x3, lsl #31
1:    
    ldp A0L, w3, [x0, #0]
    orr A0, A0, x3, lsl #31
    ldp A1L, w3, [x0, #8]
    orr A1, A1, x3, lsl #31
    ldp A2L, w3, [x0, #16]
    orr A2, A2, x3, lsl #31
    ldp A3L, w3, [x0, #24]
    orr A3, A3, x3, lsl #31
    ldp A4L, w3, [x0, #32]
    orr A4, A4, x3, lsl #31
    ldp A5L, w3, [x0, #40]
    orr A5, A5, x3, lsl #31
    ldp A6L, w3, [x0, #48]
    orr A6, A6, x3, lsl #31
Multiply8DoubleLimbs:
    /* ---- Column #0: no input carry ---- */
    FIRST A0, B0
    and   x5, x3, #0x7fffffff
    ubfx  x6, x3, #31, #31
    stp   w5, w6, [x0, #0]
    extr  x1, x4, x3, #62

    /* ---- Column #1 (2 products) ---- */
    FIRST A1, B0
    MAC   A0, B1
    CLOSE 8

    /* ---- Column #2 (3 products)  ---- */
    FIRST A2, B0
    MAC   A1, B1
    MAC   A0, B2
    CLOSE 16

    /* ---- Column #3 (4 products: the carry still fits in a word) ---- */
    FIRST A3, B0
    MAC   A2, B1
    MAC   A1, B2
    MAC   A0, B3
    CLOSE 24

    /* ---- Column #4: input carry in 1 word, output carry in 2 words ---- */
    FIRST A4, B0
    MAC   A3, B1
    MAC   A2, B2
    MAC   A1, B3
    MAC   A0, B4
    CLOSE 32
    lsr   x7, x4, #62

    /* ---- Column #5 (6 products) ---- */
    FIRST A5, B0
    MAC   A4, B1
    MAC   A3, B2
    MAC   A2, B3
    MAC   A1, B4
    MAC   A0, B5
    CLOSE 40, x7
    lsr   x7, x4, #62

    /* ---- Column #6 (7 products) ---- */
    FIRST A6, B0
    MAC   A5, B1
    MAC   A4, B2
    MAC   A3, B3
    MAC   A2, B4
    MAC   A1, B5
    MAC   A0, B6
    CLOSE 48, x7
    lsr   x7, x4, #62

    /* ---- Column #7 (8 products) ---- */
    FIRST A7, B0
    MAC   A6, B1
    MAC   A5, B2
    MAC   A4, B3
    MAC   A3, B4
    MAC   A2, B5
    MAC   A1, B6
    MAC   A0, B7
    CLOSE 56, x7
    lsr   x7, x4, #62

    /* ---- Column #8 (7 products) ---- */
    FIRST A7, B1
    MAC   A6, B2
    MAC   A5, B3
    MAC   A4, B4
    MAC   A3, B5
    MAC   A2, B6
    MAC   A1, B7
    CLOSE 64, x7
    lsr   x7, x4, #62

    /* ---- Column #9 (6 products) ---- */
    FIRST A7, B2
    MAC   A6, B3
    MAC   A5, B4
    MAC   A4, B5
    MAC   A3, B6
    MAC   A2, B7
    CLOSE 72, x7
    lsr   x7, x4, #62

    /* ---- Column #10 (5 products) ---- */
    FIRST A7, B3
    MAC   A6, B4
    MAC   A5, B5
    MAC   A4, B6
    MAC   A3, B7
    CLOSE 80, x7
    lsr   x7, x4, #62

    /* ---- Column #11: input carry in 2 words, output carry in 1 word ---- */
    FIRST A7, B4
    MAC   A6, B5
    MAC   A5, B6
    MAC   A4, B7
    CLOSE 88, x7

    /* ---- Column #12 (3 products) ---- */
    FIRST A7, B5
    MAC   A6, B6
    MAC   A5, B7
    CLOSE 96

    /* ---- Column #13 (2 products) ---- */
    FIRST A7, B6
    MAC   A6, B7
    CLOSE 104

    /* ---- Column #14 (1 product) ---- */
    FIRST A7, B7
    CLOSE 112
    tbnz w2, #0, EndMultiply8DoubleLimbs
    
    /* ---- Final carry (< 2^62): last two limbs ---- */
    FINAL_CARRY 120
    /* ---- Restore registers and go out ---- */
EndMultiply8DoubleLimbs:
    ldp   x23, x24, [sp, #32]
    ldp   x21, x22, [sp, #16]
    ldp   x19, x20, [sp], #48
    ret

mm_jump_table:
    nop
    b ClassicalMult1Limb
    b ClassicalMult2Limbs
    b ClassicalMult3Limbs
    b ClassicalMult4Limbs
    b ClassicalMult5Limbs
    b ClassicalMult6Limbs
    b ClassicalMult7Limbs
    b ClassicalMult8Limbs
    b ClassicalMult9Limbs
    b ClassicalMult10Limbs
    b ClassicalMult11Limbs
    b ClassicalMult12Limbs
    b ClassicalMult13Limbs
    b ClassicalMult14Limbs
    b ClassicalMult15Limbs
    b ClassicalMult16Limbs
