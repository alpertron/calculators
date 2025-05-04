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
 x1: 0x7FFFFFFF (mask)
 x4: x3 = Sum of products of double limbs in a column during long multiplication.
 x6: x5 = Product of double limbs in a column during long multiplication.
 x8, x9, ... = Double limbs of multplicand and multiplier.
 */

ClassicalMult:
    /* Get external values and pointers to buffers. */
    adrp x3, :got:arr
    ldr x3, [x3, #:got_lo12:arr]     /* x3 = Pointer to arr. */
    add x0, x3, x0, lsl #2           /* x0 = Pointer to multiplier. */
    mov w1, #0x7FFFFFFF
    cmp w2, #1
    bne Test2Limbs

ClassicalMult1Limb:
    ldp w5, w6, [x0]             /* w5 = Multiplicand, w6 = Multiplier */
    umull x9, w5, w6             /* Compute product */
    lsr x10, x9, #31             /* Get most significant limb */
    bic w9, w9, #0x80000000      /* Get least significant limb */
    stp w9, w10, [x0]            /* Save product */
    ret

Test2Limbs:
    cmp w2, #2
    bne Test3Limbs

ClassicalMult2Limbs:
    ldp w5, w6, [x0]
    orr x5, x5, x6, lsl #31      /* x5 = Multiplicand */
    ldp w6, w7, [x0, #8]
    orr x6, x6, x7, lsl #31      /* x6 = Multiplier */
    mul x3, x5, x6
    umulh x4, x5, x6             /* x5: x4 = Product */
    and x9, x1, x3               /* Get limb #0 of product */
    and x10, x1, x3, lsr #31     /* Get limb #1 of product */
    stp w9, w10, [x0]            /* Store product in memory */
    lsl x9, x4, #2
    orr x9, x9, x3, lsr #62
    and x9, x1, x9               /* Get limb #2 of product */
    lsr x10, x4, #29             /* Get limb #3 of product */
    stp w9, w10, [x0, #8]        /* Store product in memory */
    ret

Test3Limbs:
    cmp w2, #3
    bne MoreThan3Limbs

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
MoreThan3Limbs:
    cmp w2, #4
    bhi MoreThan4Limbs
    /* Load multiplicand and multiplier into registers */
    ldp w8, w3, [x0, #0]
    orr x8, x8, x3, lsl #31
    ldp w9, w3, [x0, #8]
    orr x9, x9, x3, lsl #31
    ldp w10, w3, [x0, #16]
    orr x10, x10, x3, lsl #31
    ldp w11, w3, [x0, #24]
    orr x11, x11, x3, lsl #31
Multiply2DoubleLimbs:
    /* Perform first multiplication */
    mul x3, x8, x10
    umulh x4, x8, x10
    and x6, x1, x3               /* Get limb #0 of product */
    and x7, x1, x3, lsr #31      /* Get limb #1 of product */
    stp w6, w7, [x0]             /* Save these limbs */
    umulh x6, x9, x10
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #1: sum of products */
    mul x5, x9, x10
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x8, x11
    adds x3, x3, x5
    umulh x6, x8, x11
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #2 of product */
    and x7, x1, x3, lsr #31    /* Get limb #3 of product */
    stp w6, w7, [x0, #8]      /* Save these limbs */
    umulh x6, x9, x11
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #2: sum of products */
    mul x5, x9, x11
    adds x3, x3, x5
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #4 of product */
    and x7, x1, x3, lsr #31    /* Get limb #5 of product */
    stp w6, w7, [x0, #16]      /* Save these limbs */
    tbnz w2, #0, EndMultiply2DoubleLimbs
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x6, x4, #2
    orr x3, x6, x3, lsr #62
    /* Store carry in product */
    and x6, x1, x3             /* Get limb #6 of product */
    and x7, x1, x3, lsr #31    /* Get limb #7 of product */
    stp w6, w7, [x0, #24]      /* Save these limbs */
EndMultiply2DoubleLimbs:
    ret
MoreThan4Limbs:
    cmp w2, #6
    bhi MoreThan6Limbs
    /* Load multiplicand and multiplier into registers */
    ldp w8, w3, [x0, #0]
    orr x8, x8, x3, lsl #31
    ldp w9, w3, [x0, #8]
    orr x9, x9, x3, lsl #31
    tbnz w2, #0, LengthIs5Limbs
    ldp w10, w3, [x0, #16]
    orr x10, x10, x3, lsl #31
    ldp w11, w3, [x0, #24]
    orr x11, x11, x3, lsl #31
    ldp w12, w3, [x0, #32]
    orr x12, x12, x3, lsl #31
    ldp w13, w3, [x0, #40]
    orr x13, x13, x3, lsl #31
    b Multiply3DoubleLimbs
LengthIs5Limbs:
    ldr w10, [x0, #16]
    ldp w11, w3, [x0, #20]
    orr x11, x11, x3, lsl #31
    ldp w12, w3, [x0, #28]
    orr x12, x12, x3, lsl #31
    ldr w13, [x0, #36]
Multiply3DoubleLimbs:
    /* Perform first multiplication */
    mul x3, x8, x11
    umulh x4, x8, x11
    and x6, x1, x3               /* Get limb #0 of product */
    and x7, x1, x3, lsr #31      /* Get limb #1 of product */
    stp w6, w7, [x0]             /* Save these limbs */
    umulh x6, x9, x11
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #1: sum of products */
    mul x5, x9, x11
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x8, x12
    adds x3, x3, x5
    umulh x6, x8, x12
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #2 of product */
    and x7, x1, x3, lsr #31    /* Get limb #3 of product */
    stp w6, w7, [x0, #8]      /* Save these limbs */
    umulh x6, x10, x11
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #2: sum of products */
    mul x5, x10, x11
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x9, x12
    adds x3, x3, x5
    umulh x6, x9, x12
    adc x4, x4, x6
    mul x5, x8, x13
    adds x3, x3, x5
    umulh x6, x8, x13
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #4 of product */
    and x7, x1, x3, lsr #31    /* Get limb #5 of product */
    stp w6, w7, [x0, #16]      /* Save these limbs */
    umulh x6, x10, x12
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #3: sum of products */
    mul x5, x10, x12
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x9, x13
    adds x3, x3, x5
    umulh x6, x9, x13
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #6 of product */
    and x7, x1, x3, lsr #31    /* Get limb #7 of product */
    stp w6, w7, [x0, #24]      /* Save these limbs */
    umulh x6, x10, x13
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #4: sum of products */
    mul x5, x10, x13
    adds x3, x3, x5
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #8 of product */
    and x7, x1, x3, lsr #31    /* Get limb #9 of product */
    stp w6, w7, [x0, #32]      /* Save these limbs */
    tbnz w2, #0, EndMultiply3DoubleLimbs
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x6, x4, #2
    orr x3, x6, x3, lsr #62
    /* Store carry in product */
    and x6, x1, x3             /* Get limb #10 of product */
    and x7, x1, x3, lsr #31    /* Get limb #11 of product */
    stp w6, w7, [x0, #40]      /* Save these limbs */
EndMultiply3DoubleLimbs:
    ret
MoreThan6Limbs:
    cmp w2, #8
    bhi MoreThan8Limbs
    /* Load multiplicand and multiplier into registers */
    ldp w8, w3, [x0, #0]
    orr x8, x8, x3, lsl #31
    ldp w9, w3, [x0, #8]
    orr x9, x9, x3, lsl #31
    ldp w10, w3, [x0, #16]
    orr x10, x10, x3, lsl #31
    tbnz w2, #0, LengthIs7Limbs
    ldp w11, w3, [x0, #24]
    orr x11, x11, x3, lsl #31
    ldp w12, w3, [x0, #32]
    orr x12, x12, x3, lsl #31
    ldp w13, w3, [x0, #40]
    orr x13, x13, x3, lsl #31
    ldp w14, w3, [x0, #48]
    orr x14, x14, x3, lsl #31
    ldp w15, w3, [x0, #56]
    orr x15, x15, x3, lsl #31
    b Multiply4DoubleLimbs
LengthIs7Limbs:
    ldr w11, [x0, #24]
    ldp w12, w3, [x0, #28]
    orr x12, x12, x3, lsl #31
    ldp w13, w3, [x0, #36]
    orr x13, x13, x3, lsl #31
    ldp w14, w3, [x0, #44]
    orr x14, x14, x3, lsl #31
    ldr w15, [x0, #52]
Multiply4DoubleLimbs:
    /* Perform first multiplication */
    mul x3, x8, x12
    umulh x4, x8, x12
    and x6, x1, x3               /* Get limb #0 of product */
    and x7, x1, x3, lsr #31      /* Get limb #1 of product */
    stp w6, w7, [x0]             /* Save these limbs */
    umulh x6, x9, x12
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #1: sum of products */
    mul x5, x9, x12
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x8, x13
    adds x3, x3, x5
    umulh x6, x8, x13
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #2 of product */
    and x7, x1, x3, lsr #31    /* Get limb #3 of product */
    stp w6, w7, [x0, #8]      /* Save these limbs */
    umulh x6, x10, x12
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #2: sum of products */
    mul x5, x10, x12
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x9, x13
    adds x3, x3, x5
    umulh x6, x9, x13
    adc x4, x4, x6
    mul x5, x8, x14
    adds x3, x3, x5
    umulh x6, x8, x14
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #4 of product */
    and x7, x1, x3, lsr #31    /* Get limb #5 of product */
    stp w6, w7, [x0, #16]      /* Save these limbs */
    umulh x6, x11, x12
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #3: sum of products */
    mul x5, x11, x12
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x10, x13
    adds x3, x3, x5
    umulh x6, x10, x13
    adc x4, x4, x6
    mul x5, x9, x14
    adds x3, x3, x5
    umulh x6, x9, x14
    adc x4, x4, x6
    mul x5, x8, x15
    adds x3, x3, x5
    umulh x6, x8, x15
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #6 of product */
    and x7, x1, x3, lsr #31    /* Get limb #7 of product */
    stp w6, w7, [x0, #24]      /* Save these limbs */
    umulh x6, x11, x13
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #4: sum of products */
    mul x5, x11, x13
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x10, x14
    adds x3, x3, x5
    umulh x6, x10, x14
    adc x4, x4, x6
    mul x5, x9, x15
    adds x3, x3, x5
    umulh x6, x9, x15
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #8 of product */
    and x7, x1, x3, lsr #31    /* Get limb #9 of product */
    stp w6, w7, [x0, #32]      /* Save these limbs */
    umulh x6, x11, x14
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #5: sum of products */
    mul x5, x11, x14
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x10, x15
    adds x3, x3, x5
    umulh x6, x10, x15
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #10 of product */
    and x7, x1, x3, lsr #31    /* Get limb #11 of product */
    stp w6, w7, [x0, #40]      /* Save these limbs */
    umulh x6, x11, x15
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #6: sum of products */
    mul x5, x11, x15
    adds x3, x3, x5
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #12 of product */
    and x7, x1, x3, lsr #31    /* Get limb #13 of product */
    stp w6, w7, [x0, #48]      /* Save these limbs */
    tbnz w2, #0, EndMultiply4DoubleLimbs
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x6, x4, #2
    orr x3, x6, x3, lsr #62
    /* Store carry in product */
    and x6, x1, x3             /* Get limb #14 of product */
    and x7, x1, x3, lsr #31    /* Get limb #15 of product */
    stp w6, w7, [x0, #56]      /* Save these limbs */
EndMultiply4DoubleLimbs:
    ret
MoreThan8Limbs:
    cmp w2, #10
    bhi MoreThan10Limbs
    /* Load multiplicand and multiplier into registers */
    ldp w8, w3, [x0, #0]
    orr x8, x8, x3, lsl #31
    ldp w9, w3, [x0, #8]
    orr x9, x9, x3, lsl #31
    ldp w10, w3, [x0, #16]
    orr x10, x10, x3, lsl #31
    ldp w11, w3, [x0, #24]
    orr x11, x11, x3, lsl #31
    tbnz w2, #0, LengthIs9Limbs
    ldp w12, w3, [x0, #32]
    orr x12, x12, x3, lsl #31
    ldp w13, w3, [x0, #40]
    orr x13, x13, x3, lsl #31
    ldp w14, w3, [x0, #48]
    orr x14, x14, x3, lsl #31
    ldp w15, w3, [x0, #56]
    orr x15, x15, x3, lsl #31
    ldp w16, w3, [x0, #64]
    orr x16, x16, x3, lsl #31
    ldp w17, w3, [x0, #72]
    orr x17, x17, x3, lsl #31
    b Multiply5DoubleLimbs
LengthIs9Limbs:
    ldr w12, [x0, #32]
    ldp w13, w3, [x0, #36]
    orr x13, x13, x3, lsl #31
    ldp w14, w3, [x0, #44]
    orr x14, x14, x3, lsl #31
    ldp w15, w3, [x0, #52]
    orr x15, x15, x3, lsl #31
    ldp w16, w3, [x0, #60]
    orr x16, x16, x3, lsl #31
    ldr w17, [x0, #68]
Multiply5DoubleLimbs:
    /* Perform first multiplication */
    mul x3, x8, x13
    umulh x4, x8, x13
    and x6, x1, x3               /* Get limb #0 of product */
    and x7, x1, x3, lsr #31      /* Get limb #1 of product */
    stp w6, w7, [x0]             /* Save these limbs */
    umulh x6, x9, x13
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #1: sum of products */
    mul x5, x9, x13
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x8, x14
    adds x3, x3, x5
    umulh x6, x8, x14
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #2 of product */
    and x7, x1, x3, lsr #31    /* Get limb #3 of product */
    stp w6, w7, [x0, #8]      /* Save these limbs */
    umulh x6, x10, x13
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #2: sum of products */
    mul x5, x10, x13
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x9, x14
    adds x3, x3, x5
    umulh x6, x9, x14
    adc x4, x4, x6
    mul x5, x8, x15
    adds x3, x3, x5
    umulh x6, x8, x15
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #4 of product */
    and x7, x1, x3, lsr #31    /* Get limb #5 of product */
    stp w6, w7, [x0, #16]      /* Save these limbs */
    umulh x6, x11, x13
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #3: sum of products */
    mul x5, x11, x13
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x10, x14
    adds x3, x3, x5
    umulh x6, x10, x14
    adc x4, x4, x6
    mul x5, x9, x15
    adds x3, x3, x5
    umulh x6, x9, x15
    adc x4, x4, x6
    mul x5, x8, x16
    adds x3, x3, x5
    umulh x6, x8, x16
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #6 of product */
    and x7, x1, x3, lsr #31    /* Get limb #7 of product */
    stp w6, w7, [x0, #24]      /* Save these limbs */
    umulh x6, x12, x13
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #4: sum of products */
    mul x5, x12, x13
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x11, x14
    adds x3, x3, x5
    umulh x6, x11, x14
    adc x4, x4, x6
    mul x5, x10, x15
    adds x3, x3, x5
    umulh x6, x10, x15
    adc x4, x4, x6
    mul x5, x9, x16
    adds x3, x3, x5
    umulh x6, x9, x16
    adc x4, x4, x6
    mul x5, x8, x17
    adds x3, x3, x5
    umulh x6, x8, x17
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #8 of product */
    and x7, x1, x3, lsr #31    /* Get limb #9 of product */
    stp w6, w7, [x0, #32]      /* Save these limbs */
    umulh x6, x12, x14
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #5: sum of products */
    mul x5, x12, x14
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x11, x15
    adds x3, x3, x5
    umulh x6, x11, x15
    adc x4, x4, x6
    mul x5, x10, x16
    adds x3, x3, x5
    umulh x6, x10, x16
    adc x4, x4, x6
    mul x5, x9, x17
    adds x3, x3, x5
    umulh x6, x9, x17
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #10 of product */
    and x7, x1, x3, lsr #31    /* Get limb #11 of product */
    stp w6, w7, [x0, #40]      /* Save these limbs */
    umulh x6, x12, x15
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #6: sum of products */
    mul x5, x12, x15
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x11, x16
    adds x3, x3, x5
    umulh x6, x11, x16
    adc x4, x4, x6
    mul x5, x10, x17
    adds x3, x3, x5
    umulh x6, x10, x17
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #12 of product */
    and x7, x1, x3, lsr #31    /* Get limb #13 of product */
    stp w6, w7, [x0, #48]      /* Save these limbs */
    umulh x6, x12, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #7: sum of products */
    mul x5, x12, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x11, x17
    adds x3, x3, x5
    umulh x6, x11, x17
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #14 of product */
    and x7, x1, x3, lsr #31    /* Get limb #15 of product */
    stp w6, w7, [x0, #56]      /* Save these limbs */
    umulh x6, x12, x17
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #8: sum of products */
    mul x5, x12, x17
    adds x3, x3, x5
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #16 of product */
    and x7, x1, x3, lsr #31    /* Get limb #17 of product */
    stp w6, w7, [x0, #64]      /* Save these limbs */
    tbnz w2, #0, EndMultiply5DoubleLimbs
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x6, x4, #2
    orr x3, x6, x3, lsr #62
    /* Store carry in product */
    and x6, x1, x3             /* Get limb #18 of product */
    and x7, x1, x3, lsr #31    /* Get limb #19 of product */
    stp w6, w7, [x0, #72]      /* Save these limbs */
EndMultiply5DoubleLimbs:
    ret
MoreThan10Limbs:
    cmp w2, #12
    bhi MoreThan12Limbs
    stp x19, x20, [sp, #-16]!
    /* Load multiplicand and multiplier into registers */
    ldp w8, w3, [x0, #0]
    orr x8, x8, x3, lsl #31
    ldp w9, w3, [x0, #8]
    orr x9, x9, x3, lsl #31
    ldp w10, w3, [x0, #16]
    orr x10, x10, x3, lsl #31
    ldp w11, w3, [x0, #24]
    orr x11, x11, x3, lsl #31
    ldp w12, w3, [x0, #32]
    orr x12, x12, x3, lsl #31
    tbnz w2, #0, LengthIs11Limbs
    ldp w13, w3, [x0, #40]
    orr x13, x13, x3, lsl #31
    ldp w14, w3, [x0, #48]
    orr x14, x14, x3, lsl #31
    ldp w15, w3, [x0, #56]
    orr x15, x15, x3, lsl #31
    ldp w16, w3, [x0, #64]
    orr x16, x16, x3, lsl #31
    ldp w17, w3, [x0, #72]
    orr x17, x17, x3, lsl #31
    ldp w19, w3, [x0, #80]
    orr x19, x19, x3, lsl #31
    ldp w20, w3, [x0, #88]
    orr x20, x20, x3, lsl #31
    b Multiply6DoubleLimbs
LengthIs11Limbs:
    ldr w13, [x0, #40]
    ldp w14, w3, [x0, #44]
    orr x14, x14, x3, lsl #31
    ldp w15, w3, [x0, #52]
    orr x15, x15, x3, lsl #31
    ldp w16, w3, [x0, #60]
    orr x16, x16, x3, lsl #31
    ldp w17, w3, [x0, #68]
    orr x17, x17, x3, lsl #31
    ldp w19, w3, [x0, #76]
    orr x19, x19, x3, lsl #31
    ldr w20, [x0, #84]
Multiply6DoubleLimbs:
    /* Perform first multiplication */
    mul x3, x8, x14
    umulh x4, x8, x14
    and x6, x1, x3               /* Get limb #0 of product */
    and x7, x1, x3, lsr #31      /* Get limb #1 of product */
    stp w6, w7, [x0]             /* Save these limbs */
    umulh x6, x9, x14
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #1: sum of products */
    mul x5, x9, x14
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x8, x15
    adds x3, x3, x5
    umulh x6, x8, x15
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #2 of product */
    and x7, x1, x3, lsr #31    /* Get limb #3 of product */
    stp w6, w7, [x0, #8]      /* Save these limbs */
    umulh x6, x10, x14
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #2: sum of products */
    mul x5, x10, x14
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x9, x15
    adds x3, x3, x5
    umulh x6, x9, x15
    adc x4, x4, x6
    mul x5, x8, x16
    adds x3, x3, x5
    umulh x6, x8, x16
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #4 of product */
    and x7, x1, x3, lsr #31    /* Get limb #5 of product */
    stp w6, w7, [x0, #16]      /* Save these limbs */
    umulh x6, x11, x14
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #3: sum of products */
    mul x5, x11, x14
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x10, x15
    adds x3, x3, x5
    umulh x6, x10, x15
    adc x4, x4, x6
    mul x5, x9, x16
    adds x3, x3, x5
    umulh x6, x9, x16
    adc x4, x4, x6
    mul x5, x8, x17
    adds x3, x3, x5
    umulh x6, x8, x17
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #6 of product */
    and x7, x1, x3, lsr #31    /* Get limb #7 of product */
    stp w6, w7, [x0, #24]      /* Save these limbs */
    umulh x6, x12, x14
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #4: sum of products */
    mul x5, x12, x14
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x11, x15
    adds x3, x3, x5
    umulh x6, x11, x15
    adc x4, x4, x6
    mul x5, x10, x16
    adds x3, x3, x5
    umulh x6, x10, x16
    adc x4, x4, x6
    mul x5, x9, x17
    adds x3, x3, x5
    umulh x6, x9, x17
    adc x4, x4, x6
    mul x5, x8, x19
    adds x3, x3, x5
    umulh x6, x8, x19
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #8 of product */
    and x7, x1, x3, lsr #31    /* Get limb #9 of product */
    stp w6, w7, [x0, #32]      /* Save these limbs */
    umulh x6, x13, x14
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #5: sum of products */
    mul x5, x13, x14
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x12, x15
    adds x3, x3, x5
    umulh x6, x12, x15
    adc x4, x4, x6
    mul x5, x11, x16
    adds x3, x3, x5
    umulh x6, x11, x16
    adc x4, x4, x6
    mul x5, x10, x17
    adds x3, x3, x5
    umulh x6, x10, x17
    adc x4, x4, x6
    mul x5, x9, x19
    adds x3, x3, x5
    umulh x6, x9, x19
    adc x4, x4, x6
    mul x5, x8, x20
    adds x3, x3, x5
    umulh x6, x8, x20
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #10 of product */
    and x7, x1, x3, lsr #31    /* Get limb #11 of product */
    stp w6, w7, [x0, #40]      /* Save these limbs */
    umulh x6, x13, x15
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #6: sum of products */
    mul x5, x13, x15
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x12, x16
    adds x3, x3, x5
    umulh x6, x12, x16
    adc x4, x4, x6
    mul x5, x11, x17
    adds x3, x3, x5
    umulh x6, x11, x17
    adc x4, x4, x6
    mul x5, x10, x19
    adds x3, x3, x5
    umulh x6, x10, x19
    adc x4, x4, x6
    mul x5, x9, x20
    adds x3, x3, x5
    umulh x6, x9, x20
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #12 of product */
    and x7, x1, x3, lsr #31    /* Get limb #13 of product */
    stp w6, w7, [x0, #48]      /* Save these limbs */
    umulh x6, x13, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #7: sum of products */
    mul x5, x13, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x12, x17
    adds x3, x3, x5
    umulh x6, x12, x17
    adc x4, x4, x6
    mul x5, x11, x19
    adds x3, x3, x5
    umulh x6, x11, x19
    adc x4, x4, x6
    mul x5, x10, x20
    adds x3, x3, x5
    umulh x6, x10, x20
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #14 of product */
    and x7, x1, x3, lsr #31    /* Get limb #15 of product */
    stp w6, w7, [x0, #56]      /* Save these limbs */
    umulh x6, x13, x17
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #8: sum of products */
    mul x5, x13, x17
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x12, x19
    adds x3, x3, x5
    umulh x6, x12, x19
    adc x4, x4, x6
    mul x5, x11, x20
    adds x3, x3, x5
    umulh x6, x11, x20
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #16 of product */
    and x7, x1, x3, lsr #31    /* Get limb #17 of product */
    stp w6, w7, [x0, #64]      /* Save these limbs */
    umulh x6, x13, x19
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #9: sum of products */
    mul x5, x13, x19
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x12, x20
    adds x3, x3, x5
    umulh x6, x12, x20
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #18 of product */
    and x7, x1, x3, lsr #31    /* Get limb #19 of product */
    stp w6, w7, [x0, #72]      /* Save these limbs */
    umulh x6, x13, x20
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #10: sum of products */
    mul x5, x13, x20
    adds x3, x3, x5
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #20 of product */
    and x7, x1, x3, lsr #31    /* Get limb #21 of product */
    stp w6, w7, [x0, #80]      /* Save these limbs */
    tbnz w2, #0, EndMultiply6DoubleLimbs
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x6, x4, #2
    orr x3, x6, x3, lsr #62
    /* Store carry in product */
    and x6, x1, x3             /* Get limb #22 of product */
    and x7, x1, x3, lsr #31    /* Get limb #23 of product */
    stp w6, w7, [x0, #88]      /* Save these limbs */
EndMultiply6DoubleLimbs:
    ldp x19, x20, [sp], #16
    ret
MoreThan12Limbs:
    cmp w2, #14
    bhi MoreThan14Limbs
    stp x19, x20, [sp, #-16]!
    stp x21, x22, [sp, #-16]!
    /* Load multiplicand and multiplier into registers */
    ldp w8, w3, [x0, #0]
    orr x8, x8, x3, lsl #31
    ldp w9, w3, [x0, #8]
    orr x9, x9, x3, lsl #31
    ldp w10, w3, [x0, #16]
    orr x10, x10, x3, lsl #31
    ldp w11, w3, [x0, #24]
    orr x11, x11, x3, lsl #31
    ldp w12, w3, [x0, #32]
    orr x12, x12, x3, lsl #31
    ldp w13, w3, [x0, #40]
    orr x13, x13, x3, lsl #31
    tbnz w2, #0, LengthIs13Limbs
    ldp w14, w3, [x0, #48]
    orr x14, x14, x3, lsl #31
    ldp w15, w3, [x0, #56]
    orr x15, x15, x3, lsl #31
    ldp w16, w3, [x0, #64]
    orr x16, x16, x3, lsl #31
    ldp w17, w3, [x0, #72]
    orr x17, x17, x3, lsl #31
    ldp w19, w3, [x0, #80]
    orr x19, x19, x3, lsl #31
    ldp w20, w3, [x0, #88]
    orr x20, x20, x3, lsl #31
    ldp w21, w3, [x0, #96]
    orr x21, x21, x3, lsl #31
    ldp w22, w3, [x0, #104]
    orr x22, x22, x3, lsl #31
    b Multiply7DoubleLimbs
LengthIs13Limbs:
    ldr w14, [x0, #48]
    ldp w15, w3, [x0, #52]
    orr x15, x15, x3, lsl #31
    ldp w16, w3, [x0, #60]
    orr x16, x16, x3, lsl #31
    ldp w17, w3, [x0, #68]
    orr x17, x17, x3, lsl #31
    ldp w19, w3, [x0, #76]
    orr x19, x19, x3, lsl #31
    ldp w20, w3, [x0, #84]
    orr x20, x20, x3, lsl #31
    ldp w21, w3, [x0, #92]
    orr x21, x21, x3, lsl #31
    ldr w22, [x0, #100]
Multiply7DoubleLimbs:
    /* Perform first multiplication */
    mul x3, x8, x15
    umulh x4, x8, x15
    and x6, x1, x3               /* Get limb #0 of product */
    and x7, x1, x3, lsr #31      /* Get limb #1 of product */
    stp w6, w7, [x0]             /* Save these limbs */
    umulh x6, x9, x15
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #1: sum of products */
    mul x5, x9, x15
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x8, x16
    adds x3, x3, x5
    umulh x6, x8, x16
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #2 of product */
    and x7, x1, x3, lsr #31    /* Get limb #3 of product */
    stp w6, w7, [x0, #8]      /* Save these limbs */
    umulh x6, x10, x15
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #2: sum of products */
    mul x5, x10, x15
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x9, x16
    adds x3, x3, x5
    umulh x6, x9, x16
    adc x4, x4, x6
    mul x5, x8, x17
    adds x3, x3, x5
    umulh x6, x8, x17
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #4 of product */
    and x7, x1, x3, lsr #31    /* Get limb #5 of product */
    stp w6, w7, [x0, #16]      /* Save these limbs */
    umulh x6, x11, x15
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #3: sum of products */
    mul x5, x11, x15
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x10, x16
    adds x3, x3, x5
    umulh x6, x10, x16
    adc x4, x4, x6
    mul x5, x9, x17
    adds x3, x3, x5
    umulh x6, x9, x17
    adc x4, x4, x6
    mul x5, x8, x19
    adds x3, x3, x5
    umulh x6, x8, x19
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #6 of product */
    and x7, x1, x3, lsr #31    /* Get limb #7 of product */
    stp w6, w7, [x0, #24]      /* Save these limbs */
    umulh x6, x12, x15
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #4: sum of products */
    mul x5, x12, x15
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x11, x16
    adds x3, x3, x5
    umulh x6, x11, x16
    adc x4, x4, x6
    mul x5, x10, x17
    adds x3, x3, x5
    umulh x6, x10, x17
    adc x4, x4, x6
    mul x5, x9, x19
    adds x3, x3, x5
    umulh x6, x9, x19
    adc x4, x4, x6
    mul x5, x8, x20
    adds x3, x3, x5
    umulh x6, x8, x20
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #8 of product */
    and x7, x1, x3, lsr #31    /* Get limb #9 of product */
    stp w6, w7, [x0, #32]      /* Save these limbs */
    umulh x6, x13, x15
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #5: sum of products */
    mul x5, x13, x15
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x12, x16
    adds x3, x3, x5
    umulh x6, x12, x16
    adc x4, x4, x6
    mul x5, x11, x17
    adds x3, x3, x5
    umulh x6, x11, x17
    adc x4, x4, x6
    mul x5, x10, x19
    adds x3, x3, x5
    umulh x6, x10, x19
    adc x4, x4, x6
    mul x5, x9, x20
    adds x3, x3, x5
    umulh x6, x9, x20
    adc x4, x4, x6
    mul x5, x8, x21
    adds x3, x3, x5
    umulh x6, x8, x21
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #10 of product */
    and x7, x1, x3, lsr #31    /* Get limb #11 of product */
    stp w6, w7, [x0, #40]      /* Save these limbs */
    umulh x6, x14, x15
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #6: sum of products */
    mul x5, x14, x15
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x13, x16
    adds x3, x3, x5
    umulh x6, x13, x16
    adc x4, x4, x6
    mul x5, x12, x17
    adds x3, x3, x5
    umulh x6, x12, x17
    adc x4, x4, x6
    mul x5, x11, x19
    adds x3, x3, x5
    umulh x6, x11, x19
    adc x4, x4, x6
    mul x5, x10, x20
    adds x3, x3, x5
    umulh x6, x10, x20
    adc x4, x4, x6
    mul x5, x9, x21
    adds x3, x3, x5
    umulh x6, x9, x21
    adc x4, x4, x6
    mul x5, x8, x22
    adds x3, x3, x5
    umulh x6, x8, x22
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #12 of product */
    and x7, x1, x3, lsr #31    /* Get limb #13 of product */
    stp w6, w7, [x0, #48]      /* Save these limbs */
    umulh x6, x14, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #7: sum of products */
    mul x5, x14, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x13, x17
    adds x3, x3, x5
    umulh x6, x13, x17
    adc x4, x4, x6
    mul x5, x12, x19
    adds x3, x3, x5
    umulh x6, x12, x19
    adc x4, x4, x6
    mul x5, x11, x20
    adds x3, x3, x5
    umulh x6, x11, x20
    adc x4, x4, x6
    mul x5, x10, x21
    adds x3, x3, x5
    umulh x6, x10, x21
    adc x4, x4, x6
    mul x5, x9, x22
    adds x3, x3, x5
    umulh x6, x9, x22
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #14 of product */
    and x7, x1, x3, lsr #31    /* Get limb #15 of product */
    stp w6, w7, [x0, #56]      /* Save these limbs */
    umulh x6, x14, x17
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #8: sum of products */
    mul x5, x14, x17
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x13, x19
    adds x3, x3, x5
    umulh x6, x13, x19
    adc x4, x4, x6
    mul x5, x12, x20
    adds x3, x3, x5
    umulh x6, x12, x20
    adc x4, x4, x6
    mul x5, x11, x21
    adds x3, x3, x5
    umulh x6, x11, x21
    adc x4, x4, x6
    mul x5, x10, x22
    adds x3, x3, x5
    umulh x6, x10, x22
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #16 of product */
    and x7, x1, x3, lsr #31    /* Get limb #17 of product */
    stp w6, w7, [x0, #64]      /* Save these limbs */
    umulh x6, x14, x19
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #9: sum of products */
    mul x5, x14, x19
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x13, x20
    adds x3, x3, x5
    umulh x6, x13, x20
    adc x4, x4, x6
    mul x5, x12, x21
    adds x3, x3, x5
    umulh x6, x12, x21
    adc x4, x4, x6
    mul x5, x11, x22
    adds x3, x3, x5
    umulh x6, x11, x22
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #18 of product */
    and x7, x1, x3, lsr #31    /* Get limb #19 of product */
    stp w6, w7, [x0, #72]      /* Save these limbs */
    umulh x6, x14, x20
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #10: sum of products */
    mul x5, x14, x20
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x13, x21
    adds x3, x3, x5
    umulh x6, x13, x21
    adc x4, x4, x6
    mul x5, x12, x22
    adds x3, x3, x5
    umulh x6, x12, x22
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #20 of product */
    and x7, x1, x3, lsr #31    /* Get limb #21 of product */
    stp w6, w7, [x0, #80]      /* Save these limbs */
    umulh x6, x14, x21
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #11: sum of products */
    mul x5, x14, x21
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x13, x22
    adds x3, x3, x5
    umulh x6, x13, x22
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #22 of product */
    and x7, x1, x3, lsr #31    /* Get limb #23 of product */
    stp w6, w7, [x0, #88]      /* Save these limbs */
    umulh x6, x14, x22
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #12: sum of products */
    mul x5, x14, x22
    adds x3, x3, x5
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #24 of product */
    and x7, x1, x3, lsr #31    /* Get limb #25 of product */
    stp w6, w7, [x0, #96]      /* Save these limbs */
    tbnz w2, #0, EndMultiply7DoubleLimbs
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x6, x4, #2
    orr x3, x6, x3, lsr #62
    /* Store carry in product */
    and x6, x1, x3             /* Get limb #26 of product */
    and x7, x1, x3, lsr #31    /* Get limb #27 of product */
    stp w6, w7, [x0, #104]      /* Save these limbs */
EndMultiply7DoubleLimbs:
    ldp x21, x22, [sp], #16
    ldp x19, x20, [sp], #16
    ret
MoreThan14Limbs:
    stp x19, x20, [sp, #-16]!
    stp x21, x22, [sp, #-16]!
    stp x23, x24, [sp, #-16]!
    /* Load multiplicand and multiplier into registers */
    ldp w8, w3, [x0, #0]
    orr x8, x8, x3, lsl #31
    ldp w9, w3, [x0, #8]
    orr x9, x9, x3, lsl #31
    ldp w10, w3, [x0, #16]
    orr x10, x10, x3, lsl #31
    ldp w11, w3, [x0, #24]
    orr x11, x11, x3, lsl #31
    ldp w12, w3, [x0, #32]
    orr x12, x12, x3, lsl #31
    ldp w13, w3, [x0, #40]
    orr x13, x13, x3, lsl #31
    ldp w14, w3, [x0, #48]
    orr x14, x14, x3, lsl #31
    tbnz w2, #0, LengthIs15Limbs
    ldp w15, w3, [x0, #56]
    orr x15, x15, x3, lsl #31
    ldp w16, w3, [x0, #64]
    orr x16, x16, x3, lsl #31
    ldp w17, w3, [x0, #72]
    orr x17, x17, x3, lsl #31
    ldp w19, w3, [x0, #80]
    orr x19, x19, x3, lsl #31
    ldp w20, w3, [x0, #88]
    orr x20, x20, x3, lsl #31
    ldp w21, w3, [x0, #96]
    orr x21, x21, x3, lsl #31
    ldp w22, w3, [x0, #104]
    orr x22, x22, x3, lsl #31
    ldp w23, w3, [x0, #112]
    orr x23, x23, x3, lsl #31
    ldp w24, w3, [x0, #120]
    orr x24, x24, x3, lsl #31
    b Multiply8DoubleLimbs
LengthIs15Limbs:
    ldr w15, [x0, #56]
    ldp w16, w3, [x0, #60]
    orr x16, x16, x3, lsl #31
    ldp w17, w3, [x0, #68]
    orr x17, x17, x3, lsl #31
    ldp w19, w3, [x0, #76]
    orr x19, x19, x3, lsl #31
    ldp w20, w3, [x0, #84]
    orr x20, x20, x3, lsl #31
    ldp w21, w3, [x0, #92]
    orr x21, x21, x3, lsl #31
    ldp w22, w3, [x0, #100]
    orr x22, x22, x3, lsl #31
    ldp w23, w3, [x0, #108]
    orr x23, x23, x3, lsl #31
    ldr w24, [x0, #116]
Multiply8DoubleLimbs:
    /* Perform first multiplication */
    mul x3, x8, x16
    umulh x4, x8, x16
    and x6, x1, x3               /* Get limb #0 of product */
    and x7, x1, x3, lsr #31      /* Get limb #1 of product */
    stp w6, w7, [x0]             /* Save these limbs */
    umulh x6, x9, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #1: sum of products */
    mul x5, x9, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x8, x17
    adds x3, x3, x5
    umulh x6, x8, x17
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #2 of product */
    and x7, x1, x3, lsr #31    /* Get limb #3 of product */
    stp w6, w7, [x0, #8]      /* Save these limbs */
    umulh x6, x10, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #2: sum of products */
    mul x5, x10, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x9, x17
    adds x3, x3, x5
    umulh x6, x9, x17
    adc x4, x4, x6
    mul x5, x8, x19
    adds x3, x3, x5
    umulh x6, x8, x19
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #4 of product */
    and x7, x1, x3, lsr #31    /* Get limb #5 of product */
    stp w6, w7, [x0, #16]      /* Save these limbs */
    umulh x6, x11, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #3: sum of products */
    mul x5, x11, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x10, x17
    adds x3, x3, x5
    umulh x6, x10, x17
    adc x4, x4, x6
    mul x5, x9, x19
    adds x3, x3, x5
    umulh x6, x9, x19
    adc x4, x4, x6
    mul x5, x8, x20
    adds x3, x3, x5
    umulh x6, x8, x20
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #6 of product */
    and x7, x1, x3, lsr #31    /* Get limb #7 of product */
    stp w6, w7, [x0, #24]      /* Save these limbs */
    umulh x6, x12, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #4: sum of products */
    mul x5, x12, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x11, x17
    adds x3, x3, x5
    umulh x6, x11, x17
    adc x4, x4, x6
    mul x5, x10, x19
    adds x3, x3, x5
    umulh x6, x10, x19
    adc x4, x4, x6
    mul x5, x9, x20
    adds x3, x3, x5
    umulh x6, x9, x20
    adc x4, x4, x6
    mul x5, x8, x21
    adds x3, x3, x5
    umulh x6, x8, x21
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #8 of product */
    and x7, x1, x3, lsr #31    /* Get limb #9 of product */
    stp w6, w7, [x0, #32]      /* Save these limbs */
    umulh x6, x13, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #5: sum of products */
    mul x5, x13, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x12, x17
    adds x3, x3, x5
    umulh x6, x12, x17
    adc x4, x4, x6
    mul x5, x11, x19
    adds x3, x3, x5
    umulh x6, x11, x19
    adc x4, x4, x6
    mul x5, x10, x20
    adds x3, x3, x5
    umulh x6, x10, x20
    adc x4, x4, x6
    mul x5, x9, x21
    adds x3, x3, x5
    umulh x6, x9, x21
    adc x4, x4, x6
    mul x5, x8, x22
    adds x3, x3, x5
    umulh x6, x8, x22
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #10 of product */
    and x7, x1, x3, lsr #31    /* Get limb #11 of product */
    stp w6, w7, [x0, #40]      /* Save these limbs */
    umulh x6, x14, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #6: sum of products */
    mul x5, x14, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x13, x17
    adds x3, x3, x5
    umulh x6, x13, x17
    adc x4, x4, x6
    mul x5, x12, x19
    adds x3, x3, x5
    umulh x6, x12, x19
    adc x4, x4, x6
    mul x5, x11, x20
    adds x3, x3, x5
    umulh x6, x11, x20
    adc x4, x4, x6
    mul x5, x10, x21
    adds x3, x3, x5
    umulh x6, x10, x21
    adc x4, x4, x6
    mul x5, x9, x22
    adds x3, x3, x5
    umulh x6, x9, x22
    adc x4, x4, x6
    mul x5, x8, x23
    adds x3, x3, x5
    umulh x6, x8, x23
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #12 of product */
    and x7, x1, x3, lsr #31    /* Get limb #13 of product */
    stp w6, w7, [x0, #48]      /* Save these limbs */
    umulh x6, x15, x16
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #7: sum of products */
    mul x5, x15, x16
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x14, x17
    adds x3, x3, x5
    umulh x6, x14, x17
    adc x4, x4, x6
    mul x5, x13, x19
    adds x3, x3, x5
    umulh x6, x13, x19
    adc x4, x4, x6
    mul x5, x12, x20
    adds x3, x3, x5
    umulh x6, x12, x20
    adc x4, x4, x6
    mul x5, x11, x21
    adds x3, x3, x5
    umulh x6, x11, x21
    adc x4, x4, x6
    mul x5, x10, x22
    adds x3, x3, x5
    umulh x6, x10, x22
    adc x4, x4, x6
    mul x5, x9, x23
    adds x3, x3, x5
    umulh x6, x9, x23
    adc x4, x4, x6
    mul x5, x8, x24
    adds x3, x3, x5
    umulh x6, x8, x24
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #14 of product */
    and x7, x1, x3, lsr #31    /* Get limb #15 of product */
    stp w6, w7, [x0, #56]      /* Save these limbs */
    umulh x6, x15, x17
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #8: sum of products */
    mul x5, x15, x17
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x14, x19
    adds x3, x3, x5
    umulh x6, x14, x19
    adc x4, x4, x6
    mul x5, x13, x20
    adds x3, x3, x5
    umulh x6, x13, x20
    adc x4, x4, x6
    mul x5, x12, x21
    adds x3, x3, x5
    umulh x6, x12, x21
    adc x4, x4, x6
    mul x5, x11, x22
    adds x3, x3, x5
    umulh x6, x11, x22
    adc x4, x4, x6
    mul x5, x10, x23
    adds x3, x3, x5
    umulh x6, x10, x23
    adc x4, x4, x6
    mul x5, x9, x24
    adds x3, x3, x5
    umulh x6, x9, x24
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #16 of product */
    and x7, x1, x3, lsr #31    /* Get limb #17 of product */
    stp w6, w7, [x0, #64]      /* Save these limbs */
    umulh x6, x15, x19
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #9: sum of products */
    mul x5, x15, x19
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x14, x20
    adds x3, x3, x5
    umulh x6, x14, x20
    adc x4, x4, x6
    mul x5, x13, x21
    adds x3, x3, x5
    umulh x6, x13, x21
    adc x4, x4, x6
    mul x5, x12, x22
    adds x3, x3, x5
    umulh x6, x12, x22
    adc x4, x4, x6
    mul x5, x11, x23
    adds x3, x3, x5
    umulh x6, x11, x23
    adc x4, x4, x6
    mul x5, x10, x24
    adds x3, x3, x5
    umulh x6, x10, x24
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #18 of product */
    and x7, x1, x3, lsr #31    /* Get limb #19 of product */
    stp w6, w7, [x0, #72]      /* Save these limbs */
    umulh x6, x15, x20
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #10: sum of products */
    mul x5, x15, x20
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x14, x21
    adds x3, x3, x5
    umulh x6, x14, x21
    adc x4, x4, x6
    mul x5, x13, x22
    adds x3, x3, x5
    umulh x6, x13, x22
    adc x4, x4, x6
    mul x5, x12, x23
    adds x3, x3, x5
    umulh x6, x12, x23
    adc x4, x4, x6
    mul x5, x11, x24
    adds x3, x3, x5
    umulh x6, x11, x24
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #20 of product */
    and x7, x1, x3, lsr #31    /* Get limb #21 of product */
    stp w6, w7, [x0, #80]      /* Save these limbs */
    umulh x6, x15, x21
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #11: sum of products */
    mul x5, x15, x21
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x14, x22
    adds x3, x3, x5
    umulh x6, x14, x22
    adc x4, x4, x6
    mul x5, x13, x23
    adds x3, x3, x5
    umulh x6, x13, x23
    adc x4, x4, x6
    mul x5, x12, x24
    adds x3, x3, x5
    umulh x6, x12, x24
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #22 of product */
    and x7, x1, x3, lsr #31    /* Get limb #23 of product */
    stp w6, w7, [x0, #88]      /* Save these limbs */
    umulh x6, x15, x22
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #12: sum of products */
    mul x5, x15, x22
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x14, x23
    adds x3, x3, x5
    umulh x6, x14, x23
    adc x4, x4, x6
    mul x5, x13, x24
    adds x3, x3, x5
    umulh x6, x13, x24
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #24 of product */
    and x7, x1, x3, lsr #31    /* Get limb #25 of product */
    stp w6, w7, [x0, #96]      /* Save these limbs */
    umulh x6, x15, x23
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #13: sum of products */
    mul x5, x15, x23
    adds x3, x3, x5
    adc x4, x4, x6
    mul x5, x14, x24
    adds x3, x3, x5
    umulh x6, x14, x24
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #26 of product */
    and x7, x1, x3, lsr #31    /* Get limb #27 of product */
    stp w6, w7, [x0, #104]      /* Save these limbs */
    umulh x6, x15, x24
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x7, x4, #2
    lsr x4, x4, #62
    orr x3, x7, x3, lsr #62
    /* Process the column #14: sum of products */
    mul x5, x15, x24
    adds x3, x3, x5
    adc x4, x4, x6
    and x6, x1, x3             /* Get limb #28 of product */
    and x7, x1, x3, lsr #31    /* Get limb #29 of product */
    stp w6, w7, [x0, #112]      /* Save these limbs */
    tbnz w2, #0, EndMultiply8DoubleLimbs
    /* Shift right 62 bits the sum of products x4: x3 */
    lsl x6, x4, #2
    orr x3, x6, x3, lsr #62
    /* Store carry in product */
    and x6, x1, x3             /* Get limb #30 of product */
    and x7, x1, x3, lsr #31    /* Get limb #31 of product */
    stp w6, w7, [x0, #120]      /* Save these limbs */
EndMultiply8DoubleLimbs:
    ldp x23, x24, [sp], #16
    ldp x21, x22, [sp], #16
    ldp x19, x20, [sp], #16
    ret
