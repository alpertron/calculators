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

    .syntax unified
    .arch armv7-a
    .thumb
    .global ClassicalMult
    .thumb_func
    .type ClassicalMult, %function
    .text
    .align 4
/*
 Multiply two groups of nbrLen limbs. The first one starts at idxFactor1
 and the second one at idxFactor2. The 2*nbrLen limb result is stored
 starting at idxFactor1. Use arrayAux as temporary storage.
 Accumulate products by result limb.

 On input:
 r0: Index to arr of multiplicand.
 r1: Index to arr of multiplier.
 r2: Number of limbs.
 Other registers used in the code:
 r3: Pointer to arr.
 r4: Pointer to arrAux.

                                 A  B | C  D | E  F
                                 H  I | J  K | L  M
                                 ======+======+=====
                                AM BM |CM DM |EM FM
                            AL |BL CL |DL EL |FL
                         ------+------+-----
                         AK BK |CK DK |EK FK
                     AJ |BJ CJ |DJ EJ |FJ
                  ------+------+------
                  AI BI |CI DI |EI FI
              AH |BH CH |DH EH |FH

  Before the outer loop, the first two lines are computed in groups of two limbs.
  The outer loop processes two lines at a time, in groups of two limbs.
  */
    .macro START_MULT_EVEN_LIMBS
    ldr r6, [r1], #4            /* Get multiplicand limbs #0 and #1 and incr ptr */
    ldr r7, [r1], #4            /* r6=M, r7=L */
    ldr r8, [r0]                /* Get multiplier limbs #0 and #1 (r8=F, r9=E) */
    ldr r9, [r0, #4]
    umull r10, r14, r8, r6      /* Compute product F*M */
    bic r12, r10, #0x80000000
    str r12, [r4]               /* Store limb #0 of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r9, r6      /* Add product of limbs of both factors (E*M) */
    umlal r10, r14, r8, r7      /* Add product of limbs of both factors (F*L) */
    bic r12, r10, #0x80000000
    str r12, [r4, #4]           /* Store limb #1 of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    .endm

    .macro MULT_FIRST_TWO_LINES, offset
    umlal r10, r14, r9, r7      /* Add product of limbs of E*L */
    ldr r8, [r0, #\offset]      /* Get limb of multiplier (r8=D) */
    ldr r9, [r0, #4+\offset]    /* Get limb of multiplier (r9=C) */
    umlal r10, r14, r8, r6      /* Add product of limbs of D*M */
    bic r12, r10, #0x80000000
    str r12, [r4, #\offset]     /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r9, r6      /* Add product of limbs of C*M */
    umlal r10, r14, r8, r7      /* Add product of limbs of D*L */
    bic r12, r10, #0x80000000
    str r12, [r4, #4+\offset]   /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    .endm

    .macro END_MULT_FIRST_TWO_LINES, offset
    umlal r10, r14, r9, r7      /* Add product of limbs of A*L */
    bic r12, r10, #0x80000000
    str r12, [r4, #\offset]
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    orr r10, r11, r10, lsr #31
    str r10, [r4, #4+\offset]
    .endm

    .macro START_MULT_OTHER_LINES
    ldr r6, [r1], #4            /* Get limb of multiplicand (r6=K) and incr ptr */
    ldr r7, [r1], #4            /* Get limb of multiplicand (r7=J) and incr ptr */
    ldr r8, [r0]                /* Get multiplier limb #0 (r8=F) */
    ldr r9, [r0, #4]            /* Get multiplier limb #1 (r9=E) */
    add r3, r4, r5, lsl #2      /* Point to limb #2 of product */
    mov r14, #0
    ldr r10, [r3]               /* Get limb from product */
    ldr r12, [r3, #4]           /* Get limb from product */
    umlal r10, r14, r8, r6      /* Add product of F*J */
    bic r11, r10, #0x80000000   /* Get least significant limb of product */
    str r11, [r3]               /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    adds r10, r10, r12          /* Sum limb of product */
    adc r14, r14, #0
    umlal r10, r14, r9, r6      /* Add product of limbs of E*K */
    umlal r10, r14, r8, r7      /* Add product of limbs of F*J */
    bic r12, r10, #0x80000000   /* Get most significant limb of product */
    str r12, [r3, #4]           /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    .endm

    .macro MULT_OTHER_LINES, offset_multiplier, offset_product
    umlal r10, r14, r9, r7      /* Add product of limbs of E*J */
    ldr r8, [r0, #\offset_multiplier]   /* Get limb of multiplier (r8=D) */
    ldr r9, [r0, #4+\offset_multiplier] /* Get limb of multiplier (r9=C) */
    umlal r10, r14, r8, r6      /* Add product of limbs of D*K */
    ldr r11, [r3, #\offset_product]  /* Get limb of product */
    ldr r12, [r3, #4+\offset_product]/* Get limb of product */
    adds r10, r10, r11          /* Add least significant limb of product */
    adc r14, r14, #0
    bic r11, r10, #0x80000000
    str r11, [r3, #\offset_product]   /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    adds r10, r10, r12          /* Add limb of product */
    adc r14, r14, #0
    umlal r10, r14, r9, r6      /* Add product of limbs of C*K */
    umlal r10, r14, r8, r7      /* Add product of limbs of D*J */
    bic r12, r10, #0x80000000
    str r12, [r3, #4+\offset_product] /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    .endm

    .macro END_MULT_OTHER_LINES_EVEN_LIMBS, offset
    umlal r10, r14, r9, r7      /* Add product of limbs of A*L */
    bic r12, r10, #0x80000000
    str r12, [r3, #\offset]
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    orr r10, r11, r10, lsr #31
    str r10, [r3, #4+\offset]
    .endm

  /*
                                 A  B | C  D | E  F  G
                                 H  I | J  K | L  M  N
                                =======+======+=======
                                AN BN |CN DN |EN FN GN
                            AM |BM CM |DM EM |FM GM
                         AL BL |CL DL |EL FL |GL
                     ---+------+------+------
                     AK |BK CK |DK EK |FK GK
                  AJ BJ |CJ DJ |EJ FJ |GJ
                  ------+------+------
              AI |BI CI |DI EI |FI GI
           AH BH |CH DH |EH FH |GH

  Before the outer loop, the first three lines are computed in groups of two limbs.
  The outer loop processes two lines at a time, in groups of two limbs.

 */

    .macro START_MULT_ODD_LIMBS
    ldr r5, [r1], #4            /* Get multiplicand limb r5=N */
    ldr r6, [r1], #4            /* Get multiplicand limb r6=M */
    ldr r7, [r1], #4            /* Get multiplicand limb r7=L */
    ldr r3, [r0]                /* Get multiplier limb #0 r3=G */
    ldr r8, [r0, #4]            /* Get multiplier limb #1 r8=F */
    ldr r9, [r0, #8]            /* Get multiplier limb #2 r9=E */
    umull r10, r14, r3, r5      /* Compute product G*N */
    bic r11, r10, #0x80000000
    str r11, [r4]               /* Store product limb #0 */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r8, r5      /* Add product of limbs of both factors (F*N) */
    umlal r10, r14, r3, r6      /* Add product of limbs of both factors (G*M) */
    bic r11, r10, #0x80000000
    str r11, [r4, #4]           /* Store limb #1 of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r9, r5      /* Add product of limbs of both factors (E*N) */
    umlal r10, r14, r8, r6      /* Add product of limbs of both factors (F*M) */
    umlal r10, r14, r3, r7      /* Add product of limbs of both factors (G*L) */
    bic r11, r10, #0x80000000
    str r11, [r4, #8]           /* Store limb #2 of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    .endm

    .macro MULT_FIRST_THREE_LINES, offset
    mov r12, r9                 /* Get most significant limb of previous pair of multiplier (E) */
    umlal r10, r14, r8, r7      /* Add product of limbs of F*L */
    ldr r8, [r0, #\offset]      /* Get limb of multiplier (r8=D) */
    ldr r9, [r0, #4+\offset]    /* Get limb of multiplier (r9=C) */
    umlal r10, r14, r8, r5      /* Add product of limbs of D*N */
    umlal r10, r14, r12, r6     /* Add product of limbs of E*M */
    bic r11, r10, #0x80000000
    str r11, [r4, #\offset]     /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r9, r5      /* Add product of limbs of C*N */
    umlal r10, r14, r8, r6      /* Add product of limbs of D*M */
    umlal r10, r14, r12, r7     /* Add product of limbs of E*L */
    bic r12, r10, #0x80000000
    str r12, [r4, #4+\offset]   /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    .endm

    .macro END_MULT_FIRST_THREE_LINES, offset
    umlal r10, r14, r9, r6      /* Add product of limbs of A*M */
    umlal r10, r14, r8, r7      /* Add product of limbs of B*L */
    bic r11, r10, #0x80000000
    str r11, [r4, #\offset]     /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r9, r7      /* Add product of limbs of A*L */
    bic r11, r10, #0x80000000
    str r11, [r4, #4+\offset]
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    orr r10, r11, r10, lsr #31
    str r10, [r4, #8+\offset]
    .endm

    .macro START_MULT_OTHER_LINES_ODD_LIMBS
    ldr r6, [r1], #4            /* Get limb of multiplicand (r6=K) and incr ptr */
    ldr r7, [r1], #4            /* Get limb of multiplicand (r7=J) and incr ptr */
    ldr r8, [r0]                /* Get multiplier limb #0 (r8=G) */
    ldr r9, [r0, #4]            /* Get multiplier limb #1 (r9=F) */
    add r3, r4, r5, lsl #2      /* Point to limb #2 of product */
    mov r14, #0
    ldr r10, [r3]               /* Get limbs from product */
    ldr r12, [r3, #4]           /* Get limbs from product */
    umlal r10, r14, r8, r6      /* Add product of G*K */
    bic r11, r10, #0x80000000
    str r11, [r3], #4           /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    adds r10, r10, r12          /* Add most significant limb of product to carry */
    adc r14, r14, #0
    umlal r10, r14, r9, r6      /* Add product of limbs of F*K */
    umlal r10, r14, r8, r7      /* Add product of limbs of G*J */
    bic r12, r10, #0x80000000
    str r12, [r3], #4           /* Store limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    .endm

    .macro END_MULT_OTHER_LINES_ODD_LIMBS, offset_multiplier, offset_product
    ldr r12, [r0, #\offset_multiplier]/* Get most significant limb of multiplier (A) */
    umlal r10, r14, r12, r6     /* Add product of limbs of A*K */
    umlal r10, r14, r9, r7      /* Add product of limbs of B*J */
    ldr r11, [r3, #\offset_product] /* Add product limb */
    adds r10, r10, r11
    adc r14, r14, #0
    bic r11, r10, #0x80000000
    str r11, [r3, #\offset_product] /* Save limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r12, r7     /* Add product of limbs of A*J */
    bic r11, r10, #0x80000000
    str r11, [r3, #4+\offset_product] /* Save limb of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    orr r10, r11, r10, lsr #31
    str r10, [r3, #8+\offset_product] /* Save limb of product */
    .endm

    .macro MULT_EVEN_LIMBS, nbrLimbs
      START_MULT_EVEN_LIMBS
      .set offset, 8
      .rept (\nbrLimbs-2)/2
        MULT_FIRST_TWO_LINES offset
        .set offset, offset + 8
      .endr
      END_MULT_FIRST_TWO_LINES offset
      /* At this moment, the first two lines are complete */
      mov r5, #2                  /* Index to multiplicand */
1:
      START_MULT_OTHER_LINES
      .set offset, 8
      .rept (\nbrLimbs-2)/2
        MULT_OTHER_LINES offset, offset
        .set offset, offset + 8
      .endr
      END_MULT_OTHER_LINES_EVEN_LIMBS offset
    add r5, r5, #2
    cmp r5, r2
    bne 1b
    .endm

    .macro MULT_ODD_LIMBS, nbrLimbs
      START_MULT_ODD_LIMBS
      .set offset, 12
      .rept (\nbrLimbs-3)/2
        MULT_FIRST_THREE_LINES offset
        .set offset, offset + 8
      .endr
      END_MULT_FIRST_THREE_LINES offset
      /* At this moment, the first three lines are complete */
      mov r5, #3                  /* Index to multiplicand */
1:
      START_MULT_OTHER_LINES
      .set offset_multiplier, 8
      .set offset_product, 8
      .rept (\nbrLimbs-3)/2
        MULT_OTHER_LINES offset_multiplier, offset_product
        .set offset_multiplier, offset_multiplier + 8
        .set offset_product, offset_product + 8
      .endr
      END_MULT_OTHER_LINES_ODD_LIMBS offset_multiplier, offset_product
    add r5, r5, #2
    cmp r5, r2
    bne 1b
    .endm

ClassicalMult:
    push {r4-r11, lr}
    /* Get external values and pointers to buffers. */
    ldr r7, locGOT
.LPIC0:
    add r7, pc, r7
    ldr r3, arr_offset
    ldr r3, [r3, r7]                 /* r3 = Address of arr */
    ldr r4, arrAux_offset
    ldr r4, [r4, r7]                 /* r4 = Address of arrAux */
    add r0, r3, r0, lsl #2           /* r0 = Pointer to multiplier. */
    add r1, r3, r1, lsl #2           /* r1 = Pointer to multiplicand. */
    cmp r2, #1
    bne Test2Limbs

ClassicalMult1Limb:
    ldr r5, [r0]                 /* r5 = Multiplicand */
    ldr r6, [r0, #4]             /* r6 = Multiplier */
    umull r10, r14, r5, r6       /* Compute product */
    bic r9, r10, #0x80000000
    str r9, [r0]                 /* Save limb #0 of product */
    lsl r11, r14, #1             /* Compute carry: R14: R10 LSR 31*/
    orr r10, r11, r10, lsr #31
    str r10, [r0, #4]            /* Save limb #1 of product */
    pop {r4-r11, pc}
locGOT:
    .word _GLOBAL_OFFSET_TABLE_-(.LPIC0+4)
arr_offset:
    .word arr(GOT)
arrAux_offset:
    .word arrAux(GOT)

Test2Limbs:
    cmp r2, #2
    bne Test3Limbs

ClassicalMult2Limbs:
    ldr r5, [r0]
    ldr r6, [r0, #4]             /* r6: r5 = Multiplicand */
    ldr r7, [r0, #8]
    ldr r8, [r0, #12]            /* r8: r7 = Multiplier */
    umull r10, r14, r5, r7
    bic r9, r10, #0x80000000
    str r9, [r0]                 /* Save limb #0 of product */
    lsl r11, r14, #1             /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r5, r8
    umlal r10, r14, r6, r7
    bic r9, r10, #0x80000000
    str r9, [r0, #4]             /* Save limb #1 of product */
    lsl r11, r14, #1             /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r6, r8
    bic r9, r10, #0x80000000
    str r9, [r0, #8]             /* Save limb #2 of product */
    lsl r11, r14, #1             /* Compute carry: R14: R10 LSR 31*/
    orr r10, r11, r10, lsr #31
    str r10, [r0, #12]           /* Save limb #3 of product */
    pop {r4-r11, pc}

Test3Limbs:
    cmp r2, #3
    bne MoreThan3Limbs

ClassicalMult3Limbs:
    ldr r5, [r0]
    ldr r6, [r0, #4]
    ldr r7, [r0, #8]            /* r7: r6: r5 = Multiplicand */
    ldr r8, [r0, #12]
    ldr r9, [r0, #16]
    ldr r12, [r0, #20]          /* r12: r9: r8 = Multiplier */
    umull r10, r14, r5, r8
    bic r11, r10, #0x80000000
    str r11, [r0]               /* Save limb #0 of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r5, r9
    umlal r10, r14, r6, r8
    bic r11, r10, #0x80000000
    str r11, [r0, #4]           /* Save limb #1 of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r12, r5
    umlal r10, r14, r6, r9
    umlal r10, r14, r7, r8
    bic r11, r10, #0x80000000
    str r11, [r0, #8]           /* Save limb #2 of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r6, r12
    umlal r10, r14, r7, r9
    bic r11, r10, #0x80000000
    str r11, [r0, #12]          /* Save limb #3 of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    lsr r14, r14, #31
    orr r10, r11, r10, lsr #31
    umlal r10, r14, r7, r12
    bic r9, r10, #0x80000000
    str r9, [r0, #16]           /* Save limb #4 of product */
    lsl r11, r14, #1            /* Compute carry: R14: R10 LSR 31*/
    orr r10, r11, r10, lsr #31
    str r10, [r0, #20]          /* Save limb #5 of product */
    pop {r4-r11, pc}

MoreThan3Limbs:
    cmp r2, #4
    bne MoreThan4Limbs
    MULT_EVEN_LIMBS 4
    b copy_temp_product_to_product

MoreThan4Limbs:
    cmp r2, #5
    bne MoreThan5Limbs
    MULT_ODD_LIMBS 5
    b copy_temp_product_to_product

MoreThan5Limbs:
    cmp r2, #6
    bne MoreThan6Limbs
    MULT_EVEN_LIMBS 6
    b copy_temp_product_to_product

MoreThan6Limbs:
    cmp r2, #7
    bne MoreThan7Limbs
    MULT_ODD_LIMBS 7
    b copy_temp_product_to_product

MoreThan7Limbs:
    cmp r2, #8
    bne MoreThan8Limbs
    MULT_EVEN_LIMBS 8
    b copy_temp_product_to_product

MoreThan8Limbs:
    cmp r2, #9
    bne MoreThan9Limbs
    MULT_ODD_LIMBS 9
    b copy_temp_product_to_product

MoreThan9Limbs:
    cmp r2, #10
    bne MoreThan10Limbs
    MULT_EVEN_LIMBS 10
    b copy_temp_product_to_product

MoreThan10Limbs:
    cmp r2, #11
    bne MoreThan11Limbs
    MULT_ODD_LIMBS 11
    b copy_temp_product_to_product

MoreThan11Limbs:
    cmp r2, #12
    bne MoreThan12Limbs
    MULT_EVEN_LIMBS 12
    b copy_temp_product_to_product

MoreThan12Limbs:
    cmp r2, #13
    bne MoreThan13Limbs
    MULT_ODD_LIMBS 13
    b copy_temp_product_to_product

MoreThan13Limbs:
    cmp r2, #14
    bne MoreThan14Limbs
    MULT_EVEN_LIMBS 14
    b copy_temp_product_to_product

MoreThan14Limbs:
    cmp r2, #15
    bne MoreThan15Limbs
    MULT_ODD_LIMBS 15
    b copy_temp_product_to_product

MoreThan15Limbs:
    MULT_EVEN_LIMBS 16

copy_temp_product_to_product:
    mov r6, r2
copy_product_loop:
    ldr r8, [r4], #4        /* Get limb of temporary product */
    ldr r9, [r4], #4        /* Get limb of temporary product */
    str r8, [r0], #4        /* Set limb of product */
    str r9, [r0], #4        /* Set limb of product */
    subs r6, r6, #1
    bne copy_product_loop
    pop {r4-r11, pc}
    .end
