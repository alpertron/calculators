/*
 This file is part of Alpertron Calculators.

 Copyright 2025-2026 Dario Alejandro Alpern

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


/*
 MontgomeryMult for AArch64, limbs of 31 bits stored in 32-bit words.
 NumberLength must be in the range 2..16 (no range check is done).

 On input:
   x0: pointer to factor1 (Montgomery notation).
   x1: pointer to factor2 (Montgomery notation).
   x2: pointer to product (may alias x0 or x1).

 Design:
   * A jump table selects a version specialised for each NumberLength.
   * Each version has the inner loop fully unrolled (immediate offsets, no
     index registers). The outer loop (over the limbs of factor1) is a loop.
     Fully unrolling both loops would be O(n^2) code (~50 KB), too big for L1I.
   * Temporary result lives on the stack, except limb 0, which is kept in w11
     (avoids a store->load forwarding stall on the critical path).
   * Carry chain per limb is a single instruction:  add x10, x14, x10, asr #31
     The products a*b[j] - m*N[j] + res[j] do not depend on the carry, so the
     CPU can overlap them across iterations.
   * Final "add N if negative" is branch-free (mask) and fused with the copy.

 Registers:
   x0  factor1 pointer (post-incremented)     w3  MontgomeryMultN
   x1  factor2 pointer                        x4  TestNbr pointer
   w5  NumberLength                           w6  outer loop counter
   w8  MontDig                                w9  limb of factor1
   x10 accumulator (carry chain)              w11 temporary result limb 0
   x12,x13 loads of factor2/TestNbr limbs     x14 product accumulator
   w15 scratch / sign mask

 Temporary result: slot k (offset 4*k) holds limb k for k = 1..n-1.
*/
    .text
    .arch armv8-a
    .global MontgomeryMult
    .type MontgomeryMult, %function
    .align 4

    .set MAX_LIMBS_MONTGOMERY, 16
    .set STACK_SIZE_RESULT_BUFFER, MAX_LIMBS_MONTGOMERY * 4   /* 64 */

/* ------------------------------------------------------------------ */
/* Version for LIMBS limbs.                                             */
/* ------------------------------------------------------------------ */
    .macro MONTMUL LIMBS
.Lmm\LIMBS:
    mov   w11, #0                      /* Temporary result limb 0. */
    mov   w6, #\LIMBS                  /* Outer loop counter. */
.Louter\LIMBS:
    ldr   w9, [x0], #4                 /* a[i]. */
    ldr   w12, [x1]                    /* b[0]. */
    ldr   w13, [x4]                    /* N[0]. */
    smaddl x10, w9, w12, x11           /* Pr = res[0] + a[i]*b[0]. */
    mul   w8, w10, w3
    and   w8, w8, #0x7fffffff          /* MontDig. */
    smsubl x10, w8, w13, x10           /* Pr -= MontDig*N[0] (low 31 bits = 0). */
    .set  mm_j, 1
    .rept \LIMBS - 1
    ldr   w12, [x1, #(mm_j*4)]         /* b[j]. */
    ldr   w13, [x4, #(mm_j*4)]         /* N[j]. */
    .if mm_j == \LIMBS - 1
    ldrsw x14, [sp, #(mm_j*4)]         /* Top limb: signed. */
    .else
    ldr   w14, [sp, #(mm_j*4)]         /* Other limbs: 31 bits, non-negative. */
    .endif
    smaddl x14, w9, w12, x14           /* res[j] + a*b[j]. */
    smsubl x14, w8, w13, x14           /* - MontDig*N[j]. */
    add   x10, x14, x10, asr #31       /* + carry (only serial dependency). */
    .if mm_j == 1
    and   w11, w10, #0x7fffffff        /* New res[0] stays in a register. */
    .else
    and   w15, w10, #0x7fffffff
    str   w15, [sp, #((mm_j-1)*4)]     /* New res[j-1]. */
    .endif
    .set  mm_j, mm_j + 1
    .endr
    asr   x10, x10, #31
    str   w10, [sp, #((\LIMBS-1)*4)]   /* New top limb. */
    subs  w6, w6, #1
    b.ne  .Louter\LIMBS

    /* Branch-free correction: if top limb < 0 add N; fused with the copy. */
    asr   w15, w10, #31                /* Mask: 0 or 0xFFFFFFFF. */
    .set  mm_j, 0
    .rept \LIMBS
    .if mm_j == 0
    mov   w12, w11
    .else
    ldr   w12, [sp, #(mm_j*4)]
    .endif
    ldr   w13, [x4, #(mm_j*4)]
    and   w13, w13, w15                /* N or 0. */
    add   w12, w12, w13
    .if mm_j > 0
    add   w12, w12, w14                /* Carry from previous limb. */
    .endif
    .if mm_j < \LIMBS - 1
    lsr   w14, w12, #31                /* New carry. */
    .endif
    and   w12, w12, #0x7fffffff
    str   w12, [x2, #(mm_j*4)]
    .set  mm_j, mm_j + 1
    .endr
    add   sp, sp, #STACK_SIZE_RESULT_BUFFER
    ret
    .endm

MontgomeryMult:
    .cfi_startproc
    adrp  x3, :got:MontgomeryMultN
    ldr   x3, [x3, #:got_lo12:MontgomeryMultN]
    ldr   w3, [x3]                     /* w3 <- MontgomeryMultN. */
    adrp  x4, :got:TestNbr
    ldr   x4, [x4, #:got_lo12:TestNbr] /* x4 <- pointer to TestNbr. */
    adrp  x5, :got:NumberLength
    ldr   x5, [x5, #:got_lo12:NumberLength]
    ldr   w5, [x5]                     /* w5 <- NumberLength (2..16). */

    sub   w12, w5, #2
    adr   x13, mm_jump_table
    sub   sp, sp, #STACK_SIZE_RESULT_BUFFER
    add   x13, x13, x12, lsl #2        /* Each table entry is one instruction. */
    stp   xzr, xzr, [sp]               /* Clear temporary result (fixed size, */
    stp   xzr, xzr, [sp, #16]          /* no loop). */
    stp   xzr, xzr, [sp, #32]
    stp   xzr, xzr, [sp, #48]
    br    x13

mm_jump_table:
    .irp size, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16
    b     .Lmm\size
    .endr

    .irp size, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16
    MONTMUL \size
    .endr
    .cfi_endproc
    .size MontgomeryMult, .-MontgomeryMult
