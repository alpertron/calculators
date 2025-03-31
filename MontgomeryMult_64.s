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
    .global MontgomeryMult
    .type MontgomeryMult, %function
    .text
    .align 4
/*
 On input:
 x0: Pointer to factor1 in Montgomery notation.
 x1: Pointer to factor2 in Montgomery notation.
 x2: Pointer to product in Montgomery notation.
 Other registers used in the code:
 w3: Value of MontgomeryMultN.
 x4: Pointer to TestNbr.
 w5: Value of NumberLength.
 w6: Index for multiplicand (0 to NumberLength - 1).
 w7: Index for inner loop (0 to NumberLength - 1).
 w9: Limb from multiplicand (index w6).
 x10: Accumulator for 32x32 bit multiplications.
 Multiply two numbers in Montgomery notation.

 For large numbers the REDC algorithm is:
 m <- ((T mod R)N') mod R
 t <- (T - mN) / R
 if t < 0 then
   return t + N
 else
   return t
 end if

 Multilimb version when multiplying S*T:
 1) P <- 0
 2) for i from 0 to n-1
 3)   t <- (S[i] * T[0]) mod R
 4)   u <- (t * N') mod R (one limb)
 5)   P <- (P + S[i] * T - u * N) / R
 6) if P < 0 then
 7)   return P - N
 8) else
 9)  return P
 10)end if
*/
    .set MAX_LIMBS_MONTGOMERY, 14
    .set STACK_SIZE_RESULT_BUFFER, (MAX_LIMBS_MONTGOMERY+3)/4 * 16
MontgomeryMult:
    /* Get external values and pointers to buffers. */
    adrp x3, :got:MontgomeryMultN
    ldr x3, [x3, #:got_lo12:MontgomeryMultN]
    ldr w3, [x3]               /* x3 <- MontgomeryMultN. */
    adrp x4, :got:TestNbr
    ldr x4, [x4, #:got_lo12:TestNbr]  /* x4 <- Pointer to TestNbr. */
    adrp x5, :got:NumberLength@GOTPAGE
    ldr x5, [x5, #:got_lo12:NumberLength]
    ldr w5, [x5]               /* x5 <- NumberLength */

    /* Adjust stack to store temporary result there */
    /* so the arguments are not overwritten with the result. */
    sub sp, sp, #STACK_SIZE_RESULT_BUFFER

    /* Step 1: Clear result. */
    mov w7, w5                 /* Clear NumberLength limbs. */
    mov w6, #0                 /* Set all limbs to zero. */
clear_result_loop:
    subs w7, w7, #1
    str w6, [sp, w7, uxtw #2]  /* Set limb of temporary result to zero. */
    bne clear_result_loop

    /* Step 2: Outer loop. */
    /* Index for multiplicand (w6) is already zero. */
outer_loop:

    /* Step 5: */
    mov w7, #1                 /* Initialize index for result. */
    ldr w9, [x0, w6, uxtw #2]  /* Get limb from multiplicand. */
    ldr w10, [sp]              /* Get the least significant limb of temp result. */
    ldr w12, [x1]              /* Get the least significant limb of multiplier. */
    smaddl x10, w9, w12, x10   /* x10 <- Pr + multiplicand[i] * multiplier[0]. */
    mul w8, w10, w3            /* w8 <- MontDig = Pr * MontgomeryMultN mod R */
    bic w8, w8, #0x80000000    /* Clear bit 31. */
    neg w8, w8                 /* r8 <- -MontDig */
    ldr w12, [x4]              /* Get least significant of TestNbr. */
    smaddl x10, w8, w12, x10   /* x10 = Pr <- Pr - MontDig * TestNbr[0] */
inner_loop:
    asr x10, x10, #31          /* Shift right Pr (x10) by 31 bits. */
    ldr w12, [x1, w7, uxtw #2] /* Get limb from multiplier */
    smaddl x10, w9, w12, x10   /* Pr <- Pr + multiplicand[i] * multiplier[j+1]. */
    ldr w12, [x4, w7, uxtw #2] /* Get limb from TestNbr */
    smaddl x10, w8, w12, x10   /* Pr <- Pr - MontDig * TestNbr[j+1]. */
    ldrsw x12, [sp, w7, uxtw #2] /* Get limb from temp result (sign extended). */
    add x10, x10, x12          /* Pr <- Pr + result[j+1]. */
    sub w7, w7, #1             /* Decrement index. */
    bic w12, w10, #0x80000000  /* Clear bit 31. */
    str w12, [sp, w7, uxtw #2] /* Store limb to temporary result[j]. */
    add w7, w7, #2
    cmp w7, w5
    bne inner_loop
    sub w7, w7, #1
    asr x10, x10, #31          /* Shift right Pr (x10) by 31 bits. */
    str w10, [sp, w7, uxtw #2] /* Store most significant limb of temp result. */
    add w6, w6, #1             /* Increment index to multiplicand. */
    cmp w6, w5
    bne outer_loop
    
    /* Step 6: Test whether the result is positive. */
    tst w10, #0x80000000       /* Is most significant limb of result positive? */
    beq copy_temp_to_result    /* If so, copy temporary result to actual result. */
    
    /* Step 7: Add TestNbr to result. */
    mov w6, #0                 /* Initialize index. */
    mov w7, #0                 /* Initialize result of previous addition. */
add_TestNbr_cycle:
    ldr w8, [sp, w6, uxtw #2]  /* Get limb from temporary result. */
    ldr w9, [x4, w6, uxtw #2]  /* Get limb from TestNbr. */
    add w7, w8, w7, lsr #31    /* Add with carry. */
    add w7, w7, w9
    bic w8, w7, #0x80000000    /* Clear bit 31. */
    str w8, [sp, w6, uxtw #2]  /* Store limb into temporary result. */
    add w6, w6, #1
    cmp w6, w5
    bne add_TestNbr_cycle
copy_temp_to_result:
copy_temp_to_result_loop:
    subs w5, w5, #1
    ldr w8, [sp, w5, uxtw #2]  /* Get limb from temporary result. */
    str w8, [x2, w5, uxtw #2]  /* Set limb to result. */
    bne copy_temp_to_result_loop
    add sp, sp, #STACK_SIZE_RESULT_BUFFER
    ret
    .end
