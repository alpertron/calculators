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
    .global MontMultGraphic
    .type MontMultGraphic, %function
    .text
    .align 4
/*
 On input:
 x0: Pointer to multiplicand in Montgomery notation.
 x1: Pointer to multiplier in Montgomery notation.
 x2: Pointer to product in Montgomery notation.
 Other registers used in the code:
 w3: Value of MontMultGraphicN.
 w6:w1: Multiplier.
 w5:w4: TestNbr.
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
MontMultGraphic:
    /* Get external values and pointers to buffers. */
    adrp x3, :got:MontMultGraphicN
    ldr x3, [x3, #:got_lo12:MontMultGraphicN]
    ldr w3, [x3]               /* x3 <- MontgomeryMultN. */
    adrp x4, :got:TestNbrGraphic
    ldr x4, [x4, #:got_lo12:TestNbrGraphic]  /* x4 <- Pointer to TestNbr. */
    
    ldr w5, [x4, #4]
    ldr w4, [x4]               /* R5:R4 = TestNbr */
    cmp w5, #0
    bne two_limbs
    ldr w8, [x0]               /* Get multiplicand */
    ldr w9, [x1]               /* Get multiplier */
    umull x7, w8, w9           /* x7 <- multiplicand * multiplier */
    udiv x6, x7, x4            /* x6 <- multiplicand * multiplier / TestNbr */
    msub x5, x6, x4, x7        /* x5 <- multiplicand * multiplier % TestNbr */
    str w5, [x2]               /* Store into result */
    ret
    
    /* Montgomery multiplication with two limbs */
    
two_limbs:
    ldr w6, [x1, #4]
    ldr w1, [x1]               /* R6:R1 = Multiplier */
    /* Processing least significant limb of multiplicand */
    ldr w9, [x0]               /* Get limb from multiplicand. */
    smull x10, w9, w1          /* x10 <- multiplicand[0] * multiplier[0]. */
    mul w8, w10, w3            /* w8 <- MontDig = Pr * MontgomeryMultN mod R */
    bic w8, w8, #0x80000000    /* Clear bit 31. */
    neg w8, w8                 /* r8 <- -MontDig */
    smaddl x10, w8, w4, x10    /* x10 = Pr <- Pr - MontDig * TestNbr[0] */
    asr x10, x10, #31          /* Shift right Pr (x10) by 31 bits. */
    smaddl x10, w9, w6, x10    /* Pr <- Pr + multiplicand[0] * multiplier[1]. */
    smaddl x10, w8, w5, x10    /* Pr <- Pr - MontDig * TestNbr[1]. */
    asr x7, x10, #31           /* Shift right Pr (x10) by 31 bits. */
                               /* and store most significant limb of temp result. */
    bic w10, w10, #0x80000000  /* R7:R10 = Temporary result. */
    
    /* Processing most significant limb of multiplicand */
    ldr w9, [x0, #4]           /* Get limb from multiplicand. */
    smaddl x10, w9, w1, x10    /* x10 <- Pr + multiplicand[1] * multiplier[0]. */
    mul w8, w10, w3            /* w8 <- MontDig = Pr * MontgomeryMultN mod R */
    bic w8, w8, #0x80000000    /* Clear bit 31. */
    neg w8, w8                 /* r8 <- -MontDig */
    smaddl x10, w8, w4, x10    /* x10 = Pr <- Pr - MontDig * TestNbr[0] */
    asr x10, x10, #31          /* Shift right Pr (x10) by 31 bits. */
    smaddl x10, w9, w6, x10    /* Pr <- Pr + multiplicand[1] * multiplier[1]. */
    smaddl x10, w8, w5, x10    /* Pr <- Pr - MontDig * TestNbr[1]. */
    sxtw x12, w7               /* Get MS limb from temp result (sign extended). */
    add x10, x10, x12          /* Pr <- Pr + result[j+1]. */
    bic w6, w10, #0x80000000   /* Clear bit 31 and store limb to temp result[0]. */
    asr x7, x10, #31           /* Shift right Pr (x10) by 31 bits */
                               /* and store most significant limb of temp result. */
    /* Test whether the result is positive. */
    tst w7, #0x80000000        /* Is most significant limb of result positive? */
    beq copy_temp_to_result_2L /* If so, copy temporary result to actual result. */
    
    /* Add TestNbr to result. */
    add w8, w4, w6             /* TestNbr[0] + temporary result[0]. */
    bic w6, w8, #0x80000000    /* Clear bit 31 and store into temp result. */   
    add w8, w5, w8, lsr #31    /* Add TestNbr[1] + temp result[1] with carry. */
    add w8, w8, w7
    bic w7, w8, #0x80000000    /* Clear bit 31 and store into temp result. */
    
copy_temp_to_result_2L:
    str w6, [x2]
    str w7, [x2, #4]
    ret
