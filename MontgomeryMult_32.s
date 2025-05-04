@
@ This file is part of Alpertron Calculators.
@
@ Copyright 2025 Dario Alejandro Alpern
@
@ Alpertron Calculators is free software: you can redistribute it and/or modify
@ it under the terms of the GNU General Public License as published by
@ the Free Software Foundation, either version 3 of the License, or
@ (at your option) any later version.
@
@ Alpertron Calculators is distributed in the hope that it will be useful,
@ but WITHOUT ANY WARRANTY; without even the implied warranty of
@ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
@ GNU General Public License for more details.
@
@ You should have received a copy of the GNU General Public License
@ along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
@
    .syntax unified
    .arch armv7-a
    .thumb
    .global MontgomeryMult
    .thumb_func
    .type MontgomeryMult, %function
    .text
@ R0: Pointer to factor1 in Montgomery notation.
@ R1: Pointer to factor2 in Montgomery notation.
@ R2: Pointer to product in Montgomery notation.
@ Other registers used in the code:
@ r3: Value of MontgomeryMultN.
@ r4: Pointer to TestNbr.
@ r5: Value of NumberLength.
@ r6: Index for multiplicand (0 to NumberLength - 1).
@ r7: Index for inner loop (0 to NumberLength - 1).
@ r9: Limb from multiplicand (index r6).
@ r11:r10: Accumulator for 32x32 bit multiplications.
@
@ Multiply two numbers in Montgomery notation.
@
@ For large numbers the REDC algorithm is:
@ m <- ((T mod R)N') mod R
@ t <- (T - mN) / R
@ if t < 0 then
@   return t + N
@ else
@   return t
@ end if
@
@ Multilimb version when multiplying S*T:
@ 1) P <- 0
@ 2) for i from 0 to n-1
@ 3)   t <- (S[i] * T[0]) mod R
@ 4)   u <- (t * N') mod R (one limb)
@ 5)   P <- (P + S[i] * T - u * N) / R
@ 6) if P < 0 then
@ 7)   return P - N
@ 8) else
@ 9)  return P
@ 10)end if
    .set MAX_LIMBS_MONTGOMERY, 14
    .set STACK_SIZE_RESULT_BUFFER, (MAX_LIMBS_MONTGOMERY+3)/4 * 16
MontgomeryMult:
    push {r4-r11, lr}
    @ Get external values and pointers to buffers.
    ldr r7, locGOT
.LPIC0:
    add r7, pc, r7
    ldr r3, MontgomeryMultN_offset
    ldr r3, [r3, r7]           @ r3 <- Address of MontgomeryMultN.
    ldr r3, [r3]               @ r3 <- MontgomeryMultN.
    ldr r4, TestNbr_offset
    ldr r4, [r4, r7]           @ r4 <- Address of MontgomeryMultN.
    ldr r5, NumberLength_offset
    ldr r5, [r5, r7]           @ r5 <- Address of NumberLength.
    ldr r5, [r5]               @ r5 <- NumberLength.
    cmp r5, #2
    bne more_than_two_limbs
    
    @ Montgomery multiplication with two limbs.
    @ Optimizations: Do not use loops. Use registers.
    
    ldr r5, [r4, #4]
    ldr r4, [r4]               @ R5:R4 = TestNbr.
    ldr r6, [r1, #4]
    ldr r1, [r1]               @ R6:R1 = Multiplier.
    /* Processing least significant limb of multiplicand */
    ldr r9, [r0]               @ Get limb from multiplicand.
    smull r10, r11, r9, r1     @ r11:r10 <- multiplicand[0] * multiplier[0].
    mul r8, r10, r3            @ r8 <- MontDig = Pr * MontgomeryMultN mod R
    bic r8, r8, #0x80000000    @ Clear bit 31.
    rsb r8, r8, #0             @ r8 <- -MontDig
    smlal r10, r11, r8, r4     @ r11:r10 = Pr <- Pr - MontDig * TestNbr[0]
                               @ Shift right Pr (r11:r10) 31 bits.
    lsl r12, r11, #1
    asr r11, r11, #31
    orr r10, r12, r10, lsr #31
    smlal r10, r11, r9, r6     @ Pr <- Pr + multiplicand[0] * multiplier[1].
    smlal r10, r11, r8, r5     @ Pr <- Pr - MontDig * TestNbr[1].
    lsl r7, r11, #1            @ Store most significant limb of result.
    orr r7, r7, r10, lsr #31   @ as shift right Pr (r11:r10) 31 bits.
    bic r10, r10, #0x80000000  @ r7:r10 = Temporary result.

    @ Processing most significant limb of multiplicand.
    ldr r9, [r0, #4]           @ Get limb from multiplicand.
    mov r11, #0
    smlal r10, r11, r9, r1     @ r11:r10 <- Pr + multiplicand[1] * multiplier[0].
    mul r8, r10, r3            @ r8 <- MontDig = Pr * MontgomeryMultN mod R
    bic r8, r8, #0x80000000    @ Clear bit 31.
    rsb r8, r8, #0             @ r8 <- -MontDig
    smlal r10, r11, r8, r4     @ r11:r10 = Pr <- Pr - MontDig * TestNbr[0]
                               @ Shift right Pr (r11:r10) 31 bits.
    lsl r12, r11, #1
    asr r11, r11, #31
    orr r10, r12, r10, lsr #31
    smlal r10, r11, r9, r6     @ Pr <- Pr + multiplicand[1] * multiplier[1].
    smlal r10, r11, r8, r5     @ Pr <- Pr - MontDig * TestNbr[1].
    adds r10, r10, r7          @ Pr <- Pr + result[1] (sign extended).
    adc r11, r11, r7, asr #31
    bic r6, r10, #0x80000000   @ Clear bit 31 and store limb of temp result[0].
                               @ Shift right Pr (r11:r10) 31 bits.
    lsl r7, r11, #1
    orrs r7, r7, r10, lsr #31 @ Store most significant limb of result.
    
    @ Test whether the result is positive.
    bpl copy_temp_to_result_2L @ If so, copy temporary result to actual result.
    
    @ Add TestNbr to result.
    add r8, r4, r6             /* TestNbr[0] + temporary result[0]. */
    bic r6, r8, #0x80000000    /* Clear bit 31 and store into temp result. */   
    add r8, r5, r8, lsr #31    /* Add TestNbr[1] + temp result[1] with carry. */
    add r8, r8, r7
    bic r7, r8, #0x80000000    /* Clear bit 31 and store into temp result. */
   
copy_temp_to_result_2L:
    str r6, [r2]
    str r7, [r2, #4]
    pop {r4-r11, pc}
    
more_than_two_limbs:
    @ Adjust stack to store temporary result there
    @ so the arguments are not overwritten with the result.
    sub sp, sp, #STACK_SIZE_RESULT_BUFFER

    @ Clear result.
    mov r7, r5
    mov r6, #0
clear_result_loop:
    subs r7, r7, #1
    str r6, [sp, r7, lsl #2]
    bne clear_result_loop
    
    @ Outer loop.
    @ Index for multiplicand (r6) is already zero.
outer_loop:
    mov r7, #1                 @ Initialize index for result.
    ldr r9, [r0, r6, lsl #2]   @ Get limb from multiplicand.
    ldr r10, [sp]              @ Get the least significant limb of result.
    mov r11, #0
    ldr r12, [r1]              @ Get the least significant limb of multiplier.
    smlal r10, r11, r9, r12    @ r11:r10 <- Pr + multiplicand[i] * multiplier[0].
    mul r8, r10, r3            @ r8 <- MontDig = Pr * MontgomeryMultN mod R
    bic r8, r8, #0x80000000    @ Clear bit 31.
    rsb r8, r8, #0             @ r8 <- -MontDig
    ldr r12, [r4]              @ Get least significant of TestNbr.
    smlal r10, r11, r8, r12    @ r11:r10 = Pr <- Pr - MontDig * TestNbr[0]
inner_loop:
                               @ Shift right Pr (r11:r10) 31 bits.
    lsl r12, r11, #1
    asr r11, r11, #31
    orr r10, r12, r10, lsr #31
    ldr r12, [r1, r7, lsl #2]  @ Get limb from multiplier.
    smlal r10, r11, r9, r12    @ Pr <- Pr + multiplicand[i] * multiplier[j+1].
    ldr r12, [r4, r7, lsl #2]  @ Get limb from TestNbr.
    smlal r10, r11, r8, r12    @ Pr <- Pr - MontDig * TestNbr[j+1].
    ldr r12, [sp, r7, lsl #2]  @ Get limb from result.
    adds r10, r10, r12         @ Pr <- Pr + result[j+1] (sign extended).
    adc r11, r11, r12, asr #31
    sub r7, r7, #1             @ Decrement index.
    bic r12, r10, #0x80000000  @ Clear bit 31.
    str r12, [sp, r7, lsl #2]  @ Store limb of result[j].
    add r7, r7, #2
    cmp r7, r5
    bne inner_loop
    sub r7, r7, #1
                               @ Shift right Pr (r11:r10) 31 bits.
    lsl r12, r11, #1
    orr r10, r12, r10, lsr #31
    str r10, [sp, r7, lsl #2]  @ Store most significant limb of result.
    add r6, r6, #1             @ Increment index to multiplicand.
    cmp r6, r5
    bne outer_loop
    
    @ Test whether the result is positive.
    tst r10, #0x80000000       @ Is most significant limb of result positive?
    beq copy_temp_to_result    @ Jump to end of routine if so.
    
    @ Add TestNbr to result.
    mov r6, #0                 @ Initialize index.
    mov r7, #0                 @ Initialize result of previous addition.
add_TestNbr_cycle:
    ldr r8, [sp, r6, lsl #2]   @ Get limb from result.
    ldr r9, [r4, r6, lsl #2]   @ Get limb from TestNbr.
    add r7, r8, r7, lsr #31    @ Compute addition with carry.
    add r7, r7, r9
    bic r8, r7, #0x80000000    @ Clear bit 31.
    str r8, [sp, r6, lsl #2]   @ Store limb into result.
    add r6, r6, #1
    cmp r6, r5
    bne add_TestNbr_cycle
copy_temp_to_result:
copy_temp_to_result_loop:
    subs r5, r5, #1
    ldr r8, [sp, r5, lsl #2]  /* Get limb from temporary result. */
    str r8, [r2, r5, lsl #2]  /* Set limb to result. */
    bne copy_temp_to_result_loop
    add sp, sp, #STACK_SIZE_RESULT_BUFFER
    pop {r4-r11, pc}
locGOT:
    .word _GLOBAL_OFFSET_TABLE_-(.LPIC0+4)
TestNbr_offset:
    .word TestNbr(GOT)
MontgomeryMultN_offset:
    .word MontgomeryMultN(GOT)
NumberLength_offset:
    .word NumberLength(GOT)
    .end
