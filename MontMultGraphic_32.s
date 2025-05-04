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
    .global MontMultGraphic
    .thumb_func
    .type MontMultGraphic, %function
    .text
@ R0: Pointer to multiplicand in Montgomery notation.
@ R1: Pointer to multiplier in Montgomery notation.
@ R2: Pointer to product in Montgomery notation.
@ Other registers used in the code:
@ r3: Value of MontgomeryMultN.
@ r6:r1: Multiplier.
@ r5:r4: TestNbr.
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
MontMultGraphic:
    push {r3-r12, lr}
    @ Get external values and pointers to buffers.
    ldr r7, locGOT
.LPIC0:
    add r7, pc, r7
    ldr r3, MontMultGraphicN_offset
    ldr r3, [r3, r7]           @ r3 <- Address of MontgomeryMultN.
    ldr r3, [r3]               @ r3 <- MontgomeryMultN.
    ldr r4, TestNbrGraphic_offset
    ldr r4, [r4, r7]           @ r4 <- Address of MontgomeryMultN.        
    ldr r5, [r4, #4]
    ldr r4, [r4]               @ R5:R4 = TestNbr.
    cmp r5, #0
    bne two_limbs
    
    @ Montgomery multiplication with one limb.
    
    ldr r8, [r0]               @ Get multiplicand.
    ldr r9, [r1]               @ Get multiplier.
    umull r10, r11, r8, r9     @ R11:R10 <- multiplicand * multiplier
    mul r8, r10, r3            @ r8 <- MontDig = Pr * MontgomeryMultN mod R
    bic r8, r8, #0x80000000    @ Clear bit 31.
    rsb r8, r8, #0             @ r8 <- -MontDig
    smlal r10, r11, r8, r4     @ r11:r10 = Pr <- Pr - MontDig * TestNbr[0]
                               @ Shift right Pr (r11:r10) 31 bits.    
    lsl r7, r11, #1
    orrs r7, r7, r10, lsr #31  @ Get result.
    @ Test whether the result is positive.
    bpl copy_temp_to_result_1L @ If so, copy temporary result to actual result.
    
    @ Add TestNbr to result.
    add r7, r4, r7             @ TestNbr[0] + temporary result[0].
copy_temp_to_result_1L:
    str r7, [r2]               @ Store into result
    str r5, [r2, #4]
    pop {r3-r12, pc}
    
    @ Montgomery multiplication with two limbs.
    
two_limbs:
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
    orrs r7, r7, r10, lsr #31  @ Store most significant limb of result.
    
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
    pop {r3-r12, pc}

locGOT:
    .word _GLOBAL_OFFSET_TABLE_-(.LPIC0+4)
TestNbrGraphic_offset:
    .word TestNbrGraphic(GOT)
MontMultGraphicN_offset:
    .word MontMultGraphicN(GOT)
    .end
