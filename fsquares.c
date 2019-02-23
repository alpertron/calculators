/*
    This file is part of Alpertron Calculators.

    Copyright 2015 Dario Alejandro Alpern

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "showtime.h"
#include "batch.h"
#define MAX_SIEVE 65536
#define SUBT 12
extern int lang;
static int primediv[256];
static int primeexp[256];
char hexadecimal;
static limb number[MAX_LEN];
static limb origNbr[MAX_LEN];
static limb p[MAX_LEN];
extern limb q[MAX_LEN];
extern limb Mult1[MAX_LEN];
static limb Mult2[MAX_LEN];
extern limb Mult3[MAX_LEN];
extern limb Mult4[MAX_LEN];
static BigInteger SquareMult1;
static BigInteger SquareMult2;
static BigInteger SquareMult3;
static BigInteger SquareMult4;
static limb Sum[MAX_LEN];
static int iMult3, iMult4;
static int Mult1Len, Mult2Len, Mult3Len, Mult4Len, power4;
static limb result[MAX_LEN];
static int nbrLimbs, origNbrLimbs, groupLength;
static int sieve[MAX_SIEVE];
static int TerminateThread, sum;
static int nbrModExp;
static int Computing3Squares;
static char tmpOutput[MAX_LEN*12];
static double logNormFirst, logNormSecond;
int app;
static char *square = "<span class=\"bigger\">²</span>";
static BigInteger biMult1, biMult2, biMult3, biMult4;
static BigInteger toProcess;
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
char texto[500];
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs);
int checkMinusOne(limb *value, int nbrLimbs);
void batchCubesCallback(char **pptrOutput);
#ifdef __EMSCRIPTEN__
void contfracText(char *input, int groupLen);
void fcubesText(char *input, int groupLen);
#endif

static void ShowStatus(void)
{
#ifdef __EMSCRIPTEN__
  char status[200];
  char *ptrStatus = status;
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  if (elapsedTime / 10 == oldTimeElapsed / 10)
  {
    return;
  }
  oldTimeElapsed = elapsedTime;
  ptrStatus = status;
  strcpy(ptrStatus, lang ? "4<p>Transcurrió " : "4<p>Time elapsed: ");
  ptrStatus += strlen(ptrStatus);
  GetDHMS(&ptrStatus, elapsedTime / 10);
  strcpy(ptrStatus, lang ? "&nbsp;&nbsp;&nbsp;Intentando suma de dos cuadrados de n &minus; " : "&nbsp;&nbsp;&nbsp;Attempting sum of two squares of <var>n</var> &minus; ");
  ptrStatus += strlen(ptrStatus);
  if (hexadecimal)
  {
    int2hex(&ptrStatus, iMult3);
  }
  else
  {
    int2dec(&ptrStatus, iMult3);
  } 
  strcpy(ptrStatus, square);
  ptrStatus += strlen(ptrStatus);
  if (!Computing3Squares)
  {
    strcpy(ptrStatus, " &minus; ");
    ptrStatus += strlen(ptrStatus);
    if (hexadecimal)
    {
      int2hex(&ptrStatus, iMult4);
    }
    else
    {
      int2dec(&ptrStatus, iMult4);
    }
    strcpy(ptrStatus, square);
  }
  databack(status);
#endif
}

 // If Mult1 < Mult2, exchange both numbers.
static void SortBigNbrs(limb *mult1, int *mult1Len, limb *mult2, int *mult2Len)
{
  limb tmp;
  int index;
  limb *ptr1, *ptr2;
  int len = *mult1Len;
  if (len > *mult2Len)
  {
    return;    // mult1 > mult2, so nothing to do.
  }
  if (len == *mult2Len)
  {
    ptr1 = &mult1[len-1];
    ptr2 = &mult2[len-1];
    for (index = *mult1Len - 1; index >= 0; index--)
    {
      if (ptr1->x > ptr2->x)
      {
        return;    // mult1 > mult2, so nothing to do.
      }
      if (ptr1->x < ptr2->x)
      {
        break;     // mult1 < mult2, so exchange them.
      }
      ptr1--;
      ptr2--;
    }
  }
    // Exchange lengths.
  len = *mult2Len;
  *mult2Len = *mult1Len;
  *mult1Len = len;
    // Exchange bytes that compose the numbers.
  ptr1 = mult1;
  ptr2 = mult2;
  for (index = 0; index < len; index++)
  {
    tmp.x = ptr1->x;
    ptr1->x = ptr2->x;
    ptr2->x = tmp.x;
    ptr1++;
    ptr2++;
  }
}

static void MultiplyComplexBy2(BigInteger *real, BigInteger *imag)
{
  BigIntMultiplyBy2(real);
  BigIntMultiplyBy2(imag);
}

static void DivideComplexBy2(BigInteger *real, BigInteger *imag)
{
  BigIntDivideBy2(real);
  BigIntDivideBy2(imag);
}
static void DivideComplexBy1PlusI(BigInteger *real, BigInteger *imag, BigInteger *temp)
{        // Multiply by 1-i and then divide by 2.
  BigIntAdd(real, imag, temp);
  BigIntSubt(real, imag, imag);
  CopyBigInt(real, temp);
  DivideComplexBy2(real, imag);
}

static void MultiplyComplexBy1PlusI(BigInteger *real, BigInteger *imag, BigInteger *temp)
{        // Multiply by 1+i.
  BigIntSubt(real, imag, temp);
  BigIntAdd(real, imag, imag);
  CopyBigInt(real, temp);
}

// Variable to split in up to four squares: number
int fsquares(void)
{
  static BigInteger biTemp;
#ifdef __EMSCRIPTEN__
  char *ptrOutput;
#endif
  int sqrtFound, r, tmp;
  int i, j, numberMod8, nbrDivisors;
  int index, nbrLimbsP, nbrLimbsQ, shRight, shRightMult3, count;
  int divisor, base, idx;
  limb carry;
  nbrLimbs = origNbrLimbs;
  memcpy(number, origNbr, nbrLimbs*sizeof(limb));
  squareRoot(number, Mult1, nbrLimbs, &Mult1Len);
  multiply(Mult1, Mult1, result, Mult1Len, NULL);
  if (memcmp(result, number, sizeof(limb)*nbrLimbs) == 0)
  {          // number is a perfect square.
    Mult2[0].x = 0;
    Mult2Len = 1;
    Mult3[0].x = 0;
    Mult3Len = 1;
    Mult4[0].x = 0;
    Mult4Len = 1;
    power4 = 0;
  }
  else
  {          // number is not a perfect square.
    nbrModExp = 0;
    for (i = 0; i < MAX_SIEVE / 2; i++)
    {
      sieve[i] = 0;
    }
    for (i = 3; i*i < MAX_SIEVE / 2; i += 2)
    {
      j = i*i - 3;
      j = (j % 2 == 0) ? j / 2 : (j + i) / 2;
      for (; j < MAX_SIEVE / 2; j += i)
      {
        sieve[j] = -1;             // Indicate number is composite.
      }
    }
#ifdef __EMSCRIPTEN__
    ptrOutput = tmpOutput;
    strcpy(ptrOutput, "1<p><var>n</var> = ");
    ptrOutput += strlen(ptrOutput);
    if (hexadecimal)
    {
      Bin2Hex(origNbr, ptrOutput, nbrLimbs, groupLen);
    }
    else
    {
      Bin2Dec(origNbr, ptrOutput, nbrLimbs, groupLen);
    }
    ptrOutput += strlen(ptrOutput);
    strcpy(ptrOutput, "</p>");
    databack(tmpOutput);
#endif
    DivideBigNbrByMaxPowerOf4(&power4, number, &nbrLimbs);
    Mult1Len = Mult2Len = 1;
    if (nbrLimbs == 1 && number[0].x < 4)
    {
      iMult3 = iMult4 = 0;
      switch (number[0].x)
      {
      case 3:
        iMult3 = 1;
        /* no break */
      case 2:
        Mult2[0].x = 1;
        /* no break */
      case 1:
        Mult1[0].x = 1;
      }
    }
    else
    {     // Number greater than 3.
      numberMod8 = number[0].x & 7;              // modulus mod 8

      // Fill sieve array

      for (i = 0; i < MAX_SIEVE / 2; i++)
      {
        if (sieve[i] >= 0)                        // If prime...
        {
          limb Divid;
          limb Rem;
          int LimbModQ;
          int Q;
          Rem.x = 0;
          Q = 2 * i + 3;                          // Prime
          LimbModQ = (int)(LIMB_RANGE % Q);
          for (j = nbrLimbs - 1; j >= 0; j--)
          {
            Divid.x = number[j].x + (Rem.x * LimbModQ);
            Rem.x = (int)((unsigned int)Divid.x % Q);
          }
          sieve[i] = (int)Rem.x;
        }
      }

      if (numberMod8 != 7)
      {              // n!=7 (mod 8) => Sum of three squares
        Computing3Squares = 1;
        iMult4 = 0;
        iMult3 = -1;
      }
      else
      {              // n==7 (mod 8) => Sum of four squares
        iMult3 = -1;
        iMult4 = 0;
        Computing3Squares = 0;
      }
      // If number is a sum of three squares, subtract a small square.
      // If number is a sum of four squares, subtract two small squares.
      // If the result is not the product of a power of 2, small primes
      // of the form 4k+1 and squares of primes of the form (4k+3)^2
      // and a (big) prime, try with other squares.
    compute_squares_loop:
      for (;;)
      {
        if (TerminateThread != 0)
        {
          return 2;
        }
        if (Computing3Squares == 0 && iMult3 >= iMult4)
        {
          iMult3 = 1;
          iMult4++;
        }
        else
        {
          iMult3++;
        }
        sum = iMult3 * iMult3 + iMult4 * iMult4;
        carry.x = number[0].x - sum;
        p[0].x = carry.x & MAX_VALUE_LIMB;
        carry.x >>= BITS_PER_GROUP;
        if (nbrLimbs > 1)
        {
          carry.x += number[1].x;
          p[1].x = carry.x & MAX_VALUE_LIMB;
          carry.x >>= BITS_PER_GROUP;
          for (index = 2; index < nbrLimbs; index++)
          {
            carry.x += number[index].x;
            p[index].x = carry.x & MAX_VALUE_LIMB;
            carry.x >>= BITS_PER_GROUP;
          }
        }
        // p should be the product of power of 2,
        // powers of small primes of form 4k+1 and
        // powers of squares of small primes of form 4k+3.
        nbrLimbsP = nbrLimbs;
        DivideBigNbrByMaxPowerOf2(&shRight, p, &nbrLimbsP);
        if ((p[0].x & 0x03) != 0x01)
        {                   // p is not congruent to 1 mod 4
          continue;         // so it cannot be a sum of two squares.
        }

        // Now compute modulus over the product of all odd divisors
        // less than maxsieve.
        if (shRight != 0)
        {
          primediv[0] = 2;                      // Store divisor.
          primeexp[0] = shRight;                // Store exponent.
          nbrDivisors = 1;
        }
        else
        {
          nbrDivisors = 0;
        }
        divisor = 3;
        for (i = 0; i < MAX_SIEVE / 2; i++)
        {
          if (sieve[i] >= 0 && (sieve[i] - sum) % divisor == 0)
          {                          // Divisor found.
            primediv[nbrDivisors] = divisor;   // Store divisor.
            primeexp[nbrDivisors] = 0;         // Store exponent.
            for (;;)
            {
              // Compute the remainder of p and divisor.
              int LimbModDivisor = (int)(LIMB_RANGE % divisor);
              carry.x = 0;
              for (index = nbrLimbsP - 1; index >= 0; index--)
              {
                carry.x = ((unsigned int)(carry.x * LimbModDivisor) + p[index].x) % divisor;
              }
              if (carry.x != 0)
              {     // Number is not a multiple of "divisor".
                break;
              }
              carry.x = 0;
              // Divide by divisor.
              DivBigNbrByInt((int *)p, divisor, (int *)p, nbrLimbsP);
              if (p[nbrLimbsP - 1].x == 0)
              {
                nbrLimbsP--;
              }
              primeexp[nbrDivisors]++;       // Increment exponent.
            }
            if ((divisor & 0x03) == 3 && (primeexp[nbrDivisors] & 1) != 0)
            {                   // Divisor of the form 4k+3 not square.
              goto compute_squares_loop;  // Discard this value.
            }
            nbrDivisors++;
          }
          divisor += 2;
        }

        if (p[0].x == 1 && nbrLimbsP == 1)
        {         // number is the product of only small primes 4k+1 and
                  // squares of primes 4k+3.
          Mult1[0].x = 1;
          Mult2[0].x = 0;
          Mult1[1].x = 0;
          Mult2[1].x = 0;
          nbrLimbs = 1;
        }
        else
        {
          // At this moment p should be prime, otherwise we must try another number.
          p[nbrLimbsP].x = 0;
          memcpy(q, p, (nbrLimbsP+1)*sizeof(p[0]));
          nbrLimbsQ = nbrLimbsP;
          q[0].x--;                     // q = p - 1 (p is odd, so there is no carry).
          memcpy(Mult3, q, (nbrLimbsQ+1)*sizeof(q[0]));
          Mult3Len = nbrLimbsP;
          DivideBigNbrByMaxPowerOf2(&shRightMult3, Mult3, &Mult3Len);
          base = 1;
          memcpy(TestNbr, p, (nbrLimbsP+1)*sizeof(p[0]));
          GetMontgomeryParms(nbrLimbsP);
          sqrtFound = 0;
          do
          {                 // Compute Mult1 = sqrt(-1) (mod p).
            ShowStatus();
            base++;
            modPowBaseInt(base, Mult3, Mult3Len, Mult1); // Mult1 = base^Mult3.
            nbrModExp++;   // Increment number of modular exponentiations.    
            for (i = 0; i < shRightMult3; i++)
            {              // Loop that squares number.
              modmult(Mult1, Mult1, Mult4);
              if (checkMinusOne(Mult4, nbrLimbsP) != 0)
              {
                sqrtFound = 1;
                break;      // Mult1^2 = -1 (mod p), so exit loop.
              }
              memcpy(Mult1, Mult4, nbrLimbsP*sizeof(limb));
            }
            // If power (Mult4) is 1, that means that number is at least PRP,
            // so continue loop trying to find square root of -1.
          } while (memcmp(Mult4, MontgomeryMultR1, nbrLimbsP*sizeof(limb)) == 0);
          if (sqrtFound == 0)
          {            // Cannot find sqrt(-1) (mod p), go to next candidate.
            continue;
          }
          // Convert Mult1 from Montgomery notation to standard number
          // by multiplying by 1 in Montgomery notation.
          memset(Mult2, 0, nbrLimbsP*sizeof(limb));
          Mult2[0].x = 1;
          modmult(Mult2, Mult1, Mult1);  // Mult1 = sqrt(-1) mod p.
          // Find the sum of two squares that equals this PRP x^2 + y^2 = p using
          // Gaussian GCD as: x + iy = gcd(s + i, p)
          // where s = sqrt(-1) mod p

          // Result is stored in biMult1 and biMult2.
          // Initialize real part to square root of (-1).
          memcpy(biMult1.limbs, Mult1, nbrLimbsP * sizeof(limb));
          biMult1.nbrLimbs = nbrLimbsP;
          biMult1.sign = SIGN_POSITIVE;
          intToBigInteger(&biMult2, 1);   // Initialize imaginary part to 1.
          while (biMult1.nbrLimbs > 1 && biMult1.limbs[biMult1.nbrLimbs - 1].x == 0)
          {
            biMult1.nbrLimbs--;
          }
          // Initialize real part to prime.
          memcpy(biMult3.limbs, TestNbr, nbrLimbsP*sizeof(limb));
          biMult3.nbrLimbs = nbrLimbsP;
          biMult3.sign = SIGN_POSITIVE;
          while (biMult3.nbrLimbs > 1 && biMult3.limbs[biMult3.nbrLimbs - 1].x == 0)
          {
            biMult3.nbrLimbs--;
          }
          intToBigInteger(&biMult4, 0);   // Initialize imaginary part to 0.
                                          // Find gcd of (biMult1 + biMult2 * i) and (biMult3 + biMult4 * i)
          GaussianGCD(&biMult1, &biMult2, &biMult3, &biMult4, &SquareMult1, &SquareMult2, &SquareMult3, &SquareMult4);
          nbrLimbs = biMult1.nbrLimbs;
          if (nbrLimbs < biMult2.nbrLimbs)
          {
            nbrLimbs = biMult2.nbrLimbs;
          }
          memset(Mult1, 0, (nbrLimbs + 1) * sizeof(limb));
          memset(Mult2, 0, (nbrLimbs + 1) * sizeof(limb));
          memcpy(Mult1, SquareMult1.limbs, SquareMult1.nbrLimbs * sizeof(limb));
          memcpy(Mult2, SquareMult2.limbs, SquareMult2.nbrLimbs * sizeof(limb));
        }
        // Use the other divisors of modulus in order to get Mult1 and Mult2.

        for (i = 0; i < nbrDivisors; i++)
        {
          divisor = primediv[i];
          for (index = primeexp[i]-1; index > 0; index -= 2)
          {
            MultBigNbrByInt((int*)Mult1, divisor, (int*)Mult1, nbrLimbs + 1);
            MultBigNbrByInt((int*)Mult2, divisor, (int*)Mult2, nbrLimbs + 1);
            if (Mult1[nbrLimbs].x != 0 || Mult2[nbrLimbs].x != 0)
            {
              nbrLimbs++;
              Mult1[nbrLimbs].x = 0;
              Mult2[nbrLimbs].x = 0;
            }
          }
          if (index == 0)
          {
            // Since the value primediv[i] is very low, it is faster to use
            // trial and error than the general method in order to find
            // the sum of two squares.
            j = 1;
            for (;;)
            {
              r = (int)sqrt((double)(divisor - j*j));
              if (r*r + j*j == divisor)
              {
                break;
              }
              j++;
            }
            // Compute Mult1 <- Previous Mult1 * j + Previous Mult2 * r
            // Compute Mult2 <- Previous Mult1 * r - Previous Mult2 * j
            // Use SquareMult1, SquareMult2, SquareMult3, SquareMult4 as temporary storage.
            MultBigNbrByInt((int*)Mult1, j, (int*)SquareMult1.limbs, nbrLimbs + 1);
            MultBigNbrByInt((int*)Mult2, r, (int*)SquareMult2.limbs, nbrLimbs + 1);
            MultBigNbrByInt((int*)Mult1, r, (int*)SquareMult3.limbs, nbrLimbs + 1);
            MultBigNbrByInt((int*)Mult2, j, (int*)SquareMult4.limbs, nbrLimbs + 1);
            AddBigNbr((int *)SquareMult1.limbs, (int *)SquareMult2.limbs, (int *)Mult1, nbrLimbs + 1);
            SubtractBigNbrB((int *)SquareMult3.limbs, (int *)SquareMult4.limbs, (int *)Mult2, nbrLimbs + 1);
            if ((unsigned int)Mult2[nbrLimbs].x >= LIMB_RANGE)
            {   // Since Mult2 is a difference of products, it can be
                // negative. In this case replace it by its absolute value.
              ChSignBigNbr((int *)Mult2, nbrLimbs+1);
            }
            if (Mult1[nbrLimbs].x != 0 || Mult2[nbrLimbs].x != 0)
            {
              nbrLimbs++;
              Mult1[nbrLimbs].x = 0;
              Mult2[nbrLimbs].x = 0;
            }
          }
        }            /* end for */
        break;
      }              /* end while */
    }
    // Shift left the number of bits that the original number was divided by 4.
    for (count = 0; count < power4; count++)
    {
      MultBigNbrByInt((int*)Mult1, 2, (int*)Mult1, nbrLimbs + 1);
      MultBigNbrByInt((int*)Mult2, 2, (int*)Mult2, nbrLimbs + 1);
      if (Mult1[nbrLimbs].x != 0 || Mult2[nbrLimbs].x != 0)
      {
        nbrLimbs++;
        Mult1[nbrLimbs].x = 0;
        Mult2[nbrLimbs].x = 0;
      }
    }
    nbrLimbsP = 1;
    Mult3[0].x = iMult3;
    Mult4[0].x = iMult4;
    Mult3[1].x = 0;
    Mult4[1].x = 0;
    for (count = 0; count < power4; count++)
    {
      MultBigNbrByInt((int*)Mult3, 2, (int*)Mult3, nbrLimbsP + 1);
      MultBigNbrByInt((int*)Mult4, 2, (int*)Mult4, nbrLimbsP + 1);
      if (Mult3[nbrLimbsP].x != 0 || Mult4[nbrLimbsP].x != 0)
      {
        nbrLimbsP++;
      }
      Mult3[nbrLimbsP].x = 0;
      Mult4[nbrLimbsP].x = 0;
    }
    Mult1Len = Mult2Len = nbrLimbs;
    while (Mult1[Mult1Len-1].x == 0 && Mult1Len > 1)
    {
      Mult1Len--;
    }
    while (Mult2[Mult2Len-1].x == 0 && Mult2Len > 1)
    {
      Mult2Len--;
    }
    Mult3Len = Mult4Len = nbrLimbsP;
    while (Mult3[Mult3Len-1].x == 0 && Mult3Len > 1)
    {
      Mult3Len--;
    }
    while (Mult1[Mult4Len-1].x == 0 && Mult4Len > 1)
    {
      Mult4Len--;
    }
    // Sort squares
    SortBigNbrs(Mult1, &Mult1Len, Mult2, &Mult2Len);
    SortBigNbrs(Mult1, &Mult1Len, Mult3, &Mult3Len);
    SortBigNbrs(Mult1, &Mult1Len, Mult4, &Mult4Len);
    SortBigNbrs(Mult2, &Mult2Len, Mult3, &Mult3Len);
    SortBigNbrs(Mult2, &Mult2Len, Mult4, &Mult4Len);
    SortBigNbrs(Mult3, &Mult3Len, Mult4, &Mult4Len);
  }
    // Validate result.
  idx = Mult1Len * 2;
  multiply(Mult1, Mult1, SquareMult1.limbs, Mult1Len, &tmp);
  SquareMult1.limbs[idx].x = 0;
  multiply(Mult2, Mult2, SquareMult2.limbs, Mult2Len, &tmp);
  memset(&SquareMult2.limbs[Mult2Len << 1], 0, ((Mult1Len - Mult2Len) << 1)*sizeof(limb));
  SquareMult2.limbs[idx].x = 0;
  multiply(Mult3, Mult3, SquareMult3.limbs, Mult3Len, &tmp);
  memset(&SquareMult3.limbs[Mult3Len << 1], 0, ((Mult1Len - Mult3Len) << 1)*sizeof(limb));
  SquareMult3.limbs[idx].x = 0;
  multiply(Mult4, Mult4, SquareMult4.limbs, Mult4Len, &tmp);
  memset(&SquareMult4.limbs[Mult4Len << 1], 0, ((Mult1Len - Mult4Len) << 1)*sizeof(limb));
  SquareMult4.limbs[idx].x = 0;
  idx++;
  AddBigInt(SquareMult1.limbs, SquareMult2.limbs, SquareMult1.limbs, idx);
  AddBigInt(SquareMult1.limbs, SquareMult3.limbs, SquareMult1.limbs, idx);
  AddBigInt(SquareMult1.limbs, SquareMult4.limbs, SquareMult1.limbs, idx);
  while (idx > 1 && SquareMult1.limbs[idx - 1].x == 0)
  {
    idx--;
  }
  if (idx != origNbrLimbs)
  {           // Invalid length.
    return 1;
  }
  for (index = 0; index < idx; index++)
  {
    if (SquareMult1.limbs[index].x != origNbr[index].x)
    {
      return 1;
    }
  }
  return 0;
}

void fsquaresText(char *input, int groupLen)
{
  char *ptrOutput;
#ifdef __EMSCRIPTEN__
  int elapsedTime;
#endif
  if (valuesProcessed == 0)
  {
    groupLength = groupLen;
  }
  BatchProcessing(input, &toProcess, &ptrOutput);
#ifdef __EMSCRIPTEN__
  strcpy(ptrOutput, lang ? "<p>Transcurrió " : "<p>Time elapsed: ");
  ptrOutput += strlen(ptrOutput);
  elapsedTime = (int)(tenths() - originalTenthSecond);
  GetDHMSt(&ptrOutput, elapsedTime);
#endif
  strcpy(ptrOutput, (lang ? "</p><p>" COPYRIGHT_SPANISH "</p>" :
    "</p><p>" COPYRIGHT_ENGLISH "</p>"));
}

#ifdef FSQUARES_APP
void batchCallback(char **pptrOutput)
{
  int result;
  char *ptrOutput;
  if (app == 1)
  {
    batchCubesCallback(pptrOutput);
    return;
  }
  ptrOutput = *pptrOutput;
  NumberLength = toProcess.nbrLimbs;
  CompressBigInteger((int *)number, &toProcess);
  origNbrLimbs = toProcess.nbrLimbs;
  memcpy(origNbr, toProcess.limbs, origNbrLimbs*sizeof(limb));
  result = fsquares();
  // Show the number to be decomposed into sum of squares.
  strcpy(ptrOutput, "<p>");
  ptrOutput += strlen(ptrOutput);
  if (hexadecimal)
  {
    BigInteger2Hex(&toProcess, ptrOutput, groupLen);
  }
  else
  {
    BigInteger2Dec(&toProcess, ptrOutput, groupLen);
  }
  ptrOutput += strlen(ptrOutput);
  if (toProcess.sign == SIGN_NEGATIVE)
  {
    *ptrOutput++ = ':';
    *ptrOutput++ = ' ';
    textError(ptrOutput, EXPR_NUMBER_TOO_LOW);
    ptrOutput += strlen(ptrOutput);
    strcpy(ptrOutput, "</p>");
    *pptrOutput = ptrOutput + strlen(ptrOutput);
    return;
  }
  switch (result)
  {
  case 1:
    strcpy(ptrOutput, (lang==0?": Internal error!\n\nPlease send the number to the author of the applet.</p>":
      ": ¡Error interno!\n\nPor favor envíe este número al autor del applet.</p>"));
    *pptrOutput = ptrOutput + strlen(ptrOutput);
    return;
  case 2:
    strcpy(ptrOutput, (lang==0?": User stopped the calculation":": El usuario detuvo el cálculo"));
    *pptrOutput = ptrOutput + strlen(ptrOutput);
    return;
  }
  // Show the decomposition.
  strcpy(ptrOutput, " = ");
  ptrOutput += strlen(ptrOutput);
  if (hexadecimal)
  {
    Bin2Hex(Mult1, ptrOutput, Mult1Len, groupLength);
  }
  else
  {
    Bin2Dec(Mult1, ptrOutput, Mult1Len, groupLength);
  }
  ptrOutput += strlen(ptrOutput);
  strcpy(ptrOutput, square);
  ptrOutput += strlen(ptrOutput);
  if (Mult2Len != 1 || Mult2[0].x != 0)
  {
    strcpy(ptrOutput, " + ");
    ptrOutput += strlen(ptrOutput);
    if (hexadecimal)
    {
      Bin2Hex(Mult2, ptrOutput, Mult2Len, groupLength);
    }
    else
    {
      Bin2Dec(Mult2, ptrOutput, Mult2Len, groupLength);
    }
    ptrOutput += strlen(ptrOutput);
    strcpy(ptrOutput, square);
    ptrOutput += strlen(ptrOutput);
  }
  if (Mult3Len != 1 || Mult3[0].x != 0)
  {
    strcpy(ptrOutput, " + ");
    ptrOutput += strlen(ptrOutput);
    if (hexadecimal)
    {
      Bin2Hex(Mult3, ptrOutput, Mult3Len, groupLength);
    }
    else
    {
      Bin2Dec(Mult3, ptrOutput, Mult3Len, groupLength);
    }
    ptrOutput += strlen(ptrOutput);
    strcpy(ptrOutput, square);
    ptrOutput += strlen(ptrOutput);
  }
  if (Mult4Len != 1 || Mult4[0].x != 0)
  {
    strcpy(ptrOutput, " + ");
    ptrOutput += strlen(ptrOutput);
    if (hexadecimal)
    {
      Bin2Hex(Mult4, ptrOutput, Mult4Len, groupLength);
    }
    else
    {
      Bin2Dec(Mult4, ptrOutput, Mult4Len, groupLength);
    }
    ptrOutput += strlen(ptrOutput);
    strcpy(ptrOutput, square);
    ptrOutput += strlen(ptrOutput);
  }
  strcpy(ptrOutput, "</p>");
  ptrOutput += strlen(ptrOutput);
  *pptrOutput = ptrOutput;
}
#endif

#ifdef __EMSCRIPTEN__
EXTERNALIZE void doWork(void)
{
  int groupLen = 0;
  char *ptrData = inputString;
#ifdef __EMSCRIPTEN__
  originalTenthSecond = tenths();
#endif
  if (*ptrData == 'C')
  {    // User pressed Continue button.
    if (app == 0)
    {
      fsquaresText(NULL, 0); // Routine does not use prameters in this case.
    }
    else
    {
      fcubesText(NULL, 0); // Routine does not use prameters in this case.
    }
#ifdef __EMSCRIPTEN__
    databack(output);
#endif
    return;
  }
  valuesProcessed = 0;
  while (*ptrData != ',')
  {
    groupLen = groupLen * 10 + (*ptrData++ - '0');
  }
  ptrData++;             // Skip comma.
  app = *ptrData - '0';
  if (*(ptrData + 1) != ',')
  {
    ptrData++;
    app = app * 10 + *ptrData - '0';
  }
  lang = app & 1;
  app >>= 1;
  if (app & 0x20)
  {
    app &= 0x1F;
    hexadecimal = 1;
  }
  else
  {
    hexadecimal = 0;
  }
  switch (app)
  {
  case 0:
    fsquaresText(ptrData+2, groupLen);
    break;
  case 1:
    fcubesText(ptrData+2, groupLen);
    break;
  case 2:
    contfracText(ptrData+2, groupLen);
    break;
  }
  databack(output);
}
#endif
