//
// This file is part of Alpertron Calculators.
//
// Copyright 2015-2024 Dario Alejandro Alpern
//
// Alpertron Calculators is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Alpertron Calculators is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
//
#include <string.h>
#include <stdlib.h>
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"
#include "polynomial.h"
#include "showtime.h"
#include "rootseq.h"
#ifdef FACTORIZATION_APP
#include "factor.h"
#endif
#ifdef __EMSCRIPTEN__
static char outputText[20000];
char* ptrPercentageOutput;
#endif

int attemptNbr;
static bool firstFactor;
static bool isCharacteristic2;
extern int poly4[1000000];
extern char* ptrOutput;
extern int grpLen;
extern bool teachMod;

static void showEqu(void)
{
  if (pretty == TEX)
  {
    showText("\\,\\equiv\\,");
  }
  else if (pretty == PRETTY_PRINT)
  {
    showText(" &equiv; ");
  }
  else
  {    // Pari-GP
    showText(" = Mod(");
  }
}

static void showEquAndPoly(const int* ptrPoly, int polyDegree)
{
  showEqu();
  if (polyDegree < 0)
  {
    showText("0");
  }
  else
  {
    showPowerVar(&ptrOutput, polyDegree, 'x');
    showMontPolynomial(&ptrOutput, ptrPoly,
      polyDegree, grpLen);
  }
  if (pretty == PARI_GP)
  {
    showText(", ");
    Bin2Dec(&ptrOutput, primeMod.limbs, primeMod.nbrLimbs, grpLen);
    showText(")");
  }
}

static void restartRecursion(const char* ptrPolyName, int degreePoly,
  const int *poly, int expon)
{
  if (teachMod)
  {
    showText("<p><var>");
    showText(ptrPolyName);
    showText("</var>(x) ");
    showText(lang ? "es una potencia perfecta de grado " :
      "is a perfect power of degree ");
    int2dec(&ptrOutput, primeMod.limbs[0].x);
    showText(lang ? ". Ejecutando el algoritmo nuevamente con:</p><p>" :
      ". Executing the algorithm again with:</p><p>");
    showVarIndex('f', expon);
    showEquAndPoly(poly, degreePoly);
    showText("</p>");
  }
}

#if 0
// Square-free factorization (from Wikipedia).
// i <- 1; g <- f';
// if g != 0 then {
//   c <- gcd(f, g);
//   w <- f / c;
//   while w != 1 do {
//     y <- gcd(w, c); z <- w / y;
//     Output(z^i); i <- i + 1;
//     w <- y; c <- c / y
//   }
//   if c != 1 then {
//     c <- c^(1/p);
//     Output(SFF(c)^p)
//   }
// }
// else {
//    f <- f^(1/p);
//    Output(SFF(f)^p)
//  }
//  end.
#endif
static void SquareFreeFactorization(int polyDegr, int* poly)
{
  int expon = 1;
  int polyDegree = polyDegr;
  int currentDegree;
  int degreeC;
  int degreeY;
  int index;
  int primeInt = primeMod.limbs[0].x;
  int nbrLimbs = primeMod.nbrLimbs + 1;
  int* ptrValue1;
  int* ptrValue2;
  int* ptrF = poly;
  // Generate derivative in poly1.
  if (teachMod)
  {
    showText(lang ? "<h3>Separación de factores con factores repetidos</h3>" :
      "<h3>Squarefree factorization</h3>");
  }
  for (;;)
  {
    int lenLimbs;
    char buf[100];
    ptrValue1 = ptrF;
    ptrValue2 = poly1;
    for (currentDegree = 1; currentDegree <= polyDegree; currentDegree++)
    {
      ptrValue1 += nbrLimbs;
      IntArray2BigInteger(ptrValue1, &operand1);
      modmultInt(operand1.limbs, currentDegree, operand1.limbs);
      BigInteger2IntArray(ptrValue2, &operand1);
      ptrValue2 += nbrLimbs;
    }
    // Check whether the derivative is zero.
    for (currentDegree = polyDegree - 1; currentDegree >= 0; currentDegree--)
    {
      ptrValue1 = &poly1[(currentDegree * nbrLimbs) + 1];
      for (index = 1; index < nbrLimbs; index++)
      {
        if (*ptrValue1 != 0)
        {         // Coefficient of derivative is not zero.
          break;
        }
        ptrValue1++;
      }
      if (index < nbrLimbs)
      {           // Coefficient of derivative is not zero.
        break;
      }
    }
    if (teachMod)
    {
      showText(lang ? "<p>La derivada de f(x) es:</p><p>" :
        "<p>The derivative of f(x) is:</p><p>");
      showText((pretty == PARI_GP) ? "deriv(f(x))" : "<var>f</var> '(x)");
      showEquAndPoly(poly1, currentDegree);
      showText("</p>");
    }
    if (currentDegree >= 0)
    {             // Derivative is not zero.
      int lenBytes;
      int degreeW;
      int i = 1;

      PolyModularGcd(ptrF, polyDegree, poly1, currentDegree, poly2, &degreeC); // poly2 = c
      if (teachMod)
      {
        showText("<p>");
        showVarIndex('c', 0);
        showText(lang ? "(x) = mcd" : "(x) = gcd");
        showText("(f(x), ");
        showText((pretty == PARI_GP)?"deriv(f(x))": "f '(x)");
        showEquAndPoly(poly2, degreeC);
        showText("</p><p>");
      }
      lenBytes = (polyDegree + 1) * nbrLimbs * (int)sizeof(int);
      (void)memcpy(poly4, ptrF, lenBytes);        // Backup poly
      SetNumberToOne(&poly2[degreeC * nbrLimbs]);
      if (degreeC == 0)
      {
        (void)memcpy(poly1, ptrF, lenBytes);           // poly1 = w
      }
      else
      {
        DividePolynomial(poly4, polyDegree, poly2, degreeC, poly1);    // poly1 = w
      }
      degreeW = polyDegree - degreeC;
      if (teachMod)
      {
        showVarIndex('w', 0);
        showText("(x) = f(x)/");
        showVarIndex('c', 0);
        showText("(x)");
        showEquAndPoly(poly1, degreeW);
        showText("</p><p>");
      }
      while (degreeW != 0)
      {
        lenBytes = (degreeW + 1) * nbrLimbs * (int)sizeof(int);
        (void)memcpy(poly4, poly1, lenBytes);        // Backup w.
        PolyModularGcd(poly2, degreeC, poly1, degreeW, poly3, &degreeY);      // poly3 = y
        if (teachMod)
        {
          showText("<p>");
          showVarIndex('w', i);
          showText(lang ? "(x) = mcd(" : "(x) = gcd(");
          showVarIndex('w', i - 1);
          showText(", ");
          showVarIndex('c', i - 1);
          showText(")");
          showEquAndPoly(poly3, degreeY);
          showText("</p><p>");
        }
        SetNumberToOne(&poly3[degreeY * nbrLimbs]);
        if (degreeW != degreeY)
        {
          int degreeZ;
          int lengthLimbs;
          struct sFactorInfo* pstFactorInfo;

          DividePolynomial(poly4, degreeW, poly3, degreeY, poly1);     // poly1 = z
          // z^i is divisor of the original polynomial.
          degreeZ = degreeW - degreeY;
          if (teachMod)
          {
            showText("<strong>");
            showVarIndex('z', i);
            showText("(x) = ");
            showVarIndex('w', i - 1);
            showText("(x) / ");
            showVarIndex('w', i);
            showText("(x)");
            showEquAndPoly(poly1, degreeZ);
            showText(lang ? " es un factor de f(x) con multiplicidad " :
              " is a factor of f(x) with multiplicity ");
            int2dec(&ptrOutput, i * expon);
            showText("</strong></p>");
          }
          for (currentDegree = 0; currentDegree < (i * expon); currentDegree++)
          {
            DividePolynomial(ptrOrigPoly, degreeOrigPoly, poly1, degreeZ, poly4);
            degreeOrigPoly -= degreeZ;
            lenBytes = (degreeOrigPoly + 1) * nbrLimbs * (int)sizeof(int);
            (void)memcpy(ptrOrigPoly, poly4, lenBytes);
          }
          pstFactorInfo = &factorInfo[nbrFactorsFound];
          nbrFactorsFound++;
          pstFactorInfo->ptr = ptrOrigPoly;
          pstFactorInfo->degree = degreeZ;
          pstFactorInfo->multiplicity = i * expon;
          pstFactorInfo->expectedDegree = 0;    // Unknown at this moment.
          lenBytes = degreeZ * nbrLimbs * (int)sizeof(int);
          (void)memcpy(ptrOrigPoly, poly1, lenBytes);
          lengthLimbs = (degreeZ + 1) * nbrLimbs;
          ptrOrigPoly += lengthLimbs;
          lenBytes = (degreeOrigPoly + 1) * nbrLimbs * (int)sizeof(int);
          (void)memcpy(ptrOrigPoly, poly4, lenBytes);
        }
        i++;
        lenBytes = (degreeY + 1) * nbrLimbs * (int)sizeof(int);
        (void)memcpy(poly1, poly3, lenBytes);  // Copy y to w.
        degreeW = getDegreePoly(poly1, degreeY);
        DividePolynomial(poly2, degreeC, poly1, degreeW, poly4); // Compute c.
        degreeC -= degreeW;
        if (teachMod)
        {
          showVarIndex('c', i - 1);
          showText("(x) = ");
          showVarIndex('c', i - 2);
          showText("(x) / ");
          showVarIndex('w', i - 1);
          showText("(x)");
          showEquAndPoly(poly4, degreeC);
          showText("</p><p>");
        }
        lenBytes = (degreeC + 1) * nbrLimbs * (int)sizeof(int);
        (void)memcpy(poly2, poly4, lenBytes);
      }
      if (degreeC == 0)
      {     // End of square free factorization.
        return;
      }
      // C is a perfect pth-power.
      lenLimbs = (degreeOrigPoly + 1) * nbrLimbs;
      ptrValue2 = poly2;
      char* ptrOutBak = ptrOutput;
      ptrOutput = buf;
      showVarIndex('c', i - 1);
      *ptrOutput = 0;
      ptrOutput = ptrOutBak;
    }
    else
    {           // Derivative of f(x) is zero.
      lenLimbs = (degreeOrigPoly + 1) * nbrLimbs;
      ptrValue2 = ptrF;
      strcpy(buf, "f");
      degreeC = polyDegree;
    }
    ptrValue1 = ptrOrigPoly + lenLimbs;
    // Copy the coefficients x^0, x^p, x^(2p), ... to
    // x^0, x^1, x^2, ... on destination.
    // This is the pth root mod p.
    for (currentDegree = 0; currentDegree <= degreeC; currentDegree += primeInt)
    {
      int lenBytes = nbrLimbs * (int)sizeof(int);
      (void)memcpy(ptrValue1, ptrValue2, lenBytes);
      ptrValue1 += nbrLimbs;
      lenLimbs = primeInt * nbrLimbs;
      ptrValue2 += lenLimbs;
    }
    lenLimbs = (degreeOrigPoly + 1) * nbrLimbs;
    ptrF = ptrOrigPoly + lenLimbs;
    polyDegree = degreeC / primeInt;
    expon *= primeInt;
    restartRecursion(buf, polyDegree, ptrF, expon);
  }
}

static void showGcdDistinctDegFact(const char* ptrExp)
{
  showText(lang ? "mcd(" : "gcd(");
  showVariable(&ptrOutput, 'f');
  showText("(x), ");
  showVariable(&ptrOutput, 'x');
  if (pretty == PRETTY_PRINT)
  {
    showText("<sup>");
  }
  else if (pretty == TEX)
  {
    showText("^{");
  }
  else
  {      // Pari-GP
    showText("^");
  }
  Bin2Dec(&ptrOutput, primeMod.limbs, primeMod.nbrLimbs, grpLen);
  if (pretty == TEX)
  {
    showText("^{");
  }
  else
  {      // Pari-GP or pretty print.
    showText("^");
  }
  strcpy(ptrOutput, ptrExp);
  ptrOutput += strlen(ptrExp);
  if (pretty == PRETTY_PRINT)
  {
    showText("</sup>");
  }
  else if (pretty == TEX)
  {
    showText("}}");
  }
  else
  {     // Nothing to display for Pari-GP.
  }
  if (pretty == PRETTY_PRINT)
  {
    showText(" &minus; ");
  }
  else
  {
    showText(" - ");
  }
  showVariable(&ptrOutput, 'x');
  showText(")");
}

// Perform distinct degree factorization
static void DistinctDegreeFactorization(int polyDeg)
{
  struct sFactorInfo* pstFactorInfo;
  struct sFactorInfo* pstNewFactorInfo;
  int polyDegree = polyDeg;
  int nbrLimbs = primeMod.nbrLimbs + 1;
  int currentDegree;
  int degreeMin;
  int degreeGcd;
  int* ptrPolyToFactor;
  int* ptrValue1;
  bool showStartOfLoop;
  firstFactor = true;
  // Set poly1 to x.
  int lenBytes = nbrLimbs * (polyDegree + 1) * (int)sizeof(int);
  if (teachMod)
  {
    showText(lang ? "<h3>Factorización de distintos grados</h3>" :
      "<h3>Distinct degree factorization</h3>");
  }
  (void)memset(poly1, 0, lenBytes);
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    poly1[currentDegree * nbrLimbs] = 1;
  }
  NumberLengthR1 = NumberLength;
  SetNumberToOne(&poly1[nbrLimbs]);
  pstFactorInfo = factorInfo;
  for (int nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    if ((pstFactorInfo->degree < 2) || (pstFactorInfo->expectedDegree != 0))
    {             // Polynomial is completely factored. Try next one.
      if (pstFactorInfo->expectedDegree == 0)
      {
        pstFactorInfo->expectedDegree = pstFactorInfo->degree;
      }
      pstFactorInfo++;
      continue;
    }
    ptrPolyToFactor = pstFactorInfo->ptr;
    polyDegree = pstFactorInfo->degree;
    if (teachMod)
    {
      if (firstFactor)
      {
        firstFactor = false;
        copyStr(&ptrOutput, "<ul>");
      }
      copyStr(&ptrOutput, lang? "<li><p>Factorizando f(x)":
        "<li><p>Factoring f(x)");
      showEquAndPoly(ptrPolyToFactor, polyDegree);
    }
    GetPolyInvParm(polyDegree, ptrPolyToFactor);
    // For each loop, raise this polynomial to the primeth power and 
    // then compute the gcd between the polynomial to be factored and
    // the computed polynomial less x. If the degree of GCD is > 0, then the
    // GCD is the product of all factors of degree indicated by currentDegree.
    currentDegree = 1;
    showStartOfLoop = true;
    while ((currentDegree * 2) <= polyDegree)
    {
      if (teachMod && showStartOfLoop)
      {
        char bufExp[100];
        char* ptrBuf;
        showStartOfLoop = false;
        showText(lang ? "<p>Para todos los grados " : "</p>For all degrees ");
        showVariable(&ptrOutput, 'd');
        showText(lang ? " entre " : " between ");
        int2dec(&ptrOutput, currentDegree);
        showText(lang? " y " : " and ");
        int2dec(&ptrOutput, polyDegree / 2);
        showText(lang ? ", el producto de todos los factores de grado " :
          ", the product of all factors of degree ");
        showVariable(&ptrOutput, 'd');
        showText(lang ? " se obtiene calculando " :
          " is found by computing ");
        ptrBuf = bufExp;
        showVariable(&ptrBuf, 'd');
        showGcdDistinctDegFact(bufExp);
        showText("</p>");
      }
#ifdef __EMSCRIPTEN__
      int elapsedTime = (int)(tenths() - originalTenthSecond);
      if ((elapsedTime / 10) != (oldTimeElapsed / 10))
      {
        char* ptrOut = outputText;
        oldTimeElapsed = elapsedTime;
        if (lang)
        {
          copyStr(&ptrOut, "1<p>Factorización de distintos grados módulo ");
          BigInteger2Dec(&ptrOut, &primeMod, 0);
          copyStr(&ptrOut, ": buscando factores de grado ");
          int2dec(&ptrOut, currentDegree);
          copyStr(&ptrOut, " (máx.  ");
          int2dec(&ptrOut, (polyDegree + 1) / 2);
          copyStr(&ptrOut, ") del factor número ");
          int2dec(&ptrOut, nbrFactor + 1);
          copyStr(&ptrOut, " de ");
        }
        else
        {
          copyStr(&ptrOut, "1<p>Distinct degree factorization mod ");
          BigInteger2Dec(&ptrOut, &primeMod, 0);
          copyStr(&ptrOut, ": searching for factors of degree ");
          int2dec(&ptrOut, currentDegree);
          copyStr(&ptrOut, " (max.  ");
          int2dec(&ptrOut, (polyDegree + 1) / 2);
          copyStr(&ptrOut, ") of factor number ");
          int2dec(&ptrOut, nbrFactor + 1);
          copyStr(&ptrOut, " of ");
        }
        int2dec(&ptrOut, nbrFactorsFound);
        copyStr(&ptrOut, lang ? ".</p><p>Transcurrió " : ".</p><p>Time elapsed: ");
        GetDHMS(&ptrOut, elapsedTime / 10);
        copyStr(&ptrOut, "</p>");
        databack(outputText);
      }
#endif
      // Copy polynomial to factor to poly3 and set leading coefficient to 1.
      ptrValue1 = &poly3[polyDegree * nbrLimbs];
      (void)memcpy(poly3, ptrPolyToFactor, (ptrValue1 - &poly3[0]) * sizeof(int));
      SetNumberToOne(ptrValue1);     // Set leading coefficient to 1.
      powerPolynomial(poly1, poly3,  // Base and polynomial modulus.
        polyDegree, &primeMod,       // Degree of polynomials and exponent.
        poly2, NULL,                 // Power and pointer to callback.
        0, 1);
      lenBytes = polyDegree * nbrLimbs * (int)sizeof(int);
      (void)memcpy(poly1, poly2, lenBytes);
      // Subtract x.
      IntArray2BigInteger(&poly2[nbrLimbs], &operand1);
      lenBytes = NumberLength * (int)sizeof(limb);
      (void)memcpy(operand2.limbs, MontgomeryMultR1, lenBytes);
      operand2.nbrLimbs = NumberLengthR1;
      SubtBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
      BigInteger2IntArray(&poly2[nbrLimbs], &operand1);
      // Perform Gcd.
      degreeMin = getDegreePoly(poly2, polyDegree - 1);
      PolyModularGcd(poly3, polyDegree, poly2, degreeMin,
        poly4, &degreeGcd);
      if (teachMod && (degreeGcd > 0))
      {
        char bufExp[100];
        char* ptrBufExp = bufExp;
        int2dec(&ptrBufExp, currentDegree);
        *ptrBufExp = 0;
        showText("<p><strong>");
        showGcdDistinctDegFact(bufExp);
        showEquAndPoly(poly4, degreeGcd);
        showText("</strong></p><p>");
        showText(lang ? "Este polinomio tiene " : "This polynomial has ");
        int2dec(&ptrOutput, degreeGcd / currentDegree);
        if (degreeGcd == currentDegree)
        {
          showText(lang ? " factor irreducible de grado " :
            " irreducible factor of degree ");
        }
        else
        {
          showText(lang ? " factores irreducibles de grado " :
            " irreducible factors of degree ");
        }
        int2dec(&ptrOutput, currentDegree);
        showText("</p>");
        if (degreeGcd != polyDegree)
        {    // New f(x) will have degree greater than zero.
          showText(lang ? "<p>El nuevo valor de " : "<p>The new value of ");
          showVariable(&ptrOutput, 'f');
          showText(lang ? "(x) es el cociente entre " : "(x) is the quotient between ");
          showVariable(&ptrOutput, 'f');
          showText(lang ? "(x) y el mcd anterior.</p>" : "(x) and the previous gcd.</p>");
        }
        showStartOfLoop = true;
      }
      if (degreeGcd == polyDegree)
      {
        pstFactorInfo->expectedDegree = currentDegree;
        polyDegree = 0;
      }
      else if (degreeGcd > 0)
      {         // Non-trivial factor of polynomial has been found.
                // Divide polynomial by GCD. Put the GCD in the first limbs
                // and the quotient in the last limbs.
        ptrValue1 = &poly4[degreeGcd * nbrLimbs];
        SetNumberToOne(ptrValue1);
        DividePolynomial(poly3, polyDegree,
          poly4, degreeGcd, poly2);
        // Quotient located in poly2.
        pstNewFactorInfo = &factorInfo[nbrFactorsFound];
        nbrFactorsFound++;
        pstNewFactorInfo->ptr = ptrPolyToFactor;
        pstNewFactorInfo->degree = degreeGcd;
        pstNewFactorInfo->multiplicity = pstFactorInfo->multiplicity;
        pstNewFactorInfo->expectedDegree = currentDegree;
        pstFactorInfo->degree = polyDegree - degreeGcd;
        pstFactorInfo->ptr = &ptrPolyToFactor[degreeGcd * nbrLimbs];
        lenBytes = degreeGcd * nbrLimbs * (int)sizeof(int);
        (void)memcpy(ptrPolyToFactor, poly4, lenBytes);
        lenBytes = (polyDegree - degreeGcd + 1) * nbrLimbs * (int)sizeof(int);
        (void)memcpy(pstFactorInfo->ptr, poly2, lenBytes);
        polyDegree -= degreeGcd;
        ptrPolyToFactor += degreeGcd * nbrLimbs;
        // Replace poly1 by poly1 mod ptrPolyToFactor
        DividePolynomial(poly1, polyDegree + degreeGcd - 1,
          ptrPolyToFactor, polyDegree, poly2);
        if (polyDegree > 0)
        {
          GetPolyInvParm(polyDegree, ptrPolyToFactor);
        }
      }
      else
      {            // Nothing to do.
      }
      currentDegree++;
    }
    if (polyDegree > 0)
    {
      if (teachMod)
      {
        showText("<p><strong>");
        showVariable(&ptrOutput, 'f');
        showText("(x)");
        showEquAndPoly(ptrPolyToFactor, polyDegree);
        showText("</strong></p><p>");
        showText(lang ? "Este polinomio es irreducible.</p>" :
          "This polynomial is irreducible.</p>");
      }
      pstFactorInfo->expectedDegree = polyDegree;
    }
  }
  if (teachMod && !firstFactor)
  {
    copyStr(&ptrOutput, "</ul>");
  }
}

static void percentageCallback(int percentage)
{
#ifdef __EMSCRIPTEN__
  int elapsedTime = (int)(tenths() - originalTenthSecond);
  char* ptrOut = ptrPercentageOutput;
  if ((elapsedTime / 10) != (oldTimeElapsed / 10))
  {
    oldTimeElapsed = elapsedTime;
    int2dec(&ptrOut, percentage);
    if (lang)
    {
      copyStr(&ptrOut, "% del ");
      int2dec(&ptrOut, attemptNbr);
      copyStr(&ptrOut, ".º intento");
    }
    else
    {
      copyStr(&ptrOut, "% of attempt #");
      int2dec(&ptrOut, attemptNbr);
    }
    copyStr(&ptrOut, lang ? ".</p><p>Transcurrió " : ".</p><p>Time elapsed: ");
    GetDHMS(&ptrOut, elapsedTime / 10);
    copyStr(&ptrOut, "</p>");
    databack(outputText);
  }
#else
  (void)percentage;
#endif
}

static void showStepByStepForEqualDegree(int *poly, int expectedDegree, int polyDegree)
{
  if (firstFactor)
  {
    firstFactor = false;
    copyStr(&ptrOutput, "<ul>");
  }
  showText(lang ? "<li><p>Factorizando en polinomios de grado " :
    "<li><p>Factoring in polynomials of degree ");
  int2dec(&ptrOutput, expectedDegree);
  showText(" ");
  showVariable(&ptrOutput, 'f');
  showText("(x)");
  showEquAndPoly(poly, polyDegree);
  showText(lang ? "<p>Eligiendo " : "<p>Choosing ");
  showVariable(&ptrOutput, 'h');
  showText(lang ? "(x) al azar, siendo " : "(x) at random, let ");
  if (isCharacteristic2)
  { // If prime is 2.
    // g = h + h ^ 2 + h ^ 4 + ... + h ^ 2 ^ (d - 1) mod f(x)
    if (expectedDegree > 3)
    {
      showVariable(&ptrOutput, 'e');
      showText(" = ");
      if (pretty == PRETTY_PRINT)
      {
        showText("2<sup>");
      }
      else if (pretty == TEX)
      {
        showText("2^{");
      }
      else
      {       // Pari-GP.
        showText("2^");
      }
      int2dec(&ptrOutput, expectedDegree - 1);
      if (pretty == PRETTY_PRINT)
      {
        showText("</sup>");
      }
      else if (pretty == TEX)
      {
        showText("}");
      }
      else
      {       // Nothing to show for Pari-GP.
      }
      showText("; ");
    }
    showVariable(&ptrOutput, 'g');
    showText(" = ");
    showVariable(&ptrOutput, 'h');
    showText(" + ");
    showPowerVar(&ptrOutput, 2, 'h');
    if (expectedDegree > 2)
    {
      showText(" + ");
      showPowerVar(&ptrOutput, 4, 'h');
      if (expectedDegree > 3)
      {
        showText(" + ... + ");
        showVariable(&ptrOutput, 'h');
        if (pretty == PRETTY_PRINT)
        {
          showText("<sup>");
          showVariable(&ptrOutput, 'e');
          showText("</sup>");
        }
        else
        {       // TEX or Pari-GP.
          showText("^e");
        }
      }
    }
  }
  else
  {        // Modulus is not power of 2.
    showVariable(&ptrOutput, 'e');
    showText(" = ");
    showText("(");
    Bin2Dec(&ptrOutput, primeMod.limbs,
      primeMod.nbrLimbs, grpLen);
    if (pretty == PRETTY_PRINT)
    {
      showText("<sup>");
    }
    else if (pretty == TEX)
    {
      showText("^{");
    }
    else
    {       // Pari-GP.
      showText("^");
    }
    int2dec(&ptrOutput, expectedDegree);
    if (pretty == PRETTY_PRINT)
    {
      showText("</sup> &minus; ");
    }
    else if (pretty == TEX)
    {
      showText("} - ");
    }
    else
    {       // Pari-GP.
      showText(" - ");
    }
    showText("1) / 2; ");
    showVariable(&ptrOutput, 'g');
    showEqu();
    showVariable(&ptrOutput, 'h');
    if (pretty == PRETTY_PRINT)
    {
      showText("<sup>");
      showVariable(&ptrOutput, 'e');
      showText("</sup> &minus; 1");
    }
    else
    {        // Tex or Pari-GP.
      showText("^e - 1");
    }
  }
  if (pretty == PRETTY_PRINT)
  {
    showText(" (mod ");
    showVariable(&ptrOutput, 'f');
  }
  else if (pretty == TEX)
  {
    showText(" (\\pmod f");
  }
  else
  {
    showText(", f");
  }
  showText(lang ? "), calcular gcd(" : "), then compute gcd(");
  showVariable(&ptrOutput, 'g');
  showText(", ");
  showVariable(&ptrOutput, 'f');
  showText(lang ? ") hasta que el gcd no sea igual a uno de sus argumentos.</p>" :
    ") until the gcd is not equal to one of its arguments.</p>");
}

static void showResultEqualDegreeFact(int degreeGcd, int degreeQuot, int expectedDegree)
{
  showText("<p>");
  if (degreeGcd == expectedDegree)
  {
    showText("<strong>");
  }
  showVariable(&ptrOutput, 'r');
  showEqu();
  showText(lang ? "mcd(" : "gcd(");
  showVariable(&ptrOutput, 'g');
  showText(", ");
  showVariable(&ptrOutput, 'f');
  showText(")");
  showEquAndPoly(poly4, degreeGcd);
  if (degreeGcd == expectedDegree)
  {
    showText("</strong>");
  }
  showText("</p><p>");
  if (degreeQuot == expectedDegree)
  {
    showText("<strong>");
  }
  showVariable(&ptrOutput, 'f');
  showText(" / ");
  showVariable(&ptrOutput, 'r');
  showEquAndPoly(poly2, degreeQuot);
  if (degreeQuot == expectedDegree)
  {
    showText("</strong>");
  }
  showText("</p>");
}

// Perform Cantor-Zassenhaus algorithm to factor polynomials of the same degree.
void SameDegreeFactorization(void)
{
  struct sFactorInfo* pstFactorInfo = factorInfo;
  struct sFactorInfo* pstNewFactorInfo;
  int* ptrValue1;
  int* ptrPolyToFactor;
  int currentDegree;
  int degreeGcd;
  int primeInt = primeMod.limbs[0].x;
  int nbrLimbs = primeMod.nbrLimbs + 1;
  int polyNbr = 1;
  isCharacteristic2 = ((primeMod.nbrLimbs == 1) && (primeMod.limbs[0].x == 2));
  firstFactor = true;
  if (teachMod)
  {
    showText(lang ? "<h3>Factorización del mismo grado</h3>" :
      "<h3>Equal degree factorization</h3>");
  }
  for (int nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    int polyDegree = pstFactorInfo->degree;
    if ((polyDegree < 2) || (polyDegree == pstFactorInfo->expectedDegree) ||
      (pstFactorInfo->expectedDegree == 0))
    {             // Polynomial is completely factored. Try next one.
      pstFactorInfo++;
      continue;
    }
    if (!isCharacteristic2)
    { // If prime is not 2,
      // Calculate operand2 <- (prime^degree-1)/2
      // Use operand1 as temporary variable to store the exponent.
      operand1.limbs[0].x = UintToInt((unsigned int)pstFactorInfo->expectedDegree & MAX_VALUE_LIMB);
      operand1.nbrLimbs = 1;
    }
    ptrPolyToFactor = pstFactorInfo->ptr;
    if (teachMod)
    {
      showStepByStepForEqualDegree(ptrPolyToFactor, pstFactorInfo->expectedDegree,
        polyDegree);
    }
    attemptNbr = 1;
    for (;;)
    {
#ifdef __EMSCRIPTEN__
      char* ptrOutputText = outputText;
      if (lang)
      {
        copyStr(&ptrOutputText, "1<p>Factorización del mismo grado: buscando ");
        int2dec(&ptrOutputText, polyDegree / pstFactorInfo->expectedDegree);
        copyStr(&ptrOutputText, " factores de grado ");
      }
      else
      {
        copyStr(&ptrOutputText, "1<p>Equal degree factorization: searching for ");
        int2dec(&ptrOutputText, polyDegree / pstFactorInfo->expectedDegree);
        copyStr(&ptrOutputText, " factors of degree ");
      }
      int2dec(&ptrOutputText, pstFactorInfo->expectedDegree);
      copyStr(&ptrOutputText, ".</p><p>");
      ptrPercentageOutput = ptrOutputText;
#endif
      // Copy polynomial to factor to poly3 and set leading coefficient to 1.
      // All operations below will be done modulo this polynomial.
      ptrValue1 = &poly3[pstFactorInfo->degree * nbrLimbs];
      (void)memcpy(poly3, ptrPolyToFactor, (ptrValue1 - &poly3[0]) * sizeof(int));
      SetNumberToOne(ptrValue1);  // Set leading coefficient to 1.
      if (attemptNbr == 1)
      {
        GetPolyInvParm(polyDegree, poly3);
      }
      // Initialize polynomial poly1 with different values of coefficients
      // in different iterations.
      ptrValue1 = poly1;
      if (nbrLimbs > 2)
      {    // Coefficient can be any number.
        *ptrValue1 = 1;
        *(ptrValue1 + 1) = polyNbr;
        ptrValue1 += nbrLimbs;
        *ptrValue1 = 1;
        *(ptrValue1 + 1) = 1;
        ptrValue1 += nbrLimbs;
        for (currentDegree = 2; currentDegree < polyDegree; currentDegree++)
        {
          *ptrValue1 = 1;
          *(ptrValue1 + 1) = 0;
          ptrValue1 += nbrLimbs;
        }
      }
      else
      {   // Coefficient range can be from 0 to prime-1.
        int polyCoeff = polyNbr;
        currentDegree = 0;
        while (polyCoeff != 0)
        {
          *ptrValue1 = 1;
          *(ptrValue1 + 1) = polyCoeff % primeInt;
          polyCoeff /= primeInt;
          ptrValue1 += nbrLimbs;
          currentDegree++;
        }
        *ptrValue1 = 1;
        *(ptrValue1 + 1) = 1;
        ptrValue1 += nbrLimbs;
        currentDegree++;
        for (; currentDegree < polyDegree; currentDegree++)
        {
          *ptrValue1 = 1;
          *(ptrValue1 + 1) = 0;
          ptrValue1 += nbrLimbs;
        }
      }
      if (!isCharacteristic2)
      { // If prime is not 2: compute base^((p^d-1)/2).
        int degreeFactor = pstFactorInfo->expectedDegree;
        CopyBigInt(&operand4, &primeMod);
        subtractdivide(&operand4, 1, 2);  // operand4 <- exponent = (p-1)/2.
        // Start by raising base to power (p-1)/2.
        powerPolynomial(poly1, poly3,   // Base and polynomial modulus.
          polyDegree, &operand4,        // Degree of polynomials and exponent.
          poly2, percentageCallback,   // Power and pointer to callback.
          0, degreeFactor);
        // Save base^((p-1)/2) on poly4.
        (void)CopyPolynomialFixedCoeffSize(poly4, poly2, polyDegree, nbrLimbs);
        for (currentDegree = 1; currentDegree < degreeFactor; currentDegree++)
        {
          // Square polynomial and multiply by base to get base^(p^(q-1)).
          multUsingInvPolynomial(poly2, poly2, // Multiplicands.
            poly2, polyDegree,                 // Product and degree of poly.
            poly3);                            // Polynomial modulus.
          multUsingInvPolynomial(poly2, poly1, // Multiplicands.
            poly1, polyDegree,                 // Product and degree of poly.
            poly3);                            // Polynomial modulus.
          CopyBigInt(&operand4, &primeMod);
          subtractdivide(&operand4, 1, 2);  // operand4 <- exponent = (p-1)/2.
          // Raise previous power to exponent (p-1)/2 to get
          // base^(p^(q-1)*(p-1)/2))
          powerPolynomial(poly1, poly3,   // Base and polynomial modulus.
            polyDegree, &operand4,        // Degree of polynomials and exponent.
            poly2, percentageCallback,    // Power and pointer to callback.
            currentDegree, degreeFactor);
          // Compute poly4 as the multiplication of base^(p^(q-1)*(p-1)/2)
          // times base^(p^(q-2)*(p-1)/2).
          multUsingInvPolynomial(poly2, poly4, // Multiplicands.
            poly4, polyDegree,                 // Product and degree of poly.
            poly3);                            // Polynomial modulus.
        }
        (void)CopyPolynomialFixedCoeffSize(poly2, poly4, polyDegree, nbrLimbs);
        // Subtract 1.
        IntArray2BigInteger(&poly2[0], &operand1);
        SubtBigNbrMod(operand1.limbs, MontgomeryMultR1, operand1.limbs);
        BigInteger2IntArray(&poly2[0], &operand1);
      }
      else
      { // If prime is 2, Compute poly2 = T+T^2+T^4+...+T^2^(d-1) mod f(x)
        // where T is the random polynomial.
        // Z <- T mod F.
        int lenBytes = polyDegree * nbrLimbs * (int)sizeof(int);
        (void)memcpy(poly2, poly1, lenBytes);
        for (currentDegree = 1; currentDegree < pstFactorInfo->expectedDegree; currentDegree++)
        {
          multPolynomialModPoly(poly1, poly1, poly1, polyDegree, poly3);
          for (int index = 0; index < polyDegree; index++)
          {
            IntArray2BigInteger(&poly1[index * nbrLimbs], &operand1);
            IntArray2BigInteger(&poly2[index * nbrLimbs], &operand2);
            AddBigNbrMod(operand1.limbs, operand2.limbs, operand1.limbs);
            BigInteger2IntArray(&poly2[index * nbrLimbs], &operand1);
          }
        }
      }
      PolyModularGcd(poly3, polyDegree,
        poly2, getDegreePoly(poly2, polyDegree - 1), poly4, &degreeGcd);
      if ((degreeGcd != 0) && (degreeGcd != polyDegree))
      {   // Non-trivial factor found.
        int lenBytes;
        ptrValue1 = &poly4[degreeGcd * nbrLimbs];
        SetNumberToOne(ptrValue1);
        DividePolynomial(poly3, polyDegree, poly4, degreeGcd, poly2);
        if (teachMod)
        {
          showResultEqualDegreeFact(degreeGcd, polyDegree - degreeGcd,
            pstFactorInfo->expectedDegree);
        }
        // Quotient located in poly2.
        pstNewFactorInfo = &factorInfo[nbrFactorsFound];
        nbrFactorsFound++;
        pstNewFactorInfo->ptr = &ptrPolyToFactor[degreeGcd * nbrLimbs];
        pstNewFactorInfo->degree = polyDegree - degreeGcd;
        pstNewFactorInfo->multiplicity = pstFactorInfo->multiplicity;
        pstNewFactorInfo->expectedDegree = pstFactorInfo->expectedDegree;
        pstFactorInfo->degree = degreeGcd;
        lenBytes = degreeGcd * nbrLimbs * (int)sizeof(int);
        (void)memcpy(ptrPolyToFactor, poly4, lenBytes);
        lenBytes = (polyDegree - degreeGcd) * nbrLimbs * (int)sizeof(int);
        (void)memcpy(pstNewFactorInfo->ptr, poly2, lenBytes);
        polyDegree = degreeGcd;
        attemptNbr = 0;
        if (pstFactorInfo->expectedDegree == pstFactorInfo->degree)
        {
          break;
        }
        GetPolyInvParm(polyDegree, ptrPolyToFactor);
      }
      attemptNbr++;
      polyNbr++;
    }
    pstFactorInfo++;
  }
  if (teachMod && !firstFactor)
  {
    copyStr(&ptrOutput, "</ul>");
  }
}

// Input: values = degree, coefficient degree 0, coefficient degree 1, etc.
// Output: factorInfo = structure that holds the factors.
int FactorModularPolynomial(bool inputMontgomery)
{
  struct sFactorInfo* ptrFactorInfo;
  enum eExprErr rc;
  const int* ptrValue1;
  int lenBytes;
  int nbrLimbsPrime = primeMod.nbrLimbs + 1; // Add 1 for length.
  degree = values[0];
  ptrValue1 = &values[1];
  for (int currentDegree = 0; currentDegree <= degree; currentDegree++)
  {
    NumberLength = numLimbs(ptrValue1);
    IntArray2BigInteger(ptrValue1, &operand1);
    lenBytes = (powerMod.nbrLimbs - NumberLength) * (int)sizeof(limb);
    (void)memset(&operand1.limbs[NumberLength], 0, lenBytes);
    NumberLength = powerMod.nbrLimbs;
    if (inputMontgomery)
    {
      // Convert from Montgomery to standard notation.
      operand2.limbs[0].x = 1;
      if (NumberLength > 1)
      {
        lenBytes = (NumberLength - 1) * (int)sizeof(limb);
        (void)memset(&operand2.limbs[1], 0, lenBytes);
      }
      modmult(operand1.limbs, operand2.limbs, operand1.limbs);
    }
    rc = BigIntRemainder(&operand1, &primeMod, &operand1);
    if (rc != EXPR_OK)
    {
      return rc;
    }
    NumberLength = primeMod.nbrLimbs;
    BigInteger2IntArray(&valuesPrime[currentDegree * nbrLimbsPrime], &operand1);
    ptrValue1 += numLimbs(ptrValue1);
    ptrValue1++;
  }
  if (BigIntIsZero(&operand1))
  {
    return EXPR_LEADING_COFF_MULTIPLE_OF_PRIME;
  }
  lenBytes = primeMod.nbrLimbs * (int)sizeof(limb);
  (void)memcpy(&TestNbr, primeMod.limbs, lenBytes);
  NumberLength = primeMod.nbrLimbs;
  TestNbr[NumberLength].x = 0;
  GetMontgomeryParms(primeMod.nbrLimbs);
  // Convert polynomial mod prime to monic (leading coefficient must be 1).
  ConvertToMonic(valuesPrime, degree);
  if (teachMod)
  {
    showText(lang ? "<p>Dividiendo el polinomio por el coeficiente principal:" :
      "<p>Dividing the polynomial by the leading coefficient:");
    showText("</p><p><var>f</var>(x)");
    showEquAndPoly(valuesPrime, degree);
    showText("</p>");
  }
  // Perform square free factorization.
  ptrOrigPoly = valuesPrime;
  degreeOrigPoly = degree;
  nbrFactorsFound = 0;
  if (degree != 0)
  {
    SquareFreeFactorization(degree, valuesPrime);
    DistinctDegreeFactorization(degree);
    if (inputMontgomery)
    {    // Do not perform same degree factorization if only counting
         // the number of modular factors of input polynomial.
      SameDegreeFactorization();
    }
  }
  if (inputMontgomery == false)
  {
    return EXPR_OK;
  }
  // Convert original polynomial from Montgomery to standard notation.
  ptrFactorInfo = factorInfo;
  for (int nbrFactor = 0; nbrFactor < nbrFactorsFound; nbrFactor++)
  {
    polyToStandardNotation(ptrFactorInfo->ptr, ptrFactorInfo->degree);
    ptrFactorInfo++;
  }
  OrigPolyFromMontgomeryToStandard();
  rc = HenselLifting(factorInfo, false);
  if (rc != EXPR_OK)
  {
    return rc;
  }
  if (teachMod)
  {
    showText(lang ? "<h3>Lista de factores</h3>": "<h3>List of factors</h3>");
  }
  SortFactors(&powerMod);
  return rc;
}
