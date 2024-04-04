//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2024 Dario Alejandro Alpern
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
#include <math.h>
#include <assert.h>
#include "rootseq.h"
#include "expression.h"

extern BigInteger Linear;
extern BigInteger Independent;
extern char* ptrOutput;

// Compute x = -c_0 / c_1
void LinearEquation(const int* polynomial, int multiplicity)
{
  const int* ptrPolynomial = polynomial;
  UncompressBigIntegerB(ptrPolynomial, &Independent);
  ptrPolynomial += numLimbs(ptrPolynomial);
  ptrPolynomial++;
  UncompressBigIntegerB(ptrPolynomial, &Linear);
  CopyBigInt(&Rat1.numerator, &Independent);
  CopyBigInt(&Rat1.denominator, &Linear);
  if (teach && !BigIntIsOne(&Linear))
  {
    showText(lang ? "<p>Dividiendo la ecuación por el coeficiente lineal:</p><p>" :
      "<p>Dividing the equation by the linear coefficient:</p><p>");
    showVariable(&ptrOutput, 'x');
    showPlusMinusRational(&Rat1);
    showText(" = 0</p>");
  }
  BigIntChSign(&Rat1.numerator);
  ForceDenominatorPositive(&Rat1);
  showX(multiplicity);
  showRationalNoParen(&Rat1);
  endShowX();
}

static void ShowRootsOfRationalNumbers(int polyDegree, int multiplicity)
{
  char rationalRoot[30000];
  char* ptrOutputBak = ptrOutput;
  bool isNegative = false;
  if (Rat1.numerator.sign == Rat1.denominator.sign)
  {
    isNegative = true;
  }
  Rat1.numerator.sign = SIGN_POSITIVE;
  Rat1.denominator.sign = SIGN_POSITIVE;
  ptrOutput = rationalRoot;
  showRational(&Rat1);
  *ptrOutput = 0;
  ptrOutput = ptrOutputBak;
  // Generate roots.
  GenerateRoots(multiplicity, rationalRoot, isNegative, polyDegree);
}

// If polynomial has the form ax^n + b, show its roots.
bool isLinearExponential(const int* ptrPolynomial, int polyDegree, int multiplicity)
{
  const int* ptrPoly = ptrPolynomial;
  // Get constant term.
  NumberLength = *ptrPoly;
  IntArray2BigInteger(ptrPoly, &Rat1.numerator);
  ptrPoly += numLimbs(ptrPoly);
  ptrPoly++;
  for (int currentDegree = 1; currentDegree < polyDegree; currentDegree++)
  {
    if ((*ptrPoly != 1) || (*(ptrPoly + 1) != 0))
    {
      return false;    // Polynomial does not have format ax^n + b.
    }
    ptrPoly += numLimbs(ptrPoly);
    ptrPoly++;
  }
  // Get leading coefficient.
  NumberLength = *ptrPoly;
  IntArray2BigInteger(ptrPoly, &Rat1.denominator);
  ShowRootsOfRationalNumbers(polyDegree, multiplicity);
  return true;
}


