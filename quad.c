//
// This file is part of Alpertron Calculators.
//
// Copyright 2017-2021 Dario Alejandro Alpern
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "commonstruc.h"
#include "quadmodLL.h"

enum eLinearSolution
{
  SOLUTION_FOUND = 0,
  NO_SOLUTIONS,
  INFINITE_SOLUTIONS,
};

enum eShowSolution
{
  ONE_SOLUTION = 0,
  TWO_SOLUTIONS,
  FIRST_SOLUTION,
  SECOND_SOLUTION,
};

enum eCallbackQuadModType
{
  CBACK_QMOD_PARABOLIC = 0,
  CBACK_QMOD_ELLIPTIC,
  CBACK_QMOD_HYPERBOLIC,
};

static enum eShowSolution showSolution;
static enum eCallbackQuadModType callbackQuadModType;
static int counters[400];
static char isDescending[400];
static int originalMultiplicities[400];
static BigInteger Aux[400];
static BigInteger ValA;
static BigInteger ValB;
static BigInteger ValC;
static BigInteger ValD;
static BigInteger ValE;
static BigInteger ValF;
static BigInteger ValH;
static BigInteger ValI;
static BigInteger ValL;
static BigInteger ValM;
static BigInteger ValN;
static BigInteger ValO;
static BigInteger ValP;
static BigInteger ValQ;
static BigInteger ValU;
static BigInteger ValV;
static BigInteger ValG;
static BigInteger ValR;
static BigInteger ValS;
static BigInteger ValJ;
static BigInteger ValK;
static BigInteger ValZ;
static BigInteger ValAlpha;
static BigInteger ValBeta;
static BigInteger ValDen;
static BigInteger ValDiv;
static BigInteger ValABak;
static BigInteger ValBBak;
static BigInteger ValCBak;
static BigInteger ValUBak;
static BigInteger ValVBak;
static BigInteger ValGcdHomog;
static BigInteger Tmp1;
static BigInteger Tmp2;
static int SolNbr;
static int showRecursiveSolution;
static BigInteger Xind;
static BigInteger Yind;
static BigInteger Xlin;
static BigInteger Ylin;
static int nbrFactors;
static bool solFound;
static bool originTranslated;
static char also;
static bool ExchXY;
static char *ptrOutput;
static char *divgcd;
static const char *varT = "t";
static const char *squareText = "&sup2;";
static BigInteger discr;
static BigInteger U1;
static BigInteger U2;
static BigInteger U3;
static BigInteger V1;
static BigInteger V2;
static BigInteger V3;
static BigInteger bigTmp;
static BigInteger startPeriodU;
static BigInteger startPeriodV;
static BigInteger coeffQuadr;
static BigInteger coeffLinear;
static BigInteger coeffIndep;
static BigInteger modulus;
static BigInteger currentFactor;
static BigInteger Xplus;
static BigInteger Xminus;
static BigInteger Yplus;
static BigInteger Yminus;
static BigInteger *Xbak;
static BigInteger *Ybak;
static int Show(const BigInteger *num, const char *str, enum eLinearSolution t);
static void Show1(const BigInteger *num, enum eLinearSolution t);
static void callbackQuadModParabolic(const BigInteger *value);
static void callbackQuadModElliptic(BigInteger *value);
static void callbackQuadModHyperbolic(BigInteger *value);
static void recursiveSolution(void);
static BigInteger Tmp[13];
static int indexEvenMultiplicity[400];
static int nbrPrimesEvenMultiplicity;
static int equationNbr;
static int contfracEqNbr;
static char positiveDenominator;
static const char *ptrVarNameX;
static const char *ptrVarNameY;
static const char *varX;
static const char *varY;
static const char *varXnoTrans;
static const char *varYnoTrans;
static bool firstSolutionX;
extern BigInteger LastModulus;
extern BigInteger prime;
extern bool teach;

#define NBR_COEFF 6
struct stValidateCoeff
{
  char *expression;
  BigInteger *bigint;
  char *textSpanish;
  char *textEnglish;
};
static struct stValidateCoeff astValidateCoeff[NBR_COEFF] =
{
  { NULL, &ValA, "Coeficiente <var>a</var>: ", "Coefficient <var>a</var>: " },
  { NULL, &ValB, "Coeficiente <var>b</var>: ", "Coefficient <var>b</var>: " },
  { NULL, &ValC, "Coeficiente <var>c</var>: ", "Coefficient <var>c</var>: " },
  { NULL, &ValD, "Coeficiente <var>d</var>: ", "Coefficient <var>d</var>: " },
  { NULL, &ValE, "Coeficiente <var>e</var>: ", "Coefficient <var>e</var>: " },
  { NULL, &ValF, "Coeficiente <var>f</var>: ", "Coefficient <var>f</var>: " },
};

static void showText(const char *text)
{
  copyStr(&ptrOutput, text);
}

static void showMinus(void)
{
  showText("&minus;");
}

static void shownbr(const BigInteger* value)
{
  BigInteger2Dec(&ptrOutput, value, groupLen);
}

static void showInt(int value)
{
  int2dec(&ptrOutput, value);
}

static void showSquare(void)
{
  showText("&sup2;");
}

static void showAlso(void)
{
  if (also && !teach)
  {
    showText(lang? "y también:<br>": "and also:<br>");
  }
  else
  {
    also = 1;
  }
}

static void startResultBox(enum eLinearSolution Ret)
{
  if (teach && (Ret != NO_SOLUTIONS))
  {
    showText("<div class=\"outerbox\"><div class=\"box\">");
  }
}

static void endResultBox(enum eLinearSolution Ret)
{
  if (teach && (Ret != NO_SOLUTIONS))
  {
    showText("</div></div>");
  }
}

static void showEqNbr(int nbr)
{
  showText(" <span class=\"eq\">(");
  showInt(nbr);
  showText(")</span>");
}

static void ShowLin(const BigInteger *coeffX, const BigInteger *coeffY, const BigInteger *coeffInd,
  const char *x, const char *y)
{
  enum eLinearSolution t;
  t = Show(coeffX, x, SOLUTION_FOUND);
  t = Show(coeffY, y, t);
  Show1(coeffInd, t);
}

static void ShowLinInd(const BigInteger *lin, const BigInteger *ind, const char *var)
{
  if (BigIntIsZero(ind) && BigIntIsZero(lin))
  {
    showText("0");
  }
  if (!BigIntIsZero(ind))
  {
    shownbr(ind);
  }
  *ptrOutput = ' ';
  ptrOutput++;
  if (lin->sign == SIGN_NEGATIVE)
  {
    showMinus();
  }
  else if (!BigIntIsZero(lin) && !BigIntIsZero(ind))
  {
    *ptrOutput = '+';
    ptrOutput++;
  }
  else
  {            // Nothing to do.
  }
  *ptrOutput = ' ';
  ptrOutput++;
  if (!BigIntIsZero(lin))
  {
    if ((lin->nbrLimbs != 1) || (lin->limbs[0].x != 1))
    {     // abs(lin) is not 1
      CopyBigInt(&Aux[0], lin);
      Aux[0].sign = SIGN_POSITIVE;   // Do not show negative sign twice.
      shownbr(&Aux[0]);
    }
    *ptrOutput = ' ';
    ptrOutput++;
    showText(var);
  }
}

static void PrintLinear(enum eLinearSolution Ret, const char *var)
{
  if (Ret == NO_SOLUTIONS)
  {
    return;
  }
  if ((var == varT) && !teach)
  {
    showAlso();
  }
  if (Ret == INFINITE_SOLUTIONS)
  {
    showText(lang?"<p>x, y: caulquier entero</p>": "<p>x, y: any integer</p>");
    return;
  }
  if (ExchXY)
  {
    // Exchange Xind and Yind
    CopyBigInt(&bigTmp, &Xind);
    CopyBigInt(&Xind, &Yind);
    CopyBigInt(&Yind, &bigTmp);
    // Exchange Xlin and Ylin
    CopyBigInt(&bigTmp, &Xlin);
    CopyBigInt(&Xlin, &Ylin);
    CopyBigInt(&Ylin, &bigTmp);
  }
  showText("<p>x = ");
  ShowLinInd(&Xlin, &Xind, var);
  showText("<br>y = ");
  ShowLinInd(&Ylin, &Yind, var);
  showText("</p>");
  return;
}

static void PrintQuad(const BigInteger *coeffT2, const BigInteger *coeffT, const BigInteger *coeffInd, 
                      const char *var1, const char *var2)
{
  if ((coeffT2->nbrLimbs == 1) && (coeffT2->limbs[0].x == 1))
  {             // abs(coeffT2) = 1
    if (coeffT2->sign == SIGN_POSITIVE)
    {           // coeffT2 = 1
      *ptrOutput = ' ';
      ptrOutput++;
    }
    else
    {           // coeffT2 = -1
      showMinus();
    }
    showText(var1);
    showSquare();
  }
  else if (!BigIntIsZero(coeffT2))
  {             // coeffT2 is not zero.
    shownbr(coeffT2);
    *ptrOutput = ' ';
    ptrOutput++;
    showText(var1);
    showSquare();
  }
  else
  {              // Nothing to do.
  }
  if (coeffT->sign == SIGN_NEGATIVE)
  {
    showText(" &minus; ");
  }
  else if (!BigIntIsZero(coeffT) && !BigIntIsZero(coeffT2))
  {
    showText(" + ");
  }
  else
  {              // Nothing to do.
  }
  if ((coeffT->nbrLimbs == 1) && (coeffT->limbs[0].x == 1))
  {     // abs(coeffT) = 1
    showText(var1);
    showText("&#8290;");
    if (var2 != NULL)
    {
      showText(var2);
    }
    showText(" ");
  }
  else if (!BigIntIsZero(coeffT))
  {
    if (coeffT->sign == SIGN_NEGATIVE)
    {
      Bin2Dec(&ptrOutput, coeffT->limbs, coeffT->nbrLimbs, groupLen);
    }
    else
    {
      shownbr(coeffT);
    }
    showText(" ");
    showText(var1);
    if (var2 != NULL)
    {
      showText("&#8290;");
      showText(var2);
    }
  }
  else
  {           // Nothing to do.
  }
  if (!BigIntIsZero(coeffInd))
  {
    if (!BigIntIsZero(coeffT) || !BigIntIsZero(coeffT2))
    {
      if (coeffInd->sign == SIGN_NEGATIVE)
      {
        showText(" &minus; ");
      }
      else
      {
        showText(" + ");
      }
    }
    else if (coeffInd->sign == SIGN_NEGATIVE)
    {
      showText(" &minus;");
    }
    else
    {           // Nothing to do.
    }
    if (var2 == NULL)
    {
      Bin2Dec(&ptrOutput, coeffInd->limbs, coeffInd->nbrLimbs, groupLen);
    }
    else if ((coeffInd->nbrLimbs > 1) || (coeffInd->limbs[0].x > 1))
    {
      Bin2Dec(&ptrOutput, coeffInd->limbs, coeffInd->nbrLimbs, groupLen);
      showText("&nbsp;&#8290;");
      showText(var2);
      showSquare();
    }
    else
    {
      showText(var2);
      showSquare();
    }
  }
}

static void showSolutionXY(const BigInteger *X, const BigInteger *Y)
{
  showText("<p>x = ");
  shownbr(ExchXY ? Y : X);
  showText("<BR>y = ");
  shownbr(ExchXY ? X : Y);
  showText("</p>");
}

void ShowXY(BigInteger *X, BigInteger *Y)
{
  if (showSolution == TWO_SOLUTIONS)
  {
    enum eSign signX;
    enum eSign signY;
    solFound = true;
    if (Xbak->nbrLimbs == 0)
    {
      CopyBigInt(Xbak, X);
      CopyBigInt(Ybak, Y);
      return;
    }
    // Use the lowest of |X| + |Y| and |Xbak| + |Ybak|
    signX = Xbak->sign;
    signY = Ybak->sign;
    Xbak->sign = SIGN_POSITIVE;
    Ybak->sign = SIGN_POSITIVE;
    BigIntAdd(Xbak, Ybak, &bigTmp);
    Xbak->sign = signX;
    Ybak->sign = signY;
    signX = X->sign;
    signY = Y->sign;
    X->sign = SIGN_POSITIVE;
    Y->sign = SIGN_POSITIVE;
    BigIntSubt(&bigTmp, X, &bigTmp);
    BigIntSubt(&bigTmp, Y, &bigTmp);
    X->sign = signX;
    Y->sign = signY;
    if (bigTmp.sign == SIGN_POSITIVE)
    {       // At this moment |x| + |y| <= |xbak| + |ybak|
      CopyBigInt(Xbak, X);
      CopyBigInt(Ybak, Y);
    }
  }
  else
  {
    // ONE_SOLUTION: Show it.
    showAlso();
    showSolutionXY(X, Y);
  }
}

static int Show(const BigInteger *num, const char *str, enum eLinearSolution t)
{
  enum eLinearSolution tOut = t;
  if (!BigIntIsZero(num))
  {     // num is not zero.
    if ((t == NO_SOLUTIONS) && (num->sign == SIGN_POSITIVE))
    {
      *ptrOutput = ' ';
      ptrOutput++;
      *ptrOutput = '+';
      ptrOutput++;
    }
    if (num->sign == SIGN_NEGATIVE)
    {
      *ptrOutput = ' ';
      ptrOutput++;
      *ptrOutput = '-';
      ptrOutput++;
    }
    if ((num->nbrLimbs != 1) || (num->limbs[0].x != 1))
    {    // num is not 1 or -1.
      *ptrOutput = ' ';
      ptrOutput++;
      Bin2Dec(&ptrOutput, num->limbs, num->nbrLimbs, groupLen);
      copyStr(&ptrOutput, "&nbsp;&#8290;");
    }
    else
    {
      copyStr(&ptrOutput, "&nbsp;");
    }
    showText(str);
    if (t == SOLUTION_FOUND)
    {
      tOut = NO_SOLUTIONS;
    }
  }
  return tOut;
}

static void Show1(const BigInteger *num, enum eLinearSolution t)
{
  int u = Show(num, "", t);
  *ptrOutput = ' ';
  ptrOutput++;
  if (((u & 1) == 0) || ((num->nbrLimbs == 1) && (num->limbs[0].x == 1)))
  {          // Show absolute value of num.
    Bin2Dec(&ptrOutput, num->limbs, num->nbrLimbs, groupLen);
  }
}

static void ShowEq(const BigInteger *coeffA, const BigInteger *coeffB, const BigInteger *coeffC, 
  const BigInteger *coeffD, const BigInteger *coeffE, const BigInteger *coeffF,
  const char *x, const char *y)
{
  char var[30];
  char *ptrVar;
  enum eLinearSolution t;
  ptrVar = var;
  copyStr(&ptrVar, x);
  copyStr(&ptrVar, squareText);
  t = Show(coeffA, var, SOLUTION_FOUND);

  ptrVar = var;
  copyStr(&ptrVar, x);
  copyStr(&ptrVar, "&#8290;");
  copyStr(&ptrVar, y);
  t = Show(coeffB, var, t);

  ptrVar = var;
  copyStr(&ptrVar, y);
  copyStr(&ptrVar, squareText);
  t = Show(coeffC, var, t);

  t = Show(coeffD, x, t);

  t = Show(coeffE, y, t);
  Show1(coeffF, t);
}

static void showFactors(const BigInteger *value)
{
  const struct sFactors *pstFactor;
  int numFactors = astFactorsMod[0].multiplicity;
  bool factorShown = false;
  shownbr(value);
  showText(" = ");
  if (value->sign == SIGN_NEGATIVE)
  {
    showMinus();
  }
  pstFactor = &astFactorsMod[1];
  for (int index = 0; index < numFactors; index++)
  {
    NumberLength = *pstFactor->ptrFactor;
    IntArray2BigInteger(pstFactor->ptrFactor, &prime);
    if (pstFactor->multiplicity == 0)
    {
      pstFactor++;
      continue;
    }
    if (factorShown)
    {
      showText(" &times; ");
    }
    shownbr(&prime);
    if (pstFactor->multiplicity != 1)
    {
      showText("<sup>");
      showInt(pstFactor->multiplicity);
      showText("</sup>");
    }
    factorShown = true;
    pstFactor++;
  }
  if (!factorShown)
  {      // No factor shown. Show 1.
    showText("1");
  }
}

static void showValue(BigInteger* value)
{
  SolNbr++;
  // If 2*value is greater than modulus, subtract modulus.
  multint(&Tmp1, value, 2);
  BigIntSubt(&modulus, &Tmp1, &Tmp1);
  if (Tmp1.sign == SIGN_NEGATIVE)
  {
    BigIntSubt(value, &modulus, value);
  }
  if (teach)
  {
    if (firstSolutionX)
    {
      showText("<ol>");
      firstSolutionX = false;
    }
    showText("<li><var>T</var> = ");
    BigInteger2Dec(&ptrOutput, value, groupLen);
    showText("</li>");
  }
}

static void SolutionX(BigInteger *value)
{
  showValue(value);
  if (callbackQuadModType == CBACK_QMOD_PARABOLIC)
  {
    callbackQuadModParabolic(value);
  }
  else if (callbackQuadModType == CBACK_QMOD_ELLIPTIC)
  {
    callbackQuadModElliptic(value);
  }
  else if (callbackQuadModType == CBACK_QMOD_HYPERBOLIC)
  {
    callbackQuadModHyperbolic(value);
  }
  else
  {         // Nothing to do.
  }
}

static void NoSolsModPrime(int expon)
{
  if (teach)
  {
    showText(lang ? "<p>No hay soluciones módulo " : "<p>There are no solutions modulo ");
    shownbr(&prime);
    if (expon != 1)
    {
      showText("<sup>");
      int2dec(&ptrOutput, expon);
      showText("</sup>");
    }
    showText(lang ? ", así que la ecuación modular no tiene soluciones.</p>" :
      ", so the modular equation does not have any solution.</p>");
  }
}

static void ShowSolutionsModPrime(int factorIndex, int expon, const BigInteger *pIncrement)
{
  bool last;
  static BigInteger primePower;
  if (!teach)
  {
    return;
  }
  (void)BigIntPowerIntExp(&prime, expon, &primePower);       // Get value of prime power.
  intToBigInteger(&ValH, 0);
  showText(lang ? "<p>Soluciones módulo " : "<p>Solutions modulo ");
  shownbr(&prime);
  if (expon != 1)
  {
    showText("<sup>");
    int2dec(&ptrOutput, expon);
    showText("</sup>");
  }
  showText(": ");
  do
  {
    bool oneSolution = BigIntEqual(&common.quad.Solution1[factorIndex],
      &common.quad.Solution2[factorIndex]);
    BigIntAdd(&ValH, pIncrement, &ValI);   // Next value.
    last = BigIntEqual(&primePower, &ValI);
    if (!BigIntIsZero(&ValH))
    {
      if (last && oneSolution)
      {
        showText(lang ? " y " : " and ");
      }
      else
      {
        showText(", ");
      }
    }
    BigIntAdd(&ValH, &common.quad.Solution1[factorIndex], &ValJ);
    shownbr(&ValJ);
    if (!oneSolution)
    {
      if (last)
      {
        showText(lang ? " y " : " and ");
      }
      else
      {
        showText(", ");
      }
      BigIntAdd(&ValH, &common.quad.Solution2[factorIndex], &ValJ);
      shownbr(&ValJ);
    }
    CopyBigInt(&ValH, &ValI);
  } while (!last);
  showText("</p>");
}

// Solve congruence an^2 + bn + c = 0 (mod n) where n is different from zero.
void SolveQuadModEquation(void)
{
  static BigInteger GcdAll;
  static BigInteger ValNn;
  static BigInteger z;

  firstSolutionX = true;
  modulus.sign = SIGN_POSITIVE;
  (void)BigIntRemainder(&coeffQuadr, &modulus, &coeffQuadr);
  if (coeffQuadr.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&coeffQuadr, &modulus, &coeffQuadr);
  }
  (void)BigIntRemainder(&coeffLinear, &modulus, &coeffLinear);
  if (coeffLinear.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&coeffLinear, &modulus, &coeffLinear);
  }
  (void)BigIntRemainder(&coeffIndep, &modulus, &coeffIndep);
  if (coeffIndep.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&coeffIndep, &modulus, &coeffIndep);
  }
  BigIntGcd(&coeffQuadr, &coeffLinear, &Tmp[0]);
  BigIntGcd(&coeffIndep, &Tmp[0], &GcdAll);
  (void)BigIntRemainder(&coeffIndep, &GcdAll, &Tmp[0]);
  if (!BigIntIsZero(&Tmp[0]))
  {  // ValC must be multiple of gcd(ValA, ValB).
     // Otherwise go out because there are no solutions.
    return;
  }
  BigIntGcd(&modulus, &GcdAll, &Tmp[0]);
  CopyBigInt(&GcdAll, &Tmp[0]);
  // Divide all coefficients by gcd(ValA, ValB).
  (void)BigIntDivide(&coeffQuadr, &GcdAll, &coeffQuadr);
  (void)BigIntDivide(&coeffLinear, &GcdAll, &coeffLinear);
  (void)BigIntDivide(&coeffIndep, &GcdAll, &coeffIndep);
  (void)BigIntDivide(&modulus, &GcdAll, &modulus);
  CopyBigInt(&ValNn, &modulus);
  if ((ValNn.nbrLimbs == 1) && (ValNn.limbs[0].x == 1))
  {     // All values from 0 to GcdAll - 1 are solutions.
    if ((GcdAll.nbrLimbs > 1) || (GcdAll.limbs[0].x > 5))
    {
      showText("<p>All values of <var>x</var> between 0 and ");
      addbigint(&GcdAll, -1);
      BigInteger2Dec(&ptrOutput, &GcdAll, groupLen);
      showText(" are solutions.</p>");
    }
    else
    {
      for (int ctr = 0; ctr < GcdAll.limbs[0].x; ctr++)
      {
        intToBigInteger(&Tmp[0], ctr);
        SolutionX(&Tmp[0]);
      }
      if (teach && !firstSolutionX)
      {
        showText("</ol>");
      }
    }
    return;
  }
  (void)BigIntRemainder(&coeffQuadr, &modulus, &Tmp[0]);
  if (BigIntIsZero(&Tmp[0]))
  {           // Linear equation.
    int lenBytes;
    BigIntGcd(&coeffLinear, &modulus, &Tmp[0]);
    if ((Tmp[0].nbrLimbs != 1) || (Tmp[0].limbs[0].x != 1))
    {         // ValB and ValN are not coprime. Go out.
      return;
    }
    // Calculate z <- -ValC / ValB (mod ValN)
    NumberLength = modulus.nbrLimbs;
    lenBytes = NumberLength * (int)sizeof(limb);
    (void)memcpy(TestNbr, modulus.limbs, lenBytes);
    TestNbr[NumberLength].x = 0;
    GetMontgomeryParms(NumberLength);
    BigIntModularDivision(&coeffIndep, &coeffLinear, &modulus, &z);
    if (!BigIntIsZero(&z))
    {
      BigIntSubt(&ValN, &z, &z);
    }
    (void)BigIntMultiply(&ValNn, &GcdAll, &Tmp[0]);
    do
    {
      SolutionX(&z);
      BigIntAdd(&z, &modulus, &z);
      BigIntSubt(&z, &Tmp[0], &Tmp[1]);
    } while (Tmp[1].sign == SIGN_NEGATIVE);
    if (teach && !firstSolutionX)
    {
      showText("</ol>");
    }
    return;
  }
  if (callbackQuadModType == CBACK_QMOD_PARABOLIC)
  {    // For elliptic case, the factorization is already done.
    char* ptrFactorDec;
    NumberLength = modulus.nbrLimbs;
    BigInteger2IntArray(nbrToFactor, &modulus);
    ptrFactorDec = tofactorDec;
    Bin2Dec(&ptrFactorDec, modulus.limbs, modulus.nbrLimbs, groupLen);
    factor(&modulus, nbrToFactor, factorsMod, astFactorsMod);
    CopyBigInt(&LastModulus, &modulus);           // Do not factor again same modulus.
    if (teach)
    {
      showText(lang ? "<p>Para resolver esta ecuación cuadrática modular debemos factorizar el módulo y hallar las soluciones módulo las potencias de los factores primos. Luego debemos combinar estas soluciones usando el teorema chino del resto.</p>" :
        "<p>To solve this quadratic modular equation we have to factor the modulus and find the solution modulo the powers of the prime factors. Then we combine them by using the Chinese Remainder Theorem.</p>");
      showFactors(&modulus);
    }
  }
  SetCallbacksForSolveEquation(SolutionX, ShowSolutionsModPrime, NoSolsModPrime);
  SolveEquation(&coeffQuadr, &coeffLinear, &coeffIndep, &modulus,
    &GcdAll, &ValNn);
  if (teach && !firstSolutionX)
  {
    showText("</ol>");
  }
}

static void paren(const BigInteger *num)
{
  if (num->sign == SIGN_NEGATIVE)
  {
    *ptrOutput = '(';
    ptrOutput++;
    shownbr(num);
    *ptrOutput = ')';
    ptrOutput++;
    *ptrOutput = 0;
  }
  else
  {
    BigInteger2Dec(&ptrOutput, num, groupLen);
  }
}

enum eLinearSolution LinearEq(BigInteger* coeffX, BigInteger* coeffY, BigInteger* coeffInd)
{
  BigInteger q;
  bool showSteps;
  int stepNbr;
  if (teach)
  {
    showText(lang ? "<p>Esta es una ecuación lineal " : "<p>This is a linear equation ");
    ShowLin(coeffX, coeffY, coeffInd, "x", "y");
    showText(" = 0</p>");
  }
  if (BigIntIsZero(coeffX))
  {
    if (BigIntIsZero(coeffY))
    {
      if (!BigIntIsZero(coeffInd))
      {
        return NO_SOLUTIONS;           // No solutions
      }
      else
      {
        return INFINITE_SOLUTIONS;     // Infinite number of solutions
      }
    }
    (void)BigIntRemainder(coeffInd, coeffY, &Aux[0]);
    if (!BigIntIsZero(&Aux[0]))
    {
      return NO_SOLUTIONS;             // No solutions
    }
    else
    {
      intToBigInteger(&Xind, 0);
      intToBigInteger(&Xlin, 1);
      (void)BigIntDivide(coeffInd, coeffY, &Yind);
      BigIntNegate(&Yind, &Yind);
      intToBigInteger(&Ylin, 0);
      return SOLUTION_FOUND;           // Solution found
    }
  }
  if (BigIntIsZero(coeffY))
  {
    (void)BigIntDivide(coeffInd, coeffX, &Xind);
    (void)BigIntMultiply(coeffX, &Xind, &Aux[0]);
    if (!BigIntEqual(coeffInd, &Aux[0]))
    {
      return NO_SOLUTIONS;             // No solutions
    }
    else
    {
      intToBigInteger(&Yind, 0);
      intToBigInteger(&Ylin, 1);
      BigIntNegate(&Xind, &Xind);
      intToBigInteger(&Xlin, 0);
      return SOLUTION_FOUND;           // Solution found
    }
  }
  BigIntGcd(coeffX, coeffY, &U1);
  if ((U1.nbrLimbs > 1) || (U1.limbs[0].x != 1))
  {                  // GCD is not 1.
    if (teach)
    {
      showText(lang ? "<p>Para resolverla, debemos hallar el máximo común divisor de los coeficientes lineales, esto es: mcd(" :
        "<p>To solve it, we first find the greatest common divisor of the linear coefficients, that is: gcd(");
      shownbr(coeffX);
      showText(", ");
      shownbr(coeffY);
      showText(") = ");
      shownbr(&U1);
      showText(".</p>");
    }
    (void)BigIntRemainder(coeffInd, &U1, &U2);
    if (!BigIntIsZero(&U2))
    {
      if (teach)
      {
        showText("<p>The independent coefficient is not multiple of gcd, so there are no solutions.</p>");
      }
      return NO_SOLUTIONS;            // No solutions
    }
    (void)BigIntDivide(coeffX, &U1, coeffX);   // Divide all coefficients by the gcd.
    (void)BigIntDivide(coeffY, &U1, coeffY);
    (void)BigIntDivide(coeffInd, &U1, coeffInd);
  }
  if (teach)
  {
    showText(divgcd);
    ShowLin(coeffX, coeffY, coeffInd, "x", "y");
    showText(" = 0</p>");
    showText(lang ? "<p>Ahora debemos aplicar el algoritmo generalizado de Euclides:</p>" :
      "<p>Now we must apply the Generalized Euclidean algorithm:</p>");
  }
  intToBigInteger(&U1, 1);    // U1 <- 1
  intToBigInteger(&U2, 0);    // U2 <- 0
  CopyBigInt(&U3, coeffX);    // U3 <- coeffX
  intToBigInteger(&V1, 0);    // V1 <- 0
  intToBigInteger(&V2, 1);    // V2 <- 1
  CopyBigInt(&V3, coeffY);    // V3 <- coeffY
  showSteps = ((coeffX->nbrLimbs == 1) && (coeffY->nbrLimbs == 1));
  stepNbr = 1;
  if (teach && showSteps)
  {
    showText("<p>");
  }
  for (;;)
  {
    if (BigIntIsZero(&V3))
    {                    // V3 equals zero, so go out of loop.
      break;
    }
    if (teach && showSteps)
    {
      showText(lang ? "Paso " : "Step ");
      showInt(stepNbr);
      showText(": ");
      paren(&U1);    // U1
      showText(" &times; ");
      paren(coeffX);
      showText(" + ");
      paren(&U2);    // U2
      showText(" &times; ");
      paren(coeffY);
      showText(" = ");
      paren(&U3);    // U3
      showText("<br>");
    }
    floordiv(&U3, &V3, &q);                 // q <- floor(U3 / V3).
    (void)BigIntMultiply(&q, &V1, &bigTmp);
    BigIntSubt(&U1, &bigTmp, &bigTmp);      // T <- U1 - q * V1
    CopyBigInt(&U1, &V1);                   // U1 <- V1
    CopyBigInt(&V1, &bigTmp);               // V1 <- T
    (void)BigIntMultiply(&q, &V2, &bigTmp);
    BigIntSubt(&U2, &bigTmp, &bigTmp);      // T <- U2 - q * V2
    CopyBigInt(&U2, &V2);                   // U2 <- V2
    CopyBigInt(&V2, &bigTmp);               // V2 <- T
    (void)BigIntMultiply(&q, &V3, &bigTmp);
    BigIntSubt(&U3, &bigTmp, &bigTmp);      // T <- U3 - q * V3
    CopyBigInt(&U3, &V3);                   // U3 <- V3
    CopyBigInt(&V3, &bigTmp);               // V3 <- T
    stepNbr++;
  }
  (void)BigIntDivide(coeffInd, &U3, &q);
  // Compute q as -coeffInd / U3
  BigIntChSign(&q);
  // Compute Xind as -U1 * coeffInd / U3
  (void)BigIntMultiply(&U1, &q, &Xind);
  // Set Xlin to coeffY
  CopyBigInt(&Xlin, coeffY);
  // Compute Yind as -U2 * coeffInd / U3
  (void)BigIntMultiply(&U2, &q, &Yind);
  // Set Ylin to -coeffX
  CopyBigInt(&Ylin, coeffX);
  BigIntChSign(&Ylin);
  if (teach)
  {
    showText(lang ? "Paso " : "Step ");
    showInt(stepNbr);
    showText(": ");
    paren(&U1);    // U1
    showText(" &times; ");
    paren(coeffX);
    showText(" + ");
    paren(&U2);    // U2
    showText(" &times; ");
    paren(coeffY);
    showText(" = ");
    paren(&U3);    // U3
    CopyBigInt(&bigTmp, coeffInd);
    BigIntChSign(&bigTmp);
    if ((q.sign != SIGN_POSITIVE) || (q.nbrLimbs != 1) || (q.limbs[0].x != 1))
    {    // Multiplier is not 1.
      showText(lang ? "</p><p>Multiplicando la última ecuación por " :
        "</p><p>Multiplying the last equation by ");
      paren(&q);
      showText(lang ? " obtenemos " : " we obtain:<br>");
      paren(&Xind);
      showText(" &times; ");
      paren(coeffX);
      showText(" + ");
      paren(&Yind);
      showText(" &times; ");
      paren(coeffY);
      showText(" = ");
      shownbr(&bigTmp);
    }
    showText(lang ? "</p><p>Sumando y restando " : "</p><p>Adding and subtracting ");
    paren(coeffX);
    showText(" &times; ");
    paren(coeffY);
    showText(" t' ");
    showText(lang ? "obtenemos" : "we obtain");
    showText(":</p><p>(");
    shownbr(&Xind);
    showText(" + ");
    paren(coeffY);
    showText(" t') &times; ");
    paren(coeffX);
    showText(" + (");
    shownbr(&Yind);
    showText(" - ");
    paren(coeffX);
    showText(" t') &times; ");
    paren(coeffY);
    showText(" = ");
    shownbr(&bigTmp);
    showText(lang ? "</p><p>Así, la solución está dada por el conjunto:</p>" :
      "</p><p>So, the solution is given by the set:</p>");
    PrintLinear(0, "t'");
    showText("</p>");
  }
  if (teach)
  {
    BigIntChSign(coeffX);
    showText(lang ? "<p>Es posible reducir los términos independientes.</p>"
      "<p>Sustituyendo t = t' + K:</p>":
      "<p>We can reduce the constant terms.</p><p>Substituting t = t' + K:</p>");
    enum eLinearSolution t;
    showText("<p>x = ");
    t = Show(coeffY, "t", SOLUTION_FOUND);
    t = Show(coeffY, "K", t);
    Show1(&Xind, t);
    showText("<br>y = ");
    t = Show(coeffX, "t", SOLUTION_FOUND);
    t = Show(coeffX, "K", t);
    Show1(&Yind, t);
    showText(lang?"</p><p>Debemos hallar el valor de K que minimice la suma de cuadrados "
                  "de los términos independientes:</p>":
                  "</p><p>We must find the value of K that minimizes the sum of squares "
                  "of the constant terms:</p>");
    showText("<p>f(K) = (");
    t = Show(coeffY, "K", SOLUTION_FOUND);
    Show1(&Xind, t);
    showText(")");
    showText(squareText);
    showText(" + (");
    t = Show(coeffX, "K", SOLUTION_FOUND);
    Show1(&Yind, t);
    showText(")");
    showText(squareText);
    showText(lang ? "</p><p>Para hallar el mínimo, la derivada de f(K) debe ser cero.":
                    "</p><p>To find the minimum, the derivative of f(K) must be zero.");
    showText("</p><p>2 &times; (");
    t = Show(coeffY, "K", SOLUTION_FOUND);
    Show1(&Xind, t);
    showText(") &times; ");
    if (coeffY->sign == SIGN_POSITIVE)
    {
      shownbr(coeffY);
    }
    else
    {
      showText("(");
      shownbr(coeffY);
      showText(")");
    }
    showText(" + 2 &times; (");
    t = Show(coeffX, "K", SOLUTION_FOUND);
    Show1(&Yind, t);
    showText(") &times; ");
    if (coeffX->sign == SIGN_POSITIVE)
    {
      shownbr(coeffX);
    }
    else
    {
      showText("(");
      shownbr(coeffX);
      showText(")");
    }
    showText(" = 0</p><p>");
    (void)BigIntMultiply(coeffY, &Xind, &U1);
    BigIntAdd(&U1, &U1, &U1);
    (void)BigIntMultiply(coeffX, &Yind, &U3);
    BigIntAdd(&U3, &U3, &U3);
    (void)BigIntMultiply(coeffY, coeffY, &U2);
    BigIntAdd(&U2, &U2, &U2);
    (void)BigIntMultiply(coeffX, coeffX, &V1);
    BigIntAdd(&V1, &V1, &V1);
    t = Show(&U2, "K", SOLUTION_FOUND);
    Show1(&U1, t);
    t = Show(&V1, "K", NO_SOLUTIONS);
    Show1(&U3, t);
    showText(" = 0</p><p>");
    BigIntAdd(&U1, &U3, &U1);
    BigIntAdd(&U2, &V1, &U2);
    t = Show(&U2, "K", SOLUTION_FOUND);
    Show1(&U1, t);
    showText(" = 0</p>");
    BigIntChSign(coeffX);
  }
  // Substitute variables so the independent coefficients can be minimized.
  // Reuse variables U1, U2, U3, V1, V2, V3.
  (void)BigIntMultiply(coeffX, coeffX, &U1);
  (void)BigIntMultiply(coeffY, coeffY, &q);
  BigIntAdd(&U1, &q, &U1);                  // U1 <- coeffX^2 + coeffY^2
  CopyBigInt(&U2, &U1);
  BigIntDivideBy2(&U2);                     // U2 <- (coeffX^2 + coeffY^2)/2
  (void)BigIntMultiply(coeffX, &Yind, &q);
  BigIntAdd(&U2, &q, &U2);
  (void)BigIntMultiply(coeffY, &Xind, &q);
  BigIntSubt(&U2, &q, &U2);
  floordiv(&U2, &U1, &U1);                  // U1 <- delta to add to t'
  if (teach)
  {
    showText(lang? "<p>Sustituyendo ": "<p>By making the substitution ");
    showText("<var>t'</var> = ");
    shownbr(&U1);
    if ((Xlin.sign == SIGN_NEGATIVE) && (Ylin.sign == SIGN_NEGATIVE))
    {    // If both coefficients are negative, change sign of transformation.
      showText(" &minus;");
    }
    else
    {
      showText(" +");
    }
    showText(" <var>t</var> ");
    showText(lang?"finalmente obtenemos:</p>": "we finally obtain:</p>");
  }
  (void)BigIntMultiply(&U1, coeffY, &q);
  BigIntAdd(&Xind, &q, &Xind);        // Xind <- Xind + coeffY * delta
  (void)BigIntMultiply(&U1, coeffX, &q);
  BigIntSubt(&Yind, &q, &Yind);       // Yind <- Yind - coeffX * delta
  if ((Xlin.sign == SIGN_NEGATIVE) && (Ylin.sign == SIGN_NEGATIVE))
  {    // If both coefficients are negative, make them positive.
    BigIntChSign(&Xlin);
    BigIntChSign(&Ylin);
  }
  return SOLUTION_FOUND;
}

static void EnsureCoeffAIsNotZero(void)
{
  ExchXY = false;
  if (BigIntIsZero(&ValA))
  {      // Next algorithm does not work if A = 0. In this case, exchange x and y.
    ExchXY = true;
    CopyBigInt(&bigTmp, &ValA);   // Exchange coefficients of x^2 and y^2.
    CopyBigInt(&ValA, &ValC);
    CopyBigInt(&ValC, &bigTmp);
    CopyBigInt(&bigTmp, &ValD);   // Exchange coefficients of x and y.
    CopyBigInt(&ValD, &ValE);
    CopyBigInt(&ValE, &bigTmp);
  }
}

static void ShowTDiscrZero(void)
{
  if (ExchXY)
  {     // Show bx + 2cy
    CopyBigInt(&ValH, &ValB);
    BigIntAdd(&ValA, &ValA, &ValI);
  }
  else
  {     // Show 2ax + by
    CopyBigInt(&ValI, &ValB);
    BigIntAdd(&ValA, &ValA, &ValH);
  }
  intToBigInteger(&ValJ, 0);
  ShowLin(&ValH, &ValI, &ValJ, "<var>x</var>", "<var>y</var>");
}

static void DiscriminantIsZero(void)
{
  EnsureCoeffAIsNotZero();
  // ax^2 + bxy + cx^2 + dx + ey + f = 0 (1)
  // Multiplying by 4a:
  // (2ax + by)^2 + 4adx + 4aey + 4af = 0
  // Let t = 2ax + by. So (1) becomes: (t + d)^2 = uy + v.
  // Compute u <- 2(bd - 2ae)
  if (teach)
  {
    showText(lang ? "<p>Multiplicando por" : "<p>Multiplying by");
    showText(ExchXY ? " 4&#8290;<var>c</var>" : " 4&#8290;<var>a</var>");
    showText("</p><p>(");
    ShowTDiscrZero();
    showText(")");
    showSquare();
    if (ExchXY)
    {
      (void)BigIntMultiply(&ValA, &ValE, &ValH);
      (void)BigIntMultiply(&ValA, &ValD, &ValI);
    }
    else
    {
      (void)BigIntMultiply(&ValA, &ValD, &ValH);
      (void)BigIntMultiply(&ValA, &ValE, &ValI);
    }
    (void)BigIntMultiply(&ValA, &ValF, &ValJ);
    multint(&ValH, &ValH, 4);
    multint(&ValI, &ValI, 4);
    multint(&ValJ, &ValJ, 4);
    if (BigIntIsZero(&ValH))
    {
      if (BigIntIsZero(&ValI))
      {
        if (ValJ.sign == SIGN_POSITIVE)
        {
          showText(" + ");
        }
      }
      else if (ValI.sign == SIGN_POSITIVE)
      {
        showText(" + ");
      }
    }
    else if (ValH.sign == SIGN_POSITIVE)
    {
      showText(" + ");
    }
    else
    {   // Nothing to do.
    }
    ShowLin(&ValH, &ValI, &ValJ, "<var>x</var>", "<var>y</var>");
    showText(" = 0</p><p>");
    showText(lang ? "Sea" : "Let");
    showText(" <var>t</var> = ");
    ShowTDiscrZero();
    intToBigInteger(&ValJ, 0);
    intToBigInteger(&ValH, 1);
    showEqNbr(1);
    showText("</p><p>(<var>t</var> + <var>");
    showText(ExchXY?"e</var>)": "d</var>)");
    showSquare();
    showText(" = (");
    ShowLin(&ValJ, &ValH, &ValD, "", "<var>t</var>");
    showText(")");
    showSquare();
    showText(" = ");
  }
  (void)BigIntMultiply(&ValB, &ValD, &Aux[0]);
  (void)BigIntMultiply(&ValA, &ValE, &Aux[1]);
  multint(&Aux[2], &Aux[1], 2);
  BigIntSubt(&Aux[0], &Aux[2], &Aux[1]);
  multint(&ValU, &Aux[1], 2);

  // Compute v <- d^2 - 4af
  (void)BigIntMultiply(&ValD, &ValD, &Aux[2]);
  (void)BigIntMultiply(&ValA, &ValF, &Aux[1]);
  multint(&Aux[3], &Aux[1], 4);
  BigIntSubt(&Aux[2], &Aux[3], &ValV);
  if (teach)
  {
    if (ExchXY)
    {
      ptrVarNameX = "<var>y</var>";
      ptrVarNameY = "<var>x</var>";
    }
    else
    {
      ptrVarNameX = "<var>x</var>";
      ptrVarNameY = "<var>y</var>";
    }
    ShowLin(&ValJ, &ValU, &ValV, ptrVarNameX, ptrVarNameY);
    if (!BigIntIsZero(&ValU) && !BigIntIsZero(&ValV))
    {
      showEqNbr(2);
    }
    showText("</p><p>");
    showText(lang ? "donde el coeficiente lineal es ": "where the linear coefficient is ");
    showText(ExchXY ? "2&#8290;(<var>b</var>&#8290;<var>e</var> &minus; 2&#8290;<var>c</var><var>d</var>)" :
      "2&#8290;(<var>b</var>&#8290;<var>d</var> &minus; 2&#8290;<var>a</var><var>e</var>)");
    showText(lang ? " y el término independiente es " : " and the constant coefficient is ");
    showText(ExchXY ? "<var>e</var>&sup2; &minus; 4&#8290;<var>c</var><var>f</var>.</p>" :
      "<var>d</var>&sup2; &minus; 4&#8290;<var>a</var><var>f</var>.</p>");
  }
  if (BigIntIsZero(&ValU))
  { // u equals zero, so (t+d)^2 = v.
    enum eLinearSolution ret;
    if (ValV.sign == SIGN_NEGATIVE)
    { // There are no solutions when v is negative.
      if (teach)
      {
        showText(lang ? "<p>Un cuadrado no puede ser igual a un número negativo, así que no hay soluciones.</p>" :
          "A square cannot be equal to a negative number, so there are no solutions.</p>");
      }
      return;
    }
    if (BigIntIsZero(&ValV))
    { // v equals zero, so (1) becomes 2ax + by + d = 0
      multint(&Aux[3], &ValA, 2);
      ret = LinearEq(&Aux[3], &ValB, &ValD);
      startResultBox(ret);
      PrintLinear(ret, "<var>t</var>");
      endResultBox(ret);
      return;
    }
    // u equals zero but v does not.
    // v must be a perfect square, otherwise there are no solutions.
    squareRoot(ValV.limbs, ValG.limbs, ValV.nbrLimbs, &ValG.nbrLimbs);
    ValG.sign = SIGN_POSITIVE;          // g <- sqrt(v).
    (void)BigIntMultiply(&ValG, &ValG, &Aux[3]);
    if (!BigIntEqual(&ValV, &Aux[3]))
    {   // v is not perfect square.
      if (teach)
      {
        showText(lang ? "<p>El término derecho no es cuadrado perfecto, así que no hay soluciones.</p>" :
          "The right hand side is not a perfect square, so there are no solutions.</p>");
      }
      return;
    }
    // The original equation is now: 2ax + by + (d +/- g) = 0
    multint(&Aux[3], &ValA, 2);
    BigIntAdd(&ValD, &ValG, &Aux[4]);
    CopyBigInt(&Aux[5], &ValB);
    if (teach)
    {
      showText("<p>");
      ShowLin(&ValJ, &ValH, &ValD, "", "<var>t</var>");
      showText(" = &pm;");
      shownbr(&ValG);
      showText(lang ? "</p><p>La ecuación representa dos rectas paralelas. La primera recta es: </p><p>" :
        "</p><p>This equation represents two parallel lines. The first line is: </p><p>");
      if (ExchXY)
      {
        ShowLin(&Aux[5], &Aux[3], &Aux[4], "<var>x</var>", "<var>y</var>");
      }
      else
      {
        ShowLin(&Aux[3], &Aux[5], &Aux[4], "<var>x</var>", "<var>y</var>");
      }
      showText(" = 0</p>");
    }
    ret = LinearEq(&Aux[3], &Aux[5], &Aux[4]);
    startResultBox(ret);
    PrintLinear(ret, "<var>t</var>");
    endResultBox(ret);
    multint(&Aux[3], &ValA, 2);
    BigIntSubt(&ValD, &ValG, &Aux[4]);
    CopyBigInt(&Aux[5], &ValB);
    if (teach)
    {
      showText(lang? "<p>La segunda recta es:</p><p>": "<p>The second line is:</p><p>");
      if (ExchXY)
      {
        ShowLin(&Aux[5], &Aux[3], &Aux[4], "<var>x</var>", "<var>y</var>");
      }
      else
      {
        ShowLin(&Aux[3], &Aux[5], &Aux[4], "<var>x</var>", "<var>y</var>");
      }
      showText(" = 0</p>");
    }
    ret = LinearEq(&Aux[3], &Aux[5], &Aux[4]);
    startResultBox(ret);
    PrintLinear(ret, "<var>t</var>");
    endResultBox(ret);
    return;
  }
  // At this moment u does not equal zero.
  // We have to solve the congruence T^2 = v (mod u) where T = t+d and t = 2ax+by.
  if (teach)
  {
    showText(lang ? "<p>Tenemos que resolver" : "<p>We have to solve");
    showText(" <var>T</var>");
    showSquare();
    showText(" = ");
    shownbr(&ValV);
    CopyBigInt(&ValH, &ValU);
    ValH.sign = SIGN_POSITIVE;
    showText(" (mod ");
    shownbr(&ValH);
    showText(")</p>");
  }
  callbackQuadModType = CBACK_QMOD_PARABOLIC;
  intToBigInteger(&coeffQuadr, 1);
  intToBigInteger(&coeffLinear, 0);
  CopyBigInt(&coeffIndep, &ValV);
  BigIntChSign(&coeffIndep);
  CopyBigInt(&modulus, &ValU);
  equationNbr = 3;
  SolveQuadModEquation();
}

// Compute coefficients of x: V1 + V2 * w + V3 * w^2
static void ComputeXDiscrZero(void)
{
  // Let m = 2be - 4cd
  (void)BigIntMultiply(&ValB, &ValE, &U3);
  (void)BigIntMultiply(&ValC, &ValD, &bigTmp);
  BigIntAdd(&bigTmp, &bigTmp, &bigTmp);
  BigIntSubt(&U3, &bigTmp, &U3);
  BigIntAdd(&U3, &U3, &U3);        // U3 <- m
                                   // Compute V1 <- (x'j - k)/2a + mx'^2
  (void)BigIntMultiply(&U2, &ValJ, &bigTmp);
  BigIntSubt(&bigTmp, &ValK, &bigTmp);
  (void)BigIntDivide(&bigTmp, &ValA, &V1);
  BigIntDivideBy2(&V1);
  (void)BigIntMultiply(&U3, &U2, &bigTmp);
  (void)BigIntMultiply(&bigTmp, &U2, &bigTmp);
  BigIntAdd(&V1, &bigTmp, &V1);
  // Compute V2 <- (j/2a + 2mx')z
  (void)BigIntDivide(&ValJ, &ValA, &V2);
  BigIntDivideBy2(&V2);
  (void)BigIntMultiply(&U3, &U2, &bigTmp);
  BigIntAdd(&bigTmp, &bigTmp, &bigTmp);
  BigIntAdd(&V2, &bigTmp, &V2);
  (void)BigIntMultiply(&V2, &ValZ, &V2);
  // Compute V3 as m*z^2
  (void)BigIntMultiply(&U3, &ValZ, &V3);
  (void)BigIntMultiply(&V3, &ValZ, &V3);
}

// Compute coefficients of y: V1 + V2 * w + V3 * w^2
static void ComputeYDiscrZero(void)
{
  // Compute V1 <- r + sx' + ux'^2
  (void)BigIntMultiply(&ValU, &U2, &V1);
  BigIntAdd(&V1, &ValS, &V1);
  (void)BigIntMultiply(&V1, &U2, &V1);
  BigIntAdd(&V1, &ValR, &V1);
  // Compute V2 <- (s + 2ux')z
  (void)BigIntMultiply(&ValU, &U2, &V2);
  BigIntAdd(&V2, &V2, &V2);
  BigIntAdd(&V2, &ValS, &V2);
  (void)BigIntMultiply(&V2, &ValZ, &V2);
  // Compute V3 <- uz^2
  (void)BigIntMultiply(&ValU, &ValZ, &V3);
  (void)BigIntMultiply(&V3, &ValZ, &V3);
}

static void callbackQuadModParabolic(const BigInteger *value)
{  // The argument of this function is T. t = T - d + uk (k arbitrary).
   // Compute ValR <- (T^2 - v)/u
  if (teach)
  {
    showText("<p><var>t</var> = <var>T</var> &minus; <var>");
    showText(ExchXY ? "e" : "d");
    showText("</var> ");
    CopyBigInt(&ValH, &ValU);
    if (ValH.sign == SIGN_POSITIVE)
    {
      showText("+ ");
    }
    else
    {
      showText("&minus; ");
    }
    if ((ValU.nbrLimbs > 1) || (ValU.limbs[0].x > 1))
    {         // Absolute value of U is not 1.
      ValH.sign = SIGN_POSITIVE;
      shownbr(&ValH);
    }
    showText(" &#8290<var>k</var> = ");
    BigIntSubt(value, &ValD, &ValH);
    ShowLinInd(&ValU, &ValH, "<var>k</var>");
    showEqNbr(equationNbr);
    showText(" (");
    showText(lang? "donde <var>k</var> puede ser cualquier entero).</p>": "where <var>k</var> is any integer).</p>");
  }
  (void)BigIntMultiply(value, value, &bigTmp);
  BigIntSubt(&bigTmp, &ValV, &bigTmp);
  (void)BigIntDivide(&bigTmp, &ValU, &ValR);
   // Compute ValS as 2*T
  BigIntAdd(value, value, &ValS);
  if (teach)
  {
    showText(lang ? "<p>Reemplazando <var>t</var> en la ecuación": "<p>Replacing <var>t</var> in equation");
    showEqNbr(2);
    showText(lang ? " y despejando " : " we can get the value of ");
    showText(ptrVarNameY);
    showText(":</p><p>");
    showText(ptrVarNameY);
    showText(" = ");
    PrintQuad(&ValU, &ValS, &ValR, "<var>k</var>", NULL);
    showEqNbr(equationNbr+1);
    showText("</p><p>");
    showText(lang ? "De" : "From");
    showEqNbr(1);
    showText(lang ? " y" : " and");
    showEqNbr(equationNbr);
    showText(":</p><p>");
    ShowTDiscrZero();
    showText(" = ");
    BigIntSubt(value, &ValD, &ValH);
    ShowLinInd(&ValU, &ValH, "<var>k</var>");
    showText("</p>");
    BigIntAdd(&ValA, &ValA, &ValK);
    if (!BigIntIsZero(&ValB))
    {
      showText(lang ? "<p>Usando" : "<p>Using");
      showEqNbr(equationNbr + 1);
      showText(lang ? " para sustituir " : " to substitute ");
      showText(ptrVarNameY);
      showText(":</p><p>");
      intToBigInteger(&ValJ, 0);
      ShowLinInd(&ValK, &ValJ, ptrVarNameX);
      showText(" = ");
      (void)BigIntMultiply(&ValB, &ValU, &ValJ);
      BigIntChSign(&ValJ);               // Quadratic coeff = -bU.
      (void)BigIntMultiply(&ValB, &ValS, &bigTmp);
      BigIntSubt(&ValU, &bigTmp, &ValI); // Linear coeff = U - bS.
      (void)BigIntMultiply(&ValB, &ValR, &bigTmp);
      BigIntSubt(&ValH, &bigTmp, &ValH); // Independent coeff = H - bR.
      PrintQuad(&ValJ, &ValI, &ValH, "<var>k</var>", NULL);
      showText("</p>");
    }
    // if independent coefficient is not multiple of GCD(I, K) then show that
    // there are no solutions. Note that J is always multiple of K.
    BigIntGcd(&ValI, &ValK, &ValJ);
    (void)BigIntRemainder(&ValH, &ValJ, &bigTmp);
    if (!BigIntIsZero(&bigTmp))
    {
      showText(lang ? "<p>El coeficiente independiente ": "<p>The independent coefficient ");
      shownbr(&ValH);
      showText(lang ? " no es múltiplo de " : " is not multiple of ");
      shownbr(&ValJ);
      showText(lang ? ", que es el máximo común divisor de los otros tres coeficientes, por lo que no hay solución.</p>" :
        ", which is the greatest common divisor of the other three coefficients, so there is no solution.</p>");
    }
  }
   // Find k from the congruence jk = K (mod z) where j = u-bs, K = d+br-T, z = 2a
   // Compute j <- u-bs
  (void)BigIntMultiply(&ValB, &ValS, &bigTmp);
  BigIntSubt(&ValU, &bigTmp, &ValJ);
   // Compute K <- d+br-T
  (void)BigIntMultiply(&ValB, &ValR, &bigTmp);
  BigIntAdd(&ValD, &bigTmp, &bigTmp);
  BigIntSubt(&bigTmp, value, &ValK);
   // Compute z <- 2a
  BigIntAdd(&ValA, &ValA, &ValZ);
  // If K is not multiple of gcd(j, z) there is no solution.
  BigIntGcd(&ValJ, &ValZ, &bigTmp);
  (void)BigIntRemainder(&ValK, &bigTmp, &U1);
  if (!BigIntIsZero(&U1))
  {
    return;
  }
  // Compute g = gcd(j, K, z), then recalculate j <- j/g, K <- K/g, z <- z/g
  BigIntGcd(&bigTmp, &ValK, &U1);
  (void)BigIntDivide(&ValJ, &U1, &U2);    // U2 <- j
  (void)BigIntDivide(&ValK, &U1, &U3);    // U3 <- K
  (void)BigIntDivide(&ValZ, &U1, &ValZ);
  ValZ.sign = SIGN_POSITIVE;        // Use positive sign for modulus.
  (void)BigIntRemainder(&U2, &ValZ, &U2);
  if (U2.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&U2, &ValZ, &U2);
  }
  (void)BigIntRemainder(&U3, &ValZ, &U3);
  if (U3.sign == SIGN_NEGATIVE)
  {
    BigIntAdd(&U3, &ValZ, &U3);
  }
  if (BigIntIsZero(&U2))
  {    // M and N equal zero.
    intToBigInteger(&ValZ, 1);     // In this case 0*k = 0 (mod z) means any k is valid.
  }
  else
  {
    BigIntGeneralModularDivision(&U2, &U3, &ValZ, &U2);   // U2 <- x'
  }
  if (ExchXY)
  {
    ComputeYDiscrZero();
  }
  else
  {
    ComputeXDiscrZero();
  }
  showAlso();
  startResultBox(SOLUTION_FOUND);
  showText("<p><var>x</var> = ");
  PrintQuad(&V3, &V2, &V1, "<var>k</var>", NULL);
  showText("<br>");
  if (ExchXY)
  {
    ComputeXDiscrZero();
  }
  else
  {
    ComputeYDiscrZero();
  }
  showText("<var>y</var> = ");
  PrintQuad(&V3, &V2, &V1, "<var>k</var>", NULL);
  if (teach)
  {
    showText("<br>");
  }
  showText("</p>");
  endResultBox(SOLUTION_FOUND);
  if (teach)
  {
    showText("<br>");
    equationNbr += 2;
  }
}

static void ShowPoint(const BigInteger *X, const BigInteger *Y)
{
  int solution = 0;
  if (teach)
  {
    showText("<p><var>X</var> = ");
    shownbr(X);
    showText(", <var>Y</var> = ");
    shownbr(Y);
    showText("</p>");
  }
  // Check first that (X+alpha) and (Y+beta) are multiple of D.
  BigIntAdd(X, &ValAlpha, &Tmp1);
  BigIntAdd(Y, &ValBeta, &Tmp2);
  if (teach && !(BigIntIsZero(&ValAlpha) && BigIntIsZero(&ValBeta)))
  {
    showText("<p><var>X</var> + <var>&alpha;</var> = ");
    shownbr(&Tmp1);
    showText(", <var>Y</var> + <var>&beta;</var> = ");
    shownbr(&Tmp2);
    showText("</p>");
  }
  (void)BigIntRemainder(&Tmp1, &ValDiv, &bigTmp);
  if (BigIntIsZero(&bigTmp))
  {
    (void)BigIntRemainder(&Tmp2, &ValDiv, &bigTmp);
    if (BigIntIsZero(&bigTmp))
    {
      if (teach && !BigIntIsOne(&ValDiv))
      {
        showText(lang ? "<p>Dividiendo estos números por" : "<p>Dividing these numbers by");
        showText(" <var>D</var> = ");
        shownbr(&ValDiv);
        showText(":</p>");
      }
      (void)BigIntDivide(&Tmp1, &ValDiv, &Tmp1);
      (void)BigIntDivide(&Tmp2, &ValDiv, &Tmp2);
      if (callbackQuadModType == CBACK_QMOD_HYPERBOLIC)
      {
        if (teach)
        {
          showSolutionXY(&Tmp1, &Tmp2);
        }
        ShowXY(&Tmp1, &Tmp2);
      }
      else
      {
        startResultBox(SOLUTION_FOUND);
        ShowXY(&Tmp1, &Tmp2);
        endResultBox(SOLUTION_FOUND);
      }
      solution = 1;
      showRecursiveSolution = 1; // Show recursive solution if it exists.
    }
  }
  if (teach && (solution == 0))
  {
    showText(lang ? "<p>Estos números no son múltiplos de" : "<p>These numbers are not multiple of");
    showText(" <var>D</var> = ");
    shownbr(&ValDiv);
    showText(".</p>");
  }
}

// Solve ax^2+bxy+cy^2 = K
// The quadratic modular equation algorithm requires that gcd(a, n) = 1.
// At this point gcd(a, b, c) = 1
// The possibilities are:
// - gcd(a, K) = 1. There is nothing to do.
// Otherwise perform the transformation x = PX + QY, y = RX + SY with PS - QR = 1.
// In particular perform: x = mX + (m-1)Y, y = X + Y
// We get: (am^2+bm+c)*X^2 + (2(am^2+bm+c) - (2am+b))*XY + ((am^2+bm+c) - (2am+b) + a)*Y^2
// Also perform: x = X + Y, y = (m-1)X + mY
// We get: (a+(m-1)*b+(m-1)^2*c)*X^2 + (2a + b(2m-1) + 2cm(m-1))*X*Y + (a+bm+cm^2)*Y^2
// The discriminant of the new formula does not change.
// Compute m=1, 2, 3,... until gcd(am^2+bm+c, K) = 1.
// When using the second formula, change sign of m so we know the formula used when
// undoing the unimodular transformation later.

// Since the algorithm discovers only primitive solutions, i.e. solutions (x,y) where
// gcd(x,y) = 1, we need to solve ax'^2+bx'y'+cy'^2 = K/R^2 where R^2 is a divisor of K.
// Then we get x = Rx', y = Ry'.
static void NonSquareDiscriminant(void)
{
  enum eSign ValKSignBak;
  int numFactors;
  struct sFactors *pstFactor;
  char* ptrFactorDec;
             // Find GCD(a,b,c)
  BigIntGcd(&ValA, &ValB, &bigTmp);
  BigIntGcd(&ValC, &bigTmp, &ValGcdHomog);
    // Divide A, B, C and K by this GCD.
  (void)BigIntDivide(&ValA, &ValGcdHomog, &ValA);
  (void)BigIntDivide(&ValB, &ValGcdHomog, &ValB);
  (void)BigIntDivide(&ValC, &ValGcdHomog, &ValC);
  (void)BigIntDivide(&ValK, &ValGcdHomog, &ValK);
  // Divide discriminant by the square of GCD.
  (void)BigIntDivide(&discr, &ValGcdHomog, &discr);
  (void)BigIntDivide(&discr, &ValGcdHomog, &discr);
  if (BigIntIsZero(&ValK))
  {          // If k=0, the only solution is (X, Y) = (0, 0)
    ShowPoint(&ValK, &ValK);
    return;
  }
  CopyBigInt(&ValABak, &ValA);
  CopyBigInt(&ValBBak, &ValB);
  CopyBigInt(&ValCBak, &ValC);
  // Factor independent term.
  ValKSignBak = ValK.sign;
  ValK.sign = SIGN_POSITIVE;
  NumberLength = ValK.nbrLimbs;
  BigInteger2IntArray(nbrToFactor, &ValK);
  ptrFactorDec = tofactorDec;
  Bin2Dec(&ptrFactorDec, ValK.limbs, ValK.nbrLimbs, groupLen);
  factor(&ValK, nbrToFactor, factorsMod, astFactorsMod);
  CopyBigInt(&LastModulus, &ValK);           // Do not factor again same modulus.
  ValK.sign = ValKSignBak;
  // Find all indexes of prime factors with even multiplicity.
  nbrPrimesEvenMultiplicity = 0;
  numFactors = astFactorsMod[0].multiplicity;
  pstFactor = &astFactorsMod[1];
  for (int factorNbr = 1; factorNbr <= numFactors; factorNbr++)
  {
    if (pstFactor->multiplicity > 1)
    {      // At least prime is squared.
      indexEvenMultiplicity[nbrPrimesEvenMultiplicity] = factorNbr;
      originalMultiplicities[nbrPrimesEvenMultiplicity] = pstFactor->multiplicity & (~1); // Convert to even.
      nbrPrimesEvenMultiplicity++;
    }
    pstFactor++;
  }
  (void)memset(counters, 0, sizeof(counters));
  (void)memset(isDescending, 0, sizeof(isDescending));
  intToBigInteger(&ValE, 1);  // Initialize multiplier to 1.
  // Loop that cycles through all square divisors of the independent term.
  equationNbr = 2;
  BigIntGcd(&ValA, &ValK, &bigTmp);
  intToBigInteger(&ValM, 0);
  varXnoTrans = "<var>X</var>";
  varYnoTrans = "<var>Y</var>";
  if ((bigTmp.nbrLimbs != 1) || (bigTmp.limbs[0].x != 1))
  {                        // gcd(a, K) is not equal to 1.
    if (teach)
    {
      showText(lang ? "<p>El algoritmo requiere que el coeficiente de <var>X</var>" :
        "<p>The algorithm requires that the coefficient of <var>X</var>");
      showSquare();
      showText(lang ? " y el término independiente sean primos entre sí. Como esto no ocurre, debemos hallar un valor de <var>m</var> tal que aplicando una de las transformaciones unimodulares siguientes</p>" : 
        " and the right hand side are coprime. This does not happen, so we have to find a value of <var>m</var> such that applying one of the unimodular transformations</p> ");
      showText("<ul><li><var>X</var> = <var>m</var>&#8290;<var>U</var> + <var>(m&minus;1)</var>&#8290;<var>V</var>, ");
      showText("<var>Y</var> = <var>U</var> + <var>V</var></li>");
      showText("<li><var>X</var> = <var>U</var> + <var>V</var>, ");
      showText("<var>Y</var> = <var>(m&minus;1)</var>&#8290;<var>U</var> + <var>m</var>&#8290;<var>V</var></li></ul>");
      showText(lang?"<p>el coeficiente de <var>U</var>": "<p>the coefficient of <var>U</var>");
      showSquare();
      showText(lang ? " y el término independiente sean primos entre sí. Este coeficiente es igual a ":
        " and the right hand side are coprime. This coefficient equals ");
      PrintQuad(&ValA, &ValB, &ValC, "<var>m</var>", NULL);
      showText(lang ? " en el primer caso y " : " in the first case and ");
      PrintQuad(&ValC, &ValB, &ValA, "(<var>m</var> &minus; 1)", NULL);
      showText(lang ? " en el segundo caso.</p>" : " in the second case.</p>");
    }
    intToBigInteger(&ValM, 0);
    do
    {                         // Compute U1 = cm^2 + bm + a and exit loop if this
                              // value is not coprime to K.
      (void)BigIntMultiply(&ValC, &ValM, &U2);
      BigIntAdd(&U2, &ValB, &U1);
      (void)BigIntMultiply(&U1, &ValM, &U1);
      BigIntAdd(&U1, &ValA, &U1);
      BigIntGcd(&U1, &ValK, &bigTmp);
      if ((bigTmp.nbrLimbs == 1) && (bigTmp.limbs[0].x == 1))
      {
        addbigint(&ValM, 1);  // Increment M.
        BigIntChSign(&ValM);  // Change sign to indicate type.
        break;
      }
      addbigint(&ValM, 1);    // Increment M.
                              // Compute U1 = am^2 + bm + c and loop while this
                              // value is not coprime to K.
      (void)BigIntMultiply(&ValA, &ValM, &U2);
      BigIntAdd(&U2, &ValB, &U1);
      (void)BigIntMultiply(&U1, &ValM, &U1);
      BigIntAdd(&U1, &ValC, &U1);
      BigIntGcd(&U1, &ValK, &bigTmp);
    } while ((bigTmp.nbrLimbs != 1) || (bigTmp.limbs[0].x != 1));
    // Compute 2am + b or 2cm + b as required.
    BigIntAdd(&U2, &U2, &U2);
    BigIntAdd(&U2, &ValB, &U2);
    if (ValM.sign == SIGN_POSITIVE)
    {
      // Compute c.
      BigIntSubt(&U1, &U2, &ValB);
      BigIntAdd(&ValB, &ValA, &ValC);
      // Compute b.
      BigIntAdd(&ValB, &U1, &ValB);
      // Compute a.
      CopyBigInt(&ValA, &U1);
    }
    else
    {
      // Compute c.
      BigIntAdd(&U1, &U2, &ValB);
      BigIntAdd(&ValB, &ValC, &ValC);
      // Compute b.
      BigIntAdd(&ValB, &U1, &ValB);
      // Compute a.
      CopyBigInt(&ValA, &U1);
    }
    if (teach)
    {
      int m = ValM.limbs[0].x;
      showText(lang ? "<p>Usaremos la ":"<p>We will use the ");
      if (ValM.sign == SIGN_POSITIVE)
      {
        showText(lang ? "primera" : "first");
      }
      else
      {
        showText(lang ? "segunda" : "second");
      }
      showText(lang? " transformación unimodular con <var>m</var> = ":
          " unimodular transformation with <var>m</var> = ");
      showInt(m);
      showText(": <var>X</var> = ");
      if (ValM.sign == SIGN_POSITIVE)
      {
        if (m == 1)
        {
          showText("<var>U</var>");
        }
        else
        {
          showInt(m);
          showText(" &#8290;<var>U</var> + ");
          if (m > 2)
          {
            showInt(m - 1);
            showText(" &#8290;");
          }
          showText("<var>V</var>");
        }
        showText(", <var>Y</var> = <var>U</var> + <var>V</var> ");
      }
      else
      {
        showText("<var>X</var> = <var>U</var> + <var>V</var>, <var>Y</var> = ");
        if (m > 2)
        {
          showInt(m - 1);
          showText("&#8290;");
        }
        if (m > 1)
        {
          showText("<var>U</var> + ");
          showInt(m);
          showText("&#8290;");
        }
        showText("<var>V</var> ");
      }
      showEqNbr(2);
      showText("</p>");
      showText(lang ? "<p>Mediante ": "<p>Using");
      showEqNbr(2);
      showText(lang ? ", la ecuación " : ", the equation ");
      showEqNbr(1);
      showText(lang ? " se convierte a:</p>" : " converts to:</p>");
      PrintQuad(&ValA, &ValB, &ValC, "<var>U</var>", "<var>V</var>");
      showText(" = ");
      shownbr(&ValK);
      showText(" ");
      showEqNbr(3);
      equationNbr = 4;
      varXnoTrans = "<var>U</var>";
      varYnoTrans = "<var>V</var>";
    }
  }
  if (teach)
  {
    enum eSign signK = ValK.sign;
    showText(lang ? "<p>A continuación hay que resolver varias ecuaciones cuadráticas modulares. Para ello debemos factorizar el módulo y hallar las soluciones módulo las potencias de los factores primos. Luego debemos combinar estas soluciones usando el teorema chino del resto.</p><p>Los diferentes módulos son divisores del término independiente, por lo que basta con factorizar dicho valor.</p>" :
      "<p>We will have to solve several quadratic modular equations. To do this we have to factor the modulus and find the solution modulo the powers of the prime factors. Then we combine them by using the Chinese Remainder Theorem.</p><p>The different moduli are divisors of the right hand side, so we only have to factor it once.</p>");
    showText("<p>");
    ValK.sign = SIGN_POSITIVE;   // Make the number to factor positive.
    showFactors(&ValK);
    ValK.sign = signK;
  }
  for (;;)
  {
    int index;
    if (teach)
    {
      showText("<p>");
      if ((ValE.nbrLimbs > 1) || (ValE.limbs[0].x > 1))
      {
        if (BigIntIsZero(&ValM))
        {       // No unimodular transforation.
          varX = "<var>X'</var>";
          varY = "<var>Y'</var>";
        }
        else
        {       // Unimodular transforation.
          varX = "<var>U'</var>";
          varY = "<var>V'</var>";
        }
        showText(lang ? "Sea " : "Let ");
        shownbr(&ValE);
        showText("&#8290;");
        showText(varX);
        showText(" = ");
        showText(varXnoTrans);
        showText(lang ? " y " : " and ");
        shownbr(&ValE);
        showText("&#8290;");
        showText(varY);
        showText(" = ");
        showText(varYnoTrans);
        showText(". ");
      }
      else
      {
        if (BigIntIsZero(&ValM))
        {       // No unimodular transforation.
          varX = "<var>X</var>";
          varY = "<var>Y</var>";
        }
        else
        {       // Unimodular transforation.
          varX = "<var>U</var>";
          varY = "<var>V</var>";
        }
      }
      showText(lang ? "Buscando soluciones " : "Searching for solutions ");
      showText(varX);
      showText(lang ? " e " : " and ");
      showText(varY);
      showText(lang ? " primos entre sí.</p>" : " coprime.</p>");
      if ((ValE.nbrLimbs > 1) || (ValE.limbs[0].x > 1))
      {
        showText(lang ? "<p>De la ecuación " : "<p>From equation ");
        showEqNbr(BigIntIsZero(&ValM)? 1: 3);
        showText(lang ? " obtenemos " : " we obtain ");
        PrintQuad(&ValA, &ValB, &ValC, varX, varY);
        showText(" = ");
        (void)BigIntMultiply(&ValK, &ValE, &ValH);
        (void)BigIntMultiply(&ValH, &ValE, &ValH);
        shownbr(&ValH);
        showText(" / ");
        shownbr(&ValE);
        showSquare();
        showText(" = ");
        shownbr(&ValK);
      }
    }
    CopyBigInt(&coeffQuadr, &ValA);
    CopyBigInt(&coeffLinear, &ValB);
    CopyBigInt(&coeffIndep, &ValC);
    CopyBigInt(&modulus, &ValK);
    modulus.sign = SIGN_POSITIVE;
    CopyBigInt(&LastModulus, &modulus);
    if (teach)
    {
      showText(lang ? "<p>Debemos resolver " : "<p>We have to solve:");
      PrintQuad(&coeffQuadr, &coeffLinear, &coeffIndep, "<var>T</var>", NULL);
      showText(" &equiv; 0 (mod ");
      showFactors(&modulus);
      showText(")<p>");
    }
    SolveQuadModEquation();
    // Adjust counters.
    for (index = 0; index < nbrPrimesEvenMultiplicity; index++)
    {            // Loop that increments counters.
      if (isDescending[index] == 0)
      {          // Ascending.
        pstFactor = &astFactorsMod[indexEvenMultiplicity[index]];
        if (counters[index] == originalMultiplicities[index])
        {
          isDescending[index] = 1;    // Next time it will be descending.
          continue;
        }
        else
        {
          NumberLength = *pstFactor->ptrFactor;
          IntArray2BigInteger(pstFactor->ptrFactor, &bigTmp);
          (void)BigIntMultiply(&bigTmp, &bigTmp, &U3);
          pstFactor->multiplicity -= 2;
          (void)BigIntDivide(&ValK, &U3, &ValK);       // Divide by square of prime.
          (void)BigIntMultiply(&ValE, &bigTmp, &ValE); // Multiply multiplier by prime.counters[index]++
          counters[index] += 2;
          break;
        }
      }
      else
      {         // Descending.
        if (counters[index] <= 1)
        {
          isDescending[index] = 0;    // Next time it will be ascending.
          continue;
        }
        else
        {
          pstFactor = &astFactorsMod[indexEvenMultiplicity[index]];
          NumberLength = *pstFactor->ptrFactor;
          IntArray2BigInteger(pstFactor->ptrFactor, &bigTmp);
          (void)BigIntMultiply(&bigTmp, &bigTmp, &U3);
          pstFactor->multiplicity += 2;
          (void)BigIntMultiply(&ValK, &U3, &ValK);     // Multiply by square of prime.
          (void)BigIntDivide(&ValE, &bigTmp, &ValE);   // Divide multiplier by prime.counters[index]++
          counters[index] -= 2;
          break;
        }
      }
    }
    CopyBigInt(&LastModulus, &ValK);                   // Do not try to factor the number again.
    if (index == nbrPrimesEvenMultiplicity)
    {               // All factors have been found. Exit loop.
      break;
    }
  }
  if (showRecursiveSolution && (callbackQuadModType == CBACK_QMOD_HYPERBOLIC))
  {   // Show recursive solution.
    recursiveSolution();
  }
}

static void NegativeDiscriminant(void)
{
  callbackQuadModType = CBACK_QMOD_ELLIPTIC;
  NonSquareDiscriminant();
}

static void showBeforeUnimodularSubstitution(void)
{
  if (teach)
  {
    showText("<p>");
    showText(varXnoTrans);
    showText(" = ");
    shownbr(&ValZ);
    showText(", ");
    showText(varYnoTrans);
    showText(" = ");
    shownbr(&ValO);
    showText("</p>");
    showText(lang ? "<p>De " : "<p>From ");
    showEqNbr(2);
    showText(": ");
  }
}

static void UnimodularSubstitution(void)
{
  if (ValM.sign == SIGN_NEGATIVE)
  {     // Perform the substitution: x = X + Y, y = (|m|-1)X + |m|Y
    showBeforeUnimodularSubstitution();
    ValM.sign = SIGN_POSITIVE;
    BigIntAdd(&ValZ, &ValO, &Tmp[0]);     // x
    (void)BigIntMultiply(&Tmp[0], &ValM, &Tmp[1]);
    BigIntSubt(&Tmp[1], &ValZ, &Tmp[1]);  // y
    ValM.sign = SIGN_NEGATIVE;
  }
  else if (BigIntIsZero(&ValM))
  {
    CopyBigInt(&Tmp[0], &ValZ);           // x
    CopyBigInt(&Tmp[1], &ValO);           // y
  }
  else
  {     // Perform the substitution: x = mX + (m-1)Y, y = X + Y
    showBeforeUnimodularSubstitution();
    BigIntAdd(&ValZ, &ValO, &Tmp[1]);     // y
    (void)BigIntMultiply(&Tmp[1], &ValM, &Tmp[0]);
    BigIntSubt(&Tmp[0], &ValO, &Tmp[0]);  // x
  }
}

// On input: ValH: value of u, ValI: value of v.
// Output: ((tu - nv)*E, u*E) and ((-tu + nv)*E, -u*E)
// If m is greater than zero, perform the substitution: x = mX + (m-1)Y, y = X + Y
// If m is less than zero, perform the substitution: x = X + Y, y = (|m|-1)X + |m|Y
// Do not substitute if m equals zero.

static void NonSquareDiscrSolution(BigInteger *value)
{
  CopyBigInt(&Tmp[12], value);            // Back up value.
  // Get value of tu - Kv
  (void)BigIntMultiply(value, &ValH, &ValZ);    // tu
  CopyBigInt(&bigTmp, &ValK);
  if (callbackQuadModType == CBACK_QMOD_HYPERBOLIC)
  {
    BigIntChSign(&bigTmp);                // Get K
    bigTmp.sign = SIGN_NEGATIVE;
  }
  (void)BigIntMultiply(&bigTmp, &ValI, &bigTmp);// Kv
  BigIntSubt(&ValZ, &bigTmp, &ValZ);      // U = tu - Kv
  if (teach && ((ValE.nbrLimbs > 1) || (ValE.limbs[0].x > 1)))
  {     // E > 1
    showText(lang? "<p>De ": "<p>From ");
    showEqNbr(equationNbr);
    showText(": ");
    showText(varX);
    showText(" = ");
    shownbr(&ValZ);
    showText(", ");
    showText(varY);
    showText(" = ");
    shownbr(&ValH);
  }
  (void)BigIntMultiply(&ValZ, &ValE, &ValZ);    // X = (tu - Kv)*E
  (void)BigIntMultiply(&ValH, &ValE, &ValO);    // Y = u*E
  Xbak = &Xplus;
  Ybak = &Yplus;
  UnimodularSubstitution();               // Undo unimodular substitution
  ShowPoint(&Tmp[0], &Tmp[1]);
  BigIntChSign(&ValZ);                    // (-tu - Kv)*E
  BigIntChSign(&ValO);                    // -u*E
  Xbak = &Xminus;
  Ybak = &Yminus;
  UnimodularSubstitution();               // Undo unimodular substitution
  ShowPoint(&Tmp[0], &Tmp[1]);
  CopyBigInt(value, &Tmp[12]);            // Restore value.
}

// Obtain next convergent of continued fraction of ValU/ValV
// Previous convergents U1/V1, U2/V2, U3/V3.
static void getNextConvergent(void)
{
  floordiv(&ValU, &ValV, &bigTmp);
  // Values of U3 and V3 are not used, so they can be overwritten now.
  // Compute new value of ValU and ValV.
  (void)BigIntMultiply(&bigTmp, &ValV, &U3);
  BigIntSubt(&ValU, &U3, &U3);
  CopyBigInt(&ValU, &ValV);
  CopyBigInt(&ValV, &U3);
  // Compute new convergents: h_n = a_n*h_{n-1} + h_{n-2}
  // and also k_n = k_n*k_{n-1} + k_{n-2}
  (void)BigIntMultiply(&bigTmp, &U1, &V3);
  BigIntAdd(&V3, &U2, &V3);
  CopyBigInt(&U3, &U2);
  CopyBigInt(&U2, &U1);
  CopyBigInt(&U1, &V3);
  (void)BigIntMultiply(&bigTmp, &V1, &bigTmp);
  BigIntAdd(&bigTmp, &V2, &bigTmp);
  CopyBigInt(&V3, &V2);
  CopyBigInt(&V2, &V1);
  CopyBigInt(&V1, &bigTmp);
}

static void ShowSolutionFromConvergent(void)
{
  if (teach)
  {
    showText(lang ? "<p>Solución de ": "Solution of ");
    showEqNbr(equationNbr + 1);
    showText(lang? " hallada mediante el convergente " : " found using the convergent ");
    showText(varY);
    showText(" / ");
    if (callbackQuadModType == CBACK_QMOD_HYPERBOLIC)
    {
      showText("(&minus;<var>k</var>) = ");
    }
    else
    {
      showText("<var>k</var> = ");
    }
    shownbr(&ValH);
    showText(" / ");
    shownbr(&ValI);
    showText(lang ? " de " : " of ");
    showEqNbr(contfracEqNbr);
    showText("</p>");
  }
}

static void showFirstSolution(const char *discrim, const char *valueP)
{
  showText(lang ? "<p>Cuando el discriminante vale ": "<p>When the discriminant equals ");
  showText(discrim);
  showText(lang ? " y" : " and");
  showText(" <var>P</var> = ");
  showText(valueP);
  showText(lang?", una solución es" : ", a solution is");
  showText(" (");
  showText(varY);
  showText(", <var>k</var>) = (");
}

static void showOtherSolution(const char *ordinal)
{
  showText(lang ? "<p>La ": "<p>The ");
  showText(ordinal);
  showText(lang?" solución es": " solution is");
  showText(" (");
  showText(varY);
  showText(", <var>k</var>) = (");
}

// Output:
// false = There are no solutions because gcd(P, Q, R) > 1
// true = gcd(P, Q, R) = 1.
static bool PerformTransformation(const BigInteger *value)
{
  // Compute P as (at^2+bt+c)/K
  (void)BigIntMultiply(&ValA, value, &ValQ);
  BigIntAdd(&ValQ, &ValB, &ValP);
  (void)BigIntMultiply(&ValP, value, &ValP);
  BigIntAdd(&ValP, &ValC, &ValP);
  (void)BigIntDivide(&ValP, &ValK, &ValP);
  // Compute Q <- -(2at + b).
  BigIntAdd(&ValQ, &ValQ, &ValQ);
  BigIntAdd(&ValQ, &ValB, &ValQ);
  BigIntChSign(&ValQ);
  // Compute R <- aK
  (void)BigIntMultiply(&ValA, &ValK, &ValR);
  if (teach)
  {
    showText(lang ? "<p>La transformación " : "<p>The transformation ");
    showText(varX);
    showText(" = ");
    CopyBigInt(&ValH, value);
    CopyBigInt(&ValI, &ValK);
    BigIntChSign(&ValI);
    intToBigInteger(&ValJ, 0);
    ShowLin(&ValH, &ValI, &ValJ, varY, "<var>k</var>");
    showText(" ");
    showEqNbr(equationNbr);
    showText(lang ? " convierte " : " converts ");
    PrintQuad(&ValA, &ValB, &ValC, varX, varY);
    showText(" = ");
    shownbr(&ValK);
    showText(lang ? " a " : " to ");
    showText("<var>P</var>&#8290;");
    showText(varY);
    showSquare();
    showText(" + <var>Q</var>&#8290;");
    showText(varY);
    showText("<var>k</var> + <var>R</var>&#8290; <var>k</var>");
    showSquare();
    showText(" = 1 ");
    showEqNbr(equationNbr + 1);
    showText("</p>");
    showText(lang ? "donde: " : "where: ");
    showText("<var>P</var> = (<var>a</var>&#8290;<var>T</var>");
    showSquare();
    showText(" + <var>b</var>&#8290;<var>T</var> + <var>c</var>) / n = ");
    shownbr(&ValP);
    showText(", <var>Q</var> = &minus;(2&#8290;<var>a</var>&#8290;<var>T</var> + <var>b</var>) = ");
    shownbr(&ValQ);
    showText(", <var>R</var> = <var>a</var>&#8290;<var>n</var> = ");
    shownbr(&ValR);
    showText("</p>");
  }
  // Compute gcd of P, Q and R.
  BigIntGcd(&ValP, &ValQ, &ValH);   // Use ValH and ValI as temporary variables.
  BigIntGcd(&ValH, &ValR, &ValI);
  if ((ValI.nbrLimbs == 1) && (ValI.limbs[0].x == 1))
  { // Gcd equals 1.
    return 1;
  }
  if (teach)
  {
    showText(lang ? "<p>No hay soluciones porque mcd(<var>P</var>, <var>Q</var>, <var>R</var>) es mayor que 1.</p>" :
      "<p>There are no solutions because gcd(<var>P</var>, <var>Q</var>, <var>R</var>) is greater than 1.</p>");
    equationNbr += 2;
  }
  return 0;
}

static void callbackQuadModElliptic(BigInteger *value)
{
  if (!PerformTransformation(value))
  {      // No solutions because gcd(P, Q, R) > 1.
    return;
  }
  if ((ValP.sign == SIGN_POSITIVE) && (ValP.nbrLimbs == 1))
  {
    int Plow = ValP.limbs[0].x;
    if (((discr.nbrLimbs > 1) || (discr.limbs[0].x > 4)) && (Plow == 1))
    {      // Discriminant is less than -4 and P equals 1.
      intToBigInteger(&ValH, 1);
      intToBigInteger(&ValI, 0);
      NonSquareDiscrSolution(value);   // (1, 0)
      equationNbr += 2;
      return;
    }
    if ((discr.nbrLimbs == 1) && (discr.limbs[0].x == 4))
    {      // Discriminant is equal to -4.
      CopyBigInt(&ValG, &ValQ);
      BigIntDivideBy2(&ValG);
      if (Plow == 1)
      {
        intToBigInteger(&ValH, 1);
        intToBigInteger(&ValI, 0);
        if (teach)
        {
          showFirstSolution("-4", "1");
          showText("1, 0)</p>");
        }
        NonSquareDiscrSolution(value);   // (1, 0)
        CopyBigInt(&ValH, &ValG);
        intToBigInteger(&ValI, -1);
        if (teach)
        {
          showOtherSolution(lang ? "segunda" : "second");
          showText("<var>Q</var>/2 = ");
          shownbr(&ValH);
          showText(", -1)</p>");
        }
        NonSquareDiscrSolution(value);   // (Q/2, -1)
        equationNbr += 2;
        return;
      }
      if (Plow == 2)
      {
        intToBigInteger(&ValI, -1);
        CopyBigInt(&ValH, &ValG);
        subtractdivide(&ValH, 1, 2);
        if (teach)
        {
          showFirstSolution("-4", "2");
          showText("(<var>Q</var>/2 &minus; 1) / 2 = ");
          shownbr(&ValH);
          showText(", -1)</p>");
        }
        NonSquareDiscrSolution(value);   // ((Q/2-1)/2, -1)
        intToBigInteger(&ValI, -1);
        CopyBigInt(&ValH, &ValG);
        subtractdivide(&ValH, -1, 2);
        if (teach)
        {
          showOtherSolution(lang ? "segunda" : "second");
          showText("(<var>Q</var>/2 + 1) / 2 = ");
          shownbr(&ValH);
          showText(", -1)</p>");
        }
        NonSquareDiscrSolution(value);   // ((Q/2+1)/2, -1)
        equationNbr += 2;
        return;
      }
    }
    if ((discr.nbrLimbs == 1) && (discr.limbs[0].x == 3))
    {      // Discriminant is equal to -3.
      if (Plow == 1)
      {
        intToBigInteger(&ValH, 1);
        intToBigInteger(&ValI, 0);
        if (teach)
        {
          showFirstSolution("-3", "1");
          showText("1, 0)</p>");
        }
        NonSquareDiscrSolution(value);   // (1, 0)
        CopyBigInt(&ValG, &ValQ);
        subtractdivide(&ValG, 1, 2);
        CopyBigInt(&ValH, &ValG);
        intToBigInteger(&ValI, -1);
        if (teach)
        {
          showOtherSolution(lang ? "segunda" : "second");
          showText("(<var>Q</var> &#8209; 1)/2 = ");
          shownbr(&ValH);
          showText(", -1)</p>");
        }
        NonSquareDiscrSolution(value);   // ((Q-1)/2, -1)
        intToBigInteger(&ValI, -1);
        BigIntSubt(&ValG, &ValI, &ValH);
        if (teach)
        {
          showOtherSolution(lang ? "tercera" : "third");
          showText("(<var>Q</var> + 1)/2 = ");
          shownbr(&ValH);
          showText(", -1)</p>");
        }
        NonSquareDiscrSolution(value);   // ((Q+1)/2, -1)
        equationNbr += 2;
        return;
      }
      if (Plow == 3)
      {
        intToBigInteger(&ValI, -1);
        CopyBigInt(&ValH, &ValQ);
        subtractdivide(&ValH, -3, 6);
        if (teach)
        {
          showFirstSolution("-3", "3");
          showText("(<var>Q</var> + 3)/6 = ");
          shownbr(&ValH);
          showText(", -1)</p>");
        }
        NonSquareDiscrSolution(value);   // ((Q+3)/6, -1)
        intToBigInteger(&ValI, -2);
        CopyBigInt(&ValH, &ValQ);
        subtractdivide(&ValH, 0, 3);
        if (teach)
        {
          showOtherSolution(lang ? "segunda" : "second");
          showText("<var>Q</var>/3 = ");
          shownbr(&ValH);
          showText(", -2)</p>");
        }
        NonSquareDiscrSolution(value);   // (Q/3, -2)
        intToBigInteger(&ValI, -1);
        CopyBigInt(&ValH, &ValQ);
        subtractdivide(&ValH, 3, 6);
        if (teach)
        {
          showOtherSolution(lang ? "tercera" : "third");
          showText("(<var>Q</var> &minus; 3)/6 = ");
          shownbr(&ValH);
          showText(", -1)</p>");
        }
        NonSquareDiscrSolution(value);   // ((Q-3)/6, -1)
        equationNbr += 2;
        return;
      }
    }
  }
  if (teach)
  {
    int coeffNbr = 0;
    contfracEqNbr = equationNbr + 2;
    showText(lang ? "<p>Para obtener las soluciones a la ecuación " : "<p>To obtain solutions to the equation ");
    showEqNbr(equationNbr + 1);
    showText(lang? " debemos calcular los convergentes de la fracción continua de ":
                   " we have to compute the convergents of the continued fraction of ");
    showText("&minus;<var>Q</var> / 2&#8290;<var>P</var> = ");
    CopyBigInt(&ValU, &ValQ);
    BigIntChSign(&ValU);
    BigIntAdd(&ValP, &ValP, &ValV);
    if (ValV.sign == SIGN_NEGATIVE)
    {    // If denominator is negative, change sign of numerator and denominator.
      BigIntChSign(&ValU);
      BigIntChSign(&ValV);
    }
    shownbr(&ValU);
    showText(" / ");
    shownbr(&ValV);
    showText("</p>");
    showText(lang?"<p>La fracción continua es: ": "The continued fraction is: ");
    for (;;)
    {
      floordiv(&ValU, &ValV, &bigTmp);
      shownbr(&bigTmp);
      // Values of U3 and V3 are not used, so they can be overwritten now.
      // Compute new value of ValU and ValV.
      (void)BigIntMultiply(&bigTmp, &ValV, &U3);
      BigIntSubt(&ValU, &U3, &U3);
      CopyBigInt(&ValU, &ValV);
      CopyBigInt(&ValV, &U3);
      if (BigIntIsZero(&ValV))
      {
        break;
      }
      if (coeffNbr == 0)
      {
        showText("+ //");
      }
      else
      {
        showText(", ");
      }
      coeffNbr++;
    }
    showText("// ");
    showEqNbr(equationNbr+2);
    showText("</p>");
  }
  // Compute bound L = sqrt(4P/(-D))
  multint(&U1, &ValP, 4);
  (void)BigIntDivide(&U1, &discr, &U1);
  BigIntChSign(&U1);               // 4P/(-D)
  squareRoot(U1.limbs, ValL.limbs, U1.nbrLimbs, &ValL.nbrLimbs);  // sqrt(4P/(-D))
  ValL.sign = SIGN_POSITIVE;
  intToBigInteger(&U1, 1);         // Initial value of last convergent: 1/0.
  intToBigInteger(&V1, 0);
  intToBigInteger(&U2, 0);         // Initial value of next to last convergent: 0/1.
  intToBigInteger(&V2, 1);
  // Compute continued fraction expansion of U/V = -Q/2P.
  CopyBigInt(&ValU, &ValQ);
  BigIntChSign(&ValU);
  BigIntAdd(&ValP, &ValP, &ValV);
  while (!BigIntIsZero(&ValV))
  {
    getNextConvergent();
    BigIntSubt(&ValL, &V1, &bigTmp);    // Check whether the denominator of convergent exceeds bound.
    if (bigTmp.sign == SIGN_NEGATIVE)
    {
      if (teach)
      {
        showText(lang? "<p>No hay convergentes para los cuales se cumpla ":
          "<p>There are no convergents for which ");
        showEqNbr(equationNbr + 1);
        showText(lang ? ".</p>" : " holds.</p>");
      }
      break;                            // Bound exceeded, so go out.
    }
    // Test whether P*U1^2 + Q*U1*V1 + R*V1^2 = 1.
    (void)BigIntMultiply(&ValP, &U1, &ValO);      // P*U1
    (void)BigIntMultiply(&ValQ, &V1, &bigTmp);    // Q*V1
    BigIntAdd(&ValO, &bigTmp, &ValO);             // P*U1 + Q*V1
    (void)BigIntMultiply(&ValO, &U1, &ValO);      // P*U1^2 + Q*U1*V1
    (void)BigIntMultiply(&ValR, &V1, &bigTmp);    // R*V1
    (void)BigIntMultiply(&bigTmp, &V1, &bigTmp);  // R*V1^2
    BigIntAdd(&ValO, &bigTmp, &ValO);       // P*U1^2 + Q*U1*V1 + R*V1^2
    if ((ValO.sign == SIGN_POSITIVE) && (ValO.nbrLimbs == 1) && (ValO.limbs[0].x == 1))
    {                                       // a*U1^2 + b*U1*V1 + c*V1^2 = 1.
      int D;
      CopyBigInt(&ValH, &U1);
      CopyBigInt(&ValI, &V1);
      ShowSolutionFromConvergent();
      NonSquareDiscrSolution(value);        // (U1, V1)
      D = discr.limbs[0].x;
      if ((discr.nbrLimbs > 1) || (D > 4))
      {                                     // Discriminant is less than -4, go out.
        break;
      }
      if ((D == 3) || (D == 4))
      {                                     // Discriminant is equal to -3 or -4.
        getNextConvergent();
        CopyBigInt(&ValH, &U1);
        CopyBigInt(&ValI, &V1);
        ShowSolutionFromConvergent();
        NonSquareDiscrSolution(value);      // (U1, V1)
        if (D == 3)
        {
          getNextConvergent();
          CopyBigInt(&ValH, &U1);
          CopyBigInt(&ValI, &V1);
          ShowSolutionFromConvergent();
          NonSquareDiscrSolution(value);    // (U1, V1)
        }
        break;
      }
    }
  }
  equationNbr += 3;
}

static void CheckSolutionSquareDiscr(void)
{
  int solutions = 0;
  (void)BigIntDivide(&ValZ, &currentFactor, &ValN);
  if (teach)
  {
    intToBigInteger(&ValJ, 0);
    showText("<li><p>");
    ShowLin(&ValH, &ValI, &ValJ, "X", "Y");
    showText(" = ");
    shownbr(&currentFactor);
    showText(",");
    ShowLin(&ValL, &ValM, &ValJ, "X", "Y");
    showText(" = ");
    shownbr(&ValN);
    showText("</p>");
  }
  // (IL - HM)X = NI - cM
  // (IL - HM)Y = cL - NH
  (void)BigIntMultiply(&ValI, &ValL, &ValO);
  (void)BigIntMultiply(&ValH, &ValM, &bigTmp);
  BigIntSubt(&ValO, &bigTmp, &ValO);           // ValO = Denominator.
  (void)BigIntMultiply(&ValN, &ValI, &ValP);
  (void)BigIntMultiply(&currentFactor, &ValM, &bigTmp);
  BigIntSubt(&ValP, &bigTmp, &ValP);           // ValP = Numerator of X.
  (void)BigIntRemainder(&ValP, &ValO, &U2);
  if (BigIntIsZero(&U2))
  {
    (void)BigIntDivide(&ValP, &ValO, &U1);     // X found.
    (void)BigIntMultiply(&currentFactor, &ValL, &ValP);
    (void)BigIntMultiply(&ValN, &ValH, &bigTmp);
    BigIntSubt(&ValP, &bigTmp, &ValP);         // ValP = Numerator of Y.
    (void)BigIntRemainder(&ValP, &ValO, &U3);
    if (BigIntIsZero(&U3))
    {
      (void)BigIntDivide(&ValP, &ValO, &U2);   // Y found.
      ShowPoint(&U1, &U2);                     // Show results.
      solutions = 1;
    }
  }
  if (teach && (solutions == 0))
  {
    showText(lang ? "<p>El sistema de dos ecuaciones no tiene soluciones <var>X</var>, <var>Y</var> enteras.</p>" :
      "<p>The system of two equations does not have integer solutions <var>X</var>, <var>Y</var>.</p>");
  }
}

static void PerfectSquareDiscriminant(void)
{
  int index;
  enum eSign signTemp;
  char* ptrFactorDec;

  if (BigIntIsZero(&ValA))
  { // Let R = gcd(b, c)
    // (bX + cY) Y = k
    BigIntGcd(&ValB, &ValC, &ValR);
    if (teach)
    {
      intToBigInteger(&ValS, 1);
      CopyBigInt(&V1, &ValB);
      CopyBigInt(&V2, &ValC);
      intToBigInteger(&ValH, 0);
      intToBigInteger(&ValI, 1);
      CopyBigInt(&ValL, &ValK);
    }
  }
  else
  {
    // Multiplying by 4a we get (2aX + (b+g)Y)(2aX + (b-g)Y) = 4ak
    // Let R = gcd(2a, b+g)
    BigIntAdd(&ValA, &ValA, &V1);
    BigIntAdd(&ValB, &ValG, &V2);
    BigIntGcd(&V1, &V2, &ValR);
    // Let S = gcd(2a, b-g)
    CopyBigInt(&ValH, &V1);
    BigIntSubt(&ValB, &ValG, &ValI);
    BigIntGcd(&ValH, &ValI, &ValS);
    // Let L = 4ak
    (void)BigIntMultiply(&ValA, &ValK, &ValL);
    multint(&ValL, &ValL, 4);
    if (teach)
    {
      showText(lang ? "<p>Multiplicando por" : "<p>Multiplying by");
      showText(" 4&#8290;<var>a</var>:</p>");
    }
  }
  if (teach)
  {
    intToBigInteger(&ValJ, 0);
    showText("<p>(");
    ShowLin(&V1, &V2, &ValJ, "X", "Y");
    showText(") &#8290;(");
    ShowLin(&ValH, &ValI, &ValJ, "X", "Y");
    showText(") = ");
    shownbr(&ValL);
    showText("</p><p>");
    (void)BigIntMultiply(&ValR, &ValS, &V3);
    if ((V3.nbrLimbs > 1) || (V3.limbs[0].x > 1))
    {
      (void)BigIntDivide(&V1, &ValR, &V1);
      (void)BigIntDivide(&V2, &ValR, &V2);
      (void)BigIntDivide(&ValH, &ValS, &ValH);
      (void)BigIntDivide(&ValI, &ValS, &ValI);
      (void)BigIntRemainder(&ValL, &V3, &bigTmp);
      if (!BigIntIsZero(&bigTmp))
      {
        shownbr(&V3);
        showText(" &#8290;");
      }
      showText("(");
      ShowLin(&V1, &V2, &ValJ, "X", "Y");
      showText(") &#8290;(");
      ShowLin(&ValH, &ValI, &ValJ, "X", "Y");
      showText(") = ");
      if (BigIntIsZero(&bigTmp))
      {
        (void)BigIntDivide(&ValL, &V3, &ValL);
      }
      shownbr(&ValL);
      showText("</p><p>");
      if (!BigIntIsZero(&bigTmp))
      {
        showText(lang ? "El término derecho no es múltiplo del número que está a la izquierda del paréntesis, por lo que no hay soluciones.</p>" :
          "The right hand side is not multiple of the number located at the left of the parentheses, so there are no solutions.</p>");
      }
    }
  }
  if (BigIntIsZero(&ValK))
  {      // k equals zero.
    enum eLinearSolution ret;
    if (teach)
    {
      showText(lang ? "Para que el producto valga cero, cualquiera de los paréntesis vale cero.</p><p>" :
        "The product is zero, so any of the values inside parenthesis equal zero.</p><p>");
    }
    if (BigIntIsZero(&ValA))
    {    // Coefficient a equals zero.
      if (teach)
      {
        ShowLin(&ValH, &ValI, &ValJ, "X", "Y");
        showText(" = 0</p>");
      }
      // Solve Dy - beta = 0
      intToBigInteger(&Aux[0], 0);
      BigIntChSign(&ValBeta);
      ret = LinearEq(&Aux[0], &discr, &ValBeta);
      BigIntChSign(&ValBeta);
      if ((ret == NO_SOLUTIONS) && teach)
      {
        showText(lang ? "Esta ecuación no tiene soluciones enteras.</p><p>" :
          "This equation does not have integer solutions.</p><p>");
      }
      startResultBox(ret);
      PrintLinear(ret, "t");
      endResultBox(ret);
      // Solve bDx + cDy - b*alpha - c*beta = 0
      (void)BigIntMultiply(&ValB, &discr, &Aux[0]);
      (void)BigIntMultiply(&ValC, &discr, &Aux[1]);
      (void)BigIntMultiply(&ValB, &ValAlpha, &Aux[2]);
      (void)BigIntMultiply(&ValC, &ValBeta, &bigTmp);
      BigIntAdd(&Aux[2], &bigTmp, &Aux[2]);
      BigIntChSign(&Aux[2]);
    }
    else
    {    // Coefficient a does not equal zero.
      // Solve 2aD x + (b+g)D y = 2a*alpha + (b+g)*beta
      (void)BigIntMultiply(&ValA, &discr, &Aux[0]);
      BigIntAdd(&Aux[0], &Aux[0], &Aux[0]);
      BigIntAdd(&ValB, &ValG, &Aux[1]);
      (void)BigIntMultiply(&Aux[1], &discr, &Aux[1]);
      (void)BigIntMultiply(&ValA, &ValAlpha, &Aux[2]);
      BigIntAdd(&Aux[2], &Aux[2], &Aux[2]);
      BigIntAdd(&ValB, &ValG, &bigTmp);
      (void)BigIntMultiply(&bigTmp, &ValBeta, &bigTmp);
      BigIntAdd(&Aux[2], &bigTmp, &Aux[2]);
      BigIntChSign(&Aux[2]);
      if (teach)
      {
        ShowLin(&ValH, &ValI, &ValJ, "X", "Y");
        showText(" = 0</p>");
      }
      ret = LinearEq(&Aux[0], &Aux[1], &Aux[2]);
      if ((ret == NO_SOLUTIONS) && teach)
      {
        showText(lang ? "Esta ecuación no tiene soluciones enteras.</p><p>" :
          "This equation does not have integer solutions.</p><p>");
      }
      startResultBox(ret);
      PrintLinear(ret, "t");
      endResultBox(ret);
      // Solve the equation 2aD x + (b-g)D y = 2a*alpha + (b-g)*beta
      (void)BigIntMultiply(&ValA, &discr, &Aux[0]);
      BigIntAdd(&Aux[0], &Aux[0], &Aux[0]);
      BigIntSubt(&ValB, &ValG, &Aux[1]);
      (void)BigIntMultiply(&Aux[1], &discr, &Aux[1]);
      (void)BigIntMultiply(&ValA, &ValAlpha, &Aux[2]);
      BigIntAdd(&Aux[2], &Aux[2], &Aux[2]);
      BigIntSubt(&ValB, &ValG, &bigTmp);
      (void)BigIntMultiply(&bigTmp, &ValBeta, &bigTmp);
      BigIntAdd(&Aux[2], &bigTmp, &Aux[2]);
      BigIntChSign(&Aux[2]);
    }
    if (teach)
    {
      ShowLin(&V1, &V2, &ValJ, "X", "Y");
      showText(" = 0</p>");
    }
    ret = LinearEq(&Aux[0], &Aux[1], &Aux[2]);
    startResultBox(ret);
    PrintLinear(ret, "t");
    endResultBox(ret);
    return;
  }
  // k does not equal zero.
  if (BigIntIsZero(&ValA))
  { // If R does not divide k, there is no solution.
    CopyBigInt(&U3, &ValK);
    CopyBigInt(&U1, &ValR);
  }
  else
  {
    // If R*S does not divide 4ak, there is no solution.
    (void)BigIntMultiply(&ValR, &ValS, &U1);
    (void)BigIntMultiply(&ValA, &ValK, &U2);
    multadd(&U3, 4, &U2, 0);
  }
  (void)BigIntRemainder(&U3, &U1, &U2);
  if (!BigIntIsZero(&U2))
  {
    return;
  }
  (void)BigIntDivide(&U3, &U1, &ValZ);
  if (teach)
  {
    showText(lang ? "Tenemos que hallar todos los factores del término derecho.</p>" :
      "We have to find all factors of the right hand side.</p>");
  }
  // Compute all factors of Z = 4ak/RS
  NumberLength = ValZ.nbrLimbs;
  signTemp = ValZ.sign;
  ValZ.sign = SIGN_POSITIVE;  // Factor positive number.
  BigInteger2IntArray(nbrToFactor, &ValZ);
  ptrFactorDec = tofactorDec;
  Bin2Dec(&ptrFactorDec, ValZ.limbs, ValZ.nbrLimbs, groupLen);
  factor(&ValZ, nbrToFactor, factorsMod, astFactorsMod);
  CopyBigInt(&LastModulus, &ValZ);           // Do not factor again same modulus.
  ValZ.sign = signTemp;       // Restore sign of Z = 4ak/RS.
  // x = (NI - JM) / D(IL - MH) and y = (JL - NH) / D(IL - MH)
  // The denominator cannot be zero here.
  // H = 2a/R, I = (b+g)/R, J = F + H * alpha + I * beta
  // L = 2a/S, M = (b-g)/S, N = Z/F + L * alpha + M * beta
  // F is any factor of Z (positive or negative).
  nbrFactors = astFactorsMod[0].multiplicity;
  if (teach)
  {
    showText("<p>");
    showFactors(&ValZ);
    showText("<ol>");
  }
  if (BigIntIsZero(&ValA))
  {
    intToBigInteger(&ValH, 0);                    // H <- 0
    intToBigInteger(&ValI, 1);                    // I <- 1
    (void)BigIntDivide(&ValB, &ValR, &ValL);      // L <- b/R
    (void)BigIntDivide(&ValC, &ValR, &ValM);      // M <- c/R
  }
  else
  {
    BigIntAdd(&ValA, &ValA, &U3);                 // 2a
    (void)BigIntDivide(&U3, &ValR, &ValH);        // H <- 2a/R
    (void)BigIntDivide(&U3, &ValS, &ValL);        // L <- 2a/S
    BigIntAdd(&ValB, &ValG, &ValI);
    (void)BigIntDivide(&ValI, &ValR, &ValI);      // I <- (b+g)/R
    BigIntSubt(&ValB, &ValG, &ValM);
    (void)BigIntDivide(&ValM, &ValS, &ValM);      // M <- (b-g)/S
  }
  (void)BigIntMultiply(&ValH, &ValAlpha, &ValK);  // H * alpha
  (void)BigIntMultiply(&ValI, &ValBeta, &bigTmp); // I * beta
  BigIntAdd(&ValK, &bigTmp, &ValK);         // K <- H * alpha + I * beta
  (void)BigIntMultiply(&ValL, &ValAlpha, &ValO);  // L * alpha
  (void)BigIntMultiply(&ValM, &ValBeta, &bigTmp); // M * beta
  BigIntAdd(&ValO, &bigTmp, &ValO);         // O <- L * alpha + M * beta
  // Compute denominator: D(IL - MH)
  (void)BigIntMultiply(&ValI, &ValL, &ValDen);    // IL
  (void)BigIntMultiply(&ValM, &ValH, &bigTmp);    // MH
  BigIntSubt(&ValDen, &bigTmp, &ValDen);    // IL - MH
  (void)BigIntMultiply(&ValDen, &discr, &ValDen); // D(IL - MH)

  // Loop that finds all factors of Z.
  // Use Gray code to use only one big number.
  // Gray code: 0->000, 1->001, 2->011, 3->010, 4->110, 5->111, 6->101, 7->100.
  // Change from zero to one means multiply, otherwise divide.
  (void)memset(counters, 0, sizeof(counters));
  (void)memset(isDescending, 0, sizeof(isDescending));
  intToBigInteger(&currentFactor, 1);
  for (;;)
  {
    const struct sFactors* pstFactor;
    CheckSolutionSquareDiscr();       // Process positive divisor.
    BigIntChSign(&currentFactor);
    CheckSolutionSquareDiscr();       // Process negative divisor.
    BigIntChSign(&currentFactor);
    pstFactor = &astFactorsMod[1];
    for (index = 0; index < nbrFactors; index++)
    {            // Loop that increments counters.
      if (isDescending[index] == 0)
      {          // Ascending.
        if (counters[index] == pstFactor->multiplicity)
        {
          isDescending[index] = 1;    // Next time it will be descending.
          pstFactor++;
          continue;
        }
        NumberLength = *pstFactor->ptrFactor;
        IntArray2BigInteger(pstFactor->ptrFactor, &prime);
        (void)BigIntMultiply(&currentFactor, &prime, &currentFactor);
        counters[index]++;
        break;
      }
               // Descending.
      if (counters[index] == 0)
      {
        isDescending[index] = 0;    // Next time it will be ascending.
        pstFactor++;
        continue;
      }
      NumberLength = *pstFactor->ptrFactor;
      IntArray2BigInteger(pstFactor->ptrFactor, &prime);
      (void)BigIntDivide(&currentFactor, &prime, &currentFactor);
      counters[index]--;
      break;
    }
    if (index == nbrFactors)
    {               // All factors have been found. Exit loop.
      if (teach)
      {
        showText("</ol>");
      }
      break;
    }
  }
}

static void PositiveDiscriminant(void)
{
  callbackQuadModType = CBACK_QMOD_HYPERBOLIC;
  NonSquareDiscriminant();
}

// First check: |u| < g.
// Second check: |u+g| > |v|
// Third check: |u-g| < |v|
// On input ValG = floor(g), g > 0.
// g is not an integer number.
static void CheckStartOfContinuedFractionPeriod(void)
{
  CopyBigInt(&bigTmp, &ValU);
  // Set bigTmp to |u|
  bigTmp.sign = SIGN_POSITIVE;
  // Compute bigTmp as floor(g) - |u|
  BigIntSubt(&ValG, &bigTmp, &bigTmp);
  if (bigTmp.sign == SIGN_POSITIVE)
  { // First check |u| < g passed.
    CopyBigInt(&Tmp1, &ValV);
    // Set Tmp1 to |v|
    Tmp1.sign = SIGN_POSITIVE;
    // Compute Tmp2 as u + floor(g) which equals floor(u+g)
    BigIntAdd(&ValU, &ValG, &Tmp2);
    if (Tmp2.sign == SIGN_NEGATIVE)
    { // Round to number nearer to zero.
      addbigint(&Tmp2, 1);
    }
    // Compute Tmp2 as floor(|u+g|)
    Tmp2.sign = SIGN_POSITIVE;
    // Compute bigTmp as floor(|u+g|) - |v|
    BigIntSubt(&Tmp2, &Tmp1, &bigTmp);
    if (bigTmp.sign == SIGN_POSITIVE)
    { // Second check |u+g| > |v| passed.
      // Conpute Tmp2 as u - floor(g)
      BigIntSubt(&ValU, &ValG, &Tmp2);
      if ((Tmp2.sign == SIGN_NEGATIVE) || (BigIntIsZero(&Tmp2)))
      { // Round down number to integer.
        addbigint(&Tmp2, -1);
      }
      // Compute Tmp2 as floor(|u-g|)
      Tmp2.sign = SIGN_POSITIVE;
      // Compute Tmp2 as |v| - floor(|u-g|)
      BigIntSubt(&Tmp1, &Tmp2, &bigTmp);
      if (bigTmp.sign == SIGN_POSITIVE)
      { // Third check |u-g| < |v| passed.
        // Save U and V to check period end.
        CopyBigInt(&startPeriodU, &ValU);
        CopyBigInt(&startPeriodV, &ValV);
      }
    }
  }
}

static void ShowArgumentContinuedFraction(void)
{
  showText(" (");
  if (positiveDenominator == 0)
  {
    showText("&minus;");
  }
  showText("<var>Q</var> + <span class = \"sqrtout\"><span class=\"sqrtin\"><var>D</var>");
  if ((ValB.limbs[0].x & 1) == 0)
  {
    showText(" / 4");
  }
  showText("</span></span>) / ");
  if ((ValB.limbs[0].x & 1) != 0)
  {
    showText(positiveDenominator?"(": "(&minus;");
    showText("2<var>R</var>) = ");
  }
  else
  {
    showText(positiveDenominator ? "<var>R</var> = " : "(&minus;<var>R</var>) = ");
  }
  if (!BigIntIsZero(&ValU))
  {
    showText("(");
    shownbr(&ValU);
    showText(" + ");
  }
  showText("<span class=\"sqrtout\"><span class=\"sqrtin\">");
  shownbr(&ValL);
  showText("</span></span>");
  if (!BigIntIsZero(&ValU))
  {
    showText(")");
  }
  showText(" / ");
  if (ValV.sign == SIGN_NEGATIVE)
  {
    showText("(");
  }
  shownbr(&ValV);
  if (ValV.sign == SIGN_NEGATIVE)
  {
    showText(")");
  }
}

//  PQa algorithm for (P+G)/Q where G = sqrt(discriminant):
//  Set U1 to 1 and U2 to 0.
//  Set V1 to 0 and V2 to 1.
//  Perform loop:
//  Compute a as floor((U + G)/V)
//  Set U3 to U2, U2 to U1 and U1 to a*U2 + U3
//  Set V3 to V2, V2 to V1 and V1 <- a*V2 + V3
//  Set U to a*V - U
//  Set V to (D - U^2)/V
//  Inside period when: 0 <= G - U < V

static enum eExprErr ContFrac(BigInteger *value, enum eShowSolution solutionNbr)
{
  int periodsToCompute;
  int index = 0;
  int periodIndex = 0;
  char isIntegerPart;
  bool isBeven = ((ValB.limbs[0].x & 1) == 0);
  if (originTranslated)
  {
    periodsToCompute = 2;
  }
  else
  {
    periodsToCompute = 1;
  }
  // If (D-U^2) is not multiple of V, exit routine.
  (void)BigIntMultiply(&ValU, &ValU, &bigTmp); // V <- (D - U^2)/V
  BigIntSubt(&ValL, &bigTmp, &bigTmp);   // D - U^2
  (void)BigIntRemainder(&bigTmp, &ValV, &bigTmp);
  if (!BigIntIsZero(&bigTmp))
  {
    return EXPR_OK;
  }
  CopyBigInt(&Tmp[11], value);   // Back up value.
  if (teach)
  {
    CopyBigInt(&U1, &ValU);      // Back up numerator and denominator.
    CopyBigInt(&V1, &ValV);
    showText(lang ? "<p>El desarrollo en fracciones continuas de " :
      "<p>The continued fraction expansion of ");
    ShowArgumentContinuedFraction();
    showText(lang ? " es:</p>" : " is:</p>");
    intToBigInteger(&startPeriodU, -1);      // Less than zero means outside period.
    index = 0;
    for (;;)
    {
      if (startPeriodU.sign == SIGN_POSITIVE)
      {               // Already inside period.
        periodIndex++;
        if (BigIntEqual(&ValU, &startPeriodU) && BigIntEqual(&ValV, &startPeriodV))
        {             // New period started.
          break;      // Go out in this case.
        }
      }
      else if (index > 0)
      {               // Check if periodic part of continued fraction has started.
        CheckStartOfContinuedFractionPeriod();
        if (startPeriodU.sign == SIGN_POSITIVE)
        {
          showText("<span class=\"bold\">");
        }
      }
      // Get continued fraction coefficient.
      BigIntAdd(&ValU, &ValG, &bigTmp);
      if (ValV.sign == SIGN_NEGATIVE)
      {   // If denominator is negative, round square root upwards.
        addbigint(&bigTmp, 1);
      }
      floordiv(&bigTmp, &ValV, &Tmp1);       // Tmp1 = Term of continued fraction.
      if (index == 1)
      {
        showText("+ // ");
      }
      else if (index > 1)
      {
        showText(", ");
      }
      index++;
      shownbr(&Tmp1);
      // Update numerator and denominator.
      (void)BigIntMultiply(&Tmp1, &ValV, &bigTmp); // U <- a*V - U
      BigIntSubt(&bigTmp, &ValU, &ValU);
      (void)BigIntMultiply(&ValU, &ValU, &bigTmp); // V <- (D - U^2)/V
      BigIntSubt(&ValL, &bigTmp, &bigTmp);
      (void)BigIntDivide(&bigTmp, &ValV, &Tmp1);
      CopyBigInt(&ValV, &Tmp1);
    }
    showText("</span>// ");
    showEqNbr(contfracEqNbr);
    showText("</p>");
    CopyBigInt(&ValU, &U1);      // Restore numerator and denominator.
    CopyBigInt(&ValV, &V1);
  }
  // Initialize variables.
  intToBigInteger(&U1, 1);
  intToBigInteger(&U2, 0);
  intToBigInteger(&V1, 0);
  intToBigInteger(&V2, 1);
  intToBigInteger(&startPeriodU, -1);      // Less than zero means outside period.
  index = 0;
  if (solutionNbr == SECOND_SOLUTION)
  {
    index++;
  }
  isIntegerPart = 1;
  for (;;)
  {
    if ((U1.nbrLimbs > MAX_LEN_MULT / 2) || (U2.nbrLimbs > MAX_LEN_MULT / 2) ||
      (V1.nbrLimbs > MAX_LEN_MULT / 2) || (V2.nbrLimbs > MAX_LEN_MULT / 2))
    {
      return EXPR_NUMBER_TOO_HIGH;
    }
    if ((ValV.nbrLimbs == 1) && (ValV.limbs[0].x == (isBeven ? 1 : 2)) &&
      ((index & 1) == ((ValK.sign == ValV.sign)? 0 : 1)))
    {         // Found solution.
      if ((discr.nbrLimbs == 1) && (discr.limbs[0].x == 5) && (ValA.sign != ValK.sign) && 
        (solutionNbr == FIRST_SOLUTION))
      {       // Discriminant is 5 and aK < 0. Use exceptional solution (U1-U2)/(V1-V2).
        BigIntSubt(&V1, &V2, &ValH);
        BigIntSubt(&U1, &U2, &ValI);
      }
      else
      {       // Discriminant is not 5 or aK > 0. Use convergent U1/V1 as solution.
        CopyBigInt(&ValH, &V1);
        CopyBigInt(&ValI, &U1);
        ShowSolutionFromConvergent();
      }
      showSolution = TWO_SOLUTIONS;
      solFound = false;
      NonSquareDiscrSolution(value);
      if (solFound)
      {
        break;                             // Solution found. Exit loop.
      }
    }
    if (startPeriodU.sign == SIGN_POSITIVE)
    {               // Already inside period.
      periodIndex++;
      if ((BigIntEqual(&ValU, &startPeriodU) && BigIntEqual(&ValV, &startPeriodV)) &&
                   // New period started.
         ((periodIndex & 1) == 0))
      {           // Two periods of period length is odd, one period if even.
        periodsToCompute--;
        if (periodsToCompute == 0)
        {
          break;      // Go out in this case.
        }
      }
    }
    else if (!isIntegerPart)
    {             // Check if periodic part of continued fraction has started.
      CheckStartOfContinuedFractionPeriod();
    }
    // Get continued fraction coefficient.
    BigIntAdd(&ValU, &ValG, &bigTmp);
    if (ValV.sign == SIGN_NEGATIVE)
    {   // If denominator is negative, round square root upwards.
      addbigint(&bigTmp, 1);
    }
    floordiv(&bigTmp, &ValV, &Tmp1);       // Tmp1 = Term of continued fraction.
    // Update convergents.
    CopyBigInt(&U3, &U2);                  // U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
    CopyBigInt(&U2, &U1);
    (void)BigIntMultiply(&Tmp1, &U2, &U1);
    BigIntAdd(&U1, &U3, &U1);
    CopyBigInt(&V3, &V2);                  // V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
    CopyBigInt(&V2, &V1);
    (void)BigIntMultiply(&Tmp1, &V2, &V1);
    BigIntAdd(&V1, &V3, &V1);
    // Update numerator and denominator.
    (void)BigIntMultiply(&Tmp1, &ValV, &bigTmp); // U <- a*V - U
    BigIntSubt(&bigTmp, &ValU, &ValU);
    (void)BigIntMultiply(&ValU, &ValU, &bigTmp); // V <- (D - U^2)/V
    BigIntSubt(&ValL, &bigTmp, &bigTmp);
    (void)BigIntDivide(&bigTmp, &ValV, &Tmp1);
    CopyBigInt(&ValV, &Tmp1);
    index++;
    isIntegerPart = 0;
  }
  CopyBigInt(value, &Tmp[11]);   // Restore value.
  return EXPR_OK;
}

static void ShowRecSol(char variable, const BigInteger *coefX,
                       const BigInteger *coefY, const BigInteger *coefInd)
{
  enum eLinearSolution t;
  *ptrOutput = variable;
  ptrOutput++;
  showText("<sub>n+1</sub> = ");
  t = Show(coefX, "x<sub>n</sub>", SOLUTION_FOUND);
  t = Show(coefY, "y<sub>n</sub>", t);
  Show1(coefInd, t);
}

static void ShowResult(const char *text, const BigInteger *value)
{
  showText(text);
  showText(" = ");
  shownbr(value);
  showText("<br>");
}

static void ShowAllRecSols(void)
{
  if ((ValP.nbrLimbs > 2) || (ValQ.nbrLimbs > 2))
  {
    if (BigIntIsZero(&ValAlpha) && BigIntIsZero(&ValBeta))
    {
      showText("x<sub>n+1</sub> = P&nbsp;&#8290;x<sub>n</sub> + Q&nbsp;&#8290;y<sub>n</sub><br>"
        "y<sub>n+1</sub> = R&nbsp;&#8290;x<sub>n</sub> + S&nbsp;&#8290;y<sub>n</sub></p><p>");
    }
    else
    {
      showText("x<sub>n+1</sub> = P&nbsp;&#8290;x<sub>n</sub> + Q&nbsp;&#8290;y<sub>n</sub> + K<br>"
        "y<sub>n+1</sub> = R&nbsp;&#8290;x<sub>n</sub> + S&nbsp;&#8290;y<sub>n</sub> + L</p><p>");
    }
    showText(lang? "donde:</p><p>": "where:</p><p>");
    ShowResult("P", &ValP);
    ShowResult("Q", &ValQ);
    if (!BigIntIsZero(&ValAlpha) || !BigIntIsZero(&ValBeta))
    {
      ShowResult("K", &ValK);
    }
    ShowResult("R", &ValR);
    ShowResult("S", &ValS);
    if (!BigIntIsZero(&ValAlpha) || !BigIntIsZero(&ValBeta))
    {
      ShowResult("L", &ValL);
    }
  }
  else
  {
    ShowRecSol('x', &ValP, &ValQ, &ValK);
    *ptrOutput = '<';
    ptrOutput++;
    *ptrOutput = 'b';
    ptrOutput++;
    *ptrOutput = 'r';
    ptrOutput++;
    *ptrOutput = '>';
    ptrOutput++;
    ShowRecSol('y', &ValR, &ValS, &ValL);
  }
  // Compute x_{n-1} from x_n and y_n
  // Compute new value of K and L as: Knew <- L*Q - K*S and Lnew <- K*R - L*P
  (void)BigIntMultiply(&ValL, &ValQ, &Tmp1);
  (void)BigIntMultiply(&ValK, &ValS, &bigTmp);
  BigIntSubt(&Tmp1, &bigTmp, &Tmp1);
  (void)BigIntMultiply(&ValK, &ValR, &Tmp2);
  (void)BigIntMultiply(&ValL, &ValP, &bigTmp);
  BigIntSubt(&Tmp2, &bigTmp, &ValL);
  CopyBigInt(&ValK, &Tmp1);
  // Compute new values of P, Q, R and S as: Pnew <- S, Qnew <- -Q, Rnew <- -R, Snew <- P
  BigIntChSign(&ValQ);
  BigIntChSign(&ValR);
  CopyBigInt(&bigTmp, &ValP);
  CopyBigInt(&ValP, &ValS);
  CopyBigInt(&ValS, &bigTmp);
  showText(lang? "<p>y también:</p>": "<p>and also:</p>");
  if ((ValP.nbrLimbs > 2) || (ValQ.nbrLimbs > 2))
  {
    showText("<p>");
    ShowResult("P", &ValP);
    ShowResult("Q", &ValQ);
    if (!BigIntIsZero(&ValAlpha) || !BigIntIsZero(&ValBeta))
    {
      ShowResult("K", &ValK);
    }
    ShowResult("R", &ValR);
    ShowResult("S", &ValS);
    if (!BigIntIsZero(&ValAlpha) || !BigIntIsZero(&ValBeta))
    {
      ShowResult("L", &ValL);
    }
  }
  else
  {
    ShowRecSol('x', &ValP, &ValQ, &ValK);
    *ptrOutput = '<';
    ptrOutput++;
    *ptrOutput = 'b';
    ptrOutput++;
    *ptrOutput = 'r';
    ptrOutput++;
    *ptrOutput = '>';
    ptrOutput++;
    ShowRecSol('y', &ValR, &ValS, &ValL);
  }
  *ptrOutput = '<';
  ptrOutput++;
  *ptrOutput = '/';
  ptrOutput++;
  *ptrOutput = 'p';
  ptrOutput++;
  *ptrOutput = '>';
  ptrOutput++;
}

static bool solutionFoundFromContFraction(bool isBeven, int limbValue)
{
  if (isBeven)
  {
    CopyBigInt(&ValQ, &ValB);
    BigIntDivideBy2(&ValQ);
    (void)BigIntMultiply(&ValQ, &V1, &bigTmp);
    BigIntSubt(&U1, &bigTmp, &ValP);   // P <- r - (b/2)s
    BigIntAdd(&U1, &bigTmp, &ValS);            // S <- r + (b/2)s
  }
  else
  {
    (void)BigIntMultiply(&ValB, &V1, &bigTmp);
    BigIntSubt(&U1, &bigTmp, &ValP);   // P <- r - bs
    BigIntAdd(&U1, &bigTmp, &ValS);            // S <- r + bs
    if (limbValue == 4)
    {
      BigIntDivideBy2(&ValP);               // P <- (r - bs)/2
      BigIntDivideBy2(&ValS);               // S <- (r + bs)/2
    }
  }
  (void)BigIntMultiply(&ValC, &V1, &ValQ);
  BigIntChSign(&ValQ);                      // Q <- -cs
  (void)BigIntMultiply(&ValA, &V1, &ValR);        // R <- as
  if (!isBeven && (limbValue == 1))
  {
    BigIntAdd(&ValQ, &ValQ, &ValQ);         // Q <- -2cs
    BigIntAdd(&ValR, &ValR, &ValR);         // R <- 2as
  }
  (void)BigIntMultiply(&ValAlpha, &ValP, &ValK);
  (void)BigIntMultiply(&ValBeta, &ValQ, &bigTmp);
  BigIntAdd(&ValK, &bigTmp, &ValK);
  (void)BigIntMultiply(&ValAlpha, &ValR, &ValL);
  (void)BigIntMultiply(&ValBeta, &ValS, &bigTmp);
  BigIntAdd(&ValL, &bigTmp, &ValL);
  // Check whether alpha - K and beta - L are multiple of discriminant.
  BigIntSubt(&ValAlpha, &ValK, &bigTmp);
  (void)BigIntRemainder(&bigTmp, &discr, &bigTmp);
  if (BigIntIsZero(&bigTmp))
  {
    BigIntSubt(&ValBeta, &ValL, &bigTmp);
    (void)BigIntRemainder(&bigTmp, &discr, &bigTmp);
    if (BigIntIsZero(&bigTmp))
    {    // Solution found.
      BigIntSubt(&ValAlpha, &ValK, &ValK);
      (void)BigIntDivide(&ValK, &discr, &ValK);
      BigIntSubt(&ValBeta, &ValL, &ValL);
      (void)BigIntDivide(&ValL, &discr, &ValL);
      ShowAllRecSols();
      return true;
    }
  }
  // Check whether alpha + K and beta + L are multiple of discriminant.
  BigIntAdd(&ValAlpha, &ValK, &bigTmp);
  (void)BigIntRemainder(&bigTmp, &discr, &bigTmp);
  if (BigIntIsZero(&bigTmp))
  {
    BigIntAdd(&ValBeta, &ValL, &bigTmp);
    (void)BigIntRemainder(&bigTmp, &discr, &bigTmp);
    if (BigIntIsZero(&bigTmp))
    {    // Solution found.
      BigIntAdd(&ValAlpha, &ValK, &ValK);
      (void)BigIntDivide(&ValK, &discr, &ValK);
      BigIntAdd(&ValBeta, &ValL, &ValL);
      (void)BigIntDivide(&ValL, &discr, &ValL);
      BigIntChSign(&ValP);
      BigIntChSign(&ValQ);
      BigIntChSign(&ValR);
      BigIntChSign(&ValS);
      ShowAllRecSols();
      return true;
    }
  }
  return false;
}

// Use continued fraction of sqrt(B^2-4AC)
// If the discriminant is 5, the method does not work: use 3, 1 and 7, 3.
// If the convergent is r/s we get:
// x(n+1) = Px(n) + Qy(n) + K
// y(n+1) = Rx(n) + Sy(n) + L
// where if b is odd: 
//        P = (r - bs)/2, Q = -cs, R = as, S = (r + bs)/2,
// if b is even:
//        P = r - (b/2)s, Q = -cs, R = as, S = r + (b/2)s,
// in any case:
//        K = (alpha*(1-P) - beta*Q) / D, L = (-alpha*R + beta*(1-S)) / D.
static void recursiveSolution(void)
{
  enum eSign sign = SIGN_POSITIVE;
  bool isBeven = ((ValB.limbs[0].x & 1) == 0);
  int periodNbr = 0;
  int periodLength;
  // Initialize variables.
  intToBigInteger(&ValU, 0);
  intToBigInteger(&ValV, 1);
  CopyBigInt(&ValH, &discr);
  if (isBeven)
  {
    subtractdivide(&ValH, 0, 4);
  }
  (void)BigIntMultiply(&discr, &ValGcdHomog, &discr);
  (void)BigIntMultiply(&discr, &ValGcdHomog, &discr);  // Obtain original discriminant.
  if ((discr.nbrLimbs == 1) && (discr.limbs[0].x == 5))
  {     // Discriminant is 5.
        // Do not use continued fraction because it does not work.
    intToBigInteger(&U1, 3);  // First solution to U1^2 - 5*V1^2 = 4
    intToBigInteger(&V1, 1);
    if (solutionFoundFromContFraction(isBeven, 4))
    {
      return;
    }
    intToBigInteger(&U1, 9);  // First solution to U1^2 - 5*V1^2 = 1
    intToBigInteger(&V1, 4);
    (void)solutionFoundFromContFraction(isBeven, 1);
    return;
  }
  squareRoot(ValH.limbs, ValG.limbs, ValH.nbrLimbs, &ValG.nbrLimbs);
  ValG.sign = SIGN_POSITIVE;          // g <- sqrt(discr).
  intToBigInteger(&U1, 1);
  intToBigInteger(&U2, 0);
  intToBigInteger(&V1, 0);
  intToBigInteger(&V2, 1);
  periodLength = 1;
  if (ValGcdHomog.limbs[0].x != 1)
  {
    periodLength = -1;
    do
    {
      BigIntAdd(&ValU, &ValG, &bigTmp);
      if (ValV.sign == SIGN_NEGATIVE)
      {   // If denominator is negative, round square root upwards.
        addbigint(&bigTmp, 1);
      }
      floordiv(&bigTmp, &ValV, &Tmp1);       // Tmp1 = Term of continued fraction.
      (void)BigIntMultiply(&Tmp1, &ValV, &bigTmp); // U <- a*V - U
      BigIntSubt(&bigTmp, &ValU, &ValU);
      (void)BigIntMultiply(&ValU, &ValU, &bigTmp); // V <- (D - U^2)/V
      BigIntSubt(&ValL, &bigTmp, &bigTmp);
      (void)BigIntDivide(&bigTmp, &ValV, &Tmp1);
      CopyBigInt(&ValV, &Tmp1);
      if (periodLength < 0)
      {
        CopyBigInt(&ValUBak, &ValU);
        CopyBigInt(&ValVBak, &ValV);
      }
      periodLength++;
    } while ((periodLength == 1) || !BigIntEqual(&ValU, &ValUBak) || !BigIntEqual(&ValV, &ValVBak));
    intToBigInteger(&ValU, 0);    // Reset values of U and V.
    intToBigInteger(&ValV, 1);
  }
  showText(lang ? "<p>Soluciones recursivas:</p><p>" : "<p>Recursive solutions:</p><p>");
  CopyBigInt(&ValA, &ValABak);
  CopyBigInt(&ValB, &ValBBak);
  CopyBigInt(&ValC, &ValCBak);
  for (;;)
  {
    if ((U1.nbrLimbs > MAX_LEN_MULT / 2) || (U2.nbrLimbs > MAX_LEN_MULT / 2) ||
      (V1.nbrLimbs > MAX_LEN_MULT / 2) || (V2.nbrLimbs > MAX_LEN_MULT / 2))
    {
      showText(lang ? "<p>Coeficientes demasiado grandes.</p>" :
        "<p>Coefficients are too large.</p>");
      return;
    }
    int limbValue;
    BigIntAdd(&ValU, &ValG, &bigTmp);
    if (ValV.sign == SIGN_NEGATIVE)
    {   // If denominator is negative, round square root upwards.
      addbigint(&bigTmp, 1);
    }
    floordiv(&bigTmp, &ValV, &Tmp1);       // Tmp1 = Term of continued fraction.
    CopyBigInt(&U3, &U2);                  // U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
    CopyBigInt(&U2, &U1);
    (void)BigIntMultiply(&Tmp1, &U2, &U1);
    BigIntAdd(&U1, &U3, &U1);
    CopyBigInt(&V3, &V2);                  // V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
    CopyBigInt(&V2, &V1);
    (void)BigIntMultiply(&Tmp1, &V2, &V1);
    BigIntAdd(&V1, &V3, &V1);
    (void)BigIntMultiply(&Tmp1, &ValV, &bigTmp); // U <- a*V - U
    BigIntSubt(&bigTmp, &ValU, &ValU);
    (void)BigIntMultiply(&ValU, &ValU, &bigTmp); // V <- (D - U^2)/V
    BigIntSubt(&ValH, &bigTmp, &bigTmp);
    (void)BigIntDivide(&bigTmp, &ValV, &Tmp1);
    CopyBigInt(&ValV, &Tmp1);
    if (sign == SIGN_POSITIVE)
    {
      sign = SIGN_NEGATIVE;
    }
    else
    {
      sign = SIGN_POSITIVE;
    }
    // Expecting denominator to be 1 (B even or odd) or 4 (B odd) with correct sign.
    if ((ValV.nbrLimbs != 1) || (ValV.sign != sign))
    {   
      continue;
    }
    limbValue = ValV.limbs[0].x;
    if ((limbValue != 1) && (isBeven || (limbValue != 4)))
    {
      continue;
    }
    periodNbr++;
    if (((periodNbr*periodLength) % ValGcdHomog.limbs[0].x) != 0)
    {
      continue;
    }
    // Found solution from continued fraction.
    if (solutionFoundFromContFraction(isBeven, limbValue))
    {
      return;
    }
  }
}

static void callbackQuadModHyperbolic(BigInteger *value)
{
  enum eExprErr rc;
  bool isBeven = ((ValB.limbs[0].x & 1) == 0);
  positiveDenominator = 1;
  if (!PerformTransformation(value))
  {      // No solutions because gcd(P, Q, R) > 1.
    return;
  }
  // Compute P as floor((2*a*theta + b)/2)
  BigIntAdd(&ValA, &ValA, &ValP);
  (void)BigIntMultiply(&ValP, value, &ValP);
  BigIntAdd(&ValP, &ValB, &ValP);
  subtractdivide(&ValP, ValP.limbs[0].x & 1, 2);
  // Compute Q = a*abs(K)
  CopyBigInt(&ValQ, &ValK);
  ValQ.sign = SIGN_POSITIVE;
  (void)BigIntMultiply(&ValQ, &ValA, &ValQ);
  // Find U, V, L so we can compute the continued fraction expansion of (U+sqrt(L))/V.
  CopyBigInt(&ValL, &discr);
  if (isBeven)
  {
    CopyBigInt(&ValU, &ValP);       // U <- P
    subtractdivide(&ValL, 0, 4);    // Argument of square root is discriminant/4.
    CopyBigInt(&ValV, &ValQ);
  }
  else
  {
    BigIntAdd(&ValP, &ValP, &ValU);
    addbigint(&ValU, 1);            // U <- 2P+1
    BigIntAdd(&ValQ, &ValQ, &ValV); // V <- 2Q
  }
  BigIntChSign(&ValU);
  // If L-U^2 is not multiple of V, there is no solution, so go out.
  (void)BigIntMultiply(&ValU, &ValU, &bigTmp);
  BigIntSubt(&ValL, &bigTmp, &bigTmp);
  (void)BigIntRemainder(&bigTmp, &ValV, &bigTmp);
  if (!BigIntIsZero(&bigTmp))
  {
    if (teach)
    {
      showText(lang ? "<p>No hay soluciones de " : "<p>There are no solutions of ");
      showEqNbr(equationNbr + 1);
      showText(lang ? " usando la fracción continua de" : " using the continued fraction of ");
      ShowArgumentContinuedFraction();
      showText(lang ? " porque " : " because ");
      showText("<var>D</var> &minus; <var>Q</var>");
      showSquare();
      showText(lang ? " no es múltiplo de " : " is not multiple of ");
      if ((ValB.limbs[0].x & 1) != 0)
      {                       // Odd discriminant.
        showText("2");
      }
      showText("<var>R</var></p>");
    }
    equationNbr += 2;
    return;
  }
  // Set G to floor(sqrt(L))
  squareRoot(ValL.limbs, ValG.limbs, ValL.nbrLimbs, &ValG.nbrLimbs);
  ValG.sign = SIGN_POSITIVE;          // g <- sqrt(discr).
  Xplus.nbrLimbs = 0;                 // Invalidate solutions.
  Xminus.nbrLimbs = 0;
  Yplus.nbrLimbs = 0;
  Yminus.nbrLimbs = 0;
  CopyBigInt(&ValUBak, &ValU);
  CopyBigInt(&ValVBak, &ValV);
  contfracEqNbr = equationNbr + 2;
  rc = ContFrac(value, FIRST_SOLUTION);    // Continued fraction of (U+G)/V
  if (rc == EXPR_NUMBER_TOO_HIGH)
  {
    showText(lang ? "<p>Solución con demasiados dígitos.</p>":
       "<p>Solution is too large to be displayed.</p>");
    return;
  }
  positiveDenominator = 0;
  CopyBigInt(&ValU, &ValUBak);
  CopyBigInt(&ValV, &ValVBak);
  BigIntChSign(&ValU);
  BigIntChSign(&ValV);
  contfracEqNbr = equationNbr + 3;
  if ((ValU.limbs[0].x == 3) && (ValV.limbs[0].x == 9))
  {
    contfracEqNbr++;
  }
  rc = ContFrac(value, SECOND_SOLUTION);   // Continued fraction of (-U+G)/(-V)
  if (rc == EXPR_NUMBER_TOO_HIGH)
  {
    showText(lang ? "<p>Solución con demasiados dígitos.</p>" :
      "<p>Solution is too large to be displayed.</p>");
    return;
  }
  showSolution = ONE_SOLUTION;
  if (Xplus.nbrLimbs != 0)
  {
    startResultBox(SOLUTION_FOUND);
    ShowXY(&Xplus, &Yplus);
    endResultBox(SOLUTION_FOUND);
  }
  if (Xminus.nbrLimbs != 0)
  {
    startResultBox(SOLUTION_FOUND);
    ShowXY(&Xminus, &Yminus);
    endResultBox(SOLUTION_FOUND);
  }
  equationNbr += 4;
}

static void PrintQuadEqConst(bool showEquationNbr)
{
  showText("<p>");
  PrintQuad(&ValA, &ValB, &ValC, "<var>X</var>", "<var>Y</var>");
  showText(" = ");
  shownbr(&ValK);
  showText(" ");
  if (showEquationNbr && (ValK.sign == SIGN_POSITIVE))
  {
    showEqNbr(1);
  }
  showText("</p>");
}

void SolveQuadEquation(void)
{
  originTranslated = false;
  also = 0;
  showSolution = ONE_SOLUTION;
  showRecursiveSolution = 0;    // Do not show recursive solution by default.
  divgcd = lang ? "<p>Dividiendo la ecuación por el máximo común divisor obtenemos:</p>" :
    "<p>Dividing the equation by the greatest common divisor we obtain:</p>";
  // Get GCD of A, B, C, D and E.
  BigIntGcd(&ValA, &ValB, &Aux[0]);
  BigIntGcd(&Aux[0], &ValC, &Aux[1]);
  BigIntGcd(&Aux[1], &ValD, &Aux[0]);
  BigIntGcd(&Aux[0], &ValE, &Aux[1]);
  (void)BigIntRemainder(&ValF, &Aux[1], &Aux[0]);
  if (!BigIntIsZero(&Aux[0]))
  {   // F is not multiple of GCD(A, B, C, D, E) so there are no solutions.
    if (teach)
    {
      showText(lang ? "<p>El término independiente no es múltiplo de " : "<p>The constant coefficient is not multiple of ");
      shownbr(&Aux[1]);
      showText(lang ? ", que es el máximo común divisor de los otros coeficientes, así que no hay soluciones.</p>" :
        ", which is the greatest common divisor of the other coefficients, so there are no solutions.</p>");
    }
    else
    {
      showText(lang ? "<p>No hay soluciones.</p>" : "<p>There are no solutions.</p>");
    }
    return;
  }
  // Divide all coefficients by GCD(A, B, C, D, E)
  (void)BigIntDivide(&ValA, &Aux[1], &ValA);
  (void)BigIntDivide(&ValB, &Aux[1], &ValB);
  (void)BigIntDivide(&ValC, &Aux[1], &ValC);
  (void)BigIntDivide(&ValD, &Aux[1], &ValD);
  (void)BigIntDivide(&ValE, &Aux[1], &ValE);
  (void)BigIntDivide(&ValF, &Aux[1], &ValF);
  // Test whether the equation is linear. A = B = C = 0.
  if (BigIntIsZero(&ValA) && BigIntIsZero(&ValB) && BigIntIsZero(&ValC))
  {
    enum eLinearSolution ret = LinearEq(&ValD, &ValE, &ValF);
    startResultBox(ret);
    PrintLinear(ret, "t");
    endResultBox(ret);
    return;
  }
  // Compute discriminant: b^2 - 4ac.
  (void)BigIntMultiply(&ValB, &ValB, &Aux[0]);
  (void)BigIntMultiply(&ValA, &ValC, &Aux[1]);
  multint(&Aux[2], &Aux[1], 4);
  BigIntSubt(&Aux[0], &Aux[2], &discr);
  if (teach)
  {
    showText(lang ? "<p>El discriminante es" : "<p>The discriminant is");
    showText(" <var>D</var> = <var>b</var>");
    showSquare();
    showText("&nbsp;&minus;&nbsp;4&#8290;<var>a</var>&#8290;<var>c</var> = ");
    shownbr(&discr);
    showText("</p>");
  }
  if (BigIntIsZero(&discr))
  {        // Discriminant is zero.
    DiscriminantIsZero();
    return;
  }
  // Compute gcd(a,b,c).
  BigIntGcd(&ValA, &ValB, &bigTmp);
  BigIntGcd(&bigTmp, &ValC, &U1);
  // Discriminant is not zero.
  if (BigIntIsZero(&ValD) && BigIntIsZero(&ValE))
  {   // Do not translate origin.
    intToBigInteger(&ValDiv, 1);
    CopyBigInt(&ValK, &ValF);
    BigIntChSign(&ValK);
    intToBigInteger(&ValAlpha, 0);
    intToBigInteger(&ValBeta, 0);
    if (teach)
    {
      PrintQuadEqConst(true);
    }
  }
  else
  {
    originTranslated = true;
    CopyBigInt(&ValDiv, &discr);
    // Translate the origin (x, y) by (alpha, beta).
    // Compute alpha = 2cd - be
    (void)BigIntMultiply(&ValC, &ValD, &ValAlpha);
    BigIntMultiplyBy2(&ValAlpha);
    (void)BigIntMultiply(&ValB, &ValE, &bigTmp);
    BigIntSubt(&ValAlpha, &bigTmp, &ValAlpha);
    // Compute beta = 2ae - bd
    (void)BigIntMultiply(&ValA, &ValE, &ValBeta);
    BigIntMultiplyBy2(&ValBeta);
    (void)BigIntMultiply(&ValB, &ValD, &bigTmp);
    BigIntSubt(&ValBeta, &bigTmp, &ValBeta);
    // We get the equation ax^2 + bxy + cy^2 = k
    // where k = -D (ae^2 - bed + cd^2 + fD)
    (void)BigIntMultiply(&ValA, &ValE, &ValK);     // ae
    (void)BigIntMultiply(&ValK, &ValE, &ValK);     // ae^2
    (void)BigIntMultiply(&ValB, &ValE, &bigTmp);   // be
    (void)BigIntMultiply(&bigTmp, &ValD, &bigTmp); // bed
    BigIntSubt(&ValK, &bigTmp, &ValK);        // ae^2 - bed
    (void)BigIntMultiply(&ValC, &ValD, &bigTmp);   // cd
    (void)BigIntMultiply(&bigTmp, &ValD, &bigTmp); // cd^2
    BigIntAdd(&ValK, &bigTmp, &ValK);              // ae^2 - bed + cd^2
    (void)BigIntMultiply(&ValF, &discr, &bigTmp);  // fD
    BigIntAdd(&ValK, &bigTmp, &ValK);              // ae^2 - bed + cd^2 + fD
    (void)BigIntMultiply(&ValK, &discr, &ValK);    // D (ae^2 - bed + cd^2 + fD)
    BigIntChSign(&ValK);                           // k
    if (teach)
    {
      showText(lang ? "<p>Aplicamos la transformación de Legendre " :
        "<p>We apply the transformation of Legendre ");
      showText("<var>D</var><var>x</var> = <var>X</var> + <var>&alpha;</var>,  <var>D</var><var>y</var> = <var>Y</var> + <var>&beta;</var>, ");
      showText(lang ? "y obtenemos" : "and we obtain");
      showText(":</p><p><var>&alpha;</var> = 2&#8290;<var>c</var>&#8290;<var>d</var> - <var>b</var>&#8290;<var>e</var> = ");
      shownbr(&ValAlpha);
      showText("</p><p><var>&beta;</var> = 2&#8290;<var>a</var>&#8290;<var>e</var> - <var>b</var>&#8290;<var>d</var> = ");
      shownbr(&ValBeta);
      showText("</p>");
      PrintQuadEqConst(BigIntIsOne(&U1));
      showText("<p>");
      showText(lang ? "donde el término derecho es igual a " : "where the right hand side equals ");
      showText("&minus;<var>D</var> (<var>a</var>&#8290;<var>e</var>");
      showSquare();
      showText(" &minus; <var>b</var>&#8290;<var>e</var>&#8290;<var>d</var> + <var>c</var>&#8290;<var>d</var>");
      showSquare();
      showText(" + <var>f</var>&#8290;<var>D</var>)</p>");
      if (!BigIntIsOne(&U1))
      {
        CopyBigInt(&ValABak, &ValA);
        CopyBigInt(&ValBBak, &ValB);
        CopyBigInt(&ValCBak, &ValC);
        CopyBigInt(&V1, &ValK);
        (void)BigIntDivide(&ValA, &U1, &ValA);
        (void)BigIntDivide(&ValB, &U1, &ValB);
        (void)BigIntDivide(&ValC, &U1, &ValC);
        (void)BigIntDivide(&ValK, &U1, &ValK);
        showText(lang ? "<p>Dividiendo ambos miembros por " : "<p>Dividing both sides by ");
        Bin2Dec(&ptrOutput, U1.limbs, U1.nbrLimbs, groupLen);
        showText(":</p>");
        PrintQuadEqConst(true);
        CopyBigInt(&ValA, &ValABak);
        CopyBigInt(&ValB, &ValBBak);
        CopyBigInt(&ValC, &ValCBak);
        CopyBigInt(&ValK, &V1);
      }  
    }
  }
  // If k is not multiple of gcd(A, B, C), there are no solutions.
  (void)BigIntRemainder(&ValK, &U1, &bigTmp);
  if (!BigIntIsZero(&bigTmp))
  {
    if (teach)
    {
      showText(lang ? "<p>El término derecho no es múltiplo de ": "<p>The right hand side is not multiple of ");
      shownbr(&U1);
      showText(lang ? ", que es el máximo común divisor de los tres coeficientes, así que no hay soluciones" :
        ", which is the greatest common divisor of all three coefficients, so there are no solutions");
    }
    return;     // There are no solutions.
  }
  if (ValK.sign == SIGN_NEGATIVE)
  {
    BigIntChSign(&ValA);
    BigIntChSign(&ValB);
    BigIntChSign(&ValC);
    BigIntChSign(&ValK);
    if (teach)
    {
      showText(lang ? "<p>El algoritmo requiere que el término independiente no sea negativo, así que multiplicamos ambos miembros por &minus;1.</p>":
        "<p>The algorithm requires the constant coefficient to be positive, so we multiply both RHS and LHS by &minus;1.</p>");
      PrintQuadEqConst(true);
    }
  }
  if (discr.sign == SIGN_NEGATIVE)
  {
    NegativeDiscriminant();
    return;
  }
  squareRoot(discr.limbs, ValG.limbs, discr.nbrLimbs, &ValG.nbrLimbs);
  ValG.sign = SIGN_POSITIVE;
  (void)BigIntMultiply(&ValG, &ValG, &bigTmp);
  if (BigIntEqual(&bigTmp, &discr))
  {   // Discriminant is a perfect square.
    PerfectSquareDiscriminant();
    return;
  }
  PositiveDiscriminant();
}

void quadText(char *coefAText, char *coefBText, char *coefCText,
              char *coefDText, char *coefEText, char *coefFText)
{
  int coeffNbr;
  struct stValidateCoeff *pstValidateCoeff = astValidateCoeff;
  astValidateCoeff[0].expression = coefAText;
  astValidateCoeff[1].expression = coefBText;
  astValidateCoeff[2].expression = coefCText;
  astValidateCoeff[3].expression = coefDText;
  astValidateCoeff[4].expression = coefEText;
  astValidateCoeff[5].expression = coefFText;
  ptrOutput = output;
  copyStr(&ptrOutput, "2<p>");
  for (coeffNbr = 0; coeffNbr < NBR_COEFF; coeffNbr++)
  {
    enum eExprErr rc = ComputeExpression(pstValidateCoeff->expression,
                           pstValidateCoeff->bigint);
    if (rc != EXPR_OK)
    {
      copyStr(&ptrOutput, lang ? pstValidateCoeff->textSpanish : pstValidateCoeff->textEnglish);
      textError(&ptrOutput, rc);
      copyStr(&ptrOutput, "</p>");
      break;
    }
    pstValidateCoeff++;
  }
  if (coeffNbr == NBR_COEFF)
  {
    const char* ptrBeginSol;
    showText("<h2>");
    ShowEq(&ValA, &ValB, &ValC, &ValD, &ValE, &ValF, "x", "y");
    showText(" = 0</h2>");
    SolNbr = 0;
    ptrBeginSol = ptrOutput;
    SolveQuadEquation();
    if (ptrBeginSol == ptrOutput)
    {
      showText(lang ? "<p>La ecuación no tiene soluciones enteras.</p>" :
        "<p>The equation does not have integer solutions.</p>");
    }
  }
  showText("<p>");
  showText(lang ? COPYRIGHT_SPANISH: COPYRIGHT_ENGLISH);
  showText("</p>");
}

#if defined(__EMSCRIPTEN__) && !defined(_MSC_VER)
#ifdef __ANDROID__
EXTERNALIZE void doWorkQuad(void)
#else
EXTERNALIZE void doWork(void)
#endif
{
  int flags;
  char *ptrData = inputString;
  char* ptrCoeffA;
  char* ptrCoeffB;
  char *ptrCoeffC;
  char* ptrCoeffD;
  char* ptrCoeffE;
  char *ptrCoeffF;
  groupLen = 0;
  while (*ptrData != ',')
  {
    groupLen = (groupLen * 10) + (*ptrData - '0');
    ptrData++;
  }
  ptrData++;                    // Skip comma.
  flags = *ptrData;
#ifndef lang  
  lang = ((flags & 1)? true: false);
#endif
  teach = ((flags & 2)? true: false);
  ptrCoeffA = ptrData + 2;  // Skip flags and comma.
  ptrCoeffB = ptrCoeffA + (int)strlen(ptrCoeffA) + 1;
  ptrCoeffC = ptrCoeffB + (int)strlen(ptrCoeffB) + 1;
  ptrCoeffD = ptrCoeffC + (int)strlen(ptrCoeffC) + 1;
  ptrCoeffE = ptrCoeffD + (int)strlen(ptrCoeffD) + 1;
  ptrCoeffF = ptrCoeffE + (int)strlen(ptrCoeffE) + 1;
  quadText(ptrCoeffA, ptrCoeffB, ptrCoeffC, ptrCoeffD, ptrCoeffE, ptrCoeffF);
  databack(output);
}
#endif
