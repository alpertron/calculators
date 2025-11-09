//
// This file is part of Alpertron Calculators.
//
// Copyright 2019-2023 Dario Alejandro Alpern
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
#include "string/strings.h"
#include "rootseq.h"
#include "expression.h"
#include "copyStr.h"

#define NBR_COEFF 6

BigInteger Quintic;
BigInteger Quartic;
BigInteger Cubic;
BigInteger Quadratic;
BigInteger Linear;
BigInteger Independent;
BigInteger discr;
BigInteger commonDenom;
BigInteger tmp0;
BigInteger tmp1;
BigInteger tmp2;
BigInteger tmp3;
BigInteger tmp4;
BigInteger tmp5;
BigInteger tmp6;
BigInteger tmp7;
BigRational RatQuartic;
BigRational RatCubic;
BigRational RatQuadratic;
BigRational RatLinear;
BigRational RatIndependent;
BigRational RatDeprCubic;
BigRational RatDeprQuadratic;
BigRational RatDeprLinear;
BigRational RatDeprIndependent;
BigRational RatDiscr;
BigRational RatDelta0;
BigRational RatDelta1;
BigRational RatD;
BigRational Rat1;
BigRational Rat2;
BigRational Rat3;
BigRational Rat4;
BigRational Rat5;
BigRational RatS;
int indexRoot;
const char* ptrPlusMinus;
const char* ptrMinus;
const char* ptrTimes;
const char* ptrSin;
const char* ptrCos;
const char* ptrACos;
const char* ptrTimesPi;
const char* ptrPi;
const char* ptrI;
static int totients[(2 * MAX_DEGREE) + 1];
extern char* ptrOutput;
static char Sine[100];
static char Cosine[100];
static char ArcCosine[100];
static char TimesPi[100];

enum toShow
{
  SHOW_REAL = 0,
  SHOW_IMAG,
};

struct cosine12
{
  char* radicands;
  int denominator;
};

static struct cosine12 aCosine12[] =
{   // S = sqrt next digit, [ = start sqrt, ] = end sqrt
  {"(S6*[5+S5]+S5-1)", 8},          // cos 12°
  {"(S6*[5-S5]+S5+1)", 8},          // cos 24°
  {"(S5+1)", 4},                    // cos 36°
  {"(S6*[5+S5]-S5+1)", 8},          // cos 48°
  {"1", 2},                         // cos 60°
  {"(S5-1)", 4},                    // cos 72°
  {"(S6*[5-S5]-S5-1)", 8},          // cos 84°
};

static struct cosine12 a2Plus2Cosine12[] =
{   // S = sqrt next digits, [ = start sqrt, ] = end sqrt
  {"7+S6*[5+S5]+S5", 4},          // 2 + 2*cos 12°
  {"9+S6*[5-S5]+S5", 4},          // 2 + 2*cos 24°
  {"10+2*S5", 4},                 // 2 + 2*cos 36°
  {"9+S6*[5+S5]-S5", 4},          // 2 + 2*cos 48°
  {"3", 2},                       // 2 + 2*cos 60°
  {"6-2*S5", 4},                  // 2 + 2*cos 72°
  {"7+S6*[5-S5]-S5", 4},          // 2 + 2*cos 84°
};

static struct cosine12 a2Minus2Cosine12[] =
{   // S = sqrt next digits, [ = start sqrt, ] = end sqrt
  {"9-S6*[5+S5]-S5", 4},          // 2 - 2*cos 12°
  {"7-S6*[5-S5]-S5", 4},          // 2 - 2*cos 24°
  {"6-2*S5", 4},                  // 2 - 2*cos 36°
  {"7-S6*[5+S5]+S5", 4},          // 2 - 2*cos 48°
  {"1", 2},                       // 2 - 2*cos 60°
  {"10-2*S5", 4},                 // 2 - 2*cos 72°
  {"9-S6*[5-S5]+S5", 4},          // 2 - 2*cos 84°
};

#if 0
-cos(Pi/17) + (1/16)*(1-sqrt(17)+sqrt(34-2*sqrt(17))+2*sqrt(17+3*sqrt(17)+sqrt(170+38*sqrt(17))))
-cos(2*Pi/17) + (1/16)*(-1+sqrt(17)+sqrt(34-2*sqrt(17))+2*sqrt(17+3*sqrt(17)-sqrt(170+38*sqrt(17))))
-cos(3*Pi/17) + (1/16)*(1+sqrt(17)+sqrt(34+2*sqrt(17))+2*sqrt(17-3*sqrt(17)-sqrt(170-38*sqrt(17))))
-cos(4*Pi/17) + (1/16)*(-1+sqrt(17)-sqrt(34-2*sqrt(17))+2*sqrt(17+3*sqrt(17)+sqrt(170+38*sqrt(17))))
-cos(5*Pi/17) + (1/16)*(1+sqrt(17)+sqrt(34+2*sqrt(17))-2*sqrt(17-3*sqrt(17)-sqrt(170-38*sqrt(17))))
-cos(6*Pi/17) + (1/16)*(-1-sqrt(17)+sqrt(34+2*sqrt(17))+2*sqrt(17-3*sqrt(17)+sqrt(170-38*sqrt(17))))
-cos(7*Pi/17) + (1/16)*(1+sqrt(17)-sqrt(34+2*sqrt(17))+2*sqrt(17-3*sqrt(17)+sqrt(170-38*sqrt(17))))
-cos(8*Pi/17) + (1/16)*(-1+sqrt(17)+sqrt(34-2*sqrt(17))-2*sqrt(17+3*sqrt(17)-sqrt(170+38*sqrt(17))))
#endif


static struct cosine12 aCosine17[] =
{   // S = sqrt next digits, [ = start sqrt, ] = end sqrt
  {"(1-S17+[34-2*S17]+2*[17+3*S17+[170+38*S17]])", 16},   // cos(pi/17)
  {"(-1+S17+[34-2*S17]+2*[17+3*S17-[170+38*S17]])", 16},  // cos(2*pi/17)
  {"(1+S17+[34+2*S17]+2*[17-3*S17-[170-38*S17]])", 16},   // cos(3*pi/17)
  {"(-1+S17-[34-2*S17]+2*[17+3*S17+[170+38*S17]])", 16},  // cos(4*pi/17)
  {"(1+S17+[34+2*S17]-2*[17-3*S17-[170-38*S17]])", 16},   // cos(5*pi/17)
  {"(-1-S17+[34+2*S17]+2*[17-3*S17+[170-38*S17]])", 16},  // cos(6*pi/17)
  {"(1+S17-[34+2*S17]+2*[17-3*S17+[170-38*S17]])", 16},   // cos(7*pi/17)
  {"(-1+S17+[34-2*S17]-2*[17+3*S17-[170+38*S17]])", 16},  // cos(8*pi/17)
};

static struct cosine12 a2PM2Cosine17[] =
{   // S = sqrt next digits, [ = start sqrt, ] = end sqrt
  {"34-2*S17+2*[34-2*S17]+4*[17+3*S17+[170+38*S17]]", 1},  // 2+2*cos(pi/17)
  {"34-2*S17-2*[34-2*S17]-4*[17+3*S17-[170+38*S17]]", 1},  // 2-2*cos(2*pi/17)
  {"34+2*S17+2*[34+2*S17]+4*[17-3*S17-[170-38*S17]]", 1},  // 2+2*cos(3*pi/17)
  {"34-2*S17+2*[34-2*S17]-4*[17+3*S17+[170+38*S17]]", 1},  // 2-2*cos(4*pi/17)
  {"34+2*S17+2*[34+2*S17]-4*[17-3*S17-[170-38*S17]]", 1},  // 2+2*cos(5*pi/17)
  {"34+2*S17-2*[34+2*S17]-4*[17-3*S17+[170-38*S17]]", 1},  // 2-2*cos(6*pi/17)
  {"34+2*S17-2*[34+2*S17]+4*[17-3*S17+[170-38*S17]]", 1},  // 2+2*cos(7*pi/17)
  {"34-2*S17-2*[34-2*S17]+4*[17+3*S17-[170+38*S17]]", 1},  // 2-2*cos(8*pi/17)
};

void startLine(void)
{
  showText("<li>");
  if (pretty == TEX)
  {
    showText("\\bullet\\,\\,");
  }
}

void endLine(void)
{
  if (pretty == TEX)
  {
    showText("\\\\");
  }
  showText("</li>");
}

void showVariable(char **pptrOutput, char letter)
{
  char* ptrOut = *pptrOutput;
  if (pretty == PRETTY_PRINT)
  {
    copyStr(&ptrOut, "<var>");
    *ptrOut = letter;
    ptrOut++;
    copyStr(&ptrOut, "</var>");
  }
  else
  {
    *ptrOut = letter;
    ptrOut++;
  }
  *pptrOutput = ptrOut;
}

void showVarIndex(char letter, int index)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<var>");
    *ptrOutput = letter;
    ptrOutput++;
    showText("</var><sub>");
    int2dec(&ptrOutput, index);
    showText("</sub>");
  }
  else if (pretty == TEX)
  {
    *ptrOutput = letter;
    ptrOutput++;
    showText("_{");
    int2dec(&ptrOutput, index);
    showText("}");
  }
  else
  {
    *ptrOutput = letter;
    ptrOutput++;
    int2dec(&ptrOutput, index);
  }
}

static void showXindex(int index)
{
  showVarIndex('x', index);
}

void showX(int multiplicity)
{
  if (teach)
  {
    showText("<div class=\"outerbox\"><div class=\"box\"><p>");
  }
  else
  {
    startLine();
  }
  if (multiplicity > 2)
  {
    char from[1000];
    char to[1000];
    char* ptrOutputBak = ptrOutput;
    showXindex(indexRoot);
    memcpy(from, ptrOutputBak, ptrOutput - ptrOutputBak);
    from[ptrOutput - ptrOutputBak] = 0;
    ptrOutput = ptrOutputBak;
    indexRoot += multiplicity;
    showXindex(indexRoot - 1);
    memcpy(to, ptrOutputBak, ptrOutput - ptrOutputBak);
    to[ptrOutput - ptrOutputBak] = 0;
    ptrOutput = ptrOutputBak;
    // $1s to $2s = 
    formatString(&ptrOutput, LITERAL_SHOW_X, from, to);
  }
  else
  {    // Multiplicity 1 or 2.
    for (int ctr = 0; ctr < multiplicity; ctr++)
    {
      showXindex(indexRoot);
      indexRoot++;
      showText(" = ");
    }
  }
}

void endShowX(void)
{
  if (teach)
  {
    showText("</p></div></div>");
  }
  else
  {
    endLine();
  }
}

void startSqrt(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<r-2><r-a>");
  }
  else if (pretty == TEX)
  {
    showText("\\sqrt{");
  }
  else
  {
    showText("(");
  }
}

void endSqrt(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("</r-a></r-2>");
  }
  else if (pretty == TEX)
  {
    showText("}");
  }
  else
  {
    showText(")^(1/2)");
  }
}

void startCbrt(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<r-t><span class=\"befrad\" aria-hidden=\"true\">3</span><span class=\"radicand3\">");
  }
  else if (pretty == TEX)
  {
    showText("\\sqrt[3]{");
  }
  else
  {
    showText("(");
  }
}

void endCbrt(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("</span></r-t>");
  }
  else if (pretty == TEX)
  {
    showText("}");
  }
  else
  {
    showText(")^(1/3)");
  }
}

void startParen(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<o-p><p-p>");
  }
  else if (pretty == TEX)
  {
    showText("\\left(");
  }
  else
  {
    *ptrOutput ='(';
    ptrOutput++;
  }
}

void endParen(void)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("</p-p></o-p>");
  }
  else if (pretty == TEX)
  {
    showText("\\right)");
  }
  else
  {
    *ptrOutput =')';
    ptrOutput++;
  }
}

void showPlusSignOn(bool condPlus, int type)
{
  if ((type & TYPE_PM_SPACE_BEFORE) != 0)
  {
    *ptrOutput = ' ';
    ptrOutput++;
  }
  if (condPlus)
  {
    *ptrOutput = '+';
    ptrOutput++;
  }
  else
  {
    showText(ptrMinus);
  }
  if ((type & TYPE_PM_SPACE_AFTER) != 0)
  {
    *ptrOutput = ' ';
    ptrOutput++;
  }
}
void showRatConstants(const char* numerator, const char* denominator)
{
  if (pretty == PRETTY_PRINT)
  {
    showText("<f-f><f-n>");
    showText(numerator);
    showText("</f-n><f-d>");
    showText(denominator);
    showText("</f-d></f-f>");
  }
  else if (pretty == TEX)
  {
    showText("\\frac{");
    showText(numerator);
    showText("}{");
    showText(denominator);
    showText("}");
  }
  else
  {
    showText("(");
    showText(numerator);
    showText("/");
    showText(denominator);
    showText(")");
  }
}

void showTimesPi(char** pptrString)
{
  char* ptrString = *pptrString;
  if (pretty == PARI_GP)
  {
    copyStr(&ptrString, "*Pi");
  }
  else if (pretty == PRETTY_PRINT)
  {
    copyStr(&ptrString, ptrTimesPi);
  }
  else
  {             // TEX
    copyStr(&ptrString, ptrTimes);
    *ptrString = ' ';
    ptrString++;
    copyStr(&ptrString, ptrPi);
  }
  *pptrString = ptrString;
}

static int gcd(int first, int second)
{
  int a = first;
  int b = second;
  while (b != 0)
  {
    int c = a % b;
    a = b;
    b = c;
  }
  return a;
}

// When gcd(first, second) = 1, find Mult1st and Mult2nd such that:
// first * Mult1st + second * Mult2nd = 1.
static void ExtendedGCD(int first, int second, int* pMult1st, int* pMult2nd)
{
  int a = first;
  int b = second;
  int mult1stOld = 1;
  int mult1st = 0;
  int mult2ndOld = 0;
  int mult2nd = 1;

  while (b != 0)
  {
    int quot = a / b;
    int temp = b;
    b = a - (quot * b);
    a = temp;
    temp = mult1st;
    mult1st = mult1stOld - (quot * mult1st);
    mult1stOld = temp;
    temp = mult2nd;
    mult2nd = mult2ndOld - (quot * mult2nd);
    mult2ndOld = temp;
  }
  *pMult1st = mult1stOld;
  *pMult2nd = mult2ndOld;
}

void showRatString(const char* num, const char* den)
{
  if (pretty != PARI_GP)
  {
    showRatConstants(num, den);
  }
  else
  {
    showText("(");
    showText(num);
    showText("/");
    showText(den);
    showText(")");
  }
}

static void ParseExpression(const char* expr)
{
  const char* ptrExpr = expr;
  while (*ptrExpr != 0)
  {
    if (*ptrExpr == 'S')
    {    // Square root of next digits.
      if (pretty != PARI_GP)
      {
        startSqrt();
      }
      for (;;)
      {
        ptrExpr++;
        char c = *ptrExpr;
        if ((c < '0') || (c > '9'))
        {
          break;
        }
        *ptrOutput = c;
        ptrOutput++;
      }
      ptrExpr--;
      if (pretty != PARI_GP)
      {
        endSqrt();
      }
      else
      {
        showText("^(1/2)");
      }
    }
    else if (*ptrExpr == '[')
    {  // Start of square root.
      startSqrt();
    }
    else if (*ptrExpr == ']')
    {  // End of square root.
      endSqrt();
    }
    else if (*ptrExpr == '-')
    {
      showText(ptrMinus);
    }
    else if (*ptrExpr == '*')
    {
      showText(ptrTimes);
    }
    else if (*ptrExpr == '(')
    {
      startParen();
    }
    else if (*ptrExpr == ')')
    {
      endParen();
    }
    else
    {
      *ptrOutput = *ptrExpr;
      ptrOutput++;
    }
    ptrExpr++;
  }
}

// Show cos(num*pi/den) as radicals.
// The output is the sign and the value of the denominator, or zero
// if the result is zero.
static int showRadicals(int numerator, int denominator, int multipl,
  int powerOf2, const char *times)
{
  int multiple = multipl;
  int power2 = powerOf2;
  int num = numerator;
  int den = denominator;
  int arraySigns[20];
  int indexSigns = 0;
  int exprDen;
  int den2 = 2 * den;
  int angle;
  const char* ptrExpr;
  int sign;
  int mult;
  int result;
  char denom[15];
  num = num % den2;  // Convert to range 0 to 360 degrees.
  if (num < 0)
  {
    num += den2;
  }
  if (num == 0)
  {
    return 1;
  }
  if (num == den)
  {
    return -1;
  }
  if (((num * 2) == den) || ((num*2) == (den*3)))
  {
    return 0;
  }
  while (((num % 2) == 0) && ((den % 2) == 0))
  {
    num /= 2;
    den /= 2;
    power2--;
  }
  if (multiple == 1)
  {
    power2 -= 2;
    multiple = 4;
    angle = num % 8;
    if ((angle == 1) || (angle == 7))
    {    // 45 or 315 degrees.
      sign = 1;
    }
    else
    {
      sign = -1;
    }
    if (power2 > 0)
    {
      if (sign > 0)
      {
        ptrExpr = "2+S2";   // 2+2*cos(Pi/4)
      }
      else
      {
        ptrExpr = "2-S2";   // 2-2*cos(Pi/4)
      }
    }
    else
    {
      ptrExpr = "S2";       // cos(Pi/4)
    }
    exprDen = 2;
  }
  else
  {  // Obtain cosine of num*pi/multiple.
    int indexExpr;
    angle = (2*num * 15 / multiple) % 60;
    if (angle < 15)
    {     // Angle between 0 and 90 degrees.
      indexExpr = (angle - 1) / 2;
      sign = 1;
    }
    else if (angle < 30)
    {     // Angle between 90 and 180 degrees.
      indexExpr = (30 - angle - 1) / 2;
      sign = -1;
    }
    else if (angle < 45)
    {     // Angle between 180 and 270 degrees.
      indexExpr = (angle - 30 - 1) / 2;
      sign = -1;
    }
    else
    {     // Angle between 270 and 360 degrees.
      indexExpr = (60 - angle - 1) / 2;
      sign = 1;
    }
    if (power2 > 0)
    {
      if (sign > 0)
      {
        ptrExpr = a2Plus2Cosine12[indexExpr].radicands;
        exprDen = a2Plus2Cosine12[indexExpr].denominator;
      }
      else
      {
        ptrExpr = a2Minus2Cosine12[indexExpr].radicands;
        exprDen = a2Minus2Cosine12[indexExpr].denominator;
      }
    }
    else
    {
      ptrExpr = aCosine12[indexExpr].radicands;
      exprDen = aCosine12[indexExpr].denominator;
    }
  }
  mult = 1;
  for (indexSigns = power2 - 1; indexSigns >= 0; indexSigns--)
  {
    if (indexSigns >= ((int)sizeof(arraySigns) / (int)sizeof(arraySigns[0])))
    {    // This cannot occur.
      return 0;
    }
    arraySigns[indexSigns] = sign;
    angle = num * 15 / multiple % (60*mult);
    if ((angle < (15*mult)) || (angle > (45*mult)))
    {     // Angle between 0 and 90 degrees.
          // or from 270 to 360 degrees.
      sign = 1;
    }
    else
    {     // Angle between 90 and 270 degrees.
      sign = -1;
    }
    mult *= 2;
  }
  for (indexSigns = 0; indexSigns < power2; indexSigns++)
  {
    if (indexSigns == (power2 - 1))
    {
      if ((indexSigns > 0) && (exprDen != 2))
      {
        char* ptrDenom = denom;
        if (exprDen == 1)
        {
          *ptrOutput = '2';
          ptrOutput++;
        }
        else
        {
          int2dec(&ptrDenom, exprDen / 2);
          *ptrDenom = 0;   // End of string.
          showRatString("1", denom);
        }
        showText(ptrTimes);
      }
      startSqrt();
      break;
    }
    startSqrt();
    *ptrOutput = '2';
    ptrOutput++;
    if (arraySigns[indexSigns] < 0)
    {
      *ptrOutput = ' ';
      ptrOutput++;
      showText(ptrMinus);
      *ptrOutput = ' ';
      ptrOutput++;
    }
    else
    {
      showText(" + ");
    }
  }
  // Interpret expression.
  if ((power2 == 0) && (strcmp(ptrExpr, "1") == 0))
  {
    return sign * exprDen;
  }
  ParseExpression(ptrExpr);
  for (indexSigns = 0; indexSigns < power2; indexSigns++)
  {
    endSqrt();
  }
  result = ((power2 > 1)? 2 : exprDen);
  if (strcmp(ptrExpr, "1") && (result != 1))
  {
    showText(times);
  }
  return sign*result;
}

// Show cos(numerator34*Pi/34)
static int showRadicals17(int numerator34)
{
  int index;
  int angle2;
  int angle = numerator34 % 68;        // Convert to range 0 to 360 degrees.
  if (angle < 0)
  {
    angle += 68;
  }
  assert((angle != 17) && (angle != 51)); // It cannot be 90 deg or 270 deg.
  if ((angle % 2) == 0)
  {
    int sign;
    angle /= 2;                       // Show cos(angle*Pi/17).
    if (angle < 9)
    {           // Range 0 to 90 degrees.
      index = angle - 1;
      sign = 1;
    }
    else if (angle < 17)
    {           // Range 90 to 180 degrees.
      index = 17 - angle - 1;
      sign = -1;
    }
    else if (angle < 26)
    {           // Range 180 to 270 degrees.
      index = angle - 17 - 1;
      sign = -1;
    }
    else
    {           // Range 270 to 360 degrees.
      index = 34 - angle - 1;
      sign = 1;
    }
    ParseExpression(aCosine17[index].radicands);
    return aCosine17[index].denominator * sign;
  }
  angle2 = angle % 34;
  if (angle2 < 9)
  {           // Range 0 to 90 degrees. Cosine positive.
    index = angle2 - 1;
  }
  else if (angle2 < 17)
  {           // Range 90 to 180 degrees. Cosine negative.
    index = 17 - angle2 - 1;
  }
  else if (angle2 < 26)
  {           // Range 180 to 270 degrees. Cosine negative.
    index = angle2 - 17 - 1;
  }
  else
  {           // Range 270 to 360 degrees. Cosine positive.
    index = 34 - angle2 - 1;
  }
  startSqrt();
  ParseExpression(a2PM2Cosine17[index].radicands);
  endSqrt();
  if ((angle < 17) || (angle > 51))
  {
    return 8;
  }
  return -8;
}

static void AdjustComponent(int denominator, char* ptrStart, enum toShow toShow,
  int isFirst, const char *realRoot)
{
  int denomin = denominator;
  char beginning[1000];
  char* ptrBeginning = beginning;
  int lenBeginning;
  *ptrBeginning = 0;
  size_t diffPtrs;
  if (denomin == 0)
  {     // Discard all output if result is zero.
    ptrOutput = ptrStart;
    *ptrOutput = 0;
    return;
  }
  if (denomin < 0)
  {
    copyStr(&ptrBeginning, (pretty == PRETTY_PRINT)? "&minus;" : "-");
    denomin = -denomin;    // Make it positive.
  }
  else if ((toShow == SHOW_IMAG) || (isFirst == 0))
  {
    copyStr(&ptrBeginning, " + ");
  }
  else
  {          // Nothing to do.
  }
  if (denomin == 1)
  {
    if (toShow == SHOW_IMAG)
    {
      if (pretty == PARI_GP)
      {
        *ptrBeginning = 'I';
        ptrBeginning++;
      }
      else
      {
        *ptrBeginning = 'i';
        ptrBeginning++;
        *ptrBeginning = ' ';
        ptrBeginning++;
      }
      copyStr(&ptrBeginning, ptrTimes);
    }
  }
  else if (pretty != PARI_GP)
  {
    copyStr(&ptrBeginning, (pretty == TEX)? "\\frac{":
      "<f-f><f-n>");
    // Start numerator.
    if (pretty == PARI_GP)
    {
      *ptrBeginning = ((toShow == SHOW_REAL)? '1' : 'I');
    }
    else
    {
      *ptrBeginning = ((toShow == SHOW_REAL)? '1' : 'i');
    }
    ptrBeginning++;
    // End numerator.
    copyStr(&ptrBeginning, (pretty == TEX)? "}{": "</f-n><f-d>");
    // Start denominator.
    int2dec(&ptrBeginning, denomin);
    // End numerator.
    copyStr(&ptrBeginning, (pretty == TEX)? "}": "</f-d></f-f> ");
    copyStr(&ptrBeginning, ptrTimes);
  }
  else
  {             // PARI-GP
    *ptrBeginning = '(';
    ptrBeginning++;
    *ptrBeginning = ((toShow == SHOW_REAL)? '1' : 'I');
    ptrBeginning++;
    *ptrBeginning = '/';
    ptrBeginning++;
    int2dec(&ptrBeginning, denomin);
    *ptrBeginning = ')';
    ptrBeginning++;
    *ptrBeginning = '*';
    ptrBeginning++;
    *ptrBeginning = 0;          // Add terminator at end of string.
  }
  copyStr(&ptrBeginning, realRoot);
  diffPtrs = ptrBeginning - &beginning[0];
  lenBeginning = (int)diffPtrs;
  if ((*realRoot != 0) && (lenBeginning != 0) && (*ptrStart != 0))
  {
    copyStr(&ptrBeginning, ptrTimes);
  }
  diffPtrs = ptrBeginning - &beginning[0];
  lenBeginning = (int)diffPtrs;
  (void)memmove(ptrStart + lenBeginning, ptrStart, strlen(ptrStart));
  (void)memcpy(ptrStart, beginning, lenBeginning);
  ptrOutput += lenBeginning;
  *ptrOutput = 0;
}

// If den is multiple of 17, compute the extended GCD of 17 and den/17.
// Let the results be A and B respectively. so num*pi/den = num*A*pi/17 +
// num*B*pi/(den/17). Use the formula cos(a+b) = cos a * cos b - sin a * sin b.
static void showComponent(int num, int den, int multipl, int power2,
  enum toShow toShow, const char *realRoot)
{
  int numer = num;
  int denom = den;
  int multiple = multipl;
  int denominCos;
  int den2 = 2 * denom;
  char * ptrStartRadicals = ptrOutput;

  *ptrOutput = 0;
  numer = numer % den2;  // Convert to range 0 to 360 degrees.
  if (numer < 0)
  {
    numer += den2;
  }
  if ((denom % 17) == 0)
  {           // Denominator is multiple of 17.
    int numerator34;
    int numeratorDen;
    int denominCos17;
    denom /= 17;
    multiple /= 17;
    ExtendedGCD(denom, 17, &numerator34, &numeratorDen);
    numerator34 *= 2*numer;
    numeratorDen *= numer;
    // The original angle is multDen*pi/den + mult17*pi/17 = A + B.
    // Show cos(A) * cos(B) - sin(A) * sin(B).
        // Show cos(A).
    denominCos = showRadicals(numeratorDen, denom, multiple, power2, ptrTimes);
    if (denominCos != 0)
    {   // Show cos(B).
      denominCos17 = showRadicals17(numerator34);
      AdjustComponent(denominCos * denominCos17, ptrStartRadicals, toShow, 1, realRoot);
    }
        // Show sin(A).
    ptrStartRadicals = ptrOutput;
    *ptrOutput = 0;
    if ((denom & 0x01) == 0x01)
    {           // denom is odd.
      denominCos = showRadicals(denom - (numeratorDen * 2), denom * 2, multiple, 1, ptrTimes);
    }
    else
    {           // denom is even.
      denominCos = showRadicals((denom / 2) - numeratorDen, denom, multiple, power2, ptrTimes);
    }
    if (denominCos != 0)
    {    // Show sin(B).
      denominCos17 = showRadicals17(17 - numerator34);
      AdjustComponent(-denominCos * denominCos17, ptrStartRadicals, toShow, 0, realRoot);
    }
    return;
  }
  denominCos = showRadicals(numer, denom, multiple, power2, "");
  AdjustComponent(denominCos, ptrStartRadicals, toShow, 1, realRoot);
}

// Try to find a radical expression for cos(num*pi/den) + 
// i*sin(num*pi/den). If it is not possible, return with no output.
static void outputRadicandsForCosSin(int num, int den, const char *realRoot)
{
  // den must be 3^e3 * 5^e5 * 17^e17 * 2^k, where the exponents can be 0 or 1.
  // At this moment we do not consider other Fermat prime numbers.
  int multiple = 1;
  int power2 = 0;
  int quot = den;
  if ((quot % 3) == 0)
  {
    multiple *= 3;
    quot /= 3;
  }
  if ((quot % 5) == 0)
  {
    multiple *= 5;
    quot /= 5;
  }
  if ((quot % 17) == 0)
  {
    multiple *= 17;
    quot /= 17;
  }
  // At this moment quot should be a power of 2.
  while ((quot % 2) == 0)
  {
    power2++;
    quot /= 2;
  }
  if (quot != 1)
  {                   // Denominator is multiple of another prime.
    return;           // Cannot express by radicals.
  }
  showText(" = ");

  showComponent(num, den, multiple, power2, SHOW_REAL, realRoot);
  if (power2 == 0)
  {
    showComponent(den-(num*2), den*2, multiple, 1, SHOW_IMAG, realRoot);
  }
  else
  {
    showComponent((den/2) - num, den, multiple, power2, SHOW_IMAG, realRoot);
  }
}

// Show multiplicand * cos(M) plus i * multiplicand * sin(M)
// where M equals realNum*pi/realDen
static void showTrig(int numerator, int denominator, const char* multiplicand)
{
  char num[300];
  char den[300];
  char* ptrNum = den;
  int2dec(&ptrNum, denominator);
  *ptrNum = 0;  // Include string terminator.
  showText(multiplicand);
  if (*multiplicand != 0)
  {
    if (pretty != PARI_GP)
    {
      *ptrOutput = ' ';
      ptrOutput++;
    }
    showText(ptrTimes);
  }
  showText(ptrCos);
  ptrNum = num;
  if (numerator != 1)
  {
    int2dec(&ptrNum, numerator);
    *ptrOutput = ' ';
    ptrOutput++;
    showTimesPi(&ptrNum);
  }
  else
  {
    copyStr(&ptrNum, ptrPi);
  }
  showRatString(num, den);
  if (pretty == TEX)
  {
    *ptrOutput = '}';
    ptrOutput++;
  }
  showText(" + ");
  showText(ptrI);
  showText(" ");
  showText(ptrTimes);
  showText(multiplicand);
  if (*multiplicand != 0)
  {
    *ptrOutput = ' ';
    ptrOutput++;
    showText(ptrTimes);
    *ptrOutput = ' ';
    ptrOutput++;
  }
  showText(ptrSin);
  if (pretty != PARI_GP)
  {
    showRatConstants(num, den);
    if (pretty == TEX)
    {
      *ptrOutput = '}';
      ptrOutput++;
    }
  }
  else
  {
    showText("(");
    showText(num);
    showText("/");
    showText(den);
    showText(")");
  }
}

static bool TestCyclotomic(const int* ptrPolynomial, int multiplicity, int polyDegree)
{
  int index;
  const int* ptrCoeff;
  int currentDegree;
  if (polyDegree <= 0)
  {
    return false;   // Constant polynomials are not cyclotomic.
  }
  // Polynomial must be palindromic of even degree
  // and the absolute value must be less than 10.
  if ((polyDegree & 0x01) == 0x01)
  {      // Polynomial has odd degree so it cannot be cyclotomic.
    return false;
  }
  ptrCoeff = ptrPolynomial;
  for (index = 0; index <= polyDegree; index++)
  {
    if ((*ptrCoeff != 1) && (*ptrCoeff != -1))
    {            // Too low ot too high.
      return false;
    }
    if (*(ptrCoeff + 1) >= 10)
    {            // Absolute value too high.
      return false;
    }
    ptrCoeff += 2;
  }
  if (totients[0] == 0)
  {    // Table not initialized.
    initializeSmallPrimes(smallPrimes);
    totients[0] = 1;          // indicate table initialized
    totients[2] = 1;
    totients[3] = 2;
    // Compute totient[index] as index*(p1-1)/p1*(p2-1)/p2*...
    // where p_n is a prime factor of index.
    for (index = 4; index <= (2*MAX_DEGREE); index++)
    {
      int prime = 2;
      int primeIndex = 0;
      int quotient = index;
      int totient = index;
      while ((quotient != 1) && ((prime*prime) <= index))
      {
        if (((quotient / prime) * prime) == quotient)
        {   // Prime dividing index was found.
          totient = totient * (prime - 1) / prime;
          do
          {
            quotient /= prime;
          } while (((quotient / prime) * prime) == quotient);
        }
        primeIndex++;
        prime = smallPrimes[primeIndex];
      }
      if (quotient > 1)
      {
        totient = totient * (quotient - 1) / quotient;
      }
      totients[index] = totient;
    }
  }
  // Degree of cyclotomic polynomial n is totient(n).
  for (index = polyDegree; index <= (2 * MAX_DEGREE); index++)
  {
    if (polyDegree != totients[index])
    {   // Degree is not correct.
      continue;
    }
    // Test whether x^degree - 1 divides this polynomial.
    int base[MAX_DEGREE];
    int prod[MAX_DEGREE];
    int lenBytes;
    int* ptrBase = base;
    ptrCoeff = ptrPolynomial;
    for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
    {
      *ptrBase = *ptrCoeff* *(ptrCoeff + 1);
      ptrBase++;
      ptrCoeff += 2;
    }
    lenBytes = polyDegree * (int)sizeof(int);
    (void)memset(prod, 0, lenBytes);
    prod[1] = 1;    // Initialize polynomial to x.
    for (int expon = 1; expon < index; expon++)
    {               // Multiply by x.
      int coeffMaxDegree = prod[polyDegree - 1];
      for (currentDegree = polyDegree - 1; currentDegree > 0; currentDegree--)
      {
        prod[currentDegree] = prod[currentDegree-1] - (coeffMaxDegree * base[currentDegree]);
      }
      prod[0] = -coeffMaxDegree * base[0];
    }
    // If result is 1, then the polynomial is cyclotomic.
    prod[0]--;
    lenBytes = polyDegree * (int)sizeof(int);
    (void)memset(base, 0, lenBytes);
    if (memcmp(base, prod, lenBytes) == 0)
    {     // The polynomial is cyclotomic.
      bool denIsOdd = ((index & 0x01) == 0x01)? true: false;
      int realDen = denIsOdd ? index : (index / 2);
      for (int numerator = 1; numerator < index; numerator++)
      {
        int realNum;
        if (gcd(numerator, index) > 1)
        {        // Numbers are not coprime. Next loop.
          continue;
        }
        showX(multiplicity);
        if (denIsOdd)
        {
          realNum = 2 * numerator;
        }
        else
        {
          realNum = numerator;
        }
        // Show cos(M) + i sin(M)
        // where M equals realNum*pi/realDen
        showTrig(realNum, realDen, "");
        outputRadicandsForCosSin(realNum, realDen, "");
        endShowX();
      }
      return true;
    }
  }
  return false;         // polynomial is not cyclotomic.
}

static bool isPalindromic(int* ptrPolynomial, int polyDegree)
{
  int currentDegree;
  int *ptrs[MAX_DEGREE+1];
  int* ptrCoeff = ptrPolynomial;
  // Get pointers to coefficients.
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    ptrs[currentDegree] = ptrCoeff;
    ptrCoeff += numLimbs(ptrCoeff);
    ptrCoeff++;
  }
  for (currentDegree = 0; (currentDegree * 2) < polyDegree; currentDegree++)
  {
    if (*ptrs[currentDegree] != *ptrs[polyDegree - currentDegree])
    {
      return false;    // Polynomial is not palindromic.
    }
    if (memcmp(ptrs[currentDegree] + 1, ptrs[polyDegree - currentDegree] + 1,
      numLimbs(ptrs[currentDegree])) != 0)
    {
      return false;    // Polynomial is not palindromic.
    }
  }
  return true;        // Polynomial is palindromic.
}

void StartRadicand(int polyDegree)
{
  if (pretty == PRETTY_PRINT)
  {
    if (polyDegree < 10)
    {       // degree has 1 digit.
      copyStr(&ptrOutput, "<r-t>");
    }
    else if (polyDegree < 100)
    {       // degree has 2 digits.
      copyStr(&ptrOutput, "<span class=\"root2dig\">");
    }
    else
    {       // degree has 3 digits.
      copyStr(&ptrOutput, "<span class=\"root3dig\">");
    }
    copyStr(&ptrOutput, "<span class=\"befrad\">");
    int2dec(&ptrOutput, polyDegree);
    copyStr(&ptrOutput, "</span><r-d>");
  }
  if (pretty == TEX)
  {
    copyStr(&ptrOutput, "\\sqrt[");
    int2dec(&ptrOutput, polyDegree);
    copyStr(&ptrOutput, "]{");
  }
}

void EndRadicand(int polyDegree)
{
  if (pretty == PRETTY_PRINT)
  {
    if (polyDegree < 10)
    {       // degree has 1 digit.
      copyStr(&ptrOutput, "</r-d></r-t>");
    }
    else
    {
      copyStr(&ptrOutput, "</span></span>");
    }
  }
  else if (pretty == TEX)
  {
    *ptrOutput = '}';
    ptrOutput++;
  }
  else
  {
    *ptrOutput = '^';
    ptrOutput++;
    *ptrOutput = '(';
    ptrOutput++;
    *ptrOutput = '1';
    ptrOutput++;
    *ptrOutput = '/';
    ptrOutput++;
    int2dec(&ptrOutput, polyDegree);
    *ptrOutput = ')';
    ptrOutput++;
  }
  *ptrOutput = 0;
}

void GenerateRoots(int multiplicity, const char* rationalRoot,
  bool isNegative, int polyDegree)
{
  for (int currentDegree = 0; currentDegree < polyDegree; currentDegree++)
  {
    int realNum;
    int realDen;
    char realRoot[30000];
    char* ptrOutputBak = ptrOutput;
    ptrOutput = realRoot;
    int numer = (2 * currentDegree) + (isNegative ? 1 : 0);
    int gcdNumDen = gcd(numer, polyDegree);
    realNum = numer / gcdNumDen;
    realDen = polyDegree / gcdNumDen;
    StartRadicand(polyDegree);
    if (pretty == PARI_GP)
    {
      *ptrOutput = '(';
      ptrOutput++;
    }
    copyStr(&ptrOutput, rationalRoot);
    if (pretty == PARI_GP)
    {
      *ptrOutput = ')';
      ptrOutput++;
    }
    EndRadicand(polyDegree);
    *ptrOutput = 0;    // Mark end of string.
    ptrOutput = ptrOutputBak;
    showX(multiplicity);
    showTrig(realNum, realDen, realRoot);
    outputRadicandsForCosSin(realNum, realDen, realRoot);
    endShowX();
  }
}

// Save factor degrees sorting by ascending order.
static void SaveFactorDegrees(int prime, int *factors, int nbrFactors)
{
  int* piFactors = factors;
  int* ptrFactors;
  int* ptrOldFactor;
  int tmpDegree;

  *piFactors = prime;
  piFactors++;
  *piFactors = nbrFactors;
  piFactors++;
  ptrFactors = piFactors;
  for (int currentFactor = 0; currentFactor < nbrFactors; currentFactor++)
  {
    int newDegree = factorInfo[currentFactor].degree;
    for (ptrOldFactor = ptrFactors; ptrOldFactor < piFactors; ptrOldFactor++)
    {
      if (newDegree < *ptrOldFactor)
      {
        for (; ptrOldFactor < piFactors; ptrOldFactor++)
        {
          tmpDegree = *ptrOldFactor;
          *ptrOldFactor = newDegree;
          newDegree = tmpDegree;
        }
        break;
      }
    }
    *ptrOldFactor = newDegree;
    piFactors++;
  }
}

static void showDegrees(const int *factors)
{
  const int* piFactors = factors;
  piFactors++;
  int nbrFactors = *piFactors;
  for (int currentFactor = 0; currentFactor < nbrFactors; currentFactor++)
  {
    if (currentFactor != 0)
    {
      if (currentFactor == (nbrFactors - 1))
      {
        // and
        showText(LITERAL_SHOW_DEGREES);
      }
      else
      {
        *ptrOutput = ',';
        ptrOutput++;
        *ptrOutput = ' ';
        ptrOutput++;
      }
    }
    piFactors++;
    int2dec(&ptrOutput, *piFactors);
  }
}

static void ShowNoSolvable(const int* firstArray, const int* secondArray,
  int numberDifferentX, int nbrFactor, int gcdDegrees)
{
  if (firstArray != secondArray)
  {
    if (secondArray == NULL)
    {
      showX(numberDifferentX);
      *(ptrOutput - 2) = ':';  // Replace equal sign by colon.
      // The roots of the polynomial $1?1+?number $1d ??cannot be expressed by radicals.
      formatString(&ptrOutput, LITERAL_SHOW_NO_SOLVABLE1, nbrFactor + 1);
      if (gcdDegrees > 1)
      {
        // We set $1v = 
        formatString(&ptrOutput, LITERAL_SYMM_OR_ALTER2, 'y');
        showPowerX(&ptrOutput, gcdDegrees);
        showText(". ");
      }
      // The degrees of the factors of polynomial modulo $1d are
      formatString(&ptrOutput, LITERAL_SHOW_NO_SOLVABLE2, *firstArray);
    }
    else
    {
      // and the degrees of the factors of polynomial modulo $1d are
      formatString(&ptrOutput, LITERAL_SHOW_NO_SOLVABLE3, *firstArray);
    }
    showDegrees(firstArray);
  }
  // (the Galois group contains a cycle of 
  showText(LITERAL_SHOW_NO_SOLVABLE4);
}

static void showExplanation(int left, const char* oper1, int middle, const char* oper2, int right)
{
  *ptrOutput = ' ';
  ptrOutput++;
  *ptrOutput = '(';
  ptrOutput++;
  int2dec(&ptrOutput, left);
  *ptrOutput = ' ';
  ptrOutput++;
  showText(oper1);
  *ptrOutput = ' ';
  ptrOutput++;
  int2dec(&ptrOutput, middle);
  *ptrOutput = ' ';
  ptrOutput++;
  showText(oper2);
  *ptrOutput = ' ';
  ptrOutput++;
  int2dec(&ptrOutput, right);
  *ptrOutput = ')';
  ptrOutput++;
}

static bool isPrime(int value)
{
  int divisor = 3;
  if (value == 2)
  {
    return true;
  }
  if ((value & 1) == 0)
  {
    return false;        // Even value different from 2: composite.
  }
  while ((divisor * divisor) <= value)
  {
    if ((value % divisor) == 0)
    {
      return false;      // Composite.
    }
    divisor += 2;
  }
  return true;           // Prime value.
}

// If polynomial is S_n or A_n, indicate that the roots are not solvable.
// Factor the polynomial modulo different primes less than 100.
// If with the same or different factorizations, there is a factor of
// degree 2 or 3 and there is no factor of degree multiple of 2 or 3 respectively
// and if there is a factor with prime degree greater than half the degree of
// the polynomial, the polynomial is not solvable with radicals.
// Discard the primes such that the factorizations give duplicated roots.
// If the factor of prime degree is n/2 < p < n-2, or if n/3 < p <= n/2 and
// the largest degree is odd and greater than n/2, then the polynomial is not
// solvable with radicals.
static bool isSymmetricOrAlternating(int nbrFactor, const int* ptrPolynomial, 
  int polyDeg, int multiplicity)
{
  int polyDegree = polyDeg;
  int factorDegreesCycle2Or3[MAX_DEGREE+2];
  int factorDegreesCycleP[MAX_DEGREE + 2];  // n/2 < p
  int factorDegreesCycleOther[MAX_DEGREE + 2]; // n/3 < p < n/2
  int cycle2Found = 0;
  int cycle3Found = 0;
  int cyclePrGtNOver2Found = 0;
  int cyclePrGtNOver2ToLess2Found = 0;
  int cycleOddGtNOver2Found = 0;
  int cyclePrGtNOver3Found = 0;
  int prime;
  int primeIndex = -1;
  int gcdDegrees = 0;
  int currentDegree;
  const int* ptrCoeff;
  int* ptrCoeffDest;
  factorDegreesCycleP[1] = 0;
  // Compute GCD of all degrees of coefficients different from zero.
  // Generate polynomial common.poly.polyNonRepeatedFactors stripping all zero
  // coefficients with degree not multiple of this GCD.
  ptrCoeff = ptrPolynomial;
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {
    if ((*ptrCoeff != 1) || (*(ptrCoeff + 1) != 0))
    {            // Coefficient is not zero.
      gcdDegrees = gcd(currentDegree, gcdDegrees);
      if (gcdDegrees == 1)
      {          // GCD of degrees is one. Further loops will not change this.
        break;
      }
    }
    ptrCoeff += numLimbs(ptrCoeff);
    ptrCoeff++;
  }
  if (gcdDegrees == 0)
  {    // No gcd computed. Go out.
    return false;
  }
  polyDegree /= gcdDegrees;
  if (polyDegree < 5)
  {     // Polynomial is solvable with radical expressions.
    showX(multiplicity * polyDegree * gcdDegrees);
    *(ptrOutput - 2) = ':';  // Replace equal sign by colon.
    // The roots of the polynomial $1?1+?number $1d ??can be expressed by radicals.
    formatString(&ptrOutput, LITERAL_SYMM_OR_ALTER1, nbrFactor + 1);
    if (gcdDegrees > 1)
    {
      // We set $1v =
      formatString(&ptrOutput, LITERAL_SYMM_OR_ALTER2, 'y');
      showPowerX(&ptrOutput, gcdDegrees);
      showText(". ");
    }
    // The polynomial has degree $1d which is less than 5.
    formatString(&ptrOutput, LITERAL_SYMM_OR_ALTER3, polyDegree);
    return true;
  }
  ptrCoeff = ptrPolynomial;
  ptrCoeffDest = common.poly.polyNonRepeatedFactors;
  *ptrCoeffDest = polyDegree;
  ptrCoeffDest++;
  for (currentDegree = 0; currentDegree <= polyDegree; currentDegree++)
  {    // Copy coefficient and skip zero coefficients.
    int nbrLen = 1 + numLimbs(ptrCoeff);
    int lenBytes = nbrLen * (int)sizeof(int);
    int offset = nbrLen + (2 * (gcdDegrees - 1));
    (void)memcpy(ptrCoeffDest, ptrCoeff, lenBytes);
    ptrCoeff += offset;
    ptrCoeffDest += nbrLen;
  }
  do
  {
    int nbrFactors;
    int cycle2FoundInThisFactor = 0;
    int cycle3FoundInThisFactor = 0;
    const struct sFactorInfo* pstFactorInfo;
    (void)memset(factorInfo, 0, sizeof(factorInfo));
    primeIndex = getNextPrimeNoDuplicatedFactors(primeIndex);
    prime = smallPrimes[primeIndex];
    FactorPolynomialModPrime(prime);
    pstFactorInfo = factorInfo;
    // Check whether there is a factor of degree 2, 3 or prime
    // greater than half the degree of the polynomial being analyzed.
    for (nbrFactors = 0; nbrFactors < MAX_DEGREE; nbrFactors++)
    {
      if (pstFactorInfo->ptr == NULL)
      {    // No more factors.
        break;
      }
      pstFactorInfo++;
    }
    pstFactorInfo = factorInfo;
    for (int factorNbr = 0; factorNbr < nbrFactors; factorNbr++)
    {
      int currDegree;
      currDegree = pstFactorInfo->degree;
      if ((currDegree % 2) == 0)
      {
        if (currDegree == 2)
        {
          cycle2FoundInThisFactor++;
        }
        else
        {
          cycle2FoundInThisFactor = 2;
        }
      }
      if ((currDegree % 3) == 0)
      {
        if (currDegree == 3)
        {
          cycle3FoundInThisFactor++;
        }
        else
        {
          cycle3FoundInThisFactor = 2;
        }
      }
      if ((currDegree > 3) && (currDegree > (degree / 2)))
      {
        if (isPrime(currDegree))
        {      // Current degree > n/2 and is prime.
          if (currDegree < (degree - 2))
          {
            SaveFactorDegrees(prime, factorDegreesCycleP, nbrFactors);
            cyclePrGtNOver2ToLess2Found = currDegree;
          }
          else if (cyclePrGtNOver2Found == 0)
          {    // If first condition holds, the polynomial is not solvable.
            SaveFactorDegrees(prime, factorDegreesCycleP, nbrFactors);
            cyclePrGtNOver2Found = currDegree;
            cycleOddGtNOver2Found = currDegree;
          }
          else
          {              // Nothing to do.
          }
        }
        else if (((currDegree & 0x01) == 0x01) && (currDegree < degree))
        {      // Current degree > n/2, it is odd and less than the polynomial degree.
          if (cycleOddGtNOver2Found == 0)
          {    // If first condition holds, the polynomial is not solvable.
               // Test that this degree is coprime to all other degrees of
               // factors of this polynomial.
            int factNbr;
            pstFactorInfo = factorInfo;
            for (factNbr = 0; factNbr < nbrFactors; factNbr++)
            {
              if ((pstFactorInfo->degree < currDegree) && 
                (gcd(pstFactorInfo->degree, currDegree) != 1))
              {
                break;
              }
              pstFactorInfo++;
            }
            if (factNbr == nbrFactors)
            {      // All degrees are coprime to current degree.
              SaveFactorDegrees(prime, factorDegreesCycleP, nbrFactors);
              cycleOddGtNOver2Found = currDegree;
            }
          }
        }
        else
        {            // Nothing to do.
        }
      }
      else if ((currDegree > 3) && ((currDegree > degree) / 3))
      {
        if (isPrime(currDegree) && ((degree % currDegree) != 0))
        {      // Current degree > n/3, it is prime, and it does not divide the degree.
               // Ensure that only this degree is multiple of itself.
          int nbrMultiples = 0;
          pstFactorInfo = factorInfo;
          for (int factNbr = 0; factNbr < nbrFactors; factNbr++)
          {
            if ((pstFactorInfo->degree % currDegree) == 0)
            {
              nbrMultiples++;
            }
            pstFactorInfo++;
          }
          if (nbrMultiples == 1)
          {           // Only one multiple of currentDegree expected.
            SaveFactorDegrees(prime, factorDegreesCycleOther, nbrFactors);
            cyclePrGtNOver3Found = currDegree;
            if ((degree & 0x01) == 0x01)
            {         // If degree is odd, the group is very transitive.
              break;
            }
          }
        }
      }
      else
      {                // Nothing to do.
      }
      pstFactorInfo++;
    }
    if ((cycle2Found == 0) && (cycle3Found == 0))
    {
      if (cycle2FoundInThisFactor == 1)
      {
        cycle2Found = 1;
        SaveFactorDegrees(prime, factorDegreesCycle2Or3, nbrFactors);
      }
      else if (cycle3FoundInThisFactor == 1)
      {
        cycle3Found = 1;
        SaveFactorDegrees(prime, factorDegreesCycle2Or3, nbrFactors);
      }
      else
      {   // No more cases.
      }
    }
    if (cyclePrGtNOver2ToLess2Found != 0)
    {           // Group is very transitive.
      break;
    }
    if ((cyclePrGtNOver3Found != 0) && 
       (((degree & 0x01) == 0x01) || (cycleOddGtNOver2Found != 0) ||
         (cyclePrGtNOver2Found != 0)))
    {           // Group is very transitive.
      break;
    }
    if (((cycle2Found != 0) || (cycle3Found != 0)) && (cyclePrGtNOver2Found != 0))
    {           // Polynomial is not solvable with radicals. Exit loop.
      break;
    }
  } while (prime < 100);
  int numberDifferentX = multiplicity * polyDegree * gcdDegrees;
  if (cyclePrGtNOver2ToLess2Found != 0)
  {      // Group is very transitive.
    ShowNoSolvable(factorDegreesCycleP, NULL, numberDifferentX,
      nbrFactor, gcdDegrees);
    // prime length greater than half the degree of polynomial
    showText(LITERAL_SYMM_OR_ALTER4);
    showExplanation(cyclePrGtNOver2ToLess2Found, "&gt;", degree, "&divide;", 2);
    // and less than the degree minus 2
    showText(LITERAL_SYMM_OR_ALTER5);
    showExplanation(cyclePrGtNOver2ToLess2Found, "&lt;", degree, "&minus;", 2);
    *ptrOutput = ')';
    ptrOutput++;
  }
  else if ((cyclePrGtNOver3Found != 0) && ((degree & 0x01) == 0x01))
  {     // Group is very transitive.
    ShowNoSolvable(factorDegreesCycleOther, NULL, numberDifferentX,
      nbrFactor, gcdDegrees);
    // prime length greater than a third of the degree of the polynomial
    showText(LITERAL_SYMM_OR_ALTER6);
    showExplanation(cyclePrGtNOver3Found, "&gt;", degree, "&divide;", 3);
    // and the latter ($1d) is odd)
    formatString(&ptrOutput, LITERAL_SYMM_OR_ALTER7, degree);
  }
  else if ((cyclePrGtNOver3Found != 0) && (cycleOddGtNOver2Found != 0))
  {     // Group is very transitive.
    ShowNoSolvable(factorDegreesCycleOther, NULL, numberDifferentX,
      nbrFactor, gcdDegrees);
    // prime length greater than a third of the degree of the polynomial
    showText(LITERAL_SYMM_OR_ALTER6);
    showExplanation(cyclePrGtNOver3Found, "&gt;", degree, "&divide;", 3);
    *ptrOutput = ')';
    ptrOutput++;
    ShowNoSolvable(factorDegreesCycleP, factorDegreesCycleOther, numberDifferentX,
      nbrFactor, gcdDegrees);
    // odd length greater than half the degree of the polynomial
    showText(LITERAL_SYMM_OR_ALTER8);
    showExplanation(cycleOddGtNOver2Found, "&gt;", degree, "&divide;", 2);
    *ptrOutput = ')';
    ptrOutput++;
  }
  else if (((cycle2Found != 0) || (cycle3Found != 0)) && (cyclePrGtNOver2Found != 0))
  {
    ShowNoSolvable(factorDegreesCycle2Or3, NULL, numberDifferentX,
      nbrFactor, gcdDegrees);
    // length 2)    length 3)
    copyStr(&ptrOutput, (cycle2Found ? LITERAL_SYMM_OR_ALTER9 : LITERAL_SYMM_OR_ALTER10));
    ShowNoSolvable(factorDegreesCycleP, factorDegreesCycle2Or3, numberDifferentX,
      nbrFactor, gcdDegrees);
    // prime length greater than half the degree of polynomial
    showText(LITERAL_SYMM_OR_ALTER4);
    showExplanation(cyclePrGtNOver2Found, "&gt;", degree, "&divide;", 2);
    *ptrOutput = ')';
    ptrOutput++;
  }
  else
  {
    return false;
  }
  return true;
}

void getRootsPolynomial(int nbrFactor, char **pptrOutput, struct sFactorInfo* pstFactorInfo, int groupLength)
{
  int nbrFactorsFoundBak = nbrFactorsFound;
  static struct sFactorInfo factorInfoIntegerBak[MAX_DEGREE];
  static int polyIntegerBak[1000000];
  int multiplicity = pstFactorInfo->multiplicity;
  if (pretty == PRETTY_PRINT)
  {
    char *ptr = Sine;
    formatString(&ptr, "<span role=\"img\" aria-label=\" $1s \">sin</span>",
      LITERAL_SINE);
    ptrSin = Sine;
    ptr = Cosine;
    formatString(&ptr, "<span role=\"img\" aria-label=\" $1s \">cos</span>",
      LITERAL_COSINE);
    ptrCos = Cosine;
    ptr = ArcCosine;
    formatString(&ptr, "<span role=\"img\" aria-label=\" $1s \">arc cos</span>",
      LITERAL_ARC_COSINE);
    ptrACos = ArcCosine;
  }
  else if (pretty == TEX)
  {
    ptrSin = "\\sin{";
    ptrCos = "\\cos{";
    ptrACos = "\\arccos";
  }
  else
  {
    ptrSin = "sin";
    ptrCos = "cos";
    ptrACos = "acos";
  }
  groupLen = groupLength;
  ptrOutput = *pptrOutput;
  if (pretty == PRETTY_PRINT)
  {
    ptrMinus = "&minus;";
    ptrTimes = "<o-t></o-t>";
    ptrTimesPi = TimesPi;
    char* pTimesPi = TimesPi;
    formatString(&pTimesPi, "<span role=\"img\" aria-label=\" $1s\">π</span>",
      LITERAL_TIMES_PI);
    ptrTimesPi = TimesPi;
    ptrPi = "<span role=\"img\" aria-label=\"pi\">π</span>";
    ptrI = "i";
  }
  else if (pretty == TEX)
  {
    ptrMinus = "-";
    ptrTimes = "";  // Do not show multiplication sign in TeX.
    ptrPi = "\\pi ";
    ptrI = "i";
  }
  else
  {
    ptrMinus = "-";
    ptrTimes = "*";
    ptrPi = "Pi";
    ptrI = "I";
  }
  if (teach)
  {
    startLine();
    // The equation to solve is:
    formatString(&ptrOutput, "<p>$1s</p><p>", LITERAL_GET_ROOTS_POLY1);
    outputPolynomialFactor(&ptrOutput, groupLength, pstFactorInfo);
    showText(" = 0</p>");
  }
  switch (pstFactorInfo->degree)
  {
  case 1:
    LinearEquation(pstFactorInfo->ptrPolyLifted, multiplicity);
    break;
  case 2:
    QuadraticEquation(pstFactorInfo->ptrPolyLifted, multiplicity);
    break;
  case 3:
    CubicEquation(pstFactorInfo->ptrPolyLifted, multiplicity);
    break;
  case 4:
    (void)memcpy(polyIntegerBak, common.poly.polyInteger, sizeof(common.poly.polyInteger));
    (void)memcpy(factorInfoIntegerBak, factorInfoInteger, sizeof(factorInfoInteger));
    QuarticEquation(pstFactorInfo->ptrPolyLifted, multiplicity);
    (void)memcpy(common.poly.polyInteger, polyIntegerBak, sizeof(common.poly.polyInteger));
    (void)memcpy(factorInfoInteger, factorInfoIntegerBak, sizeof(factorInfoInteger));
    break;
  case 5:
    if (isLinearExponential(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is ax^n+b = 0, show roots.
      break;
    }
    if (isSymmetricOrAlternating(nbrFactor, pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is S_n or A_n, indicate that the roots are not solvable.
      break;
    }
    (void)memcpy(polyIntegerBak, common.poly.polyInteger, sizeof(common.poly.polyInteger));
    (void)memcpy(factorInfoIntegerBak, factorInfoInteger, sizeof(factorInfoInteger));
    QuinticEquation(pstFactorInfo->ptrPolyLifted, pstFactorInfo->multiplicity);
    (void)memcpy(common.poly.polyInteger, polyIntegerBak, sizeof(common.poly.polyInteger));
    (void)memcpy(factorInfoInteger, factorInfoIntegerBak, sizeof(factorInfoInteger));
    break;
  default:
    if (isPalindromic(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree) &&
        TestCyclotomic(pstFactorInfo->ptrPolyLifted, multiplicity,
        pstFactorInfo->degree))
    {
      break;
    }
    
    if (isLinearExponential(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is ax^n+b = 0, show roots.
      break;
    }
    if (isQuadraticExponential(pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is ax^(2n)+bx^n+c = 0, show roots.
      break;
    }
    if (isSymmetricOrAlternating(nbrFactor, pstFactorInfo->ptrPolyLifted, pstFactorInfo->degree,
      multiplicity))
    {          // If polynomial is S_n or A_n, indicate that the roots are not solvable.
      break;
    }
    showX(multiplicity * pstFactorInfo->degree);
    *(ptrOutput - 2) = ':';
    // I cannot determine whether the roots of the polynomial 
    // $1?1+?number $1d ??can be solved using radical expressions or not.
    formatString(&ptrOutput, LITERAL_GET_ROOTS_POLY2, nbrFactor + 1);
    break;
  }
  nbrFactorsFound = nbrFactorsFoundBak;
  *pptrOutput = ptrOutput;
}
