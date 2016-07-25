// <XMP>
// Compute gaussian integer expressions
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated June 1st, 2002, See http://www.alpertron.com.ar/GAUSSIAN.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.math.*;
public final class GaussExpression {

/* Info array indices */
private static final int STACK_INDEX = 0;
private static final int EXPR_INDEX = 1;
private static final BigInteger BigInt0 = BigInteger.valueOf(0L);
private static final BigInteger BigInt1 = BigInteger.valueOf(1L);

  // Errors for next routine:
  // >0: Ok.
  // -2: Number too high (more than 1000 digits).
  // -3: Intermediate expression too high (more than 2000 digits).
  // -4: Non-integer division.
  // -5: Parenthesis mismatch.
  // -6: Syntax error
  // -7: Too many parentheses.
  // -8: Invalid parameter.
  // -100: Break.
  // Operators accepted: +, -, *, /, ^, !, F(, L(, P(.

public static int ComputeExpression(String expr, BigInteger ExpressionResult[]) {
  int retcode;
  int Info[] = new int[2];
  Info[STACK_INDEX] = 0;
  Info[EXPR_INDEX] = 0;
  int stackOperators[] = new int[400];
  BigInteger stackRealValues[] = new BigInteger[400];
  BigInteger stackImagValues[] = new BigInteger[400];
  retcode = ComputeExpr(Info, expr, stackOperators, stackRealValues,
                        stackImagValues, ExpressionResult);
  if (retcode != 0) {return retcode;}
  if (ExpressionResult[0].bitLength() > 3322 &&
      ExpressionResult[1].bitLength() > 3322) {return -2;}
  return 0;
  }

private static int ComputeExpr(int [] Info, String expr,
    int [] stackOperators, BigInteger [] stackRealValues,
    BigInteger [] stackImagValues, BigInteger [] ExpressionResult) {

  int i,j;
  char charValue;
  int exprIndexAux;
  int SubExprResult,len;
  BigInteger factorial, Tmp;
  int retcode;
  int leftNumberFlag = 0;
  int startStackIndex = Info[STACK_INDEX];

  while (Info[EXPR_INDEX] < expr.length()) {
    charValue = expr.charAt(Info[EXPR_INDEX]);
    if (charValue == '!') {           // Calculating factorial.
      if (leftNumberFlag == 0) {return -6;}
      if (stackImagValues[Info[STACK_INDEX]].signum() != 0) {return -8;}
      len = stackRealValues[Info[STACK_INDEX]].bitLength()-1;
      if (len > 16) {return -3;}
      len = stackRealValues[Info[STACK_INDEX]].intValue();
      if (len < 0 || len > 807) {return -3;}
      factorial = BigInt1;
      for (i=2; i<=len; i++) {
        factorial = factorial.multiply(BigInteger.valueOf(i));
        }
      stackRealValues[Info[STACK_INDEX]] = factorial;
      }
    else if (charValue == '#') {           // Calculating primorial.
      if (leftNumberFlag == 0) {return -6;}
      if (stackImagValues[Info[STACK_INDEX]].signum() != 0) {return -8;}
      len = stackRealValues[Info[STACK_INDEX]].bitLength()-1;
      if (len > 16) {return -3;}
      len = stackRealValues[Info[STACK_INDEX]].intValue();
      if (len < 0 || len > 4691) {return -3;}
      factorial = BigInt1;
                      // Check if number is prime
      for (i=2; i*i<=len; i++) {
        if (len/i*i==len) {return -8;}
        }
      factorial = BigInt1;
compute_primorial_loop:
      for (i=2; i<=len; i++) {
        for (j=2; j*j<=i; j++) {
          if (i/j*j==i) {continue compute_primorial_loop;}
          }
        factorial = factorial.multiply(BigInteger.valueOf(i));
        }
      stackRealValues[Info[STACK_INDEX]] = factorial;
      }
    else if ((retcode = func(Info, expr, ExpressionResult, stackOperators,
                             stackRealValues, stackImagValues,
                             "GCD", 2, leftNumberFlag)) <= 0) {
      if (retcode != 0) {return retcode;}
      retcode = ComputeGCD(Info, stackRealValues, stackImagValues);
      if (retcode != 0) {return retcode;}
      leftNumberFlag = 1;
      }
    else if ((retcode = func(Info, expr, ExpressionResult, stackOperators,
                             stackRealValues, stackImagValues,
                             "RE", 1, leftNumberFlag)) <= 0) {
      if (retcode != 0) {return retcode;}
      retcode = ComputeRe(Info, stackRealValues, stackImagValues);
      if (retcode != 0) {return retcode;}
      leftNumberFlag = 1;
      }
    else if ((retcode = func(Info, expr, ExpressionResult, stackOperators,
                             stackRealValues, stackImagValues,
                             "NORM", 1, leftNumberFlag)) <= 0) {
      if (retcode != 0) {return retcode;}
      retcode = ComputeNorm(Info, stackRealValues, stackImagValues);
      if (retcode != 0) {return retcode;}
      leftNumberFlag = 1;
      }
    else if ((retcode = func(Info, expr, ExpressionResult, stackOperators,
                             stackRealValues, stackImagValues,
                             "IM", 1, leftNumberFlag)) <= 0) {
      if (retcode != 0) {return retcode;}
      retcode = ComputeIm(Info, stackRealValues, stackImagValues);
      if (retcode != 0) {return retcode;}
      leftNumberFlag = 1;
      }
    else if ((retcode = func(Info, expr, ExpressionResult, stackOperators,
                             stackRealValues, stackImagValues,
                             "MODEXP", 3, leftNumberFlag)) <= 0) {
      if (retcode != 0) {return retcode;}
      retcode = ComputeModExp(Info, stackRealValues, stackImagValues);
      if (retcode != 0) {return retcode;}
      leftNumberFlag = 1;
      }
    else if ((retcode = func(Info, expr, ExpressionResult, stackOperators,
                             stackRealValues, stackImagValues,
                             "MODINV", 2, leftNumberFlag)) <= 0) {
      if (retcode != 0) {return retcode;}
      retcode = ComputeModInv(Info, stackRealValues, stackImagValues);
      if (retcode != 0) {return retcode;}
      leftNumberFlag = 1;
      }
    else if ((retcode = func(Info, expr, ExpressionResult, stackOperators,
                             stackRealValues, stackImagValues,
                             "F", 1, leftNumberFlag)) <= 0) {
      if (retcode != 0) {return retcode;}
      retcode = ComputeFibonacci(Info, stackRealValues, stackImagValues);
      if (retcode != 0) {return retcode;}
      leftNumberFlag = 1;
      }
    else if ((retcode = func(Info, expr, ExpressionResult, stackOperators,
                             stackRealValues, stackImagValues,
                             "L", 1, leftNumberFlag)) <= 0) {
      if (retcode != 0) {return retcode;}
      retcode = ComputeLucas(Info, stackRealValues, stackImagValues);
      if (retcode != 0) {return retcode;}
      leftNumberFlag = 1;
      }
    else if ((retcode = func(Info, expr, ExpressionResult, stackOperators,
                             stackRealValues, stackImagValues,
                             "P", 1, leftNumberFlag)) <= 0) {
      if (retcode != 0) {return retcode;}
      retcode = ComputePartition(Info, stackRealValues, stackImagValues);
      if (retcode != 0) {return retcode;}
      leftNumberFlag = 1;
      }
    else if (charValue == '+' || charValue == '-') {
      if (leftNumberFlag == 0) {      // Unary plus/minus operator
        Info[EXPR_INDEX]++;
        if (charValue == '+') {
          continue;
          }
        else {
          if (Info[STACK_INDEX] > startStackIndex && stackOperators[Info[STACK_INDEX]-1] == '_') {
            Info[STACK_INDEX]--;
            continue;
            }
          if (Info[STACK_INDEX] >= stackOperators.length) {return -7;}
          stackOperators[Info[STACK_INDEX]++] = '_'; /* Unitary minus */
          continue;
          }
        }
      if (Info[STACK_INDEX] > startStackIndex &&
          stackOperators[Info[STACK_INDEX]-1] != '(') {
        if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
          return SubExprResult;
          }
        if (Info[STACK_INDEX] > startStackIndex &&
            stackOperators[Info[STACK_INDEX]-1] != '(') {
          if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
            return SubExprResult;
            }
          if (Info[STACK_INDEX] > startStackIndex &&
              stackOperators[Info[STACK_INDEX]-1] != '(') {
            if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
              return SubExprResult;
              }
            }                         /* end if */
          }                           /* end if */
        }                             /* end if */
      stackOperators[Info[STACK_INDEX]++] = charValue;
      leftNumberFlag = 0;
      }                               /* end if */
    else if (charValue == '*' || charValue == '/' || charValue == '%') {
      if (leftNumberFlag == 0) {return -6;}
      if (Info[STACK_INDEX] > startStackIndex && (stackOperators[Info[STACK_INDEX]-1] == '^' ||
          stackOperators[Info[STACK_INDEX]-1] == '*' ||
          stackOperators[Info[STACK_INDEX]-1] == '/')) {
        if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
          return SubExprResult;
          }
        if (Info[STACK_INDEX] > startStackIndex &&
             (stackOperators[Info[STACK_INDEX]-1] == '^' ||
              stackOperators[Info[STACK_INDEX]-1] == '*' ||
              stackOperators[Info[STACK_INDEX]-1] == '/' ||
              stackOperators[Info[STACK_INDEX]-1] == '%')) {
          if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
            return SubExprResult;
            }
          }                         /* end if */
        }                           /* end if */
      stackOperators[Info[STACK_INDEX]++] = charValue;
      leftNumberFlag = 0;
      }                             
    else if (charValue == '^') {
      if (leftNumberFlag == 0) {return -6;}
      if (Info[STACK_INDEX] > 0 &&
          (stackOperators[Info[STACK_INDEX]-1] == '^')) {
        if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
          return SubExprResult;
          }
        }                         /* end if */
      stackOperators[Info[STACK_INDEX]++] = charValue;
      leftNumberFlag = 0;
      }                           /* end if */
    else if (charValue == '(') {
      if (leftNumberFlag == 1) {return -6;}
      if (Info[STACK_INDEX] >= stackOperators.length) {return -7;}
      stackOperators[Info[STACK_INDEX]++] = charValue;
      }                           
    else if (charValue == ')' || charValue == ',') {
      if (leftNumberFlag == 0) {return -6;}
      if (Info[STACK_INDEX] > startStackIndex &&
          stackOperators[Info[STACK_INDEX]-1] != '(') {
        if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
          return SubExprResult;
          }
        if (Info[STACK_INDEX] > startStackIndex &&
            stackOperators[Info[STACK_INDEX]-1] != '(') {
          if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
            return SubExprResult;
            }
          if (Info[STACK_INDEX] > startStackIndex &&
              stackOperators[Info[STACK_INDEX]-1] != '(') {
            if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
              return SubExprResult;
              }
            }
          }
        }
      if (Info[STACK_INDEX] == startStackIndex) {break;}
      if (charValue == ',') {return -5;}
      Info[STACK_INDEX]--;    /* Discard ')' */
      stackRealValues[Info[STACK_INDEX]] = stackRealValues[Info[STACK_INDEX]+1];
      stackImagValues[Info[STACK_INDEX]] = stackImagValues[Info[STACK_INDEX]+1];
      leftNumberFlag = 1;
      }
    else if (charValue == 'i' || charValue == 'I') {
      if (leftNumberFlag == 0) {
        stackRealValues[Info[STACK_INDEX]] = BigInt0;
        stackImagValues[Info[STACK_INDEX]] = BigInt1;
        leftNumberFlag = 1;
        }
      else {
        Tmp = stackImagValues[Info[STACK_INDEX]];
        stackImagValues[Info[STACK_INDEX]] = stackRealValues[Info[STACK_INDEX]];
        stackRealValues[Info[STACK_INDEX]] = Tmp.negate();
        }
      }
    else if (charValue >= '0' && charValue <= '9') {
      exprIndexAux = Info[EXPR_INDEX];
      while (exprIndexAux < expr.length()-1) {
        charValue = expr.charAt(exprIndexAux+1);
        if (charValue >= '0' && charValue <= '9') {
          exprIndexAux++;
          }
        else {
          break;
          }
        }
      stackRealValues[Info[STACK_INDEX]] = new BigInteger(expr.substring(Info[EXPR_INDEX],exprIndexAux+1));
      stackImagValues[Info[STACK_INDEX]] = BigInt0;
      leftNumberFlag = 1;
      Info[EXPR_INDEX] = exprIndexAux;
      }                            /* end if */
    Info[EXPR_INDEX]++;
    }                              /* end while */
  if (leftNumberFlag == 0) {return -6;}
  if (Info[STACK_INDEX] > startStackIndex && stackOperators[Info[STACK_INDEX]-1] != '(') {
    if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
      return SubExprResult;
      }
    if (Info[STACK_INDEX] > startStackIndex && stackOperators[Info[STACK_INDEX]-1] != '(') {
      if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
        return SubExprResult;
        }
      if (Info[STACK_INDEX] > startStackIndex && stackOperators[Info[STACK_INDEX]-1] != '(') {
        if ((SubExprResult = ComputeSubExpr(Info, stackRealValues, stackImagValues, stackOperators)) != 0) {
          return SubExprResult;
          }
        }
      }
    }
  if (Info[STACK_INDEX] != startStackIndex) {return -5;}
  ExpressionResult[0] = stackRealValues[startStackIndex];
  ExpressionResult[1] = stackImagValues[startStackIndex];
  return 0;
  }

private static int func(int [] Info, String expr, 
                        BigInteger [] ExpressionResult,
                        int [] stackOperators,
                        BigInteger [] stackRealValues,
                        BigInteger [] stackImagValues, String funcName,
                        int funcArgs, int leftNumberFlag) {
  int index, retcode;
  char compareChar;

  if (Info[EXPR_INDEX] + funcName.length() > expr.length()) {return 1;}
  if (expr.substring(Info[EXPR_INDEX],Info[EXPR_INDEX]+funcName.length()).
               toUpperCase().equals(funcName) == false) {return 1;}
  Info[EXPR_INDEX] += funcName.length();
  if (leftNumberFlag == 1) {return -6;}
  SkipSpaces(Info, expr);
  if (Info[EXPR_INDEX] == expr.length() ||
      expr.charAt(Info[EXPR_INDEX]++) != '(') {
    return -6;
    }
  for (index = 0; index < funcArgs; index++) {
    SkipSpaces(Info, expr);
    if (Info[STACK_INDEX] >= stackOperators.length) {return -7;}
    retcode = ComputeExpr(Info, expr, stackOperators, stackRealValues,
                          stackImagValues, ExpressionResult);
    if (retcode != 0) {return retcode;}
    SkipSpaces(Info, expr);
    compareChar = (index == funcArgs-1? ')': ',');
    if (Info[EXPR_INDEX] == expr.length() ||
        expr.charAt(Info[EXPR_INDEX]++) != compareChar) {
      return -6;
      }
    Info[STACK_INDEX]++;
    }
  Info[STACK_INDEX] -= funcArgs;
  Info[EXPR_INDEX]--;
  return 0;
  }

private static void SkipSpaces(int Info[], String expr) {
  while (Info[EXPR_INDEX] < expr.length()) {
    if (expr.charAt(Info[EXPR_INDEX]) > ' ') {
      break;
      }
    Info[EXPR_INDEX]++;
    }
  return;
  }

private static int ComputeSubExpr(int Info[], BigInteger [] stackRealValues, BigInteger [] stackImagValues, int [] stackOperators) {
  int i, len, val;
  double logarithm;
  int stackOper;
  BigInteger Re1, Re2, Im1, Im2, Re, Im, ReTmp, ImTmp;
  Info[STACK_INDEX]--;

  stackOper = stackOperators[Info[STACK_INDEX]];
  Re1 = stackRealValues[Info[STACK_INDEX]];
  Re2 = stackRealValues[Info[STACK_INDEX]+1];
  Im1 = stackImagValues[Info[STACK_INDEX]];
  Im2 = stackImagValues[Info[STACK_INDEX]+1];

  switch (stackOper) {
    case '+':
      stackRealValues[Info[STACK_INDEX]] = Re1.add(Re2);
      stackImagValues[Info[STACK_INDEX]] = Im1.add(Im2);
      return 0;
    case '-':
      stackRealValues[Info[STACK_INDEX]] = Re1.subtract(Re2);
      stackImagValues[Info[STACK_INDEX]] = Im1.subtract(Im2);
      return 0;
    case '_':
      stackRealValues[Info[STACK_INDEX]] = Re2.negate();
      stackImagValues[Info[STACK_INDEX]] = Im2.negate();
      return 0;
    case '/':
      BigInteger norm = Re2.multiply(Re2).add(Im2.multiply(Im2));
      if (norm.signum() == 0) {return -3;}
      Re = Re1.multiply(Re2).add(Im1.multiply(Im2));
      Im = Im1.multiply(Re2).subtract(Re1.multiply(Im2));
      if (Re.remainder(norm).signum() != 0 ||
          Im.remainder(norm).signum() != 0) {return -4;}
      stackRealValues[Info[STACK_INDEX]] = Re.divide(norm);
      stackImagValues[Info[STACK_INDEX]] = Im.divide(norm);
      return 0;
    case '%':
      BigInteger [] Result = new BigInteger[2];
      int retcode = Modulo(Re1, Im1, Re2, Im2, Result);
      if (retcode != 0) {return retcode;}
      stackRealValues[Info[STACK_INDEX]] = Result[0];
      stackImagValues[Info[STACK_INDEX]] = Result[1];
      return 0;
    case '*':
      Re = stackRealValues[Info[STACK_INDEX]] = Re1.multiply(Re2).subtract(Im1.multiply(Im2));
      Im = stackImagValues[Info[STACK_INDEX]] = Re1.multiply(Im2).add(Im1.multiply(Re2));
      if (Re.bitLength() > 6644 || Im.bitLength() > 6644) {return -3;}
      return 0;
    case '^':
      if (Im2.signum() != 0) {return -8;}
      norm = Re1.multiply(Re1).add(Im1.multiply(Im1));
      len = norm.bitLength()-1;
      if (len > 32) {
        logarithm = (double)(len-32) +
                 Math.log(norm.shiftRight(len-32).doubleValue())/Math.log(2);
        }
      else {
        logarithm = Math.log(norm.doubleValue())/Math.log(2);
        }
      if (logarithm * Re2.doubleValue() > 6644*2) {return -3;}
          /* Compute actual power */
      Re = BigInt1;
      Im = BigInt0;
      for (i=Re2.bitLength()-1; i>=0; i--) {
        ReTmp = Re.multiply(Re).subtract(Im.multiply(Im));
        ImTmp = Re.multiply(Im).shiftLeft(1);
        Re = ReTmp; Im = ImTmp;
        if (Re2.testBit(i)) {
          ReTmp = Re1.multiply(Re).subtract(Im1.multiply(Im));
          ImTmp = Re1.multiply(Im).add(Im1.multiply(Re));
          Re = ReTmp; Im = ImTmp;
          }
        }
      stackRealValues[Info[STACK_INDEX]] = Re;
      stackImagValues[Info[STACK_INDEX]] = Im;
      return 0;
    }                /* end switch */
  return 0;
  }

private static int ComputeFibonacci(int [] Info, BigInteger[] stackRealValues,
                             BigInteger[] stackImagValues) {

  BigInteger FibonPrev, FibonAct, FibonNext;
  BigInteger Re = stackRealValues[Info[STACK_INDEX]];
  BigInteger Im = stackImagValues[Info[STACK_INDEX]];
  if (Im.signum() != 0) {return -8;}
  int len = Re.bitLength()-1;
  if (len > 16) {return -3;}
  len = Re.intValue();
  if (len > 9571) {return -3;}
  if (len < 0) {return -8;}
  FibonPrev = BigInteger.valueOf(1);
  FibonAct = BigInteger.valueOf(0);
  for (int i=1; i<=len; i++) {
    FibonNext = FibonPrev.add(FibonAct);
    FibonPrev = FibonAct;
    FibonAct = FibonNext;
    }
  stackRealValues[Info[STACK_INDEX]] = FibonAct;
  stackImagValues[Info[STACK_INDEX]] = BigInteger.valueOf(0L);
  return 0;
  }

private static int ComputeLucas(int [] Info, BigInteger[] stackRealValues,
                         BigInteger[] stackImagValues) {

  BigInteger FibonPrev, FibonAct, FibonNext;
  BigInteger Re = stackRealValues[Info[STACK_INDEX]];
  BigInteger Im = stackImagValues[Info[STACK_INDEX]];

  if (Im.signum() != 0) {return -8;}
  int len = Re.bitLength()-1;
  if (len > 16) {return -3;}
  len = Re.intValue();
  if (len > 9572) {return -3;}
  if (len < 0) {return -8;}
  FibonPrev = BigInteger.valueOf(-1);
  FibonAct = BigInteger.valueOf(2);
  for (int i=1; i<=len; i++) {
    FibonNext = FibonPrev.add(FibonAct);
    FibonPrev = FibonAct;
    FibonAct = FibonNext;
    }
  stackRealValues[Info[STACK_INDEX]] = FibonAct;
  stackImagValues[Info[STACK_INDEX]] = BigInteger.valueOf(0L);
  return 0;
  }

private static int ComputePartition(int [] Info, BigInteger[] stackRealValues,
                             BigInteger[] stackImagValues) {
  int val, i, j, k, u, Tmp, indexH, indexL, indexNew, count;
  long Cy, Result;
  double Tmp1, Tmp2, Tmp3;
  long Part[];
  int Index[];
  byte Conv[];
  long DosALa63 = 1L << 63;
  BigInteger Re = stackRealValues[Info[STACK_INDEX]];
  BigInteger Im = stackImagValues[Info[STACK_INDEX]];

  if (Im.signum() != 0) {return -8;}
  int len = Re.bitLength()-1;
  if (len > 24) {return -3;}
  len = Re.intValue();
  if (len > 3520000) {return -3;}
  if (len < 0) {return -8;}
  len = 2;
  Tmp1 = 0.0578227587396094872;     // pi * sqrt(2/3) / log(2^64)
  Tmp2 = 0.9563674804631159673;     // 1 - log(4*sqrt(3)) / log(2^64)
  Tmp3 = 0.0225421100138900531;     // 1 / log(2^64)
  val = Re.intValue();
  for (i=1; i<=val; i++) {
    len += (long)(Math.floor(Tmp1 * Math.sqrt((double)i) + Tmp2 -
                Math.log((double)i) * Tmp3));
    }
  Part = new long[len];
  Index = new int[val+2];      /* Index to Part vector for start of number */
  Part[0] = DosALa63+1;        /* Initialize p(0) */
  Index[0] = 0;                /* Initialize index to p(0) */
  Index[1] = Tmp = len = 1;    /* Initialize index to p(1) */
  for (i=1; i<=val; i++) {
    len = Index[i] - Index[i-1] + 1;     /* Initialize number length */
    Tmp = Index[i] + len;
    for (k=Index[i]; k<Tmp; k++) {       /* Initialize number to zero */
      Part[k] = DosALa63;
      }
    for (k=1; (3*k-1)*k<=2*i; k++) {
      indexL = Index[i-(3*k-1)*k/2];
      indexH = Index[i-(3*k-1)*k/2+1];
      for (u=((3*k+1)*k<=2*i?0:1); u<=1; u++) {
        Cy = DosALa63;
        indexNew = Index[i];
        if (k%2 == 0) {                    /* Subtract */
          for (j=indexL; j<indexH; j++) {
            Result = Cy + Part[indexNew] - Part[j];
            Cy = (Result > Part[indexNew] || (Result == Part[indexNew] && Cy != DosALa63)? DosALa63-1:DosALa63);
            Part[indexNew] = Result;
            indexNew++;
            }
          while (indexNew < Tmp) {
            Result = Cy + Part[indexNew] + DosALa63;
            Cy = (Result > Part[indexNew] || (Result == Part[indexNew] && Cy != DosALa63)? DosALa63-1:DosALa63);
            Part[indexNew] = Result;
            indexNew++;
            }
          }
        else {                              /* Add */
          for (j=indexL; j<indexH; j++) {
            Result = Cy + Part[indexNew] + Part[j];
            Cy = (Result < Part[indexNew] || Result < Part[j]?
                  DosALa63+1: DosALa63);
            Part[indexNew] = Result;
            indexNew++;
            }   
          while (indexNew < Tmp) {
            Result = Cy + Part[indexNew] + DosALa63;
            Cy = (Result < Part[indexNew]?DosALa63+1: DosALa63);
            Part[indexNew] = Result;
            indexNew++;
            }
          }
        if (u==0) {
          indexL = Index[i-(3*k+1)*k/2];
          indexH = Index[i-(3*k+1)*k/2+1];
          }
        }              // end for u
      }
    if (Part[Tmp - 1] == DosALa63) {   /* Discard heading zeros */
      len--;
      Tmp--;
      }
    else {
      }
    Index[i+1] = Tmp;
    }                      // end for i
  count = len*8;
  Conv = new byte[count+1];
  for (i=Index[val]; i<Tmp; i++) {
    Conv[count] = (byte)(Part[i] & 0xFF);
    Conv[count-1] = (byte)(Part[i] >>> 8 & 0xFF);
    Conv[count-2] = (byte)(Part[i] >>> 16 & 0xFF);
    Conv[count-3] = (byte)(Part[i] >>> 24 & 0xFF);
    Conv[count-4] = (byte)(Part[i] >>> 32 & 0xFF);
    Conv[count-5] = (byte)(Part[i] >>> 40 & 0xFF);
    Conv[count-6] = (byte)(Part[i] >>> 48 & 0xFF);
    Conv[count-7] = (byte)((Part[i] >>> 56 & 0xFF) ^ 0x80);
    count -= 8;
    }
  Conv[0] = 0;
  stackRealValues[Info[STACK_INDEX]] = new BigInteger(Conv);
  stackImagValues[Info[STACK_INDEX]] = BigInteger.valueOf(0L);
  return 0;
  }

private static int ComputeGCD(int [] Info, BigInteger[] stackRealValues,
                              BigInteger[] stackImagValues) {

  BigInteger Re, Im;
  BigInteger Re1 = stackRealValues[Info[STACK_INDEX]];
  BigInteger Re2 = stackRealValues[Info[STACK_INDEX]+1];
  BigInteger Im1 = stackImagValues[Info[STACK_INDEX]];
  BigInteger Im2 = stackImagValues[Info[STACK_INDEX]+1];
  while (Re2.signum() != 0 || Im2.signum() != 0) {
    BigInteger norm = Re2.multiply(Re2).add(Im2.multiply(Im2));
    Re = Re1.multiply(Re2).add(Im1.multiply(Im2)).
         shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
    Im = Im1.multiply(Re2).subtract(Re1.multiply(Im2)).
         shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
    Re1 = Re1.subtract(Re.multiply(Re2)).add(Im.multiply(Im2));
    Im1 = Im1.subtract(Im.multiply(Re2)).subtract(Re.multiply(Im2));
    Re = Re1; Re1 = Re2; Re2 = Re;
    Im = Im1; Im1 = Im2; Im2 = Im;
    }
  while (Re1.compareTo(Im1.abs()) < 0) {
    Re = Re1; Re1 = Im1.negate(); Im1 = Re;  /* multiply by i */
    }
  stackRealValues[Info[STACK_INDEX]] = Re1;
  stackImagValues[Info[STACK_INDEX]] = Im1;
  return 0;
  }

private static int ComputeModExp(int [] Info, BigInteger[] stackRealValues,
                                 BigInteger[] stackImagValues) {
  BigInteger [] Result = new BigInteger[2];
  int retcode, i;
  double logarithm;
  BigInteger norm, ReTmp, ImTmp;

  BigInteger ReBase = stackRealValues[Info[STACK_INDEX]];
  BigInteger ImBase = stackImagValues[Info[STACK_INDEX]];
  BigInteger ReExp = stackRealValues[Info[STACK_INDEX]+1];
  BigInteger ImExp = stackImagValues[Info[STACK_INDEX]+1];
  BigInteger ReMod = stackRealValues[Info[STACK_INDEX]+2];
  BigInteger ImMod = stackImagValues[Info[STACK_INDEX]+2];
  if (ImExp.signum() != 0) {return -8;}
  if (ReExp.signum() < 0) {
    ReExp = ReExp.negate();
    retcode = ModInv(ReBase, ImBase, ReMod, ImMod, Result);
    if (retcode != 0) {return retcode;}
    ReBase = Result[0];
    ImBase = Result[1];
    }
  BigInteger Re = BigInt1;
  BigInteger Im = BigInteger.valueOf(0L);
  if (ReMod.signum() == 0 && ImMod.signum() == 0) {   /* Modulus is zero */
    norm = ReBase.multiply(ReBase).add(ImBase.multiply(ImBase));
    int len = norm.bitLength()-1;
    if (len > 32) {
      logarithm = (double)(len-32) +
               Math.log(norm.shiftRight(len-32).doubleValue())/Math.log(2);
      }
    else {
      logarithm = Math.log(norm.doubleValue())/Math.log(2);
      }
    if (logarithm * ReExp.doubleValue() > 6644*2) {return -3;}
        /* Compute actual power */
    for (i=ReExp.bitLength()-1; i>=0; i--) {
      ReTmp = Re.multiply(Re).subtract(Im.multiply(Im));
      ImTmp = Re.multiply(Im).shiftLeft(1);
      Re = ReTmp; Im = ImTmp;
      if (ReExp.testBit(i)) {
        ReTmp = ReBase.multiply(Re).subtract(ImBase.multiply(Im));
        ImTmp = ReBase.multiply(Im).add(ImBase.multiply(Re));
        Re = ReTmp; Im = ImTmp;
        }
      }
    stackRealValues[Info[STACK_INDEX]] = Re;
    stackImagValues[Info[STACK_INDEX]] = Im;
    }
  else {                            /* Modulus is not zero */
    norm = ReMod.multiply(ReMod).add(ImMod.multiply(ImMod));
    for (i=ReExp.bitLength()-1; i>=0; i--) {
      ReTmp = Re.multiply(Re).subtract(Im.multiply(Im));
      ImTmp = Re.multiply(Im).shiftLeft(1);
      Re = ReTmp.multiply(ReMod).add(ImTmp.multiply(ImMod)).
           shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
      Im = ImTmp.multiply(ReMod).subtract(ReTmp.multiply(ImMod)).
           shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
      ReTmp = ReTmp.subtract(Re.multiply(ReMod)).add(Im.multiply(ImMod));
      ImTmp = ImTmp.subtract(Im.multiply(ReMod)).subtract(Re.multiply(ImMod));
      Re = ReTmp; Im = ImTmp;
      if (ReExp.testBit(i)) {
        ReTmp = ReBase.multiply(Re).subtract(ImBase.multiply(Im));
        ImTmp = ReBase.multiply(Im).add(ImBase.multiply(Re));
        Re = ReTmp.multiply(ReMod).add(ImTmp.multiply(ImMod)).
             shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
        Im = ImTmp.multiply(ReMod).subtract(ReTmp.multiply(ImMod)).
            shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
        ReTmp = ReTmp.subtract(Re.multiply(ReMod)).add(Im.multiply(ImMod));
        ImTmp = ImTmp.subtract(Im.multiply(ReMod)).subtract(Re.multiply(ImMod));
        Re = ReTmp; Im = ImTmp;
        }
      }
    Modulo(Re, Im, ReMod, ImMod, Result);
    stackRealValues[Info[STACK_INDEX]] = Result[0];
    stackImagValues[Info[STACK_INDEX]] = Result[1];
    return 0;
    }
  return 0;
  }

private static int ComputeModInv(int [] Info, BigInteger[] stackRealValues,
                                 BigInteger[] stackImagValues) {
  BigInteger [] Result = new BigInteger[2];
  int retcode;

  retcode = ModInv(stackRealValues[Info[STACK_INDEX]],
                   stackImagValues[Info[STACK_INDEX]],
                   stackRealValues[Info[STACK_INDEX]+1],
                   stackImagValues[Info[STACK_INDEX]+1],
                   Result);
  if (retcode != 0) {return retcode;}
  stackRealValues[Info[STACK_INDEX]] = Result[0];
  stackImagValues[Info[STACK_INDEX]] = Result[1];
  return 0;
  }

private static int ModInv(BigInteger RealNbr, BigInteger ImagNbr,
                          BigInteger RealMod, BigInteger ImagMod,
                          BigInteger [] Result) {

  BigInteger ReG0, ReG1, ImG0, ImG1;
  BigInteger ReU0, ReU1, ImU0, ImU1;
  BigInteger ReV0, ReV1, ImV0, ImV1;
  BigInteger Re, Im, Tmp;

  if (RealMod.signum() == 0 && ImagMod.signum() == 0) {return -8;}
  ReG0 = RealNbr; ImG0 = ImagNbr;
  ReG1 = RealMod; ImG1 = ImagMod;
  ReU0 = BigInt1;
  ReU1 = ImU0 = ImU1 = BigInteger.valueOf(0L);
  while (ReG1.signum() != 0 || ImG1.signum() != 0) {
    BigInteger norm = ReG1.multiply(ReG1).add(ImG1.multiply(ImG1));
    Re = ReG0.multiply(ReG1).add(ImG0.multiply(ImG1)).
         shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
    Im = ImG0.multiply(ReG1).subtract(ReG0.multiply(ImG1)).
         shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
    ReG0 = ReG0.subtract(Re.multiply(ReG1)).add(Im.multiply(ImG1));
    ImG0 = ImG0.subtract(Im.multiply(ReG1)).subtract(Re.multiply(ImG1));
    Tmp = ReG0; ReG0 = ReG1; ReG1 = Tmp;
    Tmp = ImG0; ImG0 = ImG1; ImG1 = Tmp;
    ReU0 = ReU0.subtract(Re.multiply(ReU1)).add(Im.multiply(ImU1));
    ImU0 = ImU0.subtract(Im.multiply(ReU1)).subtract(Re.multiply(ImU1));
    Tmp = ReU0; ReU0 = ReU1; ReU1 = Tmp;
    Tmp = ImU0; ImU0 = ImU1; ImU1 = Tmp;
    }
  while (ReG0.compareTo(ImG0.abs()) < 0) {
    Tmp = ReG0; ReG0 = ImG0.negate(); ImG0 = Tmp;  /* multiply by i */
    Tmp = ReU0; ReU0 = ImU0.negate(); ImU0 = Tmp;  /* multiply by i */
    }
  if (ReG0.equals(BigInt1) == false || ImG0.signum() != 0) {return -8;}
  return Modulo(ReU0, ImU0, RealMod, ImagMod, Result);
  }

private static int ComputeRe(int [] Info, BigInteger[] stackRealValues,
                             BigInteger[] stackImagValues) {

  stackImagValues[Info[STACK_INDEX]] = BigInteger.valueOf(0L);
  return 0;
  }

private static int Modulo(BigInteger ReNum, BigInteger ImNum,
                          BigInteger ReDen, BigInteger ImDen,
                          BigInteger [] Result) {
  BigInteger Re, Im;
  if (ReDen.signum() == 0 && ImDen.signum() == 0) {
    Result[0] = ReNum;
    Result[1] = ImNum;
    return 0;
    }
  BigInteger ReMin = BigInteger.valueOf(0L);
  BigInteger ImMin = BigInteger.valueOf(0L);
  BigInteger norm = ReDen.multiply(ReDen).add(ImDen.multiply(ImDen));
  Re = ReNum.multiply(ReDen).add(ImNum.multiply(ImDen)).
       shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
  Im = ImNum.multiply(ReDen).subtract(ReNum.multiply(ImDen)).
       shiftLeft(1).divide(norm).add(BigInt1).shiftRight(1);
  ReNum = ReNum.subtract(Re.multiply(ReDen)).add(Im.multiply(ImDen));
  ImNum = ImNum.subtract(Im.multiply(ReDen)).subtract(Re.multiply(ImDen));
  BigInteger normmin = BigInteger.valueOf(-1L);
  for (int i=0; i<9; i++) {
    switch (i) {
      case 0:
        Re = ReNum; Im = ImNum;
        break;
      case 1:
        Re = ReNum.add(ReDen); Im = ImNum.add(ImDen);
        break;
      case 2:
        Re = ReNum.add(ReDen).subtract(ImDen); Im = ImNum.add(ReDen).add(ImDen);
        break;
      case 3:
        Re = ReNum.subtract(ImDen); Im = ImNum.add(ReDen);
        break;
      case 4:
        Re = ReNum.subtract(ReDen).subtract(ImDen); Im = ImNum.add(ReDen).subtract(ImDen);
        break;
      case 5:
        Re = ReNum.subtract(ReDen); Im = ImNum.subtract(ImDen);
        break;
      case 6:
        Re = ReNum.subtract(ReDen).add(ImDen); Im = ImNum.subtract(ReDen).subtract(ImDen);
        break;
      case 7:
        Re = ReNum.add(ImDen); Im = ImNum.subtract(ReDen);
        break;
      case 8:
        Re = ReNum.add(ReDen).add(ImDen); Im = ImNum.subtract(ReDen).add(ImDen);
        break;
      }
    if (Re.signum() >= 0) {
      norm = Re.multiply(Re).add(Im.multiply(Im));
      if (normmin.signum() < 0 || norm.compareTo(normmin) < 0) {
        normmin = norm;
        ReMin = Re;
        ImMin = Im;
        }
      }
    }
  Result[0] = ReMin;
  Result[1] = ImMin;
  return 0;
  }
private static int ComputeIm(int [] Info, BigInteger[] stackRealValues,
                             BigInteger[] stackImagValues) {
  stackRealValues[Info[STACK_INDEX]] = stackImagValues[Info[STACK_INDEX]];
  stackImagValues[Info[STACK_INDEX]] = BigInteger.valueOf(0L);
  return 0;
  }

private static int ComputeNorm(int [] Info, BigInteger[] stackRealValues,
                               BigInteger[] stackImagValues) {
  BigInteger Re1 = stackRealValues[Info[STACK_INDEX]];
  BigInteger Im1 = stackImagValues[Info[STACK_INDEX]];
  stackRealValues[Info[STACK_INDEX]] = Re1.multiply(Re1).add(Im1.multiply(Im1));
  stackImagValues[Info[STACK_INDEX]] = BigInteger.valueOf(0L);
  return 0;
  }

}
