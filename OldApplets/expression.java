// <XMP>
// Compute Expressions
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated May 8th, 2011, See http://www.alpertron.com.ar/ECM.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.math.*;
public final class expression {
  private static final BigInteger BigInt1 = BigInteger.valueOf(1L);
  private static final BigInteger BigInt2 = BigInteger.valueOf(2L);
  private static final BigInteger BigInt3 = BigInteger.valueOf(3L);

  // Errors for next routine:
  // >0: Ok.
  // -2: Number too high (more than 10000 digits).
  // -3: Intermediate expression too high (more than 2000 digits).
  // -4: Non-integer division.
  // -5: Parenthesis mismatch.
  // -6: Syntax error
  // -7: Too many parentheses.
  // -8: Invalid parameter.
  // -100: Break.
  // Operators accepted: +, -, *, /, ^, !, F(, L(, P(.

  public static int ComputeExpression(String expr, int type, BigInteger ExpressionResult[]) {
      BigInteger BigInt1 = BigInteger.valueOf(1L);
      int stackIndex = 0;
      int exprIndex = 0;
      int exprLength = expr.length();
      int i,j;
      char charValue;
      boolean leftNumberFlag = false;
      int exprIndexAux;
      int SubExprResult,len;
      BigInteger factorial;
      BigInteger stackValues[] = new BigInteger[400];
      int stackOperators[] = new int[400];

      while (exprIndex < exprLength) {
        charValue = expr.charAt(exprIndex);
        if (charValue == '!') {           // Calculating factorial.
          if (leftNumberFlag == false) {return -6;}
          len = stackValues[stackIndex].bitLength()-1;
          if (len > 16) {return -3;}
          len = stackValues[stackIndex].intValue();
          if (len < 0 || len > 5984) {return -3;}
          factorial = BigInt1;
          for (i=2; i<=len; i++) {
            factorial = factorial.multiply(BigInteger.valueOf(i));
          }
          stackValues[stackIndex] = factorial;
        }
        if (charValue == '#') {           // Calculating primorial.
          if (leftNumberFlag == false) {return -6;}
          len = stackValues[stackIndex].bitLength()-1;
          if (len > 16) {return -3;}
          len = stackValues[stackIndex].intValue();
          if (len < 0 || len > 46049) {return -3;}
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
          stackValues[stackIndex] = factorial;
        }
        if (charValue == 'B' || charValue == 'b' ||
            charValue == 'N' || charValue == 'n' ||
            charValue == 'F' || charValue == 'f' ||
            charValue == 'P' || charValue == 'p' ||
            charValue == 'L' || charValue == 'l') {
          if (leftNumberFlag || exprIndex == exprLength-1) {
            return -6;
          }
          exprIndex++;
          if (expr.charAt(exprIndex) != '(') {return -6;}
          if (stackIndex > 395) {return -7;}
          stackOperators[stackIndex++] = charValue & 0xDF; /* Convert to uppercase */
          charValue = '(';
        }
        if (charValue == '+' || charValue == '-') {
          if (leftNumberFlag == false) {      // Unary plus/minus operator
            exprIndex++;
            if (charValue == '+') {
              continue;
              }
            else {
              if (stackIndex > 0 && stackOperators[stackIndex-1] == '_') {
                stackIndex--;
                continue;
                }
              if (stackIndex > 395) {return -7;}
              stackOperators[stackIndex++] = '_'; /* Unitary minus */
              continue;
              }
            }
          if (stackIndex > 0 && stackOperators[stackIndex-1] != '(') {
            if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
              return SubExprResult;
            }
            if (stackIndex > 0 && stackOperators[stackIndex-1] != '(') {
              if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                return SubExprResult;
              }
              if (stackIndex > 0 && stackOperators[stackIndex-1] != '(') {
                if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                  return SubExprResult;
                }
              }                         /* end if */
            }                           /* end if */
          }                             /* end if */
          stackOperators[stackIndex++] = charValue;
          leftNumberFlag = false;
        }                               /* end if */
        else {
          if (charValue == '*' || charValue == '/' || charValue == '%') {
            if (leftNumberFlag == false) {return -6;}
            if (stackIndex > 0 && (stackOperators[stackIndex-1] == '^' ||
                  stackOperators[stackIndex-1] == '*' ||
                  stackOperators[stackIndex-1] == '/' ||
                  stackOperators[stackIndex-1] == '%' ||
                  stackOperators[stackIndex-1] == 'B' ||
                  stackOperators[stackIndex-1] == 'N' ||
                  stackOperators[stackIndex-1] == 'F' ||
                  stackOperators[stackIndex-1] == 'L' ||
                  stackOperators[stackIndex-1] == 'P')) {
              if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                return SubExprResult;
              }
              if (stackIndex > 0 && (stackOperators[stackIndex-1] == '^' ||
                    stackOperators[stackIndex-1] == '*' ||
                    stackOperators[stackIndex-1] == '/' ||
                    stackOperators[stackIndex-1] == '%' ||
                    stackOperators[stackIndex-1] == 'B' ||
                    stackOperators[stackIndex-1] == 'N' ||
                    stackOperators[stackIndex-1] == 'F' ||
                    stackOperators[stackIndex-1] == 'L' ||
                    stackOperators[stackIndex-1] == 'P')) {
                if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                  return SubExprResult;
                }
              }                         /* end if */
            }                           /* end if */
            stackOperators[stackIndex++] = charValue;
            leftNumberFlag = false;
          }                             
          else {
            if (charValue == '^') {
              if (leftNumberFlag == false) {return -6;}
              if (stackIndex > 0 && (stackOperators[stackIndex-1] == '^' ||
                    stackOperators[stackIndex-1] == 'B' ||
                    stackOperators[stackIndex-1] == 'N' ||
                    stackOperators[stackIndex-1] == 'F' ||
                    stackOperators[stackIndex-1] == 'L' ||
                    stackOperators[stackIndex-1] == 'P')) {
                if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                  return SubExprResult;
                }
              }                         /* end if */
              stackOperators[stackIndex++] = charValue;
              leftNumberFlag = false;
            }                           /* end if */
            else {
              if (charValue == '(') {
                if (leftNumberFlag == true) {return -6;}
                if (stackIndex > 395) {return -7;}
                stackOperators[stackIndex++] = charValue;
              }                           
              else {
                if (charValue == ')') {
                  if (leftNumberFlag == false) {return -6;}
                  if (stackIndex > 0 && stackOperators[stackIndex-1] != '(') {
                    if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                      return SubExprResult;
                    }
                    if (stackIndex > 0 && stackOperators[stackIndex-1] != '(') {
                      if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                        return SubExprResult;
                      }
                      if (stackIndex > 0 && stackOperators[stackIndex-1] != '(') {
                        if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                          return SubExprResult;
                        }
                      }
                    }
                  }
                  if (stackIndex == 0) {return -5;}
                  stackIndex--;             /* Discard ")" */
                  stackValues[stackIndex] = stackValues[stackIndex+1];
                  leftNumberFlag = true;
                }
                else {
                  if (charValue >= '0' && charValue <= '9') {
                    exprIndexAux = exprIndex;
                    if (charValue == '0' && exprIndexAux < exprLength-2 && 
                    		expr.charAt(exprIndexAux+1) == 'x') {  // hexadecimal
                    	exprIndexAux+=2;
                      while (exprIndexAux < exprLength-1) {
                        charValue = expr.charAt(exprIndexAux+1);
                        if ((charValue >= '0' && charValue <= '9') ||
                        		(charValue >= 'A' && charValue <= 'F') ||
                        		(charValue >= 'a' && charValue <= 'f')) {
                          exprIndexAux++;
                        }
                        else {
                          break;
                        }
                      }
                      stackValues[stackIndex] = new BigInteger(expr.substring(exprIndex+2,exprIndexAux+1),16);
                    }
                    else {                   // Decimal number.
                      while (exprIndexAux < exprLength-1) {
                        charValue = expr.charAt(exprIndexAux+1);
                        if (charValue >= '0' && charValue <= '9') {
                          exprIndexAux++;
                        }
                        else {
                          break;
                        }
                      }
                      stackValues[stackIndex] = new BigInteger(expr.substring(exprIndex,exprIndexAux+1));
                    }
                    leftNumberFlag = true;
                    exprIndex = exprIndexAux;
                  }                  /* end if number */
                }                    /* end if ) */
              }                      /* end if ( */
            }                        /* end if ^ */
          }                          /* end if *, / */
        }                            /* end if +, - */
        exprIndex++;
      }                              /* end while */
      if (leftNumberFlag == false) {return -6;}
      if (stackIndex > 0 && stackOperators[stackIndex-1] != '(') {
        if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
          return SubExprResult;
        }
        if (stackIndex > 0 && stackOperators[stackIndex-1] != '(') {
          if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
            return SubExprResult;
          }
          if (stackIndex > 0 && stackOperators[stackIndex-1] != '(') {
            if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
              return SubExprResult;
            }
          }
        }
      }
      if (stackIndex != 0) {return -5;}
      if (stackValues[0].compareTo(BigInt1) <= 0 && type==0)
      if (stackValues[0].bitLength() > 33219) {return -2;}
      ExpressionResult[0] = stackValues[0];
      return 0;
  }

  private static int modInv(int NbrMod, int currentPrime)
  {
    int QQ, T1, T3;
    int V1 = 1;
    int V3 = NbrMod;
    int U1 = 0;
    int U3 = currentPrime;
    while (V3 != 0)
    {
      if (U3 < V3+V3)
      {               // QQ = 1
        T1 = U1 - V1;
        T3 = U3 - V3;
      }
      else
      {
        QQ = U3 / V3;
        T1 = U1 - V1 * QQ;
        T3 = U3 - V3 * QQ;
      }
      U1 = V1;
      U3 = V3;
      V1 = T1;
      V3 = T3;
    }
    return U1 + (currentPrime & (U1 >> 31));
  }

   // Compute the partitions of an integer p(n)
   // The formula p(n) = p(k) = p(k - 1) + p(k - 2) - p(k - 5) -
   // - p(k - 7) + p(k - 12) + p(k - 15) - p(k - 22) - ...
   // where p(0) = 1, where the numbers have the form n(3n-1)/2 for
   // n = 1, -1, 2, -2, 3, -3, ... stopping when k-n(3n-1)/2 is negative.
   // The signs in the summation are: +, +, -, -, +, +, -, -, ...
   // Modular arithmetic is used modulo the largest primes under 2^31, and
   // the result is reconstructed.
   // The approximation p(n) = e^(pi*sqrt(2n/3))/(4n*sqrt(3)) is used
   // in order to find the number of limbs needed to compute the partition.
  static BigInteger partition(int val)
  {
    int [] partArray;
    int index, currentPrime, Q, k, n, sum, idx;
    long numerator, prodmod;
    // Compute approximate number of limbs: log(p(n))/log(2^31)
    // pi * sqrt(2/3)/log(2^31) < 0.12, so 0.12 is selected.
    int limbs = (int)(0.12*Math.sqrt(val)+1);
    partArray = new int[limbs*2+val];

    // Compute the primes which will be used for the modular arithmetic
    // operations. The primes must be ordered in ascending order.
    currentPrime = 0x7FFFFFFF;    // Greatest prime number below 2^31.
    for (index=limbs+val-1; index>=val; index--)
    {
      partArray[index] = currentPrime;
calculate_previous_prime_loop : for (;;)
      {
        currentPrime -= 2;
        if (currentPrime % 3 == 0)
        {
          continue calculate_previous_prime_loop; /* Composite */
        }
        for (Q = 5; Q <= 46341; Q += 6)  // 46341 is the square root of 2^31. 
        { /* Check if Base is prime */
          if (currentPrime % Q == 0)
          {
            continue calculate_previous_prime_loop; /* Composite */
          }
          if (currentPrime % (Q+2) == 0)
          {
            continue calculate_previous_prime_loop; /* Composite */
          }
        }
        break; /* Prime found */
      }
    }
        // Perform modular arithmetic.
    for (index=val; index<val+limbs; index++)
    {
      currentPrime = partArray[index];
      sum = 1;                 // Initialize p(0) mod currentPrime.
      for (k=1; k<=val; k++)            // Generate all partition numbers
      {                                 // up to the one wanted.
        idx = k;
        partArray[k-1] = sum;           // Store p(k-1) mod currentPrime.
        sum = 0;
        n = 1;
        for (;;)                        // Loop for n.
        {
          idx -= n+n-1;
          if (idx < 0)
          {
            break;                      // Negative index, so go out.
          }
          sum -= currentPrime - partArray[idx];
          sum += currentPrime & (sum >> 31);
          idx -= n;
          if (idx < 0)
          {
            break;                      // Negative index, so go out.
          }
          sum -= currentPrime - partArray[idx];
          sum += currentPrime & (sum >> 31);
          n++;
          idx -= n+n-1;
          if (idx < 0)
          {
            break;                      // Negative index, so go out.
          }
          sum -= partArray[idx];
          sum += currentPrime & (sum >> 31);
          idx -= n;
          if (idx < 0)
          {
            break;                      // Negative index, so go out.
          }
          sum -= partArray[idx];
          sum += currentPrime & (sum >> 31);
          n++;
        }
      }
      partArray[index+limbs] = sum;
    }
    // Reconstruct the result from p(val) mod all primes.
    // v_1 <- u_1 mod m_1
    // v_2 <- (u_2 - v_1) (m_1)^(-1) mod m_2
    // v_3 <- (u_3 - (v_1+m_1 v_2)) (m_1 m_2)^(-1) mod m_3
    // ...
    // v_r <- (u_n - (v_1+m_1(v_2+m_2(v_3+...+m_(r-2) v_(r-1)...)))*
    //        (m_1 m_2 ... m_(r-1))^(-1) mod m_r
    //
    // u = v_r m_(r-1) ..m_2 m_1 + v_3 m_2 m_1 + v_2 m_1 + v_1
    partArray[0] = partArray[val+limbs];
    for (index=1; index<limbs; index++)
    {
      currentPrime = partArray[val+index];
      prodmod = 1;
      for (k=index-1; k>=0; k--)
      {
        prodmod = prodmod * partArray[val+k] % currentPrime;
      }
      prodmod = modInv((int)prodmod, currentPrime);
      numerator = partArray[index-1];
      for (k=index-2; k>=0; k--)
      {
        numerator = (numerator*partArray[val+k] + partArray[k]) %
                    currentPrime;
      }
      sum = partArray[val+limbs+index] - (int)numerator;
      if (sum<0)
      {
        sum += currentPrime;
      }
      partArray[index] = (int)(sum*prodmod%currentPrime);
    }
    BigInteger result = BigInteger.valueOf(partArray[0]);
    BigInteger prodModulus = BigInteger.valueOf(1);
    for (index=1; index<limbs; index++)
    {
      prodModulus = prodModulus.multiply(BigInteger.valueOf(partArray[val+index-1]));
      result = result.add(prodModulus.multiply(BigInteger.valueOf(partArray[index])));
    }
    return result;
  }

  private static int ComputeSubExpr(int stackIndex, BigInteger [] stackValues, int [] stackOperators) {
      int i, j, len, val;
      double logarithm;
      BigInteger FibonPrev, FibonAct, FibonNext;
      int stackOper;

      stackOper = stackOperators[stackIndex];
      switch (stackOper) {
        case '+':
          stackValues[stackIndex] = stackValues[stackIndex].add(stackValues[stackIndex+1]);
          return 0;
        case '-':
          stackValues[stackIndex] = stackValues[stackIndex].subtract(stackValues[stackIndex+1]);
          return 0;
        case '_':
          stackValues[stackIndex] = stackValues[stackIndex+1].negate();
          return 0;
        case '/':
          if (stackValues[stackIndex + 1].signum() == 0) {return -3;}
          if (stackValues[stackIndex].remainder(stackValues[stackIndex + 1]).
            signum() != 0) {return -4;}
          stackValues[stackIndex] = stackValues[stackIndex].divide(stackValues[stackIndex+1]);
          return 0;
        case '%':
          if (stackValues[stackIndex + 1].signum() != 0) {
            stackValues[stackIndex] = stackValues[stackIndex].remainder(stackValues[stackIndex+1]);
            }
          return 0;
        case '*':
          if (stackValues[stackIndex].bitLength() + stackValues[stackIndex+1].bitLength() > 66438) {return -3;}
          stackValues[stackIndex] = stackValues[stackIndex].multiply(stackValues[stackIndex+1]);
          return 0;
        case '^':
          len = stackValues[stackIndex].bitLength()-1;
          if (len > 32) {
            logarithm = (double)(len-32) +
                Math.log(stackValues[stackIndex].shiftRight(len-32).
                doubleValue())/Math.log(2);
            }
          else {
            logarithm = Math.log(stackValues[stackIndex].
                doubleValue())/Math.log(2);
            }
          if (logarithm * stackValues[stackIndex+1].doubleValue() > 66438) {return -3;}
          stackValues[stackIndex] = stackValues[stackIndex].pow(stackValues[stackIndex+1].intValue());
          return 0;
        case 'F':
        case 'L':
          len = stackValues[stackIndex+1].bitLength()-1;
          if (len > 17) {return -3;}
          len = stackValues[stackIndex+1].intValue();
          if (len > 95662) {return -3;}
          if (len < 0) {return -8;}
          FibonPrev = BigInteger.valueOf(stackOper == 'L'?-1:1);
          FibonAct = BigInteger.valueOf(stackOper == 'L'?2:0);
          for (i=1; i<=len; i++) {
            FibonNext = FibonPrev.add(FibonAct);
            FibonPrev = FibonAct;
            FibonAct = FibonNext;
            }
          stackValues[stackIndex] = FibonAct;
          return 0;
        case 'P':
          len = stackValues[stackIndex+1].bitLength()-1;
          if (len > 24) {return -3;}
          len = stackValues[stackIndex+1].intValue();
          if (len > 3520000) {return -3;}
          if (len < 0) {return -8;}
          val = stackValues[stackIndex+1].intValue();
          stackValues[stackIndex] = partition(val);
          break;
        case 'B':
        case 'N':
          int Base, Q, baseNbr;
          BigInteger value;
          if (stackOper == 'B') {
            j = stackValues[stackIndex+1].compareTo(BigInt3);
            if (j < 0) {return -8;}
            if (j == 0) {
              stackValues[stackIndex] = BigInt2;
              return 0;
              }
            value = stackValues[stackIndex+1].subtract(BigInt2).or(BigInt1);
            }
          else {
            if (stackValues[stackIndex+1].compareTo(BigInt2) < 0) {return -8;}
            value = stackValues[stackIndex+1].add(BigInt1).or(BigInt1);
            }
outer_calculate_SPRP:
          while (true) {        /* Search for next pseudoprime */
calculate_SPRP:
            do {
              if (value.bitLength() < 16) {
                j = value.intValue();
                if (j >= 9) {
                  for (Q=3; Q*Q<=j; Q+=2) {     /* Check if Base is prime */
                    if (j%Q == 0) {
                      break calculate_SPRP;     /* Composite */
                      }  
                    }
                  }
                break outer_calculate_SPRP;     /* Prime */
                }
              Base = 3;
              for (baseNbr=100; baseNbr>0; baseNbr--) {
                if (value.mod(BigInteger.valueOf(Base)).signum() == 0) {
                  break calculate_SPRP;         /* Composite */
                  }
calculate_new_prime3:
                do {
                  Base+=2;
                  for (Q=3; Q*Q<=Base; Q+=2) {  /* Check if Base is prime */
                    if (Base%Q == 0) {
                      continue calculate_new_prime3;   /* Composite */
                      }  
                    }
                    break;                      /* Prime found */
                  } while (true);
                if (value.mod(BigInteger.valueOf(Base)).signum() == 0) {
                    break calculate_SPRP;       /* Composite */
                  }
                }
              BigInteger valuem1 = value.subtract(BigInt1);
              int exp = valuem1.getLowestSetBit();
              Base = 3;
compute_SPRP_loop:
              for (baseNbr=20; baseNbr>0; baseNbr--) {               
calculate_new_prime4:
                do {
                  Base+=2;
                  for (Q=3; Q*Q<=Base; Q+=2) {  /* Check if Base is prime */
                    if (Base%Q == 0) {
                      continue calculate_new_prime4;   /* Composite */
                      }  
                    }
                    break;                      /* Prime found */
                  } while (true);
                BigInteger bBase = BigInteger.valueOf(Base);
                BigInteger pow = bBase.modPow(valuem1.shiftRight(exp), value);
                if (pow.equals(BigInt1) || pow.equals(valuem1)) {
                  continue compute_SPRP_loop;   /* Strong pseudoprime */
                  }
                for (j=1; j<exp; j++) {
                  pow = pow.multiply(pow).mod(value);
                  if (pow.equals(valuem1)) {
                    continue compute_SPRP_loop; /* Strong pseudoprime */
                    }
                  if (pow.equals(BigInt1)) {
                    break compute_SPRP_loop;    /* Composite */
                    }
                  }
                break;                          /* Composite */
                }                               /* End for */
              if (baseNbr == 0) {
                break outer_calculate_SPRP;     /* Strong pseudoprime */
                }
              } while (false);
            if (stackOper == 'B') {
              value = value.subtract(BigInt2);
              }
            else {
              value = value.add(BigInt2);
              }
            }
          stackValues[stackIndex] = value;
        }              /* end switch */
     return 0;
     }
  }
