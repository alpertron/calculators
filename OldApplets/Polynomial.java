// <XMP>
// Polynomial class
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated January 9th, 2005. See http://www.alpertron.com.ar/POLFACT.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.math.*;
import java.util.*;

public class Polynomial implements Cloneable {
  public int degree;
  public static boolean superscript = true;
  public BigInteger modulus;
  public BigInteger [] coeff;
  private static final BigInteger BigInt0 = BigInteger.valueOf(0L);
  private static final BigInteger BigInt1 = BigInteger.valueOf(1L);
  private static final int maxDegree = 200;

  public Polynomial(int degree, BigInteger modulus) {
    this.degree = degree;
    this.modulus = modulus;
    coeff = new BigInteger[degree+1];
    }

  public static void main(String args[]) {
    Polynomial pol1 = new Polynomial(5, BigInteger.valueOf(0L));
    pol1.coeff[0] = BigInteger.valueOf(-56);
    pol1.coeff[1] = BigInteger.valueOf(117);
    pol1.coeff[2] = BigInteger.valueOf(-221);
    pol1.coeff[3] = BigInteger.valueOf(20);
    pol1.coeff[4] = BigInteger.valueOf(39);
    pol1.coeff[5] = BigInteger.valueOf(5);
    Polynomial pol2 = new Polynomial(5, BigInteger.valueOf(0L));
    pol2.coeff[0] = BigInteger.valueOf(136);
    pol2.coeff[1] = BigInteger.valueOf(-259);
    pol2.coeff[2] = BigInteger.valueOf(173);
    pol2.coeff[3] = BigInteger.valueOf(-72);
    pol2.coeff[4] = BigInteger.valueOf(7);
    pol2.coeff[5] = BigInteger.valueOf(7);
    System.out.println("pol1 = "+pol1.toString());
    System.out.println("pol2 = "+pol2.toString());
    Polynomial pol3 = pol1.gcd(pol2);
    System.out.println("pol3 = "+pol3.toString());
    }

  public static Polynomial one(BigInteger modulus) {
    Polynomial output = new Polynomial(0, modulus);
    output.coeff[0] = BigInt1;
    return output;
    }

  public static Polynomial x(BigInteger modulus) {
    Polynomial output = new Polynomial(1, modulus);
    output.coeff[1] = BigInt1;
    output.coeff[0] = BigInt0;
    return output;
    }

  public Polynomial add(Polynomial other) {
    Polynomial sum;
    int i;
    if (modulus.compareTo(other.modulus) != 0) {
      throw new ArithmeticException("Polynomials have different moduli");
      }
    if (degree > other.degree) {
      sum = new Polynomial(degree, modulus);
      for (i=degree; i>other.degree; i--) {
        sum.coeff[i] = coeff[i];
        }
      }
    else if (degree < other.degree) {
      sum = new Polynomial(other.degree, modulus);
      for (i=other.degree; i>degree; i--) {
        sum.coeff[i] = other.coeff[i];
        }
      }
    else {
      if (modulus.equals(BigInt0)) {
        for (i=degree; i>=0; i--) {
          if (coeff[i].add(other.coeff[i]).signum() != 0) break;
          }
        }
      else {
        for (i=degree; i>=0; i--) {
          if (coeff[i].add(other.coeff[i]).mod(modulus).signum() != 0) break;
          }
        }
      sum = new Polynomial(i, modulus);
      }
    if (modulus.equals(BigInt0)) {
      for (; i>=0; i--) {
        sum.coeff[i] = coeff[i].add(other.coeff[i]);
        }
      }
    else {
      for (; i>=0; i--) {
        sum.coeff[i] = coeff[i].add(other.coeff[i]).mod(modulus);
        }
      }
    return sum;
    }

  public Object clone() {
    Polynomial clone = new Polynomial(degree,modulus);
    for (int i=degree; i >= 0; i--) {
      clone.coeff[i] = coeff[i];
      }
    return clone;
    }

  public Polynomial subtract(Polynomial other) {
    Polynomial diff;
    int i;

    if (modulus.compareTo(other.modulus) != 0) {
      throw new ArithmeticException("Polynomials have different moduli");
      }
    if (degree > other.degree) {
      diff = new Polynomial(degree, modulus);
      for (i=degree; i>other.degree; i--) {
        diff.coeff[i] = coeff[i];
        }
      }
    else if (degree < other.degree) {
      diff = new Polynomial(other.degree, modulus);
      if (modulus.equals(BigInt0)) {
        for (i=other.degree; i>degree; i--) {
          diff.coeff[i] = other.coeff[i].negate();
          }
        }
      else {
        for (i=other.degree; i>degree; i--) {
          diff.coeff[i] = other.coeff[i].negate().mod(modulus);
          }
        }
      }
    else {
      if (modulus.equals(BigInt0)) {
        for (i=degree; i>=0; i--) {
          if (coeff[i].subtract(other.coeff[i]).signum() != 0) break;
          }
        }
      else {
        for (i=degree; i>=0; i--) {
          if (coeff[i].subtract(other.coeff[i]).mod(modulus).signum() != 0) break;
          }
        }
      diff = new Polynomial(i, modulus);
      }
    if (modulus.equals(BigInt0)) {
      for (; i>=0; i--) {
        diff.coeff[i] = coeff[i].subtract(other.coeff[i]);
        }
      }
    else {
      for (; i>=0; i--) {
        diff.coeff[i] = coeff[i].subtract(other.coeff[i]).mod(modulus);
        }
      }
    return diff;
    }

  public Polynomial multiply(Polynomial other) {
    Polynomial prod;
    int i,j;

    if (modulus.compareTo(other.modulus) != 0) {
      throw new ArithmeticException("Polynomials have different moduli");
      }
    if (degree < 0) return this;
    if (other.degree < 0) return other;
    prod = new Polynomial(degree + other.degree, modulus);
    for (i=prod.degree; i>=0; i--) {
      prod.coeff[i] = BigInt0;
      }
    for (i=degree; i>=0; i--) {
      for (j=other.degree; j>=0; j--) {
        prod.coeff[i+j] = prod.coeff[i+j].
                          add(coeff[i].multiply(other.coeff[j]));
        }
      }
    if (modulus.equals(BigInt0) == false) {
      for (i=degree+other.degree; i>=0; i--) {
        prod.coeff[i] = prod.coeff[i].mod(modulus);
        }
      }
    return prod;
    }

  public Polynomial multiply(BigInteger other) {
    if (modulus.equals(BigInt0) == false) other = other.mod(modulus);
    if (other.signum() == 0) {
      return new Polynomial(-1, modulus);
      }
    Polynomial output = new Polynomial(degree, modulus);
    if (modulus.equals(BigInt0)) {
      for (int i=degree; i>=0; i--) {
        output.coeff[i] = coeff[i].multiply(other);
        }
      }
    else {
      for (int i=degree; i>=0; i--) {
        output.coeff[i] = coeff[i].multiply(other).mod(modulus);
        }
      }
    return output;
    }

  public Polynomial divide(Polynomial other) {
    Polynomial quot, rem, div;
    int i,j;

    if (modulus.compareTo(other.modulus) != 0) {
      throw new ArithmeticException("Polynomials have different moduli");
      }
    if (degree < other.degree) {
      return new Polynomial(-1, modulus);
      }
    rem = new Polynomial(degree, modulus);
    quot = new Polynomial(degree-other.degree, modulus);
    div = new Polynomial(other.degree, modulus);
    BigInteger d = other.coeff[other.degree].modInverse(modulus);
    for (i=degree; i>=0; i--) {
      rem.coeff[i] = coeff[i].multiply(d).mod(modulus);
      }
    for (i=other.degree; i>=0; i--) {
      div.coeff[i] = other.coeff[i].multiply(d).mod(modulus);
      }
      /* Perform pseudo-division of rem/div */
      /* The quotient will be in quot */

    for (i=rem.degree - div.degree; i>=0; i--) {
      quot.coeff[i] = rem.coeff[i+div.degree];
      for (j=div.degree-1; j>=0; j--) {
        rem.coeff[i+j] = rem.coeff[i+j].
            subtract(quot.coeff[i].multiply(div.coeff[j])).mod(modulus);
        }
      }
    return quot;
    }
  
  public Polynomial mod(Polynomial other) {
    Polynomial rem, div;
    int i,j;
    BigInteger quot;

    if (modulus.compareTo(other.modulus) != 0) {
      throw new ArithmeticException("Polynomials have different moduli");
      }
    if (other.degree < 0) return (Polynomial)this.clone();
    if (degree < other.degree) {
      return (Polynomial)this.clone();
      }
    rem = new Polynomial(degree, modulus);
    div = new Polynomial(other.degree, modulus);
    BigInteger d = other.coeff[other.degree].modInverse(modulus);
    for (i=degree; i>=0; i--) {
      rem.coeff[i] = coeff[i].multiply(d).mod(modulus);
      }
    for (i=other.degree; i>=0; i--) {
      div.coeff[i] = other.coeff[i].multiply(d).mod(modulus);
      }
      /* Perform pseudo-division of rem/div */
      /* The remainder will be in quot */

    for (i=rem.degree - div.degree; i>=0; i--) {
      quot = rem.coeff[i+div.degree];
      for (j=div.degree-1; j>=0; j--) {
        rem.coeff[i+j] = rem.coeff[i+j].
            subtract(quot.multiply(div.coeff[j])).mod(modulus);
        }
      }

      /* Compute new degree of remainder */

    for (i=div.degree-1; i>=0; i--) {
      if (rem.coeff[i].signum() != 0) {
        break;
        }
      }
    rem.degree = i;
    d = d.modInverse(modulus);
    for (; i>=0; i--) {
      rem.coeff[i] = rem.coeff[i].multiply(d).mod(modulus);
      }
    return rem;
    }
  
  public Polynomial differential() {
    Polynomial diff;
    int i;

    if (degree <= 0) {
      return new Polynomial(-1, modulus);
      }
    if (modulus.equals(BigInt0)) {
      diff = new Polynomial(degree-1, modulus);
      for (i=degree; i>0; i--) {
        diff.coeff[i-1] = coeff[i].multiply(BigInteger.valueOf(i));
        }
      }
    else {
      for (i=degree; i>0; i--) {
        if (coeff[i].multiply(BigInteger.valueOf(i)).mod(modulus).signum() != 0) {
          break;
          }
        }
      diff = new Polynomial(i-1, modulus);
      for (; i>0; i--) {
        diff.coeff[i-1] = coeff[i].multiply(BigInteger.valueOf(i)).mod(modulus);
        }
      }
    return diff;
    }

  public BigInteger getLeadingCoeff() {
    return coeff[degree];
    }

  public Polynomial monic() {
    Polynomial monic = new Polynomial(degree, modulus);
    BigInteger d = coeff[degree].modInverse(modulus);
    for (int i=degree; i>=0; i--) {
      monic.coeff[i] = coeff[i].multiply(d).mod(modulus);
      }
    return monic;
    }

  public Polynomial gcd(Polynomial other) {
    Polynomial result;
    Polynomial first;
    Polynomial second;
    Polynomial thispr;
    Polynomial otherpr;
    Polynomial polyTemp;
    Polynomial temp;
    BigInteger d, prime, W, Z, G;
    long longprime, Q;
    int i,j;

    if (modulus.compareTo(other.modulus) != 0) {
      throw new ArithmeticException("Polynomials have different moduli");
      }
    if (modulus.equals(BigInt0)) {
      thispr = this.primitive();
      otherpr = other.primitive();
      if (thispr.degree < 0) return otherpr;
      if (otherpr.degree < 0) return thispr;
      prime = BigInteger.valueOf(longprime = 48731L);
      while (thispr.coeff[degree].mod(prime).signum() == 0 ||
             otherpr.coeff[other.degree].mod(prime).signum() == 0) {
calculate_new_prime:
        do {
          longprime += 2;
          for (Q = 3; Q*Q <= longprime; Q += 2) { // Check if it is prime
            if (longprime%Q == 0) {
              continue calculate_new_prime;       // Loop if composite
              }
            }
          break;
          } while (true);
        prime = BigInteger.valueOf(longprime);
        }
      first = new Polynomial(degree, prime);
      for (i=degree; i>=0; i--) {
        first.coeff[i] = coeff[i].mod(prime);
        }
      second = new Polynomial(other.degree, prime);
      for (i=other.degree; i>=0; i--) {
        second.coeff[i] = other.coeff[i].mod(prime);
        }
      if (first.gcd(second).degree == 0) {
        return one(BigInt0);  // no common roots
        }
      do {
        if (thispr.degree < otherpr.degree) {
          temp = thispr; thispr = otherpr; otherpr = temp;
          }
        W = thispr.coeff[thispr.degree];
        Z = otherpr.coeff[otherpr.degree];
        G = W.gcd(Z);
        W = W.divide(G); Z = Z.divide(G);
        if (W.equals(BigInt1) == false) {
          otherpr = otherpr.multiply(W);
          }
        if (Z.equals(BigInt1) == false) {
          thispr = thispr.multiply(Z);
          }
        i = thispr.degree - otherpr.degree;
        for (j=otherpr.degree-1; j>=0; j--) {
          thispr.coeff[i+j] = thispr.coeff[i+j].subtract(otherpr.coeff[j]);
          }
        for (j=thispr.degree-1; j>=0; j--) {
          if (thispr.coeff[j].signum() != 0) break;
          }
        thispr.degree = j;
        if (j < 0) break;
        thispr = thispr.primitive();
        } while (true);
      return otherpr.primitive();
      }
    else {                    // modulus not zero
      if (degree < 0) return other;
      if (other.degree < 0) return this;
      if (degree > other.degree) {
        first = new Polynomial(degree, modulus);
        for (i=degree; i>=0; i--) {
          first.coeff[i] = coeff[i];
          }
        second = new Polynomial(other.degree, modulus);
        for (i=other.degree; i>=0; i--) {
          second.coeff[i] = other.coeff[i];
          }
        }
      else {
        first = new Polynomial(other.degree, modulus);
        for (i=other.degree; i>=0; i--) {
          first.coeff[i] = other.coeff[i];
          }
        second = new Polynomial(degree, modulus);
        for (i=degree; i>=0; i--) {
          second.coeff[i] = coeff[i];
          }
        }
      while (true) {
        /* Convert second to monic polynomial */
        d = second.coeff[second.degree].modInverse(modulus);
        for (i=second.degree; i>=0; i--) {
          second.coeff[i] = second.coeff[i].multiply(d).mod(modulus);
          }
      
        /* Perform pseudo-division of first/second */
        /* The remainder will be in first */
      
        for (i=first.degree - second.degree; i>=0; i--) {
          for (j=second.degree-1; j>=0; j--) {
            first.coeff[i+j] = first.coeff[i+j].
                subtract(first.coeff[i+second.degree].multiply(second.coeff[j])).
                mod(modulus);
            }
          }

        /* Compute new degree of remainder */

        for (i=second.degree-1; i>=0; i--) {
          if (first.coeff[i].signum() != 0) {
            break;
            }
          }
        first.degree = i;
        if (i<0) break;
        if (i==0) {
          second.degree = 0;
          second.coeff[0] = BigInt1;
          break;
          }
        polyTemp = first;
        first = second;
        second = polyTemp;
        }
      first = new Polynomial(second.degree, modulus);
         /* Compute first -> primitive part of second */
      d = second.coeff[second.degree].modInverse(modulus);
      for (i=second.degree; i>=0; i--) {
        first.coeff[i] = second.coeff[i].multiply(d).mod(modulus);
        }
      return first;
      }
    }

  public Polynomial primitive() {
    int i;
    BigInteger G;

    if (degree < 0) return this;
    G = coeff[degree];
    for (i=degree-1; i>=0; i--) {
      G = G.gcd(coeff[i]);
      if (G.equals(BigInt1)) break;
      }
    for (i=degree; i>=0; i--) {
      coeff[i] = coeff[i].divide(G);
      }
    return this;
    }

  public static Polynomial random(int degree, BigInteger modulus) {
    int deg;
    BigInteger temp;
    Random random;
    Polynomial output = new Polynomial(degree, modulus);
    int length = modulus.bitLength()+1;
    random = new Random();
    for (deg=0; deg<=degree; deg++) {
      temp = new BigInteger(length, random);
      output.coeff[deg] = temp.mod(modulus);
      }
    for (deg=degree; deg>=0; deg--) {
      if (output.coeff[deg].signum() != 0) {
        break;
        }
      }
    output.degree = deg;
    return output;
    }

  public String toString() {
    int i;
    StringBuffer output = new StringBuffer();

    if (degree < 0) return "0";
    for (i=degree; i>=0; i--) {
      if (coeff[i].signum() > 0) {
        if (i<degree) {
          output = output.append(" + ");
          }
        if (coeff[i].equals(BigInt1) == false) {
          output = output.append(coeff[i].toString());
          }
        else if (i==0) {
          output = output.append("1");
          }
        if (i==1) {
          output = output.append("x");
          }
        else if (i>1) {
          if (superscript) {
            output = output.append("x<SUP>"+i+"</SUP>");
            }
          else {
            output = output.append("x^"+i);
            }
          }
        }
      else if (coeff[i].signum() < 0) {
        if (i<degree) {
          output = output.append(" - ");
          }
        else {
          output = output.append("-");
          }
        if (coeff[i].add(BigInt1).signum() != 0) {
          output = output.append(coeff[i].negate().toString());
          }
        else if (i==0) {
          output = output.append("1");
          }
        if (i==1) {
          output = output.append("x");
          }
        else if (i>1) {
          if (superscript) {
            output = output.append("x<SUP>"+i+"</SUP>");
            }
          else {
            output = output.append("x^"+i);
            }
          }
        }
      }
    return output.toString();
    }

  public static int expr(String expr, BigInteger modulus, Polynomial [] result) {
    int stackIndex = 0;
    int exprIndex = 0;
    int exprLength = expr.length();
    int i,j;
    char charValue;
    boolean leftNumberFlag = false;
    int exprIndexAux;
    int SubExprResult,len;
    Polynomial stackValues[] = new Polynomial[400];
    int stackOperators[] = new int[400];
    BigInteger temp;

    while (exprIndex < exprLength) {
      charValue = expr.charAt(exprIndex);
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
      else if (charValue == '*' || charValue == '/' ||
               charValue == 'x' || charValue == 'X' ) {
        if (leftNumberFlag == false) {
          if (charValue == '*' || charValue == '/') {
            return -6;
            }
          stackValues[stackIndex] = Polynomial.x(modulus);
          leftNumberFlag = true;
          }
        else {
          if (stackIndex > 0 && (stackOperators[stackIndex-1] == '^' ||
                stackOperators[stackIndex-1] == '*' ||
                stackOperators[stackIndex-1] == '/')) {
            if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
              return SubExprResult;
              }
            if (stackIndex > 0 && (stackOperators[stackIndex-1] == '^' ||
                  stackOperators[stackIndex-1] == '*' ||
                  stackOperators[stackIndex-1] == '/')) {
              if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                return SubExprResult;
              }
            }                         /* end if */
          }                           /* end if */
          stackOperators[stackIndex++] = (charValue == '/'? '/':'*');
          leftNumberFlag = false;
          if (charValue == 'x' || charValue == 'X' ) {
            exprIndex--;
          }
        }                             /* end if */
      }                             
      else if (charValue == '^') {
        if (leftNumberFlag == false) {return -6;}
        if (stackIndex > 0 && stackOperators[stackIndex-1] == '^') {
          if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
            return SubExprResult;
          }
        }                         /* end if */
        stackOperators[stackIndex++] = charValue;
        leftNumberFlag = false;
      }                           /* end if */
      else if (charValue == '(') {
        if (leftNumberFlag == true) {return -6;}
        if (stackIndex > 395) {return -7;}
        stackOperators[stackIndex++] = charValue;
      }                           
      else if (charValue == ')') {
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
        else if (charValue >= '0' && charValue <= '9') {
          exprIndexAux = exprIndex;
          while (exprIndexAux < exprLength-1) {
            charValue = expr.charAt(exprIndexAux+1);
            if (charValue >= '0' && charValue <= '9') {
              exprIndexAux++;
            }
            else {
              break;
            }
          }
          temp = new BigInteger(expr.substring(exprIndex,exprIndexAux+1));
          if (temp.signum() == 0) {
            stackValues[stackIndex] = new Polynomial(-1, modulus);
            }
          else {
            stackValues[stackIndex] = new Polynomial(0, modulus);
            stackValues[stackIndex].coeff[0] = temp;
            }
          leftNumberFlag = true;
          exprIndex = exprIndexAux;
        }
        else if (charValue == 'x' || charValue == 'X') {
          if (leftNumberFlag) {
            if (stackIndex > 0 && (stackOperators[stackIndex-1] == '^' ||
                stackOperators[stackIndex-1] == '*' ||
                stackOperators[stackIndex-1] == '/')) {
              if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                return SubExprResult;
              }
              if (stackIndex > 0 && (stackOperators[stackIndex-1] == '^' ||
                   stackOperators[stackIndex-1] == '*' ||
                   stackOperators[stackIndex-1] == '/')) {
                if ((SubExprResult = ComputeSubExpr(--stackIndex, stackValues, stackOperators)) != 0) {
                  return SubExprResult;
                }
              }                         /* end if */
            }                           /* end if */
            stackOperators[stackIndex++] = '*';
          }
          stackValues[stackIndex] = Polynomial.x(modulus);
          leftNumberFlag = true;
        }
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
    if (stackValues[0].degree > maxDegree/2) {return -2;}
    result[0] = stackValues[0];
    for (i=result[0].degree; i>=0; i--) {
      result[0].coeff[i] = result[0].coeff[i].mod(modulus);
      }
    for (i=result[0].degree; i>=0; i--) {
      if (result[0].coeff[i].signum() != 0) {
        break;
        }
      }
    result[0].degree = i;
    return 0;
  }

  private static int ComputeSubExpr(int stackIndex, Polynomial [] stackValues, int [] stackOperators) {
    int i, j, len;
    int stackOper;
    long Cy, Result;
    stackOper = stackOperators[stackIndex];
    BigInteger bigExp;
    switch (stackOper) {
      case '+':
        stackValues[stackIndex] = stackValues[stackIndex].add(stackValues[stackIndex+1]);
        return 0;
      case '-':
        stackValues[stackIndex] = stackValues[stackIndex].subtract(stackValues[stackIndex+1]);
        return 0;
      case '_':
        stackValues[stackIndex] = stackValues[stackIndex+1];
        for (i=stackValues[stackIndex].degree; i>=0; i--) {
          stackValues[stackIndex].coeff[i] = stackValues[stackIndex].coeff[i].negate();
        }
        return 0;
      case 'x':
        stackValues[stackIndex] = new Polynomial(stackValues[stackIndex+1].degree+1, stackValues[stackIndex+1].modulus);
        for (i=stackValues[stackIndex].degree; i>0; i--) {
          stackValues[stackIndex].coeff[i] = stackValues[stackIndex+1].coeff[i-1];
        }
        stackValues[stackIndex].coeff[0] = BigInt0;
        return 0;
      case '/':
        if (stackValues[stackIndex+1].degree < 0) return -3;
        stackValues[stackIndex] = stackValues[stackIndex].divide(stackValues[stackIndex+1]);
        return 0;
      case '*':
        if (stackValues[stackIndex].degree + stackValues[stackIndex+1].degree > maxDegree) {return -3;}
        stackValues[stackIndex] = stackValues[stackIndex].multiply(stackValues[stackIndex+1]);
        return 0;
      case '^':
        if (stackValues[stackIndex+1].degree > 0) return -9;
        if (stackValues[stackIndex+1].degree < 0) {
          stackValues[stackIndex] = Polynomial.one(stackValues[stackIndex].modulus);
          return 0;
          }
        if (stackValues[stackIndex+1].coeff[0].signum() < 0) return -9;
        if (stackValues[stackIndex].degree < 0) return 0;
        bigExp = stackValues[stackIndex+1].coeff[0];
        len = bigExp.bitLength()-1;
        if (len > 16) return -3;
        if (stackValues[stackIndex].degree * bigExp.intValue() > maxDegree) {
          return -3;
          }
        Polynomial polyZ = stackValues[stackIndex];
        Polynomial polyW = polyZ;
        for (i=bigExp.bitLength()-2; i>=0; i--) {
          polyZ = polyZ.multiply(polyZ);
          if (bigExp.testBit(i)) {
            polyZ = polyZ.multiply(polyW);
            }
          }
        stackValues[stackIndex] = polyZ;
        return 0;
    }
  return 0;
  }
}
