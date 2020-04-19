// <XMP>
// Quadratic modular equation solver
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated May 2nd, 2002, See http://www.alpertron.com.ar/QUADMOD.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.applet.*;
import java.util.*;
import java.awt.*;
import java.math.*;

public final class quadmod extends Applet implements Runnable {

private BigInteger Primes[];
private int Exponents[];
private BigInteger PrimesBak[];
private int ExponentsBak[];
private BigInteger Solution1[];
private BigInteger Solution2[];
private BigInteger Increment[];
private BigInteger Aux[];
private BigInteger ValA, ValB, ValC, ValN;
private StringBuffer info;
private String txt;
private int SolNbr;
private Frame factorWindow;
private ecm factorPanel;
private volatile Thread calcThread;
private BigInteger LastModulus = BigInteger.valueOf(0);
private int Ind;
private String textError;

void w(String texto) {
  info=info.append(texto);
  }

byte Show(BigInteger num, String str, byte t) {
  int sign;
  BigInteger BigInt1;
  sign = num.signum();
  if (t==2) {
    txt="";
    }
  if (sign!=0) {
    String str1="";
    if ((t&1)!=0 && sign>0) {
      str1+=" +";
      }
    if (sign<0) {
      str1+=" -";
      }
    if (num.abs().compareTo(BigInteger.valueOf(1L))!=0) {
      str1+=" "+num.abs();
      }
    txt=txt+str1+str;
    return (byte)(t|1);
    }
  return t;
  }

void Show1(BigInteger num, byte t) {
  byte u=Show(num,"",t);
  if ((u&1)==0 || num.abs().compareTo(BigInteger.valueOf(1L))==0) {
    txt+=num.abs();
    }
  }

BigInteger fnMod(int Nro) {
  int Val;
  BigInteger Aux;
  Val = Exponents[Nro]/2;
  Aux = Increment[Nro].multiply(BigInteger.valueOf((long)Val));
  return (Val*2 == Exponents[Nro]?Solution1[Nro]:Solution2[Nro]).add(Aux);
  }

void Solution(BigInteger value) {
  SolNbr++;
  w(SolNbr+") x = "+value.toString()+"<BR>");
  }

public void init() {
  Primes=new BigInteger[400];
  Exponents=new int[400];
  PrimesBak=new BigInteger[400];
  ExponentsBak=new int[400];
  Solution1=new BigInteger[400];
  Solution2=new BigInteger[400];
  Increment=new BigInteger[400];
  Aux=new BigInteger[400];
  }

BigInteger value(String nbr) {
  String InputField;
  int ExpressionRC;
  BigInteger [] ExpressionResult = new BigInteger[1];
  InputField = nbr.trim();
  if (InputField.equals("")) {
    return BigInteger.valueOf(0L);
    }
  try {
    ExpressionRC = expression.ComputeExpression(InputField, 1, ExpressionResult);
    } catch (OutOfMemoryError e) {
      textError = "Out of memory.";
      return BigInteger.valueOf(0L);
    }
    switch (ExpressionRC) {
      case -2:
        textError = "Number too high (more than 1000 digits).";
        break;
      case -3:
        textError = "Intermediate expression too high (more than 2000 digits).";
        break;
      case -4:
        textError = "Non-integer division.";
        break;
      case -5:
        textError = "Parentheses mismatch.";
        break;
      case -6:
        textError = "Syntax error.";
        break;
      case -7:
        textError = "Too many parentheses.";
        break;
      case -8:
        textError = "Invalid parameter.";
        break;
      default:
        textError = "";
    }
    return ExpressionResult[0];
  }


public String startCalc(String valA, String valB, String valC, String valN) {
  if (calcThread != null) {
//    TerminateThread = true;
    try {
      calcThread.join();        /* Wait until the solving thread dies */
      } catch (InterruptedException ie) {};
    }
  calcThread = new Thread(this);  /* Start solving thread */
  try {
    ValA = value(valA);
    if (textError.length() > 0) {return "Parameter a: "+textError;}
    ValB = value(valB);
    if (textError.length() > 0) {return "Parameter b: "+textError;}
    ValC = value(valC);
    if (textError.length() > 0) {return "Parameter c: "+textError;}
    ValN = value(valN);
    if (textError.length() > 0) {return "Parameter n: "+textError;}
    } catch (Exception e) {return "Invalid data entered";}
  calcThread.start();
  return "";
  }

public String resultCalc() {
  if (calcThread == null) {
    return info.toString();
    }
  return "";
  }

public void run() {
  String parms1[];
  StringBuffer parms2[];
  parms1=new String[0];
  parms2=new StringBuffer[1];
  info=new StringBuffer();
  w("<TITLE>Quadratic modular equation</TITLE><CENTER><FONT COLOR=Red><B>");
  byte u=Show(ValA," x^2",(byte)2);
  u=Show(ValB," x",u);
  Show1(ValC,u);
  txt+=" = 0 (mod "+ValN.toString()+")";
  w(txt+"</B></FONT></CENTER><P><I>by Dario Alejandro Alpern</I><P>");
  Date OldDate=new Date();
  long Old=OldDate.getTime();
  SolNbr = 0;
  SolveEquation();
  if (SolNbr==0) {
    w("There are no solutions");
    }
  Date NewDate=new Date();
  long New=NewDate.getTime();
  w("<P>Calculation time: ");
  int t=(int)(((New-Old)/1000+86400)%86400);
  w(t/3600+"h "+((t%3600)/60)+"m "+(t%60)+"s");
  w("<P>If you found any mistake or any solution is missing, please send me an <A HREF=\"mailto:dario@alpern.com.ar?subject=Error solving: "+txt+"\">e-mail</A><P>");
  calcThread = null;
  }

void SolveEquation() {
  BigInteger BigInt0, BigInt1, BigInt2;
  BigInteger BigInt3, BigInt4, BigInt5;
  int Nro, N, T1, E, E1, Expo;
  BigInteger GcdAll;
  BigInteger ValNn;
  BigInteger z, Mult;
  BigInteger D, T, P, K1, Y1, V1, W1, Z, Den;
  BigInteger U, V, W, X, Y, K, Q, R, R1, Sol1, L;
  BigInt0 = BigInteger.valueOf(0L);
  BigInt1 = BigInteger.valueOf(1L);
  BigInt2 = BigInteger.valueOf(2L);
  BigInt3 = BigInteger.valueOf(3L);
  BigInt4 = BigInteger.valueOf(4L);
  BigInt5 = BigInteger.valueOf(5L);
  if (ValN.signum() == 0) {      /* Mod 0 => Equation in integer numbers */
    if (ValA.signum() == 0) {    /* Not Quadratic Equation */
      if (ValB.signum() == 0) {  /* Constant Equation */
        if (ValC.signum() == 0) {
          w("The equation is satisfied by any integer <B>x</B>.");
          }
        else {
          return;
          }
        }
      else {                     /* Linear Equation */
        if (ValC.remainder(ValB).signum() == 0) {
          Solution(ValC.divide(ValB).negate());
          }
        else {
          return;
          }
        }
      }
    else {                       /* Quadratic Equation */
      U = ValB.multiply(ValB).subtract(ValA.multiply(ValC).shiftLeft(2));
      if (U.signum() < 0) {return;}
      V = BigInt0;
      if (U.signum() != 0) {
        V = V.setBit((U.bitLength()+1)/2);
        W = U.divide(V).add(V).shiftRight(1);
        while (V.compareTo(W) > 0) {
          V = W;
          W = U.divide(V).add(V).shiftRight(1);
          }
        if (V.multiply(V).compareTo(U) != 0) {return;}
        }
      if (V.subtract(ValB).remainder(ValA.shiftLeft(1)).signum() == 0) {
        Solution (V.subtract(ValB).divide(ValA.shiftLeft(1)));
        }
      if (V.negate().subtract(ValB).remainder(ValA.shiftLeft(1)).signum() == 0) {
        Solution (V.negate().subtract(ValB).divide(ValA.shiftLeft(1)));
        }
      }
    SolNbr = 1;
    return;
    }
  ValN = ValN.abs();
  GcdAll = ValA.gcd(ValB.gcd(ValN));
  if (ValC.mod(GcdAll).signum() != 0) {return;}
  ValA = ValA.divide(GcdAll);
  ValB = ValB.divide(GcdAll);
  ValC = ValC.divide(GcdAll);
  ValN = ValN.divide(GcdAll);
  ValNn = ValN;
  if (ValNn.compareTo(BigInt1) == 0) {
    for (z=BigInt0; z.compareTo(ValNn.multiply(GcdAll)) < 0; z=z.add(BigInt1)) {
      Solution(z);
      }
    return;
    }
  if (ValA.mod(ValN).signum() == 0) {
    if (ValB.gcd(ValN).compareTo(BigInt1) != 0) {return;}
    z = ValC.multiply(ValB.modInverse(ValN)).negate().mod(ValN);
    do {
      Solution(z);
      z = z.add(ValN);
      } while (z.compareTo(ValNn.multiply(GcdAll)) < 0);
    return;
    }
  if (LastModulus.equals(ValN) == false) {
    LastModulus = ValN;
    if (LastModulus.bitLength() < 32) {
      long modulus = LastModulus.longValue();
      int Exp = 0;
      Ind = 0;
      while (modulus % 2 == 0) {
        Exp++;
        modulus /= 2;
        }
      if (Exp > 0) {
        PrimesBak[0] = BigInt2;
        ExponentsBak[0] = Exp;
        Ind++;
        }
      long Div = 3;
      while (Div*Div <= modulus) {
        Exp = 0;
        while (modulus % Div == 0) {
          Exp++;
          modulus /= Div;
          }
        if (Exp > 0) {
          PrimesBak[Ind] = BigInteger.valueOf(Div);
          ExponentsBak[Ind] = Exp;
          Ind++;
          }
        Div += 2;
        }
      if (modulus > 1) {
        PrimesBak[Ind] = BigInteger.valueOf(modulus);
        ExponentsBak[Ind] = 1;
        Ind++;
        }
      }
    else {
      factorPanel = new ecm();
      factorWindow = new Frame("Modulus factorization");
      factorWindow.setResizable(false);
      factorWindow.add(factorPanel);
      factorPanel.setSize(600, 390);
      Insets in = factorWindow.getInsets();
      factorWindow.setSize(600+in.right+in.left, 390+in.top+in.bottom);
      factorWindow.show();
      Ind = factorPanel.getFactors(ValN, PrimesBak, ExponentsBak);
      factorWindow.remove(factorPanel);
      factorWindow.dispose();
      }
    }
  for (Nro=0; Nro<Ind; Nro++) {
    Primes[Nro] = PrimesBak[Nro];
    Exponents[Nro] = ExponentsBak[Nro];
    }
  Q = BigInt0;
  for (Nro=0; Nro<Ind; Nro++) {
    P = Primes[Nro];
    N = Exponents[Nro];
    if (ValA.mod(P).signum() == 0) {     /* If linear equation mod p */
      V = P.pow(N);
      if (ValA.mod(V).signum() == 0) {
        U = ValB.gcd(V);
        Q = V.divide(U);
        Solution1[Nro] = ValC.divide(U).multiply(ValB.divide(U).negate().modInverse(Q)).mod(Q);
        }
      else {
        Q = ValA.gcd(V);
        Sol1 = ValC.multiply(ValB.negate().modInverse(Q));
        T1=0; V=P;
        while (ValA.mod(V).signum() == 0) {
          V = V.multiply(P);
          T1++;
          }
        do {
          L = ValA.multiply(Sol1).add(ValB).multiply(Sol1).add(ValC).divide(Q);
          K1 = ValA.multiply(Sol1).shiftLeft(1).add(ValB);
          K = L.gcd(K1);
          R = (L.divide(K).negate()).multiply(K1.divide(K).modInverse(Q)).mod(Q);
          Sol1 = Sol1.add(R.multiply(Q));
          T1 = T1+T1;
          Q=Q.multiply(Q);
          } while (T1 < N);
        Q = P.pow(N);
        T = Sol1.mod(Q);Mult=BigInt1;
        Solution1[Nro] = T;
        }
      Solution2[Nro] = Solution1[Nro];
      }
    else {                   /* If quadratic equation mod p */
      Expo = 0;
      if (P.equals(BigInt2)) {         /* Prime p is 2 */
        BigInteger coefX2 = ValA;
        BigInteger coefX1 = ValB;
        BigInteger coefX0 = ValC;
        BigInteger increm = BigInt0;
        BigInteger deltaincrem = BigInt1;
        Expo = 0;
        while (Expo < N-1 && coefX1.testBit(0) == false) {
          if (coefX0.testBit(0) == false) {   /* constant term is even */
            if (coefX0.testBit(1)) {          /* constant term = 2 (mod 4) */
              return;         /* No solutions */
              }
            coefX1 = coefX1.shiftRight(1);
            coefX0 = coefX0.shiftRight(2);
            Expo += 2;
            deltaincrem = deltaincrem.shiftLeft(1);
            }
          else {                              /* constant term is odd */
            coefX0 = coefX2.add(coefX1).add(coefX0);
            coefX1 = coefX2.shiftLeft(1).add(coefX1);
            increm = increm.add(deltaincrem);
            }
          }
        if (Expo == N-1) {            /* equation modulo 2 */
          if (coefX1.testBit(0) == false) {
            Q = BigInt1.shiftLeft((N+1)/2);
            Solution1[Nro] = Solution2[Nro] =
                      (coefX0.testBit(0)?BigInt1: BigInt0).shiftLeft(Expo/2).add(increm);
            }
          else {
            if (coefX0.testBit(0)) {    /* All three coefficients odd */
              return;                   /* No solution */
              }
            Q = BigInt1.shiftLeft(Expo/2);
            Solution1[Nro] = Solution2[Nro] = increm;
            }
          }
        else {
          if (Expo == N) {            /* equation modulo 1 */
            Q = BigInt1.shiftLeft(Expo/2);
            Solution1[Nro] = Solution2[Nro] = increm;
            }
          else {        /* coefX1 odd -> 1 solution odd, 1 solution even */
            U = coefX1.multiply(coefX1).subtract(coefX2.multiply(coefX0).shiftLeft(2));
            if (coefX0.testBit(0)) {
              return;          /* The three numbers cannot be all odd */
              }  
            Sol1 = BigInt1; T1 = 1; Q = BigInt2;
            if (N != 1) {
              do {
                L = ValA.multiply(Sol1).add(coefX1).multiply(Sol1).add(coefX0).divide(Q);
                R = L.negate().multiply(Sol1.shiftLeft(1).multiply(coefX2).add(coefX1).modInverse(Q)).mod(Q);
                Sol1 = R.multiply(Q).add(Sol1);
                T1 = T1+T1;
                Q = Q.multiply(Q);
                } while (T1 < N-Expo/2);
              }
            Q = BigInt1.shiftLeft(N-Expo);
            Solution1[Nro] = Sol1.mod(Q).shiftLeft(Expo/2).add(increm);
            Solution2[Nro] = coefX1.multiply(coefX2.modInverse(Q)).add(Sol1).
                             negate().mod(Q).shiftLeft(Expo/2).add(increm);
            if (Solution1[Nro].compareTo(Solution2[Nro]) > 0) {
              Sol1 = Solution1[Nro];
              Solution1[Nro] = Solution2[Nro];
              Solution2[Nro] = Sol1;
              }
            Q = Q.shiftLeft(Expo/2);
            }
          }
        }
      else {                        /* Prime is not 2 */
        U = ValB.multiply(ValB).subtract(ValA.multiply(ValC).shiftLeft(2));
        for (;Expo < N; Expo++) {   /* Find exponent */
          if (U.mod(P).signum() == 0) {U = U.divide(P);}
          else {break;}
          }
        if (Expo == N) {
          Q = P.pow((Expo+1)/2);
          Den = (ValA.shiftLeft(1)).modInverse(Q);
          Solution1[Nro] = Solution2[Nro] = ValB.negate().multiply(Den).mod(Q);
          }
        else {
          if ((Expo & 1) == 1) {return;}     /* No solutions if odd expo */
          Mult = P.pow(Expo/2);
          if (U.modPow(P.subtract(BigInt1).shiftRight(1),P).compareTo(BigInt1) != 0) {return;}
            if (P.mod(BigInt4).compareTo(BigInt3) == 0) {
              T = U.modPow(P.add(BigInt1).shiftRight(2),P);
              }
            else {
              if (P.mod(BigInt4.add(BigInt4)).compareTo(BigInt5) == 0) {
              V = (U.shiftLeft(1)).modPow(P.subtract(BigInt5).shiftRight(3),P);
              T = U.shiftLeft(1).multiply(V).multiply(V).subtract(BigInt1).multiply(U).multiply(V).mod(P);
              }
            else {          /* p = 1 (mod 8) */
              Q = P.subtract(BigInt1);
              E = Q.getLowestSetBit();   /* E >= 3 */
              Q = Q.shiftRight(E);
              X = BigInt1;
              do {
                X = X.add(BigInt1);
                Z = X.modPow(Q, P);
                } while (Z.modPow(BigInt0.setBit(E-1), P).compareTo(BigInt1) == 0);
              Y=Z;
              X = U.modPow(Q.subtract(BigInt1).shiftRight(1), P);
              V = U.multiply(X).mod(P);
              W = V.multiply(X).mod(P);
              while (W.compareTo(BigInt1) != 0) {
                T1=0; D=W;
                while (D.compareTo(BigInt1) != 0) {
                  D = D.multiply(D).mod(P);
                  T1++;
                  }
                D = Y.modPow(BigInt0.setBit(E-T1-1),P);
                Y1 = D.multiply(D).mod(P);
                E = T1;
                V1 = V.multiply(D).mod(P);
                W1 = W.multiply(Y1).mod(P);
                Y = Y1; V = V1; W = W1;
                }    /* end while */
              T=V;
              }      /* end if */
            }        /* end if */
          Sol1 = T; T1 = 1; Q = P;
          if (N != 1) {
            do {
              L = Sol1.multiply(Sol1).subtract(U).divide(Q);
              K = L.gcd(Sol1.shiftLeft(1));
              R = L.negate().divide(K).multiply((Sol1.shiftLeft(1).divide(K)).modInverse(Q)).mod(Q);
              Sol1 = R.multiply(Q).add(Sol1);
              T1 = T1+T1;
              Q = Q.multiply(Q);
              } while (T1 < N-Expo/2);
            Q = P.pow(N-Expo/2);
            }
          T = Sol1.mod(Q).multiply(Mult);
          Den = ValA.shiftLeft(1).modInverse(Q);
          Solution1[Nro] = ValB.negate().subtract(T).multiply(Den).mod(Q);
          Solution2[Nro] = ValB.negate().add(T).multiply(Den).mod(Q);
          if (Solution1[Nro].compareTo(Solution2[Nro]) > 0) {
            Sol1 = Solution1[Nro];
            Solution1[Nro] = Solution2[Nro];
            Solution2[Nro] = Sol1;
            }
          }
        }
      }
    Increment[Nro] = Q;
    Primes[Nro] = P.pow(N);
    Exponents[Nro] = 0;
    }   /* end for */
  do {
    Aux[0] = fnMod(0);
    U = Aux[0];
    Mult = Primes[0];
    for (T1 = 1; T1<Ind; T1++) {
      Aux[T1] = fnMod(T1);
      for (E = 0; E<T1; E++) {
        Aux[T1] = Aux[T1].subtract(Aux[E]).multiply(Primes[E].modInverse(Primes[T1])).mod(Primes[T1]);
        }
      Aux[T1] = Aux[T1].mod(Primes[T1]);
      U = Aux[T1].multiply(Mult).add(U);
      Mult = Mult.multiply(Primes[T1]);
      }   /* end for */
    for (T = BigInt0; T.compareTo(GcdAll) < 0; T=T.add(BigInt1)) {
      Solution(T.multiply(ValNn).add(U));
      }
    for (T1 = Ind-1; T1 >= 0; T1--) {
      if (Solution1[T1].compareTo(Solution2[T1]) == 0) {
        Exponents[T1] += 2;
        }
      else {
        Exponents[T1]++;
        }
      if (BigInteger.valueOf((long)Exponents[T1]).multiply(Increment[T1]).compareTo(Primes[T1].shiftLeft(1)) < 0) {break;}
      Exponents[T1] = 0;
      }   /* end for */
    } while (T1 >= 0);
  }
}
