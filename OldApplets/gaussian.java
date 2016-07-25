// <XMP>
// Gaussian Integer factorization applet
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated May 31st, 2002, See http://www.alpertron.com.ar/GAUSSIAN.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.applet.*;
import java.util.*;
import java.awt.*;
import java.math.*;

public final class gaussian extends Applet implements Runnable {

private BigInteger Primes[];
private int Exponents[];
private BigInteger Solution1[];
private BigInteger Solution2[];
private BigInteger Increment[];
private BigInteger Aux[];
private BigInteger ValA, ValB;
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

void Solution(BigInteger value) {
  SolNbr++;
  w(SolNbr+") x = "+value.toString()+"<BR>");
  }

public void init() {
  Primes=new BigInteger[400];
  Exponents=new int[400];
  Solution1=new BigInteger[400];
  Solution2=new BigInteger[400];
  Increment=new BigInteger[400];
  Aux=new BigInteger[400];
  }

/* type = 0: factor expression, type = 1: just compute expression */
public String startCalc(String expr, int type) {
  String InputField;
  if (calcThread != null) {
//    TerminateThread = true;
    try {
      calcThread.join();        /* Wait until the solving thread dies */
      } catch (InterruptedException ie) {};
    }
  calcThread = new Thread(this);  /* Start solving thread */
  try {
    int ExpressionRC;
    BigInteger [] ExpressionResult = new BigInteger[2];
    InputField = expr.trim();
    if (InputField.equals("")) {
      ExpressionRC = -6;
      }
    try {
      ExpressionRC = GaussExpression.ComputeExpression(InputField, ExpressionResult);
      } catch (OutOfMemoryError e) {ExpressionRC = -1;}
    switch (ExpressionRC) {
      case -1:
        textError = "Not enough memory.";
        break;
      case -2:
        textError = "Number too big (more than 500 digits).";
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
        if (ExpressionResult[0].compareTo(BigInteger.valueOf(10L).pow(500))>=0) {
          textError = "Number too big (more than 500 digits).";
          }
        else {
          textError = "";
          }
      }
    if (textError.length() > 0) {return textError;}
    ValA = ExpressionResult[0];
    ValB = ExpressionResult[1];
    } catch (Exception e) {return "Invalid data entered";}
  if (type == 1) {     /* Only computing expression */
    info=new StringBuffer();
    w("<HTML><HEAD><TITLE>Gaussian integer calculator</TITLE></HEAD><BODY>");
    w("<CENTER><H3><FONT COLOR=Red>Gaussian integer calculator</FONT></H3></CENTER>");
    w("<P><I>by Dario Alejandro Alpern</I><P><DL><DT>Input expression:<DD>");
    w(InputField);
    w("<P><DT>Value:<DD>");
    w(ValA.toString()+(ValB.signum()>=0?" + ":" - ")+ValB.abs().toString()+" i</DL>");
    w("<P>If you found any mistake, please send me an <A HREF=\"mailto:dario@alpern.com.ar?subject=Gaussian Factorization applet\">e-mail</A><P>");
    w("</BODY></HTML>");
    return info.toString();
    }
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
  w("<TITLE>Gaussian factorization applet</TITLE><CENTER><FONT COLOR=Red><B>");
  w("Factors of "+ValA.toString()+(ValB.signum()>=0?" + ":" - ")+ValB.abs().toString()+" i</B></FONT></CENTER><P><I>by Dario Alejandro Alpern</I><P>");
  Date OldDate=new Date();
  long Old=OldDate.getTime();
  SolNbr = 0;
  GaussianFactorization();
  Date NewDate=new Date();
  long New=NewDate.getTime();
  w("<P>Calculation time: ");
  int t=(int)(((New-Old)/1000+86400)%86400);
  w(t/3600+"h "+((t%3600)/60)+"m "+(t%60)+"s");
  w("<P>If you found any mistake, please send me an <A HREF=\"mailto:alpertron@hotmail.com?subject=Gaussian Factorization applet\">e-mail</A><P>");
  w("<A HREF=\"GBOOK.HTM\">Click here</A> to view or sign my Guestbook.</BODY></HTML>");
  calcThread = null;
  }

void GaussianFactorization() {
  BigInteger BigInt0, BigInt1, BigInt2;
  BigInt0 = BigInteger.valueOf(0L);
  BigInt1 = BigInteger.valueOf(1L);
  BigInt2 = BigInteger.valueOf(2L);
  BigInteger K, Mult1, Mult2, p, q, M1, M2, Tmp;
  int index, index2;

  BigInteger norm = ValA.multiply(ValA).add(ValB.multiply(ValB));
  Ind = 0;
  if (norm.signum() == 0) {
    w("Any gaussian prime divides this number");
    return;
    }
  w("<UL>");
  if (norm.compareTo(BigInt1) > 0) {
    if (norm.bitLength() < 32) {
      long modulus = norm.longValue();
      int Exp = 0;
      Ind = 0;
      while (modulus % 2 == 0) {
        Exp++;
        modulus /= 2;
        }
      if (Exp > 0) {
        Primes[0] = BigInt2;
        Exponents[0] = Exp;
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
          Primes[Ind] = BigInteger.valueOf(Div);
          Exponents[Ind] = Exp;
          Ind++;
          }
        Div += 2;
        }
      if (modulus > 1) {
        Primes[Ind] = BigInteger.valueOf(modulus);
        Exponents[Ind] = 1;
        Ind++;
        }
      }
    else {     /* Factor norm > 2^32 */
      factorPanel = new ecm();
      factorWindow = new Frame("Norm factorization");
      factorWindow.setResizable(false);
      factorWindow.add(factorPanel);
      factorPanel.setSize(600, 390);
      Insets in = factorWindow.getInsets();
      factorWindow.setSize(600+in.right+in.left, 390+in.top+in.bottom);
      factorWindow.show();
      Ind = factorPanel.getFactors(norm, Primes, Exponents);
      factorWindow.remove(factorPanel);
      factorWindow.dispose();
      }
    for (index = 0; index < Ind; index++) {
      p = Primes[index];
      if (p.equals(BigInt2)) {
        for (index2 = 0; index2 < Exponents[index]; index2++) {
          DivideGaussian(BigInt1, BigInt1);           /* Divide by 1+i */
          DivideGaussian(BigInt1, BigInt1.negate());  /* Divide by 1-i */
          }
        }
      if (p.testBit(1) == false) {    /* if p = 1 (mod 4) */
        q = p.subtract(BigInt1);      /* q = p-1 */
        K = BigInt1;
        do {                 // Compute Mult1 = sqrt(-1) mod p
          K = K.add(BigInt1);
          Mult1 = K.modPow(q.shiftRight(2),p);
          } while (Mult1.equals(BigInt1) || Mult1.equals(q));
        Mult2 = BigInt1;
        while (true) {
          K = Mult1.multiply(Mult1).add(Mult2.multiply(Mult2)).divide(p);
          if (K.equals(BigInt1)) {
            break;
            }
          M1 = Mult1.mod(K);
          M2 = Mult2.mod(K);
          if (M1.compareTo(K.shiftRight(1)) > 0) {M1 = M1.subtract(K);}
          if (M2.compareTo(K.shiftRight(1)) > 0) {M2 = M2.subtract(K);}
          Tmp = Mult1.multiply(M1).add(Mult2.multiply(M2)).divide(K);
          Mult2 = Mult1.multiply(M2).subtract(Mult2.multiply(M1)).divide(K);
          Mult1 = Tmp;
          }            /* end while */
        if (Mult1.abs().compareTo(Mult2.abs()) < 0) {
          Tmp = Mult1;
          Mult1 = Mult2;
          Mult2 = Tmp;
          }
        for (index2 = 0; index2 < Exponents[index]; index2++) {
          DivideGaussian(Mult1,Mult2);
          DivideGaussian(Mult1,Mult2.negate());
          }
        }              /* end p = 1 (mod 4) */
      else {           /* if p = 3 (mod 4) */
        for (index2 = 0; index2 < Exponents[index]; index2++) {
          DivideGaussian(Primes[index],BigInt0);
          }            /* end p = 3 (mod 4) */
        }
      }
    }
  if (ValA.equals(BigInt1)) {
    if (Ind == 0) {
      w("No gaussian prime divides this number");
      }
    }
  else if (ValA.add(BigInt1).signum() == 0) {
    w("<LI>-1");
    }
  else if (ValB.equals(BigInt1)) {
    w("<LI>i");
    }
  else {
    w("<LI>-i");
    }
  w("</UL>");
  }

private void DivideGaussian(BigInteger real, BigInteger imag) {
  real = real.abs();
  BigInteger norm = real.multiply(real).add(imag.multiply(imag));
  BigInteger realNum = ValA.multiply(real).add(ValB.multiply(imag));
  BigInteger imagNum = ValB.multiply(real).subtract(ValA.multiply(imag));
  if (realNum.mod(norm).signum() == 0 &&
      imagNum.mod(norm).signum() == 0) {
    ValA = realNum.divide(norm);
    ValB = imagNum.divide(norm);
    w("<LI>");
    w(real.toString()+(imag.signum()>=0?" + ":" - ")+imag.abs().toString()+" i");
    }
  }
}
