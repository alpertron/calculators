// <XMP>
// Continued fraction calculator
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated April 28th, 2002, See http://www.alpertron.com.ar/CONTFRAC.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.applet.*;
import java.util.*;
import java.awt.*;
import java.math.*;

public final class fraccont extends Applet {

private BigInteger num, den, delta;
private StringBuffer info;
private String txt;
private BigInteger BigInt0 = BigInteger.valueOf(0L);
private BigInteger BigInt1 = BigInteger.valueOf(1L);
private String textError;

void w(String texto) {
  info=info.append(texto);
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
      textError = "No hay memoria disponible.";
      return BigInteger.valueOf(0L);
    }
    switch (ExpressionRC) {
      case -2:
        textError = "Número muy grande (más de 1000 dígitos).";
        break;
      case -3:
        textError = "Expresión intermedia con más de 2000 dígitos.";
        break;
      case -4:
        textError = "División no entera.";
        break;
      case -5:
        textError = "Error de paréntesis.";
        break;
      case -6:
        textError = "Error de sintaxis.";
        break;
      case -7:
        textError = "Demasiados paréntesis.";
        break;
      case -8:
        textError = "Parámetro inválido.";
        break;
      default:
        textError = "";
    }
    return ExpressionResult[0];
  }

public String startCalc(String num, String delta, String den) {
  try {
    this.num = value(num);
    if (textError.length() > 0) {return "Numerador: "+textError;}
    this.delta = value(delta);
    if (textError.length() > 0) {return "Discriminante: "+textError;}
    this.den = value(den);
    if (textError.length() > 0) {return "Denominador: "+textError;}
    } catch (Exception e) {return "Se entraron datos inválidos.";}
  info=new StringBuffer();
  w("<TITLE>Calculadora de fracciones continuas</TITLE><CENTER><FONT COLOR=Red><B>");
  w("Calculadora de fracciones continuas</B></FONT></CENTER><P><I>por Dario Alejandro Alpern</I><P>");
  Date OldDate=new Date();
  long Old=OldDate.getTime();
  ContFrac();
  Date NewDate=new Date();
  long New=NewDate.getTime();
  w("<P>Tiempo de cálculo: ");
  int t=(int)(((New-Old)/1000+86400)%86400);
  w(t/3600+"h "+((t%3600)/60)+"m "+(t%60)+"s");
  w("<P>Si encontró algún error, por favor envíeme un <A HREF=\"mailto:alpertron@hotmail.com?subject=Calculadora de fracciones continuas\">e-mail</A><P>");
  w("<A HREF=\"GBOOK.HTM\">Apriete aquí</A> para ver o firmar mi libro de invitados.<P>");
  return info.toString();
  }

private void ContFrac() {
  BigInteger GcdAll;
  BigInteger sqrtDelta;
  BigInteger sqrtDeltaNext;
  BigInteger biK, biP, biM;
  BigInteger K, L, M, P, Z;

  if (den.signum() == 0) {
    w("Error: El denominador es cero.");
    return;
    }
  if (delta.signum() < 0) {   /* Complex number */
    w("El número no es real, por lo que no se puede desarrollar en fracciones continuas.");
    return;
    }
  w("x = ");
  if (delta.signum() == 0) {   /* Rational number */
    ShowRational(num, den);
    return;
    }
  if (delta.subtract(num.multiply(num)).remainder(den).signum() != 0) {
    delta = delta.multiply(den).multiply(den);
    num = num.multiply(den.abs());
    den = den.multiply(den.abs());
    }
  sqrtDelta = BigInt0.setBit((delta.bitLength()+1)/2);
  sqrtDeltaNext = delta.divide(sqrtDelta).add(sqrtDelta).shiftRight(1);
  while (sqrtDelta.compareTo(sqrtDeltaNext) > 0) {
    sqrtDelta = sqrtDeltaNext;
    sqrtDeltaNext = delta.divide(sqrtDelta).add(sqrtDelta).shiftRight(1);
    }
  if (sqrtDelta.multiply(sqrtDelta).equals(delta)) {
    ShowRational(num.add(sqrtDelta), den);
    return;
    }
  biP = den;
  biK = biP.signum()>0? sqrtDelta: sqrtDelta.add(BigInt1);
  biK = biK.add(num);
  if (biK.signum() > 0) {
    if (den.signum() > 0) {
      biM = biK.divide(den);
      }
    else {
      biM = den.add(BigInt1).subtract(biK).divide(den.negate());
      }
    }
  else {
    if (den.signum() > 0) {
      biM = biK.add(BigInt1).subtract(den).divide(den);
      }
    else {
      biM = biK.negate().divide(den.negate());
      }
    }

  w(biM.toString());   /* Show integer part */
  biM = biM.multiply(den).subtract(num);
  String sep="+ //";
  int cont=-1;
  K = P = L = M = BigInteger.valueOf(-1L);
  while (cont<0 || K.equals(P) == false || L.equals(M) == false) {
    w(sep);
    sep=", ";
    if (cont<0 &&
        biP.signum()>0 && biP.compareTo(sqrtDelta.add(biM)) <= 0 &&
        biM.signum()>0 && biM.compareTo(sqrtDelta) <= 0) {
      w("<B>");
      K = P = biP;
      L = M = biM;
      cont=0;
      }
    if (cont >= 0) {
                    /* both numerator and denominator are positive */
      P = delta.subtract(M.multiply(M)).divide(P);
      Z = sqrtDelta.add(M).divide(P);
      M = Z.multiply(P).subtract(M);
      cont++;
      }
    else {
      biP = delta.subtract(biM.multiply(biM)).divide(biP);
      if (biP.signum() > 0) {
        Z = sqrtDelta.add(biM).divide(biP);
        }
      else {
        Z = sqrtDelta.add(BigInt1).add(biM).divide(biP);
        }
      biM = Z.multiply(biP).subtract(biM);
      }
    w(Z.toString());      /* Show convergent */
    if (cont > 100000) {  /* Too many convergents */
      w(", ... //</B><BR>donde la parte periódica (truncada después de 100000 convergentes) está marcada en negrita.<P>");
      return;
      }
    }
  w("//</B><BR>donde la parte periódica está marcada en negrita");
  if (cont>1) {
    w(" (el período tiene "+cont+" coeficientes)");
    }
  w(".<P>");
  }

void ShowRational(BigInteger num, BigInteger den) {
  BigInteger Temp;

  BigInteger GcdAll = num.gcd(den);
  num = num.divide(GcdAll);
  den = den.divide(GcdAll);
  if (den.signum() < 0) {
    num = num.negate();
    den = den.negate();
    }
  w(num.subtract(num.mod(den)).divide(den).toString());
  num = num.mod(den);
  String sep = "+ //";
  while (num.signum() > 0) {
    if (den.signum() > 0) {
      w(sep + den.divide(num).toString());
      sep = ", ";
      }
    Temp = den.mod(num);
    den = num;
    num = Temp;
    }
  if (sep.equals(", ")) {
    w("//");
    }
  w("<P>");
  }
}

