// <XMP>
// Polynomial factorization class
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated January 9th, 2005. See http://www.alpertron.com.ar/POLFACT.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.math.*;
import java.applet.*;
import java.util.*;

public class polfact extends Applet implements Runnable {
  static final BigInteger BigInt0 = BigInteger.valueOf(0L);
  static final BigInteger BigInt1 = BigInteger.valueOf(1L);
  static final BigInteger BigInt2 = BigInteger.valueOf(2L);
  private StringBuffer info;
  private volatile Thread calcThread;
  private String textError;
  Polynomial polyU;
  Polynomial polyX;
  Polynomial polyZ;
  Polynomial poly;
  Polynomial [] polyArray;
  int indexPolyArray;
  BigInteger modulus;
  int deg;
  Vector answer;
  int [] expon;
  boolean english;
  private int modexponent;

  public static void main(String args[]) {
    int i,j,k,m;
    String result;
    polfact polfact1 = new polfact();
    polfact1.english = true;
    polfact1.modulus = BigInteger.valueOf(2);
    polfact1.startCalc("x^60+x^6+2","5^2",0);
    do {
      try {
        Thread.sleep(200);
        } catch (Exception e) {}
      result = polfact1.resultCalc();
      } while (result.equals(""));
    System.out.println(result);
    }

  public void init() {
    english = (getParameter("lang").equals("ENGLISH"));
    }

  void w(String texto) {
    info=info.append(texto);
    }

  void factor(Polynomial inputPoly) {
    int i;
    modulus = inputPoly.modulus;
    answer = new Vector();
    if (inputPoly.degree <= 0) {
      answer.addElement(inputPoly);
      expon = new int[1];
      expon[0] = 1;
      return;
      }
    Polynomial leading = new Polynomial(0, modulus);
    polyU = (Polynomial)inputPoly.clone();
    if (polyU.getLeadingCoeff().equals(BigInt1) == false) {
      leading.coeff[0] = polyU.getLeadingCoeff();
      answer.addElement(leading);
      }
    polyArray = new Polynomial[inputPoly.degree+1];
    indexPolyArray = 0;
    polyU = polyU.monic();
    Polynomial polyMF = removeMultiFactor();
    makeZpowerPminusZ();
    deg = 1;
    while (true) {
      polyX = polyU.gcd(polyZ).monic();
      if (polyX.degree != 0) {
        irreducibleFactor(polyX);
        polyU = polyU.divide(polyX);
        }
      if (polyU.degree <= 2*deg+1) {
        break;
        }
      raiseZpowerPminusZ();
      deg++;
      }
    if (polyU.degree > 0) {
      answer.addElement(polyU);
      }

    /* Get multiple factors */

    expon = new int[answer.size()];
    i = (((Polynomial)(answer.elementAt(0))).degree == 0? 1: 0);
    if (i > 0) {
      expon[0] = 1;
      }
    for ( ; i<answer.size(); i++) {
      polyZ = (Polynomial)(answer.elementAt(i));
      expon[i] = 1;
      while (polyMF.mod(polyZ).degree < 0) {
        polyMF = polyMF.divide(polyZ);
        expon[i]++;
        }
      }
    }

  Polynomial polyPower(Polynomial poly, BigInteger exp) {
    Polynomial output;
    int i;

    if (exp.equals(BigInt0)) {
      output = Polynomial.one(modulus);
      }
    else {
      output = (Polynomial)poly.clone();
      for (i=exp.bitLength()-2; i>=0; i--) {
        output = output.multiply(output).mod(polyU);
        if (exp.testBit(i)) {
          output = output.multiply(poly).mod(polyU);
          }
        }
      }
    return output;
    }

  Polynomial removeMultiFactor() {
    int i,j;

    Polynomial polyW, polyG;
    Polynomial polyMF = (Polynomial)polyU.clone();
    while (true) {
      Polynomial polyDU = polyU.differential();
      if (polyDU.degree == 0) break;
      if (polyDU.degree < 0) {
        polyW = (Polynomial)polyU.clone();
        polyU = Polynomial.one(modulus);
        }
      else {
        polyG = polyU.gcd(polyDU).monic();
        if (polyG.degree == 0) break;   /* No multiple factor */
        polyU = polyU.divide(polyG);
        polyW = polyG;
        do {                            /* Extract p-th power factor */
          polyG = polyU.gcd(polyW);
          polyW = polyW.divide(polyG);
          } while (polyG.degree != 0);
        if (polyW.degree <= 0) break;   /* No p-th power */
        }
      j = modulus.intValue();
      polyG = new Polynomial(polyW.degree/j, modulus);
      for (i=polyG.degree; i>=0; i--) {
        polyG.coeff[i] = polyW.coeff[i*j];
        }
      polyU = polyU.multiply(polyG);
      }
    return polyMF.divide(polyU);
    }

  /* Get irreducible factors of degree deg */
  void irreducibleFactor(Polynomial polyF) {
    Polynomial polyW, polyG, polyU, polyZ;
    int i, newDeg;
    Vector temp = new Vector();
    temp.addElement((Polynomial)polyF.clone());
    if (polyF.degree == deg) {
      answer.addElement((Polynomial)polyF.clone());
      return;
      }
    polyU = polyF;
    do {
      newDeg = 2+(int)(Math.floor(Math.random()*Math.random()*(2*deg-1)));
      polyW = Polynomial.random(newDeg, modulus);
      if (modulus.compareTo(BigInt2) == 0) {
        polyZ = (Polynomial)polyW.clone();
        for (i=deg-1; i>0; i--) {
          polyZ = polyZ.multiply(polyZ).mod(polyU);
          polyW = polyZ.add(polyW);
          }
        }
      else {
        polyW = makeWPD(polyW).add(Polynomial.one(modulus));
        }
      polyU = Polynomial.one(modulus);
      for (i=0; i<temp.size(); i++) {
        polyF = (Polynomial)temp.elementAt(i);
        if (polyF.degree < 0) continue;
        polyG = polyF.gcd(polyW).monic();
        if (polyG.degree > 0 && polyG.degree < polyF.degree) {
          polyF = polyF.divide(polyG);
          if (polyF.degree == deg) {
            answer.addElement(polyF);
            polyF = new Polynomial(-1, modulus);
            }
          else {
            polyU = polyU.multiply(polyF);
            }
          temp.setElementAt(polyF, i);
          if (polyG.degree == deg) {
            answer.addElement(polyG);
            }
          else {
            temp.addElement(polyG);
            polyU = polyU.multiply(polyG);
            }
          }
        else {
          polyU = polyU.multiply(polyF);
          }
        }
      } while (polyU.degree != 0);
    return;
    }

  void makeZpowerPminusZ() {
    int i;

    Polynomial polyW;
    if (BigInteger.valueOf(polyU.degree).compareTo(modulus) > 0) {
      polyZ = new Polynomial(modulus.intValue(), modulus);
      polyZ.coeff[polyZ.degree] = BigInt1;
      for (i=polyZ.degree-1; i>=0; i--) {
        polyZ.coeff[i] = BigInt0;
        }
      }
    else {
      polyZ = Polynomial.one(modulus);        // Z <- 1
      polyW = Polynomial.x(modulus);          // W <- X
      for (i=0; i<=modulus.bitLength()-2; i++) {
        if (modulus.testBit(i)) {
          polyZ = polyZ.multiply(polyW).mod(polyU);
          }
        polyW = polyW.multiply(polyW).mod(polyU);
        }
      polyZ = polyZ.multiply(polyW).mod(polyU);
      }
    polyW = polyArray[1] = polyZ;
    for (i=2; i<=polyU.degree; i++) {
      polyW = polyW.multiply(polyZ).mod(polyU);
      polyArray[i] = polyW;
      }
    polyZ = polyZ.subtract(Polynomial.x(modulus)).mod(polyU);
    }

  void raiseZpowerPminusZ() {
    int i;

    Polynomial polyW = polyZ.add(Polynomial.x(modulus));
    polyZ = new Polynomial(0, modulus);
    polyZ.coeff[0] = polyW.coeff[0];
    for (i=1; i<=polyW.degree; i++) {
      polyZ = polyZ.add(polyArray[i].multiply(polyW.coeff[i]));
      }
    polyZ = polyZ.subtract(Polynomial.x(modulus)).mod(polyU);
    }

             // Make PolyW^((P^deg-1)/2)
  Polynomial makeWPD(Polynomial polyW) {
    int i,j;

    polyW = polyPower(polyW, modulus.shiftRight(1));
    Polynomial polyV = polyW;
    for (i=1; i<deg; i++) {
      Polynomial polyX = new Polynomial(0, modulus);
      polyX.coeff[0] = polyV.coeff[0];
      for (j=1; j<=polyV.degree; j++) {
        polyX = polyX.add(polyArray[j].multiply(polyV.coeff[j]));
        }
      polyV = polyX;
      polyW = polyW.multiply(polyV).mod(polyU);
      }
    return polyW;
    }

/* type = 0: factor expression, type = 1: just compute expression */
public String startCalc(String expr, String strmod, int type) {
  String InputField;
  int caretLocation;
  if ((type & 8) != 0) {
    type -= 8;
    Polynomial.superscript = true;
    }
  else {
    Polynomial.superscript = false;
    }
  if (calcThread != null) {
//    TerminateThread = true;
    try {
      calcThread.join();          // Wait until the solving thread dies
      } catch (InterruptedException ie) {};
    }
  calcThread = new Thread(this);  // Start solving thread
  try {
    int ExpressionRC;
    Polynomial [] result = new Polynomial[1];
    InputField = strmod.trim();
    caretLocation = strmod.indexOf("^");
    if (caretLocation >= 0) {
      modulus = new BigInteger(InputField.substring(0, caretLocation).trim());
      modexponent = Integer.parseInt(InputField.substring(caretLocation+1).trim());
      if (modexponent <= 0)
      return english? "Exponent must be greater than zero": "El exponente debe ser mayor que cero";
      }
    else {
      modulus = new BigInteger(InputField);
      modexponent = 1;
      }
    InputField = expr.trim();
    if (modulus.compareTo(BigInt2) < 0 ||
        modulus.isProbablePrime(5) == false) {
      return english? "Modulus not prime": "El módulo no es primo";
      }
    if (InputField.equals("")) {
      ExpressionRC = -6;
      }
    try {
      ExpressionRC = Polynomial.expr(InputField, modulus.pow(modexponent), result);
      } catch (OutOfMemoryError e) {ExpressionRC = -1;}
    switch (ExpressionRC) {
      case -1:
        textError = english?"Not enough memory.": "No hay suficiente memoria.";
        break;
      case -2:
        textError = english?"Polynomial too big (degree more than 100).": "Polinomio muy grande (grado mayor que 100).";
        break;
      case -3:
        textError = english?"Intermediate polynomial too big (degree more than 200).": "Polinomio intermedio muy grande (grado mayor que 200).";
        break;
      case -4:
        textError = english?"Non-integer division.":"División no entera.";
        break;
      case -5:
        textError = english?"Parentheses mismatch.":"Error de paréntesis.";
        break;
      case -6:
        textError = english?"Syntax error.":"Error de sintaxis.";
        break;
      case -7:
        textError = english?"Too many parentheses.":"Demasiados paréntesis.";
        break;
      case -8:
        textError = english?"Invalid parameter.":"Parámetro inválido.";
        break;
      case -9:
        textError = english?"Exponent must be constant and non-negative.":"El exponente debe ser constante y no negativo.";
        break;
      default:
        textError = "";
      }
    if (textError.length() > 0) {return textError;}
    poly = result[0];
    } catch (Exception e) {return "Invalid data entered";}
  if (type == 1) {     // Only computing expression
    info=new StringBuffer();
    w("<HTML><HEAD><TITLE>");
    w(english?"Polynomial calculator":"Calculadora de polinomios");
    w("</TITLE></HEAD><BODY><CENTER><H3><FONT COLOR=Red>");
    w(english?"Polynomial calculator":"Calculadora de polinomios");
    w("</FONT></H3></CENTER><P><I>");
    w(english?"by Dario Alejandro Alpern":"por Darío Alejandro Alpern");
    w("</I><P><DL><DT>");
    w(english?"Input expression":"Expresión de entrada");
    w(":<DD>");
    w(InputField);
    w("<P><DT>");
    w(english?"Value":"Valor");
    w(":<DD>");
    w(poly.toString());
    w(english?"<P>If you found any mistake or you have suggestions, please fill in the <A HREF=\"FORM.HTM?Polynomial factorization applet\" TARGET=\"_blank\">form</A>":"<P>Si encontró algún error o si tiene sugerencias, por favor llene el <A HREF=\"FORMULAR.HTM?Applet de factorizacion de polinomios\" TARGET=\"_blank\">formulario</A>");
    w(english?"<P><A HREF=\"GBOOK.HTM\">Click here</A> to view or sign my Guestbook.</BODY></HTML>":"<P><A HREF=\"EGBOOK.HTM\">Apriete aquí</A> para ver o firmar mi libro de invitados.<BODY></HTML>");
    return info.toString();
    }
  if (poly.getLeadingCoeff().mod(modulus).signum() == 0 &&
      modexponent > 1) {
    info=new StringBuffer();
    w("<HTML><HEAD><TITLE>");
    w(english?"Polynomial calculator":"Calculadora de polinomios");
    w("</TITLE></HEAD><BODY><CENTER><H3><FONT COLOR=Red>");
    w(english?"Polynomial calculator":"Calculadora de polinomios");
    w("</FONT></H3></CENTER><P><I>");
    w(english?"by Dario Alejandro Alpern":"por Darío Alejandro Alpern");
    w("</I><P>");
    w(poly.toString());
    w(english?"<P>Cannot factor because the leading coefficient is multiple of the modulus base.":
      "<P>No se puede factorizar porque el coeficiente m&aacute;s significativo es m&uacute;ltiplo de la base del m&oacute;dulo.");
    w(english?"<P>If you found any mistake or you have suggestions, please fill in the <A HREF=\"FORM.HTM?Polynomial factorization applet\" TARGET=\"_blank\">form</A>":"<P>Si encontró algún error o si tiene sugerencias, por favor llene el <A HREF=\"FORMULAR.HTM?Applet de factorizacion de polinomios\" TARGET=\"_blank\">formulario</A>");
    w(english?"<P><A HREF=\"GBOOK.HTM\">Click here</A> to view or sign my Guestbook.</BODY></HTML>":"<P><A HREF=\"EGBOOK.HTM\">Apriete aquí</A> para ver o firmar mi libro de invitados.<BODY></HTML>");
    return info.toString();
    }
  calcThread.start();
  return "";
  }

public void run() {
  int i, j, actualexp, nbrfactors;
  String polstr, exponstr;
  Polynomial polynomial, product, PolQ, factor;
  Polynomial PolU1, PolU3, PolV1, PolV3, PolQt, PolT;
  Polynomial [] PolPf, PolQf;
  BigInteger oldmodulus, newmodulus;     
  info=new StringBuffer();
  w("<TITLE>");
  w(english?"Polynomial factorization applet":"Applet de factorizacion de polinomios");
  w("</TITLE><CENTER><FONT COLOR=Red><B>");
  w(english?"Factors of ":"Factores de ");
  exponstr = (modexponent>1? (Polynomial.superscript?"<SUP>":"^")+modexponent+(Polynomial.superscript?"</SUP>":""):"");
  w(poly.toString()+" (mod "+modulus.toString()+exponstr+")</B></FONT></CENTER><P><I>");
  w(english?"by Dario":"por Darío");
  w(" Alejandro Alpern</I><P>");
  Date OldDate=new Date();
  long Old=OldDate.getTime();
  polynomial = (Polynomial)(poly.clone());
  poly.modulus = modulus;
  factor(poly);
  actualexp = 1;
  nbrfactors = answer.size();
  oldmodulus = modulus;
  if (modexponent == 1) {
    }
  else if (answer.size() == 1) {           // Irreducible polynomial.
    answer.setElementAt(polynomial, 0);
    }
  else {
    for (i=0; i<nbrfactors; i++) {
      if (expon[i] != 1) {
        w(english?"Cannot perform Hensel lift because the polynomial has duplicate factors modulo "+modulus.toString()+":":
                  "No se puede efectuar el algoritmo de Hensel porque el polinomio tiene factores duplicados m&oacute;dulo "+modulus.toString()+":");
        break;
        }
      }
    if (i == nbrfactors) {
      PolPf = new Polynomial[nbrfactors];
      PolQf = new Polynomial[nbrfactors];
      for (i=0; i<nbrfactors; i++) {
        polynomial.modulus = oldmodulus;
        PolV3 = PolPf[i] = (Polynomial)(answer.elementAt(i));
                 // Perform Extended Euclidean algorithm
                 // between factor and polynomial / factor.
        PolU1 = Polynomial.one(oldmodulus);        // PolU1 = 1.
        PolU3 = polynomial.divide(PolV3);
        PolV1 = new Polynomial(-1, oldmodulus);    // PolV1 = 0.
        while (PolV3.degree > 0) {
          PolQt = PolU3.divide(PolV3);
          PolT = PolU1.subtract(PolV1.multiply(PolQt));
          PolU1 = PolV1; PolV1 = PolT;
          PolT = PolU3.subtract(PolV3.multiply(PolQt));
          PolU3 = PolV3; PolV3 = PolT;
          }
        PolQf[i] = PolV1.mod(PolPf[i]).
                   multiply(PolV3.coeff[0].modInverse(modulus));
        }
      while (actualexp < modexponent) {
        newmodulus = oldmodulus.multiply(oldmodulus);
        product = Polynomial.one(newmodulus);
        for (i=0; i<nbrfactors; i++) {
          PolPf[i].modulus = newmodulus;
          product = product.multiply(PolPf[i]);
          }
        polynomial.modulus = newmodulus;
        PolQ = polynomial.subtract(product);
        polynomial.modulus = oldmodulus;
        for (i=PolQ.degree; i>=0; i--) {
          PolQ.coeff[i] = PolQ.coeff[i].divide(oldmodulus);
          }
        PolQ.modulus = oldmodulus;
        for (i=0; i<nbrfactors; i++) {
          PolPf[i].modulus = oldmodulus;
          PolV1 = PolQ.mod(PolPf[i]).multiply(PolQf[i]).mod(PolPf[i]);
          PolPf[i].modulus = newmodulus;
          PolV1.modulus = newmodulus;
          PolPf[i] = PolV1.multiply(oldmodulus).add(PolPf[i]);
          }
        actualexp *= 2;
        if (actualexp == modexponent) {
          break;
          }
        PolQ = Polynomial.one(newmodulus);
        for (i=0; i<nbrfactors; i++) {
          product = Polynomial.one(newmodulus);
          for (j=0; j<nbrfactors; j++) {
            if (i!=j) {
              PolPf[j].modulus = newmodulus;
              product = product.multiply(PolPf[j]);
              }
            else {
              PolQf[j].modulus = newmodulus;
              product = product.multiply(PolQf[j]);
              }
            }
          PolQ = PolQ.subtract(product);
          }
        for (j=PolQ.degree; j>=0; j--) {
          PolQ.coeff[j] = PolQ.coeff[j].divide(oldmodulus);
          }
        PolQ.modulus = oldmodulus;
        for (i=0; i<nbrfactors; i++) {
          factor = PolPf[i];
          factor.modulus = oldmodulus;
          PolQf[i].modulus = oldmodulus;
          PolV1 = PolQ.mod(factor).multiply(PolQf[i]).mod(factor);
          PolQf[i].modulus = newmodulus;
          PolV1.modulus = newmodulus;
          PolQf[i] = PolV1.multiply(oldmodulus).add(PolQf[i]);
          }
        oldmodulus = newmodulus;
        }              // End while
      newmodulus = modulus.pow(modexponent);
      for (i=0; i<nbrfactors; i++) {
        factor = PolPf[i];
        factor.modulus = newmodulus;
        for (j=factor.degree; j>=0; j--) {
          factor.coeff[j] = factor.coeff[j].mod(newmodulus);
          }
        answer.setElementAt(factor, i);
        }
      }              // End if
    }                // End if
  w("<UL>");
  for (i=0; i<nbrfactors; i++) {
    PolV1 = (Polynomial)(answer.elementAt(i));
    polstr = PolV1.toString();
    if (expon[i] == 1) {
      w("<LI>"+polstr.toString());
      }
    else {
      if (PolV1.degree == 1 &&
          PolV1.coeff[0].equals(BigInt0) &&
          PolV1.coeff[1].equals(BigInt1)) {
        if (Polynomial.superscript) {
          w("<LI>"+polstr.toString()+"<SUP>"+expon[i]+"</SUP>");
          }
        else {
          w("<LI>"+polstr.toString()+"^"+expon[i]);
          }
        }
      else {
        if (Polynomial.superscript) {
          w("<LI>("+polstr.toString()+")<SUP>"+expon[i]+"</SUP>");
          }
        else {
          w("<LI>("+polstr.toString()+")^"+expon[i]);
          }
        }
      }
    }
  w("</UL>");
                     // Validate factors
  product = Polynomial.one(((Polynomial)(answer.elementAt(nbrfactors-1))).modulus);
  polynomial.modulus = product.modulus;
  for (i=0; i<nbrfactors; i++) {
    for (j=expon[i]; j>0; j--) {
      product = product.multiply((Polynomial)(answer.elementAt(i)));
      }
    }
  product = product.subtract(polynomial);
  if (product.degree >= 0) {
    w(english?"The product of factors does not equal input polynomial. Please send this polynomial and the modulus to the author.":
              "El producto de los factores no es igual al polinomio entrado. Por favor env&iacute;e el polinomio y el m&oacute;dulo al autor.");
    }
  Date NewDate=new Date();
  long New=NewDate.getTime();
  w(english?"<P>Calculation time: ":"<P>Tiempo de cálculo: ");
  int t=(int)(((New-Old)/1000+86400)%86400);
  w(t/3600+"h "+((t%3600)/60)+"m "+(t%60)+"s");
  w(english?"<P>If you found any mistake or you have suggestions, please fill in the <A HREF=\"FORM.HTM?Polynomial factorization applet\" TARGET=\"_blank\">form</A>":"<P>Si encontró algún error o si tiene sugerencias, por favor llene el <A HREF=\"FORMULAR.HTM?Applet de factorizacion de polinomios\" TARGET=\"_blank\">formulario</A>");
  w(english?"<P><A HREF=\"GBOOK.HTM\">Click here</A> to view or sign my Guestbook.</BODY></HTML>":"<P><A HREF=\"EGBOOK.HTM\">Apriete aquí</A> para ver o firmar mi libro de invitados.<BODY></HTML>");
  calcThread = null;
  }

public String resultCalc() {
  if (calcThread == null) {
    return info.toString();
    }
  return "";
  }
}
