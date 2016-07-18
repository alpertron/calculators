// <XMP>
// Decomposition of a number in a sum of four cubes
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated March 27th, 2005.
// See http://www.alpertron.com.ar/FCUBES.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.applet.*;
import java.awt.*;
import java.math.*;

public final class sumcubos extends Applet implements Runnable {

  private TextField textNumber;
  private TextArea lowerTextArea;
  private boolean TerminateThread;
  private Thread calcThread;
  private String StringToLabel;
  private StringBuffer textAreaContents;
  private int digitsInGroup = 6;
  private String paneFooter;
  private BigInteger Nbr, mult1, mult2, mult3, mult4;
  private int[] sums =
  {
    6, 0, 1, -1, -1, 0, -1, 0, 1, 1,
    6, 3, 1, 0, -1, 4, 2, -5, -2, 4,
    18, 1, 2, 14, -2, -23, -3, -26, 3, 30,
    18, 7, 1, 2, 6, -1, 8, -2, -9, 2,
    18, 8, 1, -5, -1, 14, -3, 29, 3, -30,
    18, 10, 1, 6, -1, -15, -3, -32, 3, 33,
    18, 11, 1, -1, 6, 7, 8, 10, -9, -11,
    18, 17, 2, -12, -2, 21, -3, 23, 3, -27,
    54, 20, 3, -11, -3, 10, 1, 2, -1, 7,
    72, 56, -9, 4, 1, 4, 6, -2, 8, -4,
    108, 2, -1, -22, 1, 4, -3, -41, 3, 43,
    216, 92, 3, -164, -3, 160, 1, -35, -1, 71,
    270, 146, -60, 91, -3, 13, 22, -37, 59, -89,
    270, 200, 3, 259, -3, -254, 1, 62, -1, -107,
    270, 218, -3, -56, 3, 31, -5, -69, 5, 78,
    432, 380, -3, 64, 3, -80, 2, -29, -2, 65,
    540, 38, 5, -285, -5, 267, 3, -140, -3, 190,
    810, 56, 5, -755, -5, 836, 9, -1445, -9, 1420,
    1080, 380, -1, -1438, 1, 1258, -3, -4037, 3, 4057,
    1620, 1334, -5, -3269, 5, 3107, -9, -5714, 9, 5764,
    1620, 1352, -5, 434, 5, -353, 9, -722, -9, 697,
    2160, 362, -5, -180, 5, 108, -6, -149, 6, 199,
    6480, 794, -5, -83, 5, 11, -6, -35, 6, 85, 
  };

  public void init() {
    Label labelTop, labelBottom;
    setBackground(Color.lightGray);
    setLayout(null);
    labelTop = new Label("Digite un número o expresión numérica y presione Enter");
    labelTop.reshape(10, 10, 570, 14);
    labelTop.setFont(new Font("Courier", Font.PLAIN, 12));
    labelTop.setAlignment(Label.CENTER);
    add(labelTop);
    textNumber = new TextField(1000);
    textNumber.reshape(10, 30, 570, 30);
    textNumber.setEditable(true);
    add(textNumber);
    lowerTextArea = new TextArea("", 6,75, TextArea.SCROLLBARS_VERTICAL_ONLY);
    lowerTextArea.reshape(10, 70, 570, 200);
    lowerTextArea.setEditable(false); 
    lowerTextArea.setFont(new Font("Courier", Font.PLAIN, 12));
    add(lowerTextArea);
    labelBottom = new Label("Escrito por Dario Alejandro Alpern. Modificado el 27 de marzo de 2005");
    labelBottom.reshape(10, 280, 570, 14);
    labelBottom.setFont(new Font("Courier", Font.PLAIN, 12));
    labelBottom.setAlignment(Label.CENTER);
    add(labelBottom);
    textNumber.requestFocus();
    validate();
    }

  public void destroy() {           /* Applet end */
    TerminateThread = true;
  }

public boolean handleEvent(Event e) {
  if (e.id == Event.KEY_PRESS && e.key == Event.ENTER) {
    if (calcThread != null) {
      TerminateThread = true;
      try {
        calcThread.join();        /* Wait until the solving thread dies */
        } catch (InterruptedException ie) {};
      }
    calcThread = new Thread(this);  /* Start solving thread */
    calcThread.start();
    return true;
    }
  return false;
  }

  public void run() {
    BigInteger nbr, k, x, Mult1, Mult2, Mult3, Mult4;
    BigInteger mod = BigInteger.valueOf(0);
    BigInteger [] ExpressionResult = new BigInteger[1];
    int ExpressionRC;
    long StartTime;
    int i, j, mask;
    boolean converted = false;

    StartTime = System.currentTimeMillis();
    TerminateThread = false;
    mult1 = mult2 = mult3 = mult4 = BigInteger.valueOf(0);
    lowerTextArea.setText("Calculando expresión de entrada...");
    try {
      ExpressionRC = expression.ComputeExpression(textNumber.getText().trim(),1,ExpressionResult);
    } catch (OutOfMemoryError e) {
      lowerTextArea.setText("Sin memoria.");
      return;
    } catch (ArithmeticException e) {
      return;
    }
    Nbr = nbr = ExpressionResult[0];
    if (ExpressionRC != 0) {
      switch (ExpressionRC) {
        case -2:
          lowerTextArea.setText("Número muy grande (más de 1000 dígitos).");
          break;
        case -3:
          lowerTextArea.setText("Expresión intermedia con más de 2000 dígitos.");
          break;
        case -4:
          lowerTextArea.setText("División no entera.");
          break;
        case -5:
          lowerTextArea.setText("Error de paréntesis.");
          break;
        case -6:
          lowerTextArea.setText("Error de sintaxis.");
          break;
        case -7:
          lowerTextArea.setText("Demasiados paréntesis.");
          break;
        case -8:
          lowerTextArea.setText("Parámetro inválido.");
          break;
        }
      return;
      }
    int mod18 = nbr.mod(BigInteger.valueOf(18)).intValue();
    if (mod18 == 4 || mod18 == 5 || mod18 == 13 || mod18 == 14) {
      lowerTextArea.setText("Este applet no funciona si el número es\ncongruente a 4 ó 5 (mod 9).");
      return;
      }
    if (mod18 == 16) {
      converted = true;
      nbr = nbr.negate();
      }
    for (i=0; i<sums.length; i+=10) {
      mod = BigInteger.valueOf(sums[i]);
      if (nbr.mod(mod).intValue() == sums[i+1]) {
        break;
        }
      }
    if (i<sums.length) {
      x = nbr.subtract(BigInteger.valueOf(sums[i+1])).divide(mod);
      Mult1 = BigInteger.valueOf(sums[i+2]).multiply(x).
              add(BigInteger.valueOf(sums[i+3]));
      Mult2 = BigInteger.valueOf(sums[i+4]).multiply(x).
              add(BigInteger.valueOf(sums[i+5]));
      Mult3 = BigInteger.valueOf(sums[i+6]).multiply(x).
              add(BigInteger.valueOf(sums[i+7]));
      Mult4 = BigInteger.valueOf(sums[i+8]).multiply(x).
              add(BigInteger.valueOf(sums[i+9]));
      }
    else {
      mod = BigInteger.valueOf(54);
      if (nbr.mod(mod).intValue() == 2) {
        x = nbr.subtract(BigInteger.valueOf(2)).divide(mod);
        Mult1 = BigInteger.valueOf(29484).multiply(x).
                add(BigInteger.valueOf(2211)).multiply(x).
                add(BigInteger.valueOf(43));
        Mult2 = BigInteger.valueOf(-29484).multiply(x).
                add(BigInteger.valueOf(-2157)).multiply(x).
                add(BigInteger.valueOf(-41));
        Mult3 = BigInteger.valueOf(9828).multiply(x).
                add(BigInteger.valueOf(485)).multiply(x).
                add(BigInteger.valueOf(4));
        Mult4 = BigInteger.valueOf(-9828).multiply(x).
                add(BigInteger.valueOf(-971)).multiply(x).
                add(BigInteger.valueOf(-22));
        }
      else {
        mod = BigInteger.valueOf(83*108);
        if (nbr.mod(mod).intValue() == 83*46) {
          x = nbr.subtract(BigInteger.valueOf(83*46)).divide(mod);
          Mult1 = BigInteger.valueOf(29484).multiply(x).
                  add(BigInteger.valueOf(25143)).multiply(x).
                  add(BigInteger.valueOf(5371));
          Mult2 = BigInteger.valueOf(-29484).multiply(x).
                  add(BigInteger.valueOf(-25089)).multiply(x).
                  add(BigInteger.valueOf(-5348));
          Mult3 = BigInteger.valueOf(9828).multiply(x).
                  add(BigInteger.valueOf(8129)).multiply(x).
                  add(BigInteger.valueOf(1682));
          Mult4 = BigInteger.valueOf(-9828).multiply(x).
                  add(BigInteger.valueOf(-8615)).multiply(x).
                  add(BigInteger.valueOf(-1889));
          }
        else {
          lowerTextArea.setText("Calculando soluciones...");
          // Perform Pell solution of Demjanenko's theorem
            // Using these values of P, Q, R and S, a and b will be
            // always one and zero (mod 6) respectively.
          BigInteger P = new BigInteger("-112488782561");
          BigInteger Q = new BigInteger("-6578430178320");
          BigInteger R = new BigInteger("-1923517596");
          BigInteger S = P;
          BigInteger tmpA;
          int mod83 = nbr.mod(BigInteger.valueOf(83)).intValue();
          int pow = 71;
          int exp = 0;
          while (pow != mod83) {
            exp++;
            pow = (pow*50)%83;
            }
          if (exp > 82/2) {
            exp = 82 - exp;
            Q = Q.negate();
            R = R.negate();
            }
          BigInteger P1 = BigInteger.valueOf(1);
          BigInteger Q1 = BigInteger.valueOf(0);
          BigInteger R1 = Q1;
          BigInteger S1 = P1;
          mask = 32;
          while (mask > 0) {
            BigInteger tmpP1 = P1.pow(2).add(Q1.multiply(R1));
            BigInteger tmpQ1 = P1.add(S1).multiply(Q1);
            BigInteger tmpR1 = P1.add(S1).multiply(R1);
            BigInteger tmpS1 = S1.pow(2).add(Q1.multiply(R1));
            P1 = tmpP1;
            Q1 = tmpQ1;
            R1 = tmpR1;
            S1 = tmpS1;
            if ((exp & mask) != 0) {
              tmpP1 = P.multiply(P1).add(Q.multiply(R1));
              tmpQ1 = P.multiply(Q1).add(Q.multiply(S1));
              tmpR1 = R.multiply(P1).add(S.multiply(R1));
              tmpS1 = R.multiply(Q1).add(S.multiply(S1));
              P1 = tmpP1;
              Q1 = tmpQ1;
              R1 = tmpR1;
              S1 = tmpS1;
              }
            mask /= 2;
            }
          BigInteger a = P1.multiply(BigInteger.valueOf(-3041)).
                         add(Q1.multiply(BigInteger.valueOf(-52)));
          BigInteger b = R1.multiply(BigInteger.valueOf(-3041)).
                         add(S1.multiply(BigInteger.valueOf(-52)));

          Mult1 = a.multiply(BigInteger.valueOf(27)).
                  add(b.multiply(BigInteger.valueOf(-928)));
          Mult2 = a.multiply(BigInteger.valueOf(-9)).
                  add(b.multiply(BigInteger.valueOf(-602)));
          Mult3 = a.multiply(BigInteger.valueOf(25)).
                  add(b.multiply(BigInteger.valueOf(-2937)));
          Mult4 = a.multiply(BigInteger.valueOf(-19)).
                  add(b.multiply(BigInteger.valueOf(2746)));
          k = nbr.subtract(Mult1.pow(3)).subtract(Mult2.pow(3)).
                  subtract(Mult3.pow(3)).subtract(Mult4.pow(3)).
                  divide(BigInteger.valueOf(18*83));
          Mult1 = Mult1.add(k.multiply(BigInteger.valueOf(10)));
          Mult2 = Mult2.add(k.multiply(BigInteger.valueOf(-19)));
          Mult3 = Mult3.add(k.multiply(BigInteger.valueOf(-24)));
          Mult4 = Mult4.add(k.multiply(BigInteger.valueOf(27)));
          }
        }
      }

    if (converted) {
      Mult1 = Mult1.negate();
      Mult2 = Mult2.negate();
      Mult3 = Mult3.negate();
      Mult4 = Mult4.negate();
      nbr = nbr.negate();
      }
    // Validate

    lowerTextArea.setText("Validando...");
    if (Mult1.pow(3).add(Mult2.pow(3)).add(Mult3.pow(3)).add(Mult4.pow(3)).
        equals(nbr) == false)
    {
      lowerTextArea.setText("¡Error interno!\n\nPor favor envíe el número al autor del applet");
      return;
    }

    lowerTextArea.setText("Formateando la salida...");
    mult1 = Mult1;
    mult2 = Mult2;
    mult3 = Mult3;
    mult4 = Mult4;
    Nbr = nbr;
    long t=(System.currentTimeMillis()-StartTime)/1000;
    String timeElapsed = "Tiempo transcurrido: "+t/86400+"d "+(t%86400)/3600+"h "+((t%3600)/60)+"m "+(t%60)+"s";
    paneFooter = "\n\n"+timeElapsed;
    ShowPane();
  }

  private void ShowPane()
  {
    textAreaContents = new StringBuffer();
    StringToLabel = "n = ";
    insertBigNbr(Nbr);
    textAreaContents = textAreaContents.append(StringToLabel + "\n\n").
                       append("n = a^3 + b^3 + c^3 + d^3\n");
    StringToLabel = "a = ";            // Start new line.
    insertBigNbr(mult1);
    textAreaContents = textAreaContents.append(StringToLabel + "\n");
    StringToLabel = "b = ";            // Start new line.
    insertBigNbr(mult2);
    textAreaContents = textAreaContents.append(StringToLabel + "\n");
    StringToLabel = "c = ";            // Start new line.
    insertBigNbr(mult3);
    textAreaContents = textAreaContents.append(StringToLabel + "\n");
    StringToLabel = "d = ";            // Start new line.
    insertBigNbr(mult4);
    textAreaContents = textAreaContents.append(StringToLabel + paneFooter);
    lowerTextArea.setText(textAreaContents.toString());
  }

  private void addStringToLabel(String value) {
    if (value.length() + 1 + StringToLabel.length() >= 79) {
      textAreaContents = textAreaContents.append(StringToLabel + "\n");
      StringToLabel = value+" ";
    }
    else {
      StringToLabel += value+" ";
    }
  }

  private void insertBigNbr(BigInteger N)
  {
    int i,dig,len;
    dig = digitsInGroup & 0x3FF;
    String value = N.toString();
    len = value.length();
    i = (len + dig - 1)%dig + 1;
    addStringToLabel(value.substring(0,i));
    while (i<value.length())
    {
      addStringToLabel(value.substring(i,i+dig));
      i += dig;
    }
    if (value.charAt(0) == '-')
    {
      len--;
    }
    if (len > 100)
    {
      addStringToLabel("("+len+" dígitos)");
    }
  }

  public String getTextFromPane()
  {
    if ((calcThread == null || calcThread.isAlive() == false) &&
         lowerTextArea.getText() != "")
    {
      String text="<HTML><TITLE>Suma de cuatro cubos</TITLE><BODY><H2>Suma de cuatro cubos</H2><P><I>Escrito por Dario Alejandro Alpern (alpertron@hotmail.com)</I><P><PRE>"
               + lowerTextArea.getText() + "</PRE><BODY></HTML>";
      return text;
    }
    return "";
  }

  public void setDigits(int digits)
  {
    digitsInGroup = digits;
    if ((calcThread == null || calcThread.isAlive() == false) &&
         lowerTextArea.getText() != "")
    {
      ShowPane();
    }
  }
}       // End applet class

// </XMP>
