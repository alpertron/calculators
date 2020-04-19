// <XMP>
// Decomposition of a number in a sum of up to four squares
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated December 20th, 2002.
// See http://www.alpertron.com.ar/FSQUARES.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.applet.*;
import java.util.*;
import java.awt.*;
import java.math.*;

public final class fsquares extends Applet implements Runnable {

  private static final int maxsieve = 65536;
  private TextField textNumber;
  private TextArea lowerTextArea;
  private static final BigInteger BigInt0 = BigInteger.valueOf(0L);
  private static final BigInteger BigInt1 = BigInteger.valueOf(1L);
  private static final BigInteger BigInt2 = BigInteger.valueOf(2L);
  private static final BigInteger BigInt3 = BigInteger.valueOf(3L);
  private boolean TerminateThread;
  private Thread calcThread;
  private String StringToLabel;
  private String textAreaContents;
  private int digitsInGroup = 6;
  private int [] sieve = new int [maxsieve/2];
  private int [] primediv = new int [maxsieve/4];
  private int [] primeexp = new int [maxsieve/4];
  private long [] TestNbr = new long [1200];
  private BigInteger mult1, mult2, mult3, mult4;
  private BigInteger Nbr = BigInt0;
  private String paneFooter;

  public static void main(String args[])
  {
    BigInteger nbr = new BigInteger(args[0]);
    System.out.println(pollardBrentRho(nbr,100000));
  }

  public void init() {
    Label labelTop, labelBottom;
    int i,j;

    // Determine primes (for primes sieve[i] >= 0)

    for (i=0; i<maxsieve/2; i++)
    {
      sieve[i] = 0;
    }
    for (i=3; i*i<maxsieve/2; i+=2)
    {
      j = i*i-3;
      j = (j%2 == 0)? j/2: (j+i)/2;
      for (; j<maxsieve/2; j+=i)
      {
        sieve[j] = -1;             // Indicate number is composite.
      }
    }
    setBackground(Color.lightGray);
    setLayout(null);
    labelTop = new Label("Type number or numerical expression here and press Return");
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
    labelBottom = new Label("Written by Dario Alejandro Alpern. Last updated December 20th, 2002");
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
    BigInteger nbr, modulus, Mult1, Mult2, Mult3, Mult4, p, q, K, M1, M2, Tmp;
    BigInteger [] ExpressionResult = new BigInteger[1];
    int power4, r, s, shRight, ExpressionRC;
    int modularexp;
    boolean Computing3Squares;
    long StartTime, remainder;
    BigInteger Divisor;
    int i, j, sum;
    int nbrDivisors = 0;
    int nbrModExp = 0;
    int iMult3, iMult4;
    int NumberLength;
    byte Result[];
    long L, P, Q;
    long [] TestNbr = this.TestNbr;

    StartTime = System.currentTimeMillis();
    TerminateThread = false;

    lowerTextArea.setText("Computing input expression...");
    try {
      ExpressionRC = expression.ComputeExpression(textNumber.getText().trim(),0,ExpressionResult);
    } catch (OutOfMemoryError e) {
      lowerTextArea.setText("Out of memory.");
      return;
    } catch (ArithmeticException e) {
      return;
    }
    modulus = nbr = ExpressionResult[0];
    if (ExpressionRC != 0) {
      switch (ExpressionRC) {
        case -1:
          lowerTextArea.setText("Number too low (less than 2).");
          break;
        case -2:
          lowerTextArea.setText("Number too high (more than 1000 digits).");
          break;
        case -3:
          lowerTextArea.setText("Intermediate expression too high (more than 2000 digits).");
          break;
        case -4:
          lowerTextArea.setText("Non-integer division.");
          break;
        case -5:
          lowerTextArea.setText("Parentheses mismatch.");
          break;
        case -6:
          lowerTextArea.setText("Syntax error.");
          break;
        case -7:
          lowerTextArea.setText("Too many parentheses.");
          break;
        case -8:
          lowerTextArea.setText("Invalid parameter.");
          break;
        }
      return;
      }
    // If modulus is a perfect square, show the square root.
    Mult1 = squareRoot(modulus);
    if (Mult1.multiply(Mult1).equals(modulus))
    {
      Mult2 = BigInt0;
      iMult3 = iMult4 = 0;
      power4 = 0;
    }
    else
    {
      power4 = modulus.getLowestSetBit() / 2;
      modulus = modulus.shiftRight(2*power4);    // modulus is not multiple of 4.
      if (modulus.compareTo(BigInt3) <= 0)
      {
        Mult1 = Mult2 = BigInt0;
        iMult3 = iMult4 = 0;
        switch (modulus.intValue())
        {
          case 3:
            iMult3 = 1;
          case 2:
            Mult2 = BigInt1;
          case 1:
            Mult1 = BigInt1;
        }
      }
      else
      {

        // Store in TestNbr the value of modulus.
        Result = modulus.toByteArray();
        NumberLength = (Result.length+3)/4;
        j = 0;
        L = 1;
        P = 0;
        for (i=Result.length-1; i>=0; i--)
        {
          P += L*(long)(Result[i]>=0?Result[i]:Result[i]+256);
          L *= 0x100;
          if (L == 0x100000000L)
          {
            TestNbr[j] = P;
            j++;
            L = 1;
            P = 0;
          }
        }
        TestNbr[j] = P;
      
        r = (int)TestNbr[0] & 7;                    // modulus mod 8
    
         // Fill sieve array with r%n

        for (i=0; i<maxsieve/2; i++)
        {
          if (sieve[i] >= 0)                        // If prime...
          {
            Q = 2*i+3;                              // Prime
            long Divid, Rem = 0;

            for (j=NumberLength-1; j>=0; j--)
            {
              Divid = TestNbr[j] + (Rem << 32);
              Rem = Divid % Q;
            }
            sieve[i] = (int)Rem;
          }
        }

        if (r != 7)   // n!=7 (mod 8) => Sum of three squares
        {
          lowerTextArea.setText("Computing sum of three squares...");
          Computing3Squares = true;
          iMult4 = 0;
          iMult3 = -1;
        }
        else
        {
          iMult3 = -1;
          iMult4 = 0;
          lowerTextArea.setText("Computing sum of four squares...");
          Computing3Squares = false;
        }
compute_squares_loop:
        while (true)
        {
          if (TerminateThread)
          {
            return;                           
          }
          if (Computing3Squares == false && iMult3 >= iMult4)
          {
            iMult3 = 1;
            iMult4++;
          }
          else
          {
            iMult3++;
          }
          sum = iMult3 * iMult3 + iMult4 * iMult4;
          p = modulus.subtract(BigInteger.valueOf(sum));
          shRight = p.getLowestSetBit();
          p = p.shiftRight(shRight);
          if (p.testBit(0) == false || p.testBit(1) == true)
          {
            continue;         /* p is not congruent to 1 mod 4 */
          }
          
          // Now compute modulus over the product of all odd divisors
          // less than maxsieve.

          nbrDivisors = 0;
          for (i = 0; i<maxsieve/2; i++)
          {
            if (sieve[i] >= 0 && (sieve[i] - sum) % (i*2+3) == 0)
            {                          // Divisor found.
              primediv[nbrDivisors] = i*2+3;   // Store divisor.
              primeexp[nbrDivisors] = 0;       // Store exponent.
              Divisor = BigInteger.valueOf(i*2+3);
              while (p.remainder(Divisor).signum() == 0)
              {
                p = p.divide(Divisor);
                primeexp[nbrDivisors]++;       // Increment exponent.
              }
              if ((i*2+3)%4 == 3 && primeexp[nbrDivisors]%2 != 0)
              {                   // Divisor of the form 4k+3 not square.
                continue compute_squares_loop;  // Discard this value.
              }
              nbrDivisors++;
            }
          }

          if (p.equals(BigInt1))
          {
            Mult1 = BigInt1;
            Mult2 = BigInt0;
          }
          else
          {
            // At this moment p must be prime, otherwise we must discard it.

            q = p.subtract(BigInt1);              // q = p - 1.
            K = BigInt1;
            j = q.getLowestSetBit();
            Mult3 = q.shiftRight(j);
            Mult4 = BigInt0;
            do
            {                 // Compute Mult1 = sqrt(-1) (mod p).
              K = K.add(BigInt1);
              Mult1 = K.modPow(Mult3,p);
              nbrModExp++;
              for (i=0; i<j; i++)
              {
                Mult4 = Mult1.multiply(Mult1).mod(p);
                if (Mult4.equals(q))
                {
                  break;      // Mult1^2 = -1 (mod p), so exit loop.
                }
                Mult1 = Mult4;
              }
            } while (Mult4.equals(BigInt1));  // Loop if number is at least PRP.
            if (Mult4.equals(q) == false)
            {            // Cannot find sqrt(-1) (mod p), go to next candidate.
              continue; 
            }
            Mult2 = BigInt1;
            while (true)
            {
              K = Mult1.multiply(Mult1).add(Mult2.multiply(Mult2)).divide(p);
              if (K.equals(BigInt1))
              {                                   // Calculation finished
                break;
              }
              if (p.mod(K).signum() == 0)
              {                                   // The number is not prime
                continue; 
              }
              M1 = Mult1.mod(K);
              M2 = Mult2.mod(K);
              if (M1.compareTo(K.shiftRight(1)) > 0) {M1 = M1.subtract(K);}
              if (M2.compareTo(K.shiftRight(1)) > 0) {M2 = M2.subtract(K);}
              Tmp = Mult1.multiply(M1).add(Mult2.multiply(M2)).divide(K);
              Mult2 = Mult1.multiply(M2).subtract(Mult2.multiply(M1)).divide(K);
              Mult1 = Tmp;
            }      /* end while */
          }
    
          // Use the other divisors of modulus in order to get Mult1 and Mult2.

          for (i=0; i<nbrDivisors; i++)
          {
            if (primeexp[i] >= 2)
            {
              Tmp = BigInteger.valueOf(primediv[i]).pow(primeexp[i]/2);
              Mult1 = Mult1.multiply(Tmp);
              Mult2 = Mult2.multiply(Tmp);
            }
            if (primeexp[i]%2 == 1)
            {
              // Since the value primediv[i] is very low it is faster to use
              // trial and error than the general method in order to find
              // the sum of two squares.

              s = primediv[i];
              j = 1;
              while (true)
              {
                r = (int)Math.sqrt(s-j*j);
                if (r*r+j*j == s)
                {
                  break;
                }
                j++;
              }

              M1 = BigInteger.valueOf(j);
              M2 = BigInteger.valueOf(r);
              Tmp = Mult1.multiply(M1).add(Mult2.multiply(M2));
              Mult2 = Mult1.multiply(M2).subtract(Mult2.multiply(M1));
              Mult1 = Tmp;
            }
          }            /* end for */
          break;
        }              /* end while */

        Mult1 = Mult1.shiftLeft(shRight/2);
        Mult2 = Mult2.shiftLeft(shRight/2);
        if (shRight%2 == 1)
        {
          Tmp = Mult1.add(Mult2);
          Mult2 = Mult1.subtract(Mult2).abs();
          Mult1 = Tmp;
        }
      }
    }
    Mult1 = Mult1.shiftLeft(power4).abs();
    Mult2 = Mult2.shiftLeft(power4).abs();
    Mult3 = BigInteger.valueOf(iMult3).shiftLeft(power4);
    Mult4 = BigInteger.valueOf(iMult4).shiftLeft(power4);
    // Sort squares
    if (Mult1.compareTo(Mult2)<0) {Tmp=Mult1; Mult1=Mult2; Mult2=Tmp;}
    if (Mult1.compareTo(Mult3)<0) {Tmp=Mult1; Mult1=Mult3; Mult3=Tmp;}
    if (Mult2.compareTo(Mult3)<0) {Tmp=Mult2; Mult2=Mult3; Mult3=Tmp;}
    if (Mult1.compareTo(Mult4)<0) {Tmp=Mult1; Mult1=Mult4; Mult4=Tmp;}
    if (Mult2.compareTo(Mult4)<0) {Tmp=Mult2; Mult2=Mult4; Mult4=Tmp;}
    if (Mult3.compareTo(Mult4)<0) {Tmp=Mult3; Mult3=Mult4; Mult4=Tmp;}

    // Validate

    if (Mult1.multiply(Mult1).add(Mult2.multiply(Mult2)).
        add(Mult3.multiply(Mult3)).add(Mult4.multiply(Mult4)).equals(nbr)
        == false)
    {
      lowerTextArea.setText("Internal error!\n\nPlease send the number to the author of the applet");
      Nbr = BigInt0;
      return;
    }

    mult1 = Mult1;
    mult2 = Mult2;
    mult3 = Mult3;
    mult4 = Mult4;
    Nbr = nbr;
    long t=(System.currentTimeMillis()-StartTime)/1000;
    String timeElapsed = "Time elapsed: "+t/86400+"d "+(t%86400)/3600+"h "+((t%3600)/60)+"m "+(t%60)+"s";
    paneFooter = "\n\n"+timeElapsed+
                 "\n\nNumber of modular exponentiations: "+nbrModExp;
    ShowPane();
  }

  private void ShowPane()
  {
    textAreaContents = "";
    StringToLabel = "n = ";
    insertBigNbr(Nbr);
    textAreaContents += StringToLabel + "\n\n";

    if (mult4.signum() == 0)
    {
      if (mult3.signum() == 0)
      {
        if (mult2.signum() == 0)
        {
          textAreaContents += "n = a^2\n";
          StringToLabel = "a = ";            // Start new line.
          insertBigNbr(mult1);
        }
        else
        {
          textAreaContents += "n = a^2 + b^2\n";
          StringToLabel = "a = ";            // Start new line.
          insertBigNbr(mult1);
          textAreaContents += StringToLabel + "\n";
          StringToLabel = "b = ";            // Start new line.
          insertBigNbr(mult2);
        }
      }
      else
      {
        textAreaContents += "n = a^2 + b^2 + c^2\n";
        StringToLabel = "a = ";            // Start new line.
        insertBigNbr(mult1);
        textAreaContents += StringToLabel + "\n";
        StringToLabel = "b = ";            // Start new line.
        insertBigNbr(mult2);
        textAreaContents += StringToLabel + "\n";
        StringToLabel = "c = ";            // Start new line.
        insertBigNbr(mult3);
      }
    }
    else
    {
      textAreaContents += "n = a^2 + b^2 + c^2 + d^2\n";
      StringToLabel = "a = ";            // Start new line.
      insertBigNbr(mult1);
      textAreaContents += StringToLabel + "\n";
      StringToLabel = "b = ";            // Start new line.
      insertBigNbr(mult2);
      textAreaContents += StringToLabel + "\n";
      StringToLabel = "c = ";            // Start new line.
      insertBigNbr(mult3);
      textAreaContents += StringToLabel + "\n";
      StringToLabel = "d = ";            // Start new line.
      insertBigNbr(mult4);
    }
    lowerTextArea.setText(textAreaContents + StringToLabel + paneFooter);
  }

  private void addStringToLabel(String value) {
    if (value.length() + 1 + StringToLabel.length() >= 79) {
      textAreaContents += StringToLabel + "\n";
      StringToLabel = value+" ";
    }
    else {
      StringToLabel += value+" ";
    }
  }

  private void insertBigNbr(BigInteger N)
  {
    int i,dig;
    dig = digitsInGroup & 0x3FF;
    String value = N.toString();
    i = (value.length() + dig - 1)%dig + 1;
    addStringToLabel(value.substring(0,i));
    while (i<value.length())
    {
      addStringToLabel(value.substring(i,i+dig));
      i += dig;
    }
  }

  // Compute square root using Newton's algorithm. 
  private BigInteger squareRoot(BigInteger value)
  {
    BigInteger sqrt = BigInt0;
    BigInteger newsqrt = BigInt1.shiftLeft((value.bitLength()+1)/2);
    do
    {
      sqrt = newsqrt;
      newsqrt = sqrt.add(value.divide(sqrt)).shiftRight(1);
    } while (sqrt.compareTo(newsqrt) > 0);
    return sqrt;
  }

  private static BigInteger pollardBrentRho(BigInteger value, int nbrLoops)
  {
    int r = 1;
    int i, k;
    int m = 20;
    BigInteger x;
    BigInteger ys = BigInt0;
    BigInteger G = value;
    BigInteger y = BigInt1;
    BigInteger q = BigInt1;
    BigInteger BigInt3 = BigInteger.valueOf(3L);
    do
    {
      x = y;
      for (i=1; i<=r; i++)
      {
        y = y.multiply(y).add(BigInt3).mod(value);
      }
      k = 0;
      do
      {
        ys = y;
        i = (m<r-k? m:r-k);
        for (; i>=1; i--)
        {
          y = y.multiply(y).add(BigInt3).mod(value);
          q = q.multiply(x.subtract(y)).mod(value);
        }
        G = q.gcd(value);
        k += m;
        if (k > nbrLoops)
        {
          return value;      /* Factor not found */
        }
      } while (k<r && G.compareTo(BigInt1)<=0);
      r *= 2;
    } while (G.compareTo(BigInt1)<=0);
    if (G.equals(value))
    {
      do
      {
        ys = ys.multiply(ys).add(BigInt3).mod(value);
        G = x.subtract(ys).abs().gcd(value);
      } while (G.compareTo(BigInt1)<=0);
    }
    return G;
  }

  public String getTextFromPane()
  {
    if ((calcThread == null || calcThread.isAlive() == false) &&
         lowerTextArea.getText() != "")
    {
      String text="<HTML><TITLE>Sum of squares</TITLE><BODY><H2>Sum of squares</H2><P><I>Written by Dario Alejandro Alpern (alpertron@hotmail.com)</I><P><PRE>"
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
