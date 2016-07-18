// <XMP>
// Discrete logarithm calculator
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated April 15th, 2005, See http://www.alpertron.com.ar/DILOG.HTM
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
// 
import java.applet.*;
import java.util.*;
import java.awt.*;
import java.math.*;

public final class dilog extends ecm {

  private TextField textBase, textPower, textMod, textExp, textPeriod;
  private Button btnSolve;
  private BigInteger base, power, modulus;
  private BigInteger LastModulus = BigInt0;
  private Frame factorWindow;
  private ecm factorPanel;
  private BigInteger PrimesMod[];    // Mod = Modulus.
  private int ExponentsMod[];
  private BigInteger PrimesModBak[];
  private int ExponentsModBak[];
  private BigInteger PrimesGO[];     // GO = Group Order.
  private int ExponentsGO[];
  private BigInteger PrimesGOBak[];
  private int ExponentsGOBak[];
  private BigInteger nbrV[];
  private int NumberLengthOther;
  private int NbrFactorsMod=0;
  private long TestNbrOther[] = new long[NLen];
  private long nbrA[] = new long[NLen];
  private long nbrA2[] = new long[NLen];
  private long nbrB[] = new long[NLen];
  private long nbrB2[] = new long[NLen];
  private long nbrK[] = new long[NLen];
  private long nbrR[] = new long[NLen];
  private long nbrROther[] = new long[NLen];
  private long nbrR2[] = new long[NLen];
  private long nbrPower[] = new long[NLen];
  private long nbrBase[] = new long[NLen];
  private double dNOther;

  public void init() {
    Label labelBase, labelPower, labelMod, labelExp, labelPeriod, labelBottom;

    PrimesMod=new BigInteger[400];
    ExponentsMod=new int[400];
    PrimesModBak=new BigInteger[400];
    ExponentsModBak=new int[400];
    PrimesGO=new BigInteger[400];
    ExponentsGO=new int[400];
    nbrV = new BigInteger[400];
    setLayout(null);
    labelBase = new Label("Base");
    labelBase.reshape(10, 15, 40, 14);
    labelBase.setFont(new Font("Courier", Font.PLAIN, 12));
    labelBase.setAlignment(Label.CENTER);
    add(labelBase);
    textBase = new TextField(64);
    textBase.reshape(60, 10, 520, 30);
    textBase.setEditable(true);
    add(textBase);
    labelPower = new Label("Power");
    labelPower.reshape(10, 65, 40, 14);
    labelPower.setFont(new Font("Courier", Font.PLAIN, 12));
    labelPower.setAlignment(Label.CENTER);
    add(labelPower);
    textPower = new TextField(64);
    textPower.reshape(60, 60, 520, 30);
    textPower.setEditable(true);
    add(textPower);
    labelMod = new Label("Mod");
    labelMod.reshape(10, 115, 40, 14);
    labelMod.setFont(new Font("Courier", Font.PLAIN, 12));
    labelMod.setAlignment(Label.CENTER);
    add(labelMod);
    textMod = new TextField(64);
    textMod.reshape(60, 110, 520, 30);
    textMod.setEditable(true);
    add(textMod);
    btnSolve = new Button("Find discrete logarithm");
    btnSolve.reshape(210, 160, 180, 30);
    btnSolve.setFont(new Font("Courier", Font.PLAIN, 12));
    add(btnSolve);
    labelStatus = new Label("");
    labelStatus.reshape(30, 205, 520, 14);
    labelStatus.setFont(new Font("Courier", Font.PLAIN, 12));
    labelStatus.setAlignment(Label.CENTER);
    add(labelStatus);
    textExp = new TextField(64);
    labelExp = new Label("Exp");
    labelExp.reshape(10, 235, 40, 14);
    labelExp.setFont(new Font("Courier", Font.PLAIN, 12));
    labelExp.setAlignment(Label.CENTER);
    add(labelExp);
    textExp = new TextField(64);
    textExp.reshape(60, 230, 520, 30);
    textExp.setEditable(false);
    add(textExp);
    textPeriod = new TextField(64);
    labelPeriod = new Label("Period");
    labelPeriod.reshape(10, 285, 40, 14);
    labelPeriod.setFont(new Font("Courier", Font.PLAIN, 12));
    labelPeriod.setAlignment(Label.CENTER);
    add(labelPeriod);
    textPeriod = new TextField(64);
    textPeriod.reshape(60, 280, 520, 30);
    textPeriod.setEditable(false);
    add(textPeriod);
    labelBottom = new Label("Written by Dario Alejandro Alpern. Last updated April 15th, 2005");
    labelBottom.reshape(10, 330, 570, 14);
    labelBottom.setFont(new Font("Courier", Font.PLAIN, 12));
    labelBottom.setAlignment(Label.CENTER);
    add(labelBottom);
    validate();
    textBase.requestFocus();
    }

public boolean handleEvent(Event e) {
  if (e.id == Event.ACTION_EVENT && e.target == btnSolve ||
      e.id == Event.KEY_PRESS && e.key == 13) {
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

  if (e.id == Event.KEY_PRESS && e.key == Event.TAB) {
    if (e.target == textBase) {textPower.requestFocus();}
    if (e.target == textPower) {textMod.requestFocus();}
    if (e.target == textMod) {textBase.requestFocus();}
    return true;
    }
  return super.handleEvent(e);
  }

  public void run() {
    BigInteger groupOrder, subGroupOrder, range, powSubGroupOrder;
    BigInteger GroupOrder = null;
    BigInteger Exponent, runningExp, baseExp, mod;
    BigInteger basePH, powerPH;
    BigInteger logar=null, logarMult;
    BigInteger divideExp, nbrD;
    BigInteger bigNbrA2, bigNbrB2;
    BigInteger DiscreteLog, DiscreteLogPeriod;
    int d, j, indexBase, indexExp;
    int i, index;
    long addA, addB, addA2, addB2;
    long mult1, mult2;
    long magnitude;
    int mostSignificantDword, leastSignificantDword;
    long firstLimit, secondLimit;
    long brentK, brentR;
    int ExpressionRC;
    BigInteger ExpressionResult[] = new BigInteger[1];

    Old=System.currentTimeMillis();
    OldTimeElapsed = 0;
    lModularMult = 0;
    TerminateThread = false;
    try {

      try {
        ExpressionRC = expression.ComputeExpression(textMod.getText().trim(),1,ExpressionResult);
      } catch (OutOfMemoryError e) {
        textExp.setText("Out of memory.");
        return;
      } catch (ArithmeticException e) {
        return;
      }
      modulus= ExpressionResult[0].abs();
      if (ExpressionRC != 0) {
        textExp.setText("Modulus: "+expressionText[-1-ExpressionRC]);
        return;
      }

      try {
        ExpressionRC = expression.ComputeExpression(textBase.getText().trim(),1,ExpressionResult);
      } catch (OutOfMemoryError e) {
        textExp.setText("Out of memory.");
        return;
      } catch (ArithmeticException e) {
        return;
      }
      base = ExpressionResult[0].mod(modulus);
      if (ExpressionRC != 0) {
        textExp.setText("Base: "+expressionText[-1-ExpressionRC]);
        return;
      }

      try {
        ExpressionRC = expression.ComputeExpression(textPower.getText().trim(),1,ExpressionResult);
      } catch (OutOfMemoryError e) {
        textExp.setText("Out of memory.");
        return;
      } catch (ArithmeticException e) {
        return;
      }
      power = ExpressionResult[0].abs();
      if (ExpressionRC != 0) {
        textExp.setText("Power: "+expressionText[-1-ExpressionRC]);
        return;
      }

      } catch (Exception e) {
        textExp.setText("Invalid data entered");
        return;
        }
    if (modulus.gcd(base).equals(BigInt1) == false) {
      textExp.setText("Modulus and base are not relatively prime.");
      return;
      }
    if (modulus.gcd(power).equals(BigInt1) == false) {
      textExp.setText("Modulus and power are not relatively prime.");
      return;
      }
    if (modulus.compareTo(BigInt1) <= 0) {
      textExp.setText("Modulus must be greater than 1.");
      return;
      }
    if (LastModulus.equals(modulus) == false) {
      if (GetSmallFactors(modulus, PrimesModBak, ExponentsModBak, Typ, 0) != 1) {
        LastModulus = modulus;
        factorPanel = new ecm();
        factorWindow = new Frame("Modulus factorization");
        factorWindow.setResizable(false);
        factorWindow.add(factorPanel);
        factorPanel.setSize(600, 390);
        Insets in = factorWindow.getInsets();
        factorWindow.setSize(600+in.right+in.left, 390+in.top+in.bottom);
        factorWindow.show();
        NbrFactorsMod = factorPanel.getFactors(modulus, PrimesModBak, ExponentsModBak);
        factorWindow.remove(factorPanel);
        factorWindow.dispose();
        }
      else {
        NbrFactorsMod = NbrFactors;
        }
      for (j=0; j<NbrFactorsMod; j++) {
        PrimesMod[j] = PrimesModBak[j];
        ExponentsMod[j] = ExponentsModBak[j];
        }
      }
    DiscreteLog = BigInt0;
    DiscreteLogPeriod = BigInt1;
    for (index = 0; index < NbrFactorsMod; index++) {
      groupOrder = PrimesMod[index].subtract(BigInt1);
      textExp.setText("Computing discrete logarithm...");
      if (GetSmallFactors(groupOrder, PrimesGO, ExponentsGO, Typ, 0) != 1) {
        LastModulus = modulus;
        factorPanel = new ecm();
        factorWindow = new Frame("Group order factorization");
        factorWindow.setResizable(false);
        factorWindow.add(factorPanel);
        factorPanel.setSize(600, 390);
        Insets in = factorWindow.getInsets();
        factorWindow.setSize(600+in.right+in.left, 390+in.top+in.bottom);
        factorWindow.show();
        NbrFactors = factorPanel.getFactors(groupOrder, PrimesGO, ExponentsGO);
        factorWindow.remove(factorPanel);
        factorWindow.dispose();
        }
      mod = PrimesMod[index];
      range = mod.divide(BigInt3);
      logar = BigInt0;
      logarMult = BigInt1;
      BigNbrToBigInt(mod);
      yieldFreq = 1000000/(NumberLength*NumberLength);
      GetMontgomeryParms();
      if (NumberLength == 1) {
        mostSignificantDword = NLen - 1;
        leastSignificantDword = NumberLength - 1;
        }
      else {
        mostSignificantDword = NumberLength - 1;
        leastSignificantDword = NumberLength - 2;
        }
      secondLimit = ((TestNbr[mostSignificantDword] << 31) +
                      TestNbr[leastSignificantDword]) * 2 / 3;
      firstLimit = secondLimit / 2;

      for (indexBase = 0; indexBase < NbrFactors; indexBase++) {
        textExp.setText("Computing discrete logarithm in subgroup of "+PrimesGO[indexBase]+" elements.");
        runningExp = BigInt0;
        subGroupOrder = PrimesGO[indexBase];
        ExchangeMods();                    // TestNbr <- subGroupOrder
        BigNbrToBigInt(subGroupOrder);
        dN = (double)TestNbr[NumberLength-1];
        if (NumberLength > 1) {
          dN += (double)TestNbr[NumberLength-2]/dDosALa31;
        }
        if (NumberLength > 2) {
          dN += (double)TestNbr[NumberLength-3]/dDosALa62;
        }
        ExchangeMods();                    // TestNbr <- modulus
        baseExp = groupOrder;
        powSubGroupOrder = BigInt1;
        for (indexExp = 0; indexExp < ExponentsGO[indexBase]; indexExp++) {

        /* PH below comes from Pohlig-Hellman algorithm */

          divideExp = BigInt1;
          j = ExponentsGO[indexBase] - indexExp;
          do {
            divideExp = divideExp.multiply(subGroupOrder);
            basePH = base.modPow(groupOrder.divide(divideExp), mod);
            j--;
            } while (j>0 && basePH.equals(BigInt1));
          powerPH = power.multiply(base.modPow(runningExp.negate(), mod)).
                          modPow(baseExp.divide(divideExp),mod);
          baseExp = baseExp.divide(subGroupOrder);
          if (subGroupOrder.compareTo(BigInteger.valueOf(20L)) < 0) {
            for (j=subGroupOrder.intValue()-1; j>=0; j--) {
              if (basePH.modPow(BigInteger.valueOf(j), mod).equals(powerPH)) {
                break;
                }
              }
            if (j<0) {
              textExp.setText("Cannot compute discrete logarithm    subgroup="+PrimesGO[indexBase]+"   exponent="+indexExp);
              return;
              }
            Exponent = BigInteger.valueOf(j);
            }
          else {        // Use Pollard's rho method with Brent's modification
            BigNbrToMont(powerPH, nbrPower);
            BigNbrToMont(basePH, nbrBase);
        
            for (j=0; j<NumberLength; j++) {
              nbrR2[j] = nbrBase[j];
              nbrA2[j] = nbrB2[j] = 0;
              }
            nbrB2[0] = 1;
            addA2 = addB2 = 0;
            mult2 = 1;
            brentR = 1;
            brentK = 0;
PollardBrentRho:
            do {
              for (j=0; j<NumberLength; j++) {
                nbrR[j] = nbrR2[j];
                nbrA[j] = nbrA2[j];
                nbrB[j] = nbrB2[j];
                }
              addA = addA2;
              addB = addB2;
              mult1 = mult2;
              brentR *= 2;
              do {
                if (TerminateThread) {
                  return;                /* End thread */
                }
                brentK++;
                magnitude = (nbrR2[mostSignificantDword] << 31) +
                             nbrR2[leastSignificantDword];
                if (magnitude < firstLimit) {
                  MontgomeryMult(nbrR2, nbrPower, nbrROther);
                  addA2++;
                  }
                else {
                  if (magnitude < secondLimit) {
                    MontgomeryMult(nbrR2, nbrR2, nbrROther);
                    mult2 *= 2;
                    addA2 *= 2;
                    addB2 *= 2;
                    }
                  else {
                    MontgomeryMult(nbrR2, nbrBase, nbrROther);
                    addB2++;
                    }
                  }
                long [] nbrTemp = nbrR2;
                nbrR2 = nbrROther;
                nbrROther = nbrTemp;
                if (mult2 > 10000000) {
                  ExchangeMods();              // TestNbr <- subGroupOrder
                  AdjustExponent(nbrA2, mult2, addA2);
                  AdjustExponent(nbrB2, mult2, addB2);
                  ExchangeMods();              // TestNbr <- modulus
                  mult2 = 1;
                  addA2 = addB2 = 0;
                  }
                for (j=0; j<NumberLength; j++) {
                  if (nbrR[j] != nbrR2[j]) {break;}
                  }
                if (j==NumberLength) {break PollardBrentRho;}
                } while (brentK < brentR);
              } while (true);
            ExchangeMods();                  // TestNbr <- subGroupOrder
            AdjustExponent(nbrA, mult1, addA);
            AdjustExponent(nbrB, mult1, addB);
            AdjustExponent(nbrA2, mult2, addA2);
            AdjustExponent(nbrB2, mult2, addB2);
            SubtractBigNbrModN(nbrA, nbrA2, nbrA);
            SubtractBigNbrModN(nbrB2, nbrB, nbrB);
            if (BigNbrIsZero(nbrA)) {
              for (j=0; j<NumberLength; j++) {
                nbrK[j] = 0;
                }
              bigNbrA2 = BigInt0;
              bigNbrB2 = BigIntToBigNbr(nbrB);
              d = subGroupOrder.intValue();
              }
            else {                           // gcd is 1.
              ModInvBigNbr(nbrA, nbrB2, TestNbr);
              MultBigNbrModN(nbrB2, nbrB, nbrK);
              bigNbrA2 = BigInt1;
              bigNbrB2 = basePH.modPow(BigIntToBigNbr(nbrK), mod);
              d = 1;
              }
            ExchangeMods();                  // TestNbr <- modulus
            nbrD = BigInteger.valueOf(d);
            for (j=0; j<d; j++) {
              if (powerPH.equals(bigNbrB2)) {break;}
              bigNbrB2 = bigNbrB2.multiply(bigNbrA2).mod(mod);
              }
            if (j==d) {
              textExp.setText("Cannot compute discrete logarithm    subgroup="+PrimesGO[indexBase]+"   exponent="+indexExp);
              return;
              }
            Exponent = BigIntToBigNbr(nbrK).add(subGroupOrder.divide(nbrD).
                       multiply(BigInteger.valueOf(j)).mod(subGroupOrder));
            }
          runningExp = runningExp.add(Exponent.multiply(powSubGroupOrder));
          powSubGroupOrder = powSubGroupOrder.multiply(subGroupOrder);
          }
        nbrV[indexBase] = runningExp;
        for (indexExp=0; indexExp < indexBase; indexExp++) {
          nbrV[indexBase] = nbrV[indexBase].subtract(nbrV[indexExp]).
                            multiply(PrimesGO[indexExp].pow(ExponentsGO[indexExp]).modInverse(powSubGroupOrder));
          }
        nbrV[indexBase] = nbrV[indexBase].mod(powSubGroupOrder);
        logar = logar.add(nbrV[indexBase].multiply(logarMult));
        logarMult = logarMult.multiply(powSubGroupOrder);
        }

      // Compute group order.
      BigInteger F = PrimesModBak[index];
      BigInteger G = F.subtract(BigInt1);
      BigInteger P = F;
      for (indexBase = 0; indexBase < NbrFactors; indexBase++) {
        for (indexExp = 0; indexExp < ExponentsGO[indexBase]; indexExp++) {
          if (base.modPow(G.divide(PrimesGO[indexBase]),F).equals(BigInt1)) {
            G = G.divide(PrimesGO[indexBase]);
            }
          }
        }
      logar = logar.mod(G);
      GroupOrder = G;
      for (i=2; i<ExponentsModBak[index]*2; i*=2) {
        BigInteger H = F.pow(2);
        BigInteger K0 = base.modPow(logar, H);
        BigInteger K1 = base.modPow(logar.add(G), H);
        BigInteger Num = power.subtract(K0).divide(F).mod(F);
        BigInteger Den = K1.subtract(K0).divide(F).mod(F);
        if (Den.equals(BigInt0)) {
          if (Num.equals(BigInt0)) {
            textExp.setText("Indetermination 0/0 -> Cannot compute discrete logarithm");
            }
          else {
            textExp.setText("There is no discrete logarithm");
            }
          return;
          }
        logar = logar.add(Num.multiply(Den.modInverse(F)).
                      mod(F).multiply(G));
        G = G.multiply(F);
        F = H;
        }
      GroupOrder = GroupOrder.
               multiply(PrimesMod[index].pow(ExponentsMod[index]-1));
      logar = logar.mod(GroupOrder);
      if (DiscreteLogPeriod.equals(BigInt1)) {
        DiscreteLogPeriod = GroupOrder;
        DiscreteLog = logar;
        }
      else {
        BigInteger gcd = GroupOrder.gcd(DiscreteLogPeriod);
        if (logar.mod(gcd).equals(DiscreteLog.mod(gcd)) == false) {
          textExp.setText("Cannot compute discrete logarithm");
          return;
          }
        if (GroupOrder.divide(gcd).gcd(gcd).equals(BigInt1) == false) {
          BigInteger tmp = GroupOrder;
          GroupOrder = DiscreteLogPeriod;
          DiscreteLogPeriod = tmp;
          tmp = logar;
          logar = DiscreteLog;
          DiscreteLog = tmp;
          }
        GroupOrder = GroupOrder.divide(gcd);
        logar = logar.mod(GroupOrder);
        DiscreteLog = (logar.subtract(DiscreteLog)).
                multiply(DiscreteLogPeriod.modInverse(GroupOrder)).
                mod(GroupOrder).multiply(DiscreteLogPeriod).add(DiscreteLog);
        DiscreteLogPeriod = DiscreteLogPeriod.multiply(GroupOrder);
        }
      }
    textExp.setText(DiscreteLog.toString());
    textPeriod.setText(DiscreteLogPeriod.toString());
    long t=OldTimeElapsed/1000;
    labelStatus.setText("Time elapsed: "+
           t/86400+"d "+(t%86400)/3600+"h "+((t%3600)/60)+"m "+(t%60)+
           "s    mod mult: "+lModularMult);
    }

  private void BigNbrToMont(BigInteger bigNbr, long [] nbr) {
    int i;
    long tmp;
    int NumberLengthBak = NumberLength;
    System.arraycopy(TestNbr, 0, CalcAuxGcdU, 0, NLen);
    BigNbrToBigInt(bigNbr);
    for (i=NumberLength; i<NumberLengthBak; i++) {
      TestNbr[i] = 0;
      }
    NumberLength = NumberLengthBak;
    for (i=0; i<NumberLength; i++) {
      tmp = TestNbr[i];
      TestNbr[i] = CalcAuxGcdU[i];
      CalcAuxGcdU[i] = tmp;
      }
    MultBigNbrModN(CalcAuxGcdU, MontgomeryMultR1, nbr);
    }

  // nbr = (nbr * mult + add) % TestNbr
  private void AdjustExponent(long [] nbr, long mult, long add) {
    long Pr;
    int j;

    Pr = add * DosALa31;
    for (j=0; j<NumberLength; j++) {
      Pr = (Pr >>> 31) + mult*nbr[j];
      nbr[j] = Pr & 0x7FFFFFFFl;
    }
    nbr[j] = (Pr >>> 31);
    AdjustModN(nbr);
    }

  private void ExchangeMods() {
    long [] Tmp;
    int Temp;
    double dTemp;

    Tmp = TestNbr;
    TestNbr = TestNbrOther;
    TestNbrOther = Tmp;
    Temp = NumberLength;
    NumberLength = NumberLengthOther;
    NumberLengthOther = Temp;
    dTemp = dN;
    dN = dNOther;
    dNOther = dTemp;
    }

}       // End applet class

// </XMP>
