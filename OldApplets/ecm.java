// Elliptic Curve Method (ECM) Prime Factorization
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// Last updated September 10th, 2011. See http://www.alpertron.com.ar/ECM.HTM
//
// Based in Yuji Kida's implementation for UBASIC interpreter
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely
// except that you have to write somewhere in the code this header.
//
import java.applet.*;
import java.awt.*;
import java.math.*;
import java.awt.event.*;
import java.net.*;
import java.io.*;

public class ecm extends Applet implements Runnable, factorApplet
{
  static final long serialVersionUID = 20L;
  static final boolean KARATSUBA_ENABLED = false;
  static final boolean debugging_main = false;
  boolean onlyFactoring = true;
  static final int ACTION_TEXT_NBR = 1;
  static final int ACTION_TEXT_CURVE = 2;
  static final int ACTION_BTN_CURVE = 3;
  static final int ACTION_TEXT_FACTOR = 4;
  static final int ACTION_BTN_FACTOR = 5;

  static final int TYP_AURIF = 100000000;
  static final int TYP_TABLE = 150000000;
  static final int TYP_SIQS = 200000000;
  static final int TYP_LEHMAN = 250000000;
  static final int TYP_RABIN = 300000000;
  static final int TYP_EC = 350000000;

  int digitsInGroup = 6;
  int numberThreads = 1;
  boolean forcedECM = false;
  final BigInteger SS[] = new BigInteger[4000]; // For intermediate factors 
  final BigInteger PD[] = new BigInteger[4000]; // and prime factors
  boolean numberIsNegative;
  private BigInteger probableFactor;
  private final int Exp[] = new int[4000];
  private final int Typ[] = new int[4000];
  final BigInteger PD1[] = new BigInteger[4000];
  final int Exp1[] = new int[4000];
  final int Typ1[] = new int[4000];
  String inputStr;
  boolean foundByLehman;
  boolean performLehman;
  static final BigInteger BigInt0 = BigInteger.valueOf(0L);
  static final BigInteger BigInt1 = BigInteger.valueOf(1L);
  static final BigInteger BigInt2 = BigInteger.valueOf(2L);
  static final BigInteger BigInt3 = BigInteger.valueOf(3L);
  static final int PWmax = 32, Qmax = 30241, LEVELmax = 11;
  static final int NLen = 1200;
  final int aiIndx[] = new int[Qmax];
  final int aiF[] = new int[Qmax];
  static final int aiP[] = { 2, 3, 5, 7, 11, 13 };
  static final int aiQ[] =
    {
      2,
      3,
      5,
      7,
      13,
      11,
      31,
      61,
      19,
      37,
      181,
      29,
      43,
      71,
      127,
      211,
      421,
      631,
      41,
      73,
      281,
      2521,
      17,
      113,
      241,
      337,
      1009,
      109,
      271,
      379,
      433,
      541,
      757,
      2161,
      7561,
      15121,
      23,
      67,
      89,
      199,
      331,
      397,
      463,
      617,
      661,
      881,
      991,
      1321,
      2311,
      2377,
      2971,
      3697,
      4159,
      4621,
      8317,
      9241,
      16633,
      18481,
      23761,
      101,
      151,
      401,
      601,
      701,
      1051,
      1201,
      1801,
      2801,
      3301,
      3851,
      4201,
      4951,
      6301,
      9901,
      11551,
      12601,
      14851,
      15401,
      19801,
      97,
      353,
      673,
      2017,
      3169,
      3361,
      5281,
      7393,
      21601,
      30241,
      53,
      79,
      131,
      157,
      313,
      521,
      547,
      859,
      911,
      937,
      1093,
      1171,
      1249,
      1301,
      1873,
      1951,
      2003,
      2081,
      41,
      2731,
      2861,
      3121,
      3433,
      3511,
      5851,
      6007,
      6553,
      7151,
      7723,
      8009,
      8191,
      8581,
      8737,
      9829,
      11701,
      13729,
      14561,
      15601,
      16381,
      17551,
      20021,
      20593,
      21841,
      24571,
      25741,
      26209,
      28081 };
  static final int aiG[] =
    {
      1,
      2,
      2,
      3,
      2,
      2,
      3,
      2,
      2,
      2,
      2,
      2,
      3,
      7,
      3,
      2,
      2,
      3,
      6,
      5,
      3,
      17,
      3,
      3,
      7,
      10,
      11,
      6,
      6,
      2,
      5,
      2,
      2,
      23,
      13,
      11,
      5,
      2,
      3,
      3,
      3,
      5,
      3,
      3,
      2,
      3,
      6,
      13,
      3,
      5,
      10,
      5,
      3,
      2,
      6,
      13,
      15,
      13,
      7,
      2,
      6,
      3,
      7,
      2,
      7,
      11,
      11,
      3,
      6,
      2,
      11,
      6,
      10,
      2,
      7,
      11,
      2,
      6,
      13,
      5,
      3,
      5,
      5,
      7,
      22,
      7,
      5,
      7,
      11,
      2,
      3,
      2,
      5,
      10,
      3,
      2,
      2,
      17,
      5,
      5,
      2,
      7,
      2,
      10,
      3,
      5,
      3,
      7,
      3,
      2,
      7,
      5,
      7,
      2,
      3,
      10,
      7,
      3,
      3,
      17,
      6,
      5,
      10,
      6,
      23,
      6,
      23,
      2,
      3,
      3,
      5,
      11,
      7,
      6,
      11,
      19 };
  final int aiInv[] = new int[PWmax];
  static final int aiNP[] = { 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 };
  static final int aiNQ[] = { 5, 8, 11, 18, 22, 27, 36, 59, 79, 89, 136 };
  static final int aiT[] =
    {
      2 * 2 * 3,
      2 * 2 * 3 * 5,
      2 * 2 * 3 * 3 * 5,
      2 * 2 * 3 * 3 * 5 * 7,
      2 * 2 * 2 * 3 * 3 * 5 * 7,
      2 * 2 * 2 * 2 * 3 * 3 * 5 * 7,
      2 * 2 * 2 * 2 * 3 * 3 * 3 * 5 * 7,
      2 * 2 * 2 * 2 * 3 * 3 * 3 * 5 * 7 * 11,
      2 * 2 * 2 * 2 * 3 * 3 * 3 * 5 * 5 * 7 * 11,
      2 * 2 * 2 * 2 * 2 * 3 * 3 * 3 * 5 * 5 * 7 * 11,
      2 * 2 * 2 * 2 * 2 * 3 * 3 * 3 * 5 * 5 * 7 * 11 * 13 };
  final int biTmp[] = new int[NLen];
  final int biExp[] = new int[NLen];
  final int biN[] = new int[NLen];
  final int biR[] = new int[NLen];
  final int biS[] = new int[NLen];
  final int biT[] = new int[NLen];
  final int biU[] = new int[NLen]; /* Temp */
  final int biV[] = new int[NLen]; /* Temp */
  final int biW[] = new int[NLen]; /* Temp */
  final int aiJS[][] = new int[PWmax][NLen];
  final int aiJW[][] = new int[PWmax][NLen];
  final int aiJX[][] = new int[PWmax][NLen];
  final int aiJ0[][] = new int[PWmax][NLen];
  final int aiJ1[][] = new int[PWmax][NLen];
  final int aiJ2[][] = new int[PWmax][NLen];
  final int aiJ00[][] = new int[PWmax][NLen];
  final int aiJ01[][] = new int[PWmax][NLen];
  int NumberLength; /* Length of multiple precision nbrs */
  int NbrFactors, NbrFactors1;
  int EC; /* Elliptic Curve number */
  /* Used inside GCD calculations in multiple precision numbers */
  final int CalcAuxGcdU[] = new int[NLen];
  final int CalcAuxGcdV[] = new int[NLen];
  final int CalcAuxGcdT[] = new int[NLen];
  final int CalcBigNbr[] = new int[NLen];
  final long CalcAuxModInvA[] = new long[NLen];
  final long CalcAuxModInvB[] = new long[NLen];
  final long CalcAuxModInvBB[] = new long[NLen];
  final long CalcAuxModInvMu[] = new long[NLen];
  final long CalcAuxModInvGamma[] = new long[NLen];
  int TestNbr[] = new int[NLen];
  int KaratsubaTestNbr[] = new int[NLen];
  final long TestNbr2[] = new long[NLen];
  final int GcdAccumulated[] = new int[NLen];
  final long Gamma[] = new long[386];
  final long Delta[] = new long[386];
  final long AurifQ[] = new long[386];
  final BigInteger Factores[] = new BigInteger[200];
  int[] fieldAA, fieldTX, fieldTZ, fieldUX, fieldUZ;
  int[] fieldAux1, fieldAux2, fieldAux3, fieldAux4;
  int NroFact, DegreeAurif;
  static final long DosALa32 = (long) 1 << 32;
  static final long DosALa31 = (long) 1 << 31;
  static final long DosALa31_1 = DosALa31 - 1;
  static final long DosALa63 = DosALa32 * DosALa31;
  static final double dDosALa31 = (double) DosALa31;
  static final double dDosALa62 = dDosALa31 * dDosALa31;
  static final long Mi = 1000000000;
  static long timePrimalityTests;
  static long timeSIQS;
  static long timeECM;
  static int nbrPrimalityTests;
  static int nbrSIQS;
  static int nbrECM;
  double dN;
  final int BigNbr1[] = new int[NLen];
  final int SmallPrime[] = new int[670]; /* Primes < 5000 */
  final int MontgomeryMultR1[] = new int[NLen];
  final int MontgomeryMultR2[] = new int[NLen];
  final int MontgomeryMultAfterInv[] = new int[NLen];
  long MontgomeryMultN;
  final double ProbArray[] = new double[6];
  TextArea upperTextArea;
  TextArea lowerTextArea;
  Label labelStatus;
  Label labelBottom;
  Label labelTop;
  TextField textNumber;
  TextField textCurve;
  TextField textFactor;
  Button btnCurve;
  Button btnFactor;
  BigInteger NumberToFactor;
  String StringToLabel;
  StringBuffer outputStr;
  boolean batchFinished = true;
  boolean batchPrime = false;
  volatile Thread calcThread;
  String textAreaContents;
  long OldTimeElapsed;
  long Old;
  long yieldFreq;
  boolean TerminateThread = true;
  private int NextEC = -1;
  private BigInteger InputFactor = BigInt0;
  private int StepECM;
  private int indexM, maxIndexM;
  private int nbrPrimes, indexPrimes;
  private BigInteger Quad1, Quad2, Quad3, Quad4;
  private boolean Computing3Squares;
  long polynomialsSieved;
  long trialDivisions;
  long smoothsFound;
  long totalPartials;
  long partialsFound;
  long ValuesSieved;
  private long primeModMult;
  private long lModularMult;
  private final int[] arrNbr = new int[4*NLen];
  private final int[] arrNbrM = new int[4*NLen];
  private final int[] arrNbrAux = new int[4*NLen];
  private final int[] montgKaratsubaArr = new int[2*NLen];
  private static final int KARATSUBA_CUTOFF = 64;
  private int karatLength;

  /* ECM limits for 30, 35, ..., 85 digits */
  static final int limits[] = { 5, 8, 15, 25, 25, 27, 32, 43, 70, 150, 300, 350, 600 };
  static final String[] expressionText =
  {
    "Number too high (more than 10000 digits).",
    "Intermediate expression too high (more than 20000 digits).",
    "Non-integer division.",
    "Parentheses mismatch.",
    "Syntax error.",
    "Too many parentheses.",
    "Invalid parameter."
  };

  public static void main(final String args[])
  {
    if (debugging_main)
    {
      if (true)
      {
        int [] queue = new int[200];
        int i, a=0;
        long value = Long.parseLong(args[0]);
        long time = System.currentTimeMillis();
        for (i=0; i<100000; i++)
        {
          a = Siqs.SQUFOF(value,queue);
        }
        System.out.println("factor="+a+", time="+(System.currentTimeMillis()-time)+" msec");
        if (Siqs.isProbablePrime(value))
        {
          System.out.println(value + " is pseudoprime.");
        }
        else
        {
          System.out.println(value + " is composite.");
        }
      }
      /*
      final ecm ecm1 = new ecm();
      ecm1.init();
      ecm1.StartFactorExprBatch("-991", 1);
      ecm1.StartFactorExprBatch("10^59+213", 0);
      ecm1.StartFactorExprBatch("10^36+68", 0);
      BigInteger big1, big3;
      final int exp=344;
      big3 = new BigInteger("10").pow(exp).divide(BigInteger.valueOf(9));
      NumberLength = ecm1.BigNbrToBigInt(big3, TestNbr);
      System.out.println("Initial number = 10^"+exp+" / 9");
      System.out.println("NumberLength = "+ecm1.NumberLength);
      System.out.println("TestNbr[0] = "+ecm1.TestNbr[0]);
      System.out.println("TestNbr[1] = "+ecm1.TestNbr[1]);
      System.out.println(BigIntToBigNbr(ecm1.TestNbr, ecm1.NumberLength));
      ecm1.TerminateThread = false;
//    ecm1.GetMontgomeryParms();
//    big1 = BigIntToBigNbr(ecm1.MontgomeryMultR2, ecm1.NumberLength);
//    System.out.println("ANTES ModInvBigNbr: "+big1.toString());
//    ModInvBigNbr(ecm1.MontgomeryMultR2, ecm1.MontgomeryMultR2, ecm1.TestNbr);
//    big2 = BigIntToBigNbr(ecm1.MontgomeryMultR2, ecm1.NumberLength);
//    System.out.println("DESPUES ModInvBigNbr: "+big2.toString());
//    System.out.println(big1.multiply(big2).mod(big3).toString());
      MultBigNbrByLong(ecm1.TestNbr, -7, ecm1.TestNbr, NumberLength);
      MultBigNbrByLong(ecm1.TestNbr, -7, ecm1.TestNbr, NumberLength);
      DivBigNbrByLong(ecm1.TestNbr, -7, ecm1.TestNbr, NumberLength);
      DivBigNbrByLong(ecm1.TestNbr, -7, ecm1.TestNbr, NumberLength);
      big1 = BigIntToBigNbr(ecm1.TestNbr, ecm1.NumberLength);
      System.out.println(big1.toString());
      NumberLength = ecm1.BigNbrToBigInt(big3, TestNbr);
      LongToBigNbr(-7, ecm1.MontgomeryMultR2, NumberLength);
      MultBigNbr(ecm1.MontgomeryMultR2, ecm1.TestNbr, ecm1.MontgomeryMultR1,
                 ecm1.NumberLength);
      MultBigNbr(ecm1.MontgomeryMultR1, ecm1.MontgomeryMultR2, ecm1.TestNbr,
                 ecm1.NumberLength);
      big1 = BigIntToBigNbr(ecm1.TestNbr, int NumberLength);
      System.out.println(big1.toString());
    */
    }
    else
    {
      Frame frame = new Frame("Integer factorization using ECM/SIQS");
      ecm ecm1 = new ecm();
      frame.addWindowListener(new PanelWindowListener());
      frame.setLayout(new BorderLayout());
      frame.add("Center", ecm1);
      frame.setSize(620, 420);
      ecm1.init();
      ecm1.start();
      frame.setVisible(true);
    }
  }

  public void init()
  {
    if (onlyFactoring)
    {
      layout(
        "Type number or numerical expression to factor here and press Return:",
        true);
    }
  }

  public int getFactors(
    final BigInteger NbrToFactor,
    final BigInteger[] Primes,
    final int[] Exponents)
  {
    layout("Number to factor:", false);
    textNumber.setText(NbrToFactor.toString());
    startNewFactorization(true);       // Request complete factorization
    System.arraycopy(PD, 0, Primes, 0, NbrFactors);
    System.arraycopy(Exp, 0, Exponents, 0, NbrFactors);
    return NbrFactors;
  }

  void layout(final String caption, final boolean editable)
  {
    setBackground(Color.lightGray);
    setLayout(null);
    labelTop = new Label(caption);
    labelTop.setBounds(new Rectangle(10, 10, 570, 14));
    labelTop.setFont(new Font("Courier", Font.PLAIN, 12));
    labelTop.setAlignment(Label.CENTER);
    add(labelTop);
    textNumber = new TextField(1000);
    textNumber.setBounds(new Rectangle(10, 30, 570, 30));
    textNumber.setEditable(editable);
    add(textNumber);
    upperTextArea = new TextArea("", 7, 75, TextArea.SCROLLBARS_VERTICAL_ONLY);
    upperTextArea.setBounds(new Rectangle(10, 70, 570, 125));
    upperTextArea.setEditable(false);
    upperTextArea.setFont(new Font("Courier", Font.PLAIN, 12));
    add(upperTextArea);
    lowerTextArea = new TextArea("", 6, 75, TextArea.SCROLLBARS_VERTICAL_ONLY);
    lowerTextArea.setBounds(new Rectangle(10, 205, 570, 100));
    lowerTextArea.setEditable(false);
    lowerTextArea.setFont(new Font("Courier", Font.PLAIN, 12));
    add(lowerTextArea);
    labelStatus = new Label("");
    labelStatus.setBounds(new Rectangle(10, 315, 570, 14));
    labelStatus.setFont(new Font("Courier", Font.PLAIN, 12));
    add(labelStatus);
    textCurve = new TextField(10);
    textCurve.setBounds(new Rectangle(10, 335, 80, 30));
    add(textCurve);
    btnCurve = new Button("New curve");
    btnCurve.setBounds(new Rectangle(100, 335, 80, 30));
    btnCurve.setFont(new Font("Courier", Font.PLAIN, 12));
    add(btnCurve);
    textFactor = new TextField(200);
    textFactor.setBounds(new Rectangle(250, 335, 240, 30));
    add(textFactor);
    btnFactor = new Button("Factor");
    btnFactor.setBounds(new Rectangle(500, 335, 70, 30));
    btnFactor.setFont(new Font("Courier", Font.PLAIN, 12));
    add(btnFactor);
    labelBottom =
      new Label("Written by Dario Alejandro Alpern. Last updated January 21st, 2014");
    labelBottom.setBounds(new Rectangle(10, 370, 570, 14));
    labelBottom.setFont(new Font("Courier", Font.PLAIN, 12));
    labelBottom.setAlignment(Label.CENTER);
    add(labelBottom);
    if (onlyFactoring)
    {
      textNumber.addActionListener(new Command(ACTION_TEXT_NBR, this));
      textCurve.addActionListener(new Command(ACTION_TEXT_CURVE, this));
      btnCurve.addActionListener(new Command(ACTION_BTN_CURVE, this));
      textFactor.addActionListener(new Command(ACTION_TEXT_FACTOR, this));
      btnFactor.addActionListener(new Command(ACTION_BTN_FACTOR, this));
    }
    validate();
    textNumber.requestFocus();
    ProbArray[0] = 5E4;
    ProbArray[1] = 9.9E5;
    ProbArray[2] = 1.5E7;
    ProbArray[3] = 1.75E8;
    ProbArray[4] = 1.8E9;
    ProbArray[5] = 1.53E10;
  }

  public void destroy()
  { /* Applet end */
    TerminateThread = true;
  }

  public int setState(final String state)
  {
    if (onlyFactoring)
    {
      int index = 0;
      int indexNew;
      final int indexEnd = state.length();
      int indexComma = state.indexOf(',', index);
      int indexPlus, indexPlus2;
      int indexTimes;
      int indexExp;
      int indexParen;
      int i;

      if (indexComma == -1)
      {
        indexComma = indexEnd;
      }
      NumberToFactor = BigInt1;
      NbrFactors = 0;
      for (i = 0; i < 400; i++)
      {
        Exp[i] = Typ[i] = 0;
      }
      while (index < indexComma)
      {
        indexTimes = state.indexOf('x', index);
        if (indexTimes == -1)
        {
          indexTimes = indexComma;
        }
        indexExp = state.indexOf('^', index);
        indexParen = state.indexOf('(', index);
        if (indexParen > indexTimes || indexParen == -1)
        {
          indexParen = indexTimes;
        }
        if (indexExp > indexParen || indexExp == -1)
        {
          indexExp = indexParen;
        }
        indexNew = indexTimes + 1;
        switch (state.charAt(indexTimes - 1))
        {
          case ')' :
            Typ[NbrFactors] =
              Integer.parseInt(state.substring(indexParen + 1, indexTimes - 1));
            break;
          default :
            Typ[NbrFactors] = 0; /* Prime */
        }
        PD[NbrFactors] = new BigInteger(state.substring(index, indexExp));
        if (indexExp == indexParen)
        {
          Exp[NbrFactors] = 1;
        }
        else
        {
          Exp[NbrFactors] =
            Integer.parseInt(state.substring(indexExp + 1, indexParen));
        }
        NumberToFactor =
          NumberToFactor.multiply(PD[NbrFactors].pow(Exp[NbrFactors]));
        NbrFactors++;
        index = indexNew;
      } /* end while */
      lModularMult = -1;
      if (indexComma < indexEnd)
      {
        indexPlus = state.indexOf('+', indexComma);
        if (indexPlus > 0)
        {
          OldTimeElapsed =
            Long.parseLong(state.substring(indexComma + 1, indexPlus));
          indexPlus2 = state.indexOf('+', indexPlus + 1);
          if (indexPlus2 > 0)
          {
            lModularMult =
              Long.parseLong(state.substring(indexPlus + 1, indexPlus2));
            digitsInGroup = Integer.parseInt(state.substring(indexPlus2 + 1));
          }
          else
          {
            lModularMult = Long.parseLong(state.substring(indexPlus + 1));
          }
        }
        else
        {
          OldTimeElapsed = Long.parseLong(state.substring(indexComma + 1));
        }
      }
      else
      {
        OldTimeElapsed = -1; /* Old cookie format */
      }
      if (NumberToFactor.compareTo(BigInt1) > 0)
      {
        textNumber.setText(NumberToFactor.toString());
        calcThread = new Thread(this);
        calcThread.start();
      }
      return digitsInGroup; /* Number of digits in a group */
    }
    else
    {
      return 0;
    }
  }

  public void setDigits(final int nbrDigits)
  {
    if (onlyFactoring)
    {
      digitsInGroup = (digitsInGroup & 0xFC00) | nbrDigits;
      final String Tmp1 = textAreaContents;
      final String Tmp2 = StringToLabel;
      ShowUpperPane();
      textAreaContents = Tmp1;
      StringToLabel = Tmp2;
    }
  }

  public void setThreads(final int nbrThreads)
  {
    numberThreads = nbrThreads;
  }

  public void switchSIQS(final int newSwitchSIQS)
  {
    if (newSwitchSIQS == 1)
    {
      digitsInGroup &= 0xFBFF;
    }
    else
    {
      digitsInGroup |= 0x0400;
    }
  }

  public void useCunnTable(final int newUseCunnTable)
  {
    if (newUseCunnTable == 1)
    {
      digitsInGroup &= 0xEFFF;
    }
    else
    {
      digitsInGroup |= 0x1000;
    }
  }

  public void Verbose(final int verboseType)
  {
    if (verboseType == 0)
    {
      digitsInGroup &= 0xF7FF;
    }
    else
    {
      digitsInGroup |= 0x0800;
    }
    final String Tmp1 = textAreaContents;
    final String Tmp2 = StringToLabel;
    ShowUpperPane();
    textAreaContents = Tmp1;
    StringToLabel = Tmp2;
  }

  public String getState()
  {
    if (onlyFactoring)
    {
      int i;
      StringBuffer state = new StringBuffer();

      if (calcThread != null)
      {
        TerminateThread = true;
        try
        {
          calcThread.join(); /* Wait until the factorization thread dies */
        }
        catch (InterruptedException ie)
        {
        }
      }
      for (i = 0; i < (Computing3Squares ? NbrFactors1 : NbrFactors); i++)
      {
        if (i != 0)
        {
          state.append('x');
        }
        state.append(PD[i].toString());
        if (Exp[i] != 1)
        {
          state.append('^');
          state.append(Exp[i]);
        }
        if (Typ[i] != 0)
        {
          state.append('(');
          state.append(Typ[i]);
          state.append(')');
        }
      }
      if (OldTimeElapsed >= 0)
      {
        state.append(',');
        state.append(OldTimeElapsed);
      }
      if (lModularMult >= 0)
      {
        state.append('+');
        state.append(lModularMult);
      }
      state.append('+');
      state.append(digitsInGroup);
      return state.toString();
    }
    else
    {
      return null;
    }
  }

  public String getStringsFromBothPanes()
  {
    if (onlyFactoring)
    {
      final String lower = lowerTextArea.getText();
      StringBuffer text = new StringBuffer(
        "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>" +
        "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">" + 
        "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">" +
        "<head><title>Elliptic Curve Method factorization</title></head>" +
        "<body><h2>Elliptic Curve Method factorization</h2><p><i>Written by Dario Alejandro Alpern (dario@alpern.com.ar)</i></p><pre>");
      text.append(upperTextArea.getText());
      text.append("</pre><hr /><pre>");
      text.append(lower);
      if (lower.indexOf("complete") < 0)
      {
        final long New = System.currentTimeMillis();
        if (OldTimeElapsed >= 0)
        {
          OldTimeElapsed += New - Old;
          Old = New;
          text.append("</pre><hr /><pre>Elapsed time: ");
          text.append(GetDHMS(OldTimeElapsed / 1000));
          if (lModularMult >= 0)
          {
            text.append('\n');
            text.append(lModularMult);
            text.append(" modular multiplications have been done");
          }
        }
      }
      text.append("</pre></body></html>");
      return text.toString();
    }
    else
    {
      return null;
    }
  }

  String GetDHMS(final long time)
  {
    return time / 86400 + "d " + (time % 86400) / 3600 + "h " +
           ((time % 3600) / 60) + "m " + (time % 60) + "s";
  }

  String GetDHMSd(final long time)
  {
    return time / 864000 + "d " + (time % 864000) / 36000 + "h " +
           ((time % 36000) / 600) + "m " + ((time % 600)/10)+"."+
           (time%10) + "s";
  }

  void addStringToLabel(String value)
  {
    if (value.length() + 1 + StringToLabel.length() >= 79)
    {
      textAreaContents += StringToLabel + "\n";
      StringToLabel = value + " ";
    }
    else
    {
      StringToLabel += value + " ";
    }
  }

  void ShowUpperPane()
  {
    int i;

    textAreaContents = "";
    StringToLabel = "";
    if (numberIsNegative)
    {
      addStringToLabel("-");         // Indicate number is negative.
    }
    insertBigNbr(NumberToFactor);
    if (NbrFactors == 1 && Exp[0] == 1)
    {
      if (Typ[0] > 0)
      {
        addStringToLabel("is composite");
      }
      else
      {
        if (Typ[0] < 0)
        {
          addStringToLabel("is unknown");
        }
        else
        {
          addStringToLabel("is prime");
        }
      }
    }
    else
    {
      if (numberIsNegative)
      {
        addStringToLabel("= -1 *");
      }
      else
      {
        addStringToLabel("=");
      }
      for (i = 0; i < NbrFactors; i++)
      {
        if (i != 0)
        {
          addStringToLabel("x");
        }
        insertBigNbr(PD[i]);
        if (Exp[i] != 1)
        {
          addStringToLabel("^");
          addStringToLabel(String.valueOf(Exp[i]));
        }
        if (Typ[i] > 0)
        {
          int type = Typ[i];
          int compositeType = Typ[i] / 50000000 * 50000000;
          if ((digitsInGroup & 0x800) == 0x800)
          {
            if (type == TYP_AURIF)
            {
              addStringToLabel("(Aurifeuille)");
            }
            else if (compositeType == TYP_AURIF)
            {
              addStringToLabel("(Aurif - Composite)");
            }
            else if (type == TYP_TABLE)
            {
              addStringToLabel("(Table)");
            }
            else if (compositeType == TYP_TABLE)
            {
              addStringToLabel("(Table - Composite)");
            }
            else if (type == TYP_SIQS)
            {
              addStringToLabel("(SIQS)");
            }
            else if (compositeType == TYP_SIQS)
            {
              addStringToLabel("(SIQS - Composite)");
            }
            else if (type == TYP_LEHMAN)
            {
              addStringToLabel("(Lehman)");
            }
            else if (compositeType == TYP_LEHMAN)
            {
              addStringToLabel("(Lehman - Composite)");
            }
            else if (type == TYP_RABIN)
            {
              addStringToLabel("(Miller & Rabin)");
            }
            else if (compositeType == TYP_RABIN)
            {
              addStringToLabel("(Miller & Rabin - Composite)");
            }
            else if (type > TYP_EC)
            {
              addStringToLabel("(Curve " + (Typ[i] - TYP_EC) + ")");
            }
            else
            {
              addStringToLabel("(Composite)");
            }
          }
          else if (
            Typ[i] < TYP_EC
              && Typ[i] != TYP_LEHMAN
              && Typ[i] != TYP_SIQS
              && Typ[i] != TYP_AURIF
              && Typ[i] != TYP_RABIN
              && Typ[i] != TYP_TABLE)
          {
            addStringToLabel("(Composite)");
          }
        }
        else
        {
          if (Typ[i] < 0)
          {
            addStringToLabel("(Unknown)");
          }
        }
      }
    }
    upperTextArea.setText(textAreaContents + StringToLabel);
  }

  void insertBigNbr(final BigInteger N)
  {
    int i, dig;
    dig = digitsInGroup & 0x3FF;
    final String value = N.toString();
    i = (value.length() + dig - 1) % dig + 1;
    addStringToLabel(value.substring(0, i));
    while (i < value.length())
    {
      addStringToLabel(value.substring(i, i + dig));
      i += dig;
    }
  }

  void InsertNewFactor(final BigInteger InputFactor)
  {
    int g,exp;

    /* Insert input factor */
    for (g = NbrFactors - 1; g >= 0; g--)
    {
      PD[NbrFactors] = PD[g].gcd(InputFactor);
      if (PD[NbrFactors].equals(BigInt1) || PD[NbrFactors].equals(PD[g]))
      {
        continue;
      }
      for (exp=0; PD[g].remainder(PD[NbrFactors]).signum() == 0; exp++)
      {
        PD[g] = PD[g].divide(PD[NbrFactors]);
      }
      Exp[NbrFactors] = Exp[g] * exp;
      if (Typ[g] < 100000000)
      {
        Typ[g] = -EC;
        Typ[NbrFactors] = -TYP_EC - EC;
      }
      else if (Typ[g] < 150000000)
      {
        Typ[NbrFactors] = -Typ[g];
        Typ[g] = TYP_AURIF - Typ[g];
      }
      else if (Typ[g] < 200000000)
      {
        Typ[NbrFactors] = -Typ[g];
        Typ[g] = TYP_TABLE - Typ[g];
      }
      else if (Typ[g] < 250000000)
      {
        Typ[NbrFactors] = -Typ[g];
        Typ[g] = TYP_SIQS - Typ[g];
      }
      else if (Typ[g] < 300000000)
      {
        Typ[NbrFactors] = -Typ[g];
        Typ[g] = TYP_LEHMAN - Typ[g];
      }
      else
      {
        Typ[NbrFactors] = -Typ[g];
        Typ[g] = TYP_RABIN - Typ[g];
      }
      NbrFactors++;
    }
    SortFactorsInputNbr();
  }

  void SortFactorsInputNbr()
  {
    int g, i, j;
    BigInteger Nbr1;

    for (g = 0; g < NbrFactors - 1; g++)
    {
      for (j = g + 1; j < NbrFactors; j++)
      {
        if (PD[g].compareTo(PD[j]) > 0)
        {
          Nbr1 = PD[g];
          PD[g] = PD[j];
          PD[j] = Nbr1;
          i = Exp[g];
          Exp[g] = Exp[j];
          Exp[j] = i;
          i = Typ[g];
          Typ[g] = Typ[j];
          Typ[j] = i;
        }
      }
    }
  }

  void startNewFactorization(final boolean completefactorization)
  {
    if (calcThread != null && batchFinished)
    {
      TerminateThread = true;
      try
      {
        calcThread.join(); /* Wait until the factorization thread dies */
      }
      catch (InterruptedException ie)
      {
      }
    }
    onlyFactoring = completefactorization ^ true;
    lModularMult = OldTimeElapsed = NbrFactors = EC = 0;
    if (completefactorization)
    {
      factorize();
    }
    else
    {
      calcThread = new Thread(Thread.currentThread().getThreadGroup(),this); /* Start factorization thread */
      calcThread.start();
    }
  }

  public void run()
  {
    if (!batchFinished)
    {
      BatchThread();
      return;
    }
    polynomialsSieved = 0;
    trialDivisions = 0;
    smoothsFound = 0;
    totalPartials = 0;
    partialsFound = 0;
    factorize();
  }

  void factorize()
  {
    BigInteger NN;
    long TestComp, New;
    BigInteger N1, N2, Tmp;
    int i, j;
    int ExpressionRC;
    final BigInteger ExpressionResult[] = new BigInteger[1];

    timePrimalityTests = 0;  // Reset to zero all timings.
    timeSIQS = 0;
    timeECM = 0;
    nbrPrimalityTests = 0;   // Reset to zero all counters.
    nbrSIQS = 0;
    nbrECM = 0;
    StepECM = 0;
    primeModMult = 0;
    Computing3Squares = false;
    TerminateThread = false;
    Old = System.currentTimeMillis();
    if (onlyFactoring)
    {
      if (NbrFactors == 0)
      {
        lowerTextArea.setText("Computing input expression...");
        try
        {
          ExpressionRC =
            expression.ComputeExpression(
              textNumber.getText().trim(),
              0,
              ExpressionResult);
        }
        catch (OutOfMemoryError e)
        {
          lowerTextArea.setText("Out of memory.");
          return;
        }
        catch (ArithmeticException e)
        {
          return;
        }
        NumberToFactor = ExpressionResult[0];
        if (ExpressionRC < -1)
        {
          lowerTextArea.setText(expressionText[-2 - ExpressionRC]);
          return;
        }
      }
    }
    else
    {
      if (NbrFactors == 0)
      {
        NumberToFactor = new BigInteger(textNumber.getText().trim());
      }
    }
    if (NumberToFactor.abs().compareTo(BigInt1) <= 0)
    {                    // Factor number -1, 0 or 1.
      lowerTextArea.setText(NumberToFactor.toString()+" = "+
                            NumberToFactor.toString());
      return;
    }
    else if (NumberToFactor.signum() > 0)
    {                    // Factor positive number.
      numberIsNegative = false;
    }
    else
    {                    // Factor negative number.
      numberIsNegative = true;
      NumberToFactor = NumberToFactor.negate();  // Convert to positive.
    }
    probableFactor = null;   // No factor to insert now
    BigNbr1[0] = 1;
    for (i = 1; i < NLen; i++)
    {
      BigNbr1[i] = 0;
    }
    try
    {
      if (NbrFactors == 0)
      {
        lowerTextArea.setText(
          "Searching for small factors (less than 131072).");
        TestComp = GetSmallFactors(NumberToFactor, PD, Exp, Typ, 0);
        if (TestComp != 1)
        {                      // There are factors greater than 131071.
          PD[NbrFactors] = BigIntToBigNbr(TestNbr, NumberLength);
          Exp[NbrFactors] = 1;
          Typ[NbrFactors] = -1; /* Unknown */
          NbrFactors++;
          ShowUpperPane();
          if (batchFinished || !batchPrime)
          {
            lowerTextArea.setText("Searching for perfect power plus/minus 1.");
            PowerPM1Check(); /* Find algebraic and Aurifeuillian factors */
            lowerTextArea.setText("Searching for Lucas number.");
            LucasCheck();
            lowerTextArea.setText("Searching for Fibonacci number.");
            FibonacciCheck();
          }
        }
        else                   // No more factors
        {
          if (batchPrime && NbrFactors == 1 && Typ[0] == 0)
          {
            NbrFactors = 0;    // Indicate number is prime.
            return;
          }
        }
      }
      performLehman = true;
      factor_loop : for (;;)
      {
        ShowUpperPane();
        for (i = 0; i < NbrFactors; i++)
        {
          if (Typ[i] < 0)
          { /* Unknown */
            lowerTextArea.setText("Searching for perfect power.");
            if (PowerCheck(i) != 0)
            {
              SortFactorsInputNbr();
              continue factor_loop;
            }
            if (PD[i].bitLength() <= 33)
            {
              j = 0;
            }
            else
            {
              lowerTextArea.setText("Before calling prime check routine.");
              final long oldModularMult = lModularMult;
              timePrimalityTests -= System.currentTimeMillis();
              j = AprtCle(PD[i]);
              timePrimalityTests += System.currentTimeMillis();
              nbrPrimalityTests++;
              primeModMult += lModularMult - oldModularMult;
              if (!batchFinished && batchPrime)
              {
                if (NbrFactors < 2)
                {            // If no factors found in trial factoring...
                  NbrFactors = j;
                }
                return;
              }
            }
            if (j == 0)
            {
              if (Typ[i] < -TYP_EC)
              {
                Typ[i] = -Typ[i]; /* Prime */
              }
              else if (Typ[i] < -TYP_RABIN)
              {
                Typ[i] = TYP_RABIN; /* Prime */
              }
              else if (Typ[i] < -TYP_LEHMAN)
              {
                Typ[i] = TYP_LEHMAN; /* Prime */
              }
              else if (Typ[i] < -TYP_SIQS)
              {
                Typ[i] = TYP_SIQS; /* Prime */
              }
              else if (Typ[i] < -TYP_AURIF)
              {
                Typ[i] = TYP_AURIF; /* Prime */
              }
              else
              {
                Typ[i] = 0; /* Prime */
              }
            }
            else
            {
              if (Typ[i] < -TYP_EC)
              {
                Typ[i] = -TYP_EC - Typ[i]; /* Composite */
              }
              else
              {
                Typ[i] = -Typ[i]; /* Composite */
              }
              if (probableFactor != null)
              { 
                Typ[i] = TYP_RABIN+1;
                InsertNewFactor(probableFactor);
                probableFactor = null;
              } 
            }
            continue factor_loop;
          }
        }
        for (i = 0; i < NbrFactors; i++)
        {
          EC = Typ[i];
          if (EC > 0 && EC < TYP_EC && EC != TYP_AURIF
              && EC != TYP_SIQS && EC != TYP_LEHMAN && EC != TYP_RABIN)
          { /* Composite */
            EC %= 50000000;
            timeECM -= System.currentTimeMillis();
            NN = fnECM(PD[i], i);
            timeECM += System.currentTimeMillis();
            nbrECM++;
            if (NN.equals(BigInt1))
            {
              timeSIQS -= System.currentTimeMillis();
              Siqs SIQS = new Siqs(this, onlyFactoring, NumberLength,
                                   numberThreads);
              NN = SIQS.FactoringSIQS(PD[i]);
              SIQS = null;
              timeSIQS += System.currentTimeMillis();
              nbrSIQS++;
            }
            if (foundByLehman)
            {              // Factor found using Lehman method
              Typ[i] = TYP_LEHMAN + EC + 1;
            }
            else
            {
              Typ[i] = EC;
            }
            InsertNewFactor(NN);
            continue factor_loop;
          }
        }
        break;
      }
      if (onlyFactoring)
      {
        if (!numberIsNegative)
        {
          textAreaContents = upperTextArea.getText() + "\n\n";
          StringToLabel = ""; // Start new line.
          addStringToLabel("Number of divisors: ");
          N1 = BigInt1;
          for (i = 0; i < NbrFactors; i++)
          {
            N1 = N1.multiply(BigInteger.valueOf(Exp[i] + 1));
          }
          insertBigNbr(N1); // Show number of divisors.
          textAreaContents += StringToLabel + "\n\n";
          StringToLabel = ""; // Start new line.
          addStringToLabel("Sum of divisors: ");
          N1 = BigInt1;
          for (i = 0; i < NbrFactors; i++)
          {
            N1 =
              N1.multiply(PD[i].pow(Exp[i] + 1).subtract(BigInt1)).divide(
                PD[i].subtract(BigInt1));
          }
          insertBigNbr(N1); // Show sum of divisors.
          textAreaContents += StringToLabel + "\n\n";
          StringToLabel = ""; // Start new line.
          addStringToLabel("Euler's Totient: ");
          N1 = NumberToFactor;
          for (i = 0; i < NbrFactors; i++)
          {
            N1 = N1.multiply(PD[i].subtract(BigInt1)).divide(PD[i]);
          }
          insertBigNbr(N1);
          j = 1; // Compute Moebius
          for (i = 0; i < NbrFactors; i++)
          {
            if (Exp[i] == 1)
            {
              j = -j;
            }
            else
            {
              j = 0;
            }
          } // Show Euler's totient and Moebius
          textAreaContents += StringToLabel + "\n\nMoebius: " + j;
          StringToLabel = "\n\nSum of squares: "; // Start new line.
          ComputeFourSquares(PD, Exp); // Quad1^2 + Quad2^2 + Quad3^2 + Quad4^2
          NbrFactors1 = NbrFactors;
          if (Quad4.signum() != 0)
          { // Check if four squares are really needed.
            j = NumberToFactor.getLowestSetBit();
            if (j % 2 != 0
              || !NumberToFactor.shiftRight(j).and(BigInteger.valueOf(7)).equals(
                BigInteger.valueOf(7)))
            {
              /* Only three squares are required here */
  
              Computing3Squares = true;
              j = j / 2;
              lowerTextArea.setText("Computing sum of three squares...");
              for (N1 = BigInt1.shiftLeft(j);; N1 = N1.add(BigInt1.shiftLeft(j)))
              {
                if (TerminateThread)
                {
                  throw new ArithmeticException();
                }
                New = System.currentTimeMillis();
                if (OldTimeElapsed >= 0
                  && OldTimeElapsed / 1000 != (OldTimeElapsed + New - Old) / 1000)
                {
                  OldTimeElapsed += New - Old;
                  Old = New;
                  labelStatus.setText(
                    "Time elapsed: "
                      + GetDHMS(OldTimeElapsed / 1000)
                      + "    mod mult: "
                      + (lModularMult >= 0 ? String.valueOf(lModularMult) : "I don't know"));
                }
                N2 = NumberToFactor.subtract(N1.multiply(N1));
                TestComp = GetSmallFactors(N2, PD1, Exp1, Typ1, 1);
                if (TestComp >= 0)
                {
                  if (TestComp == 1)
                  { // Number has all factors < 2^17
                    ComputeFourSquares(PD1, Exp1); // Quad1^2 + Quad2^2
                    break;
                  }
                  if (TestNbr[0] % 4 == 3)
                  {
                    continue;
                  } // This value of c does not work
                  PD1[NbrFactors] = BigIntToBigNbr(TestNbr, NumberLength);
                  Exp1[NbrFactors] = 1;
                  NbrFactors++;
                  if (ComputeFourSquares(PD1, Exp1))
                  { // Quad1^2 + Quad2^2
                    break;
                  }
                }
              } /* end for */
              Quad3 = N1;
              // Sort squares (only Quad3 can be out of order).
              if (Quad1.compareTo(Quad3) < 0)
              {
                Tmp = Quad1;
                Quad1 = Quad3;
                Quad3 = Tmp;
              }
              if (Quad2.compareTo(Quad3) < 0)
              {
                Tmp = Quad2;
                Quad2 = Quad3;
                Quad3 = Tmp;
              }
              Computing3Squares = false;
            }
          }
          NbrFactors = NbrFactors1;
          if (Quad4.signum() == 0)
          {
            if (Quad3.signum() == 0)
            {
              if (Quad2.signum() == 0)
              {
                insertBigNbr(Quad1);
                addStringToLabel("^2");
              }
              else
              {
                textAreaContents += StringToLabel + "a^2 + b^2\n";
                StringToLabel = "a = "; // Start new line.
                insertBigNbr(Quad1);
                textAreaContents += StringToLabel + "\n";
                StringToLabel = "b = "; // Start new line.
                insertBigNbr(Quad2);
              }
            }
            else
            {
              textAreaContents += StringToLabel + "a^2 + b^2 + c^2\n";
              StringToLabel = "a = "; // Start new line.
              insertBigNbr(Quad1);
              textAreaContents += StringToLabel + "\n";
              StringToLabel = "b = "; // Start new line.
              insertBigNbr(Quad2);
              textAreaContents += StringToLabel + "\n";
              StringToLabel = "c = "; // Start new line.
              insertBigNbr(Quad3);
            }
          }
          else
          {
            textAreaContents += StringToLabel + "a^2 + b^2 + c^2 + d^2\n";
            StringToLabel = "a = "; // Start new line.
            insertBigNbr(Quad1);
            textAreaContents += StringToLabel + "\n";
            StringToLabel = "b = "; // Start new line.
            insertBigNbr(Quad2);
            textAreaContents += StringToLabel + "\n";
            StringToLabel = "c = "; // Start new line.
            insertBigNbr(Quad3);
            textAreaContents += StringToLabel + "\n";
            StringToLabel = "d = "; // Start new line.
            insertBigNbr(Quad4);
          }
        }
        upperTextArea.setText(textAreaContents + StringToLabel);
        New = System.currentTimeMillis();
        if (OldTimeElapsed >= 0)
        {
          OldTimeElapsed += New - Old;
          String timeElapsed = GetDHMS(OldTimeElapsed / 1000);
          String textAreaInfo = "Factorization complete in " + timeElapsed;
          if (lModularMult >= 0)
          {
            textAreaInfo += "\nECM: "
              + (lModularMult - primeModMult)
              + " modular multiplications\nPrime checking: "
              + primeModMult
              + " modular multiplications";
          }
          if (smoothsFound > 0)
          {
            textAreaInfo += "\nSIQS: "
              + polynomialsSieved
              + " polynomials sieved\n      "
              + trialDivisions
              + " sets of trial divisions\n      "
              + smoothsFound
              + " smooth congruences found (1 out of every "
              + ValuesSieved/smoothsFound+" values)\n      "
              + totalPartials
              + " partial congruences found (1 out of every "
              + ValuesSieved/totalPartials+" values)\n      "
              + partialsFound
              + " useful partial congruences";
          }
          if (nbrSIQS > 0 || nbrECM > 0 || nbrPrimalityTests > 0)
          {
            textAreaInfo += "\n\nTimings:";
            if (nbrPrimalityTests > 0)
            {
              textAreaInfo += "\nPrimality test of "+nbrPrimalityTests+
                              " number"+(nbrPrimalityTests!=1?"s":"")+
                              ": " + GetDHMSd(timePrimalityTests/100);
            }
            if (nbrECM > 0)
            {
              textAreaInfo += "\nFactoring "+nbrECM+
                              " number"+(nbrECM!=1?"s":"")+
                              " using ECM: " + GetDHMSd(timeECM/100);
            }
            if (nbrSIQS > 0)
            {
              textAreaInfo += "\nFactoring "+nbrSIQS+
                              " number"+(nbrSIQS!=1?"s":"")+
                              " using SIQS: " + GetDHMSd(timeSIQS/100);
            }
          }
          lowerTextArea.setText(textAreaInfo);
          labelStatus.setText(
            "Time elapsed: "
              + timeElapsed
              + "    mod mult: "
              + (lModularMult >= 0 ? String.valueOf(lModularMult) :
                 "I don't know"));
        }
        else
        {
          lowerTextArea.setText("Factorization complete");
        }
        NextEC = -1; /* First curve of new number should be 1 */
      }
    }
    catch (ArithmeticException e)
    {
      New = System.currentTimeMillis();
      if (OldTimeElapsed >= 0)
      {
        OldTimeElapsed += New - Old;
      }
      System.gc();
      return;
    }
    System.gc();
  }

  long GetSmallFactors(
    final BigInteger NumberToFactor,
    BigInteger PD[],
    int Exp[],
    int Typ[],
    final int Type)
  {

    long Div, TestComp;
    int i;
    boolean checkExpParity = false;

    NumberLength = BigNbrToBigInt(NumberToFactor, TestNbr);
    NbrFactors = 0;
    for (i = 0; i < 400; i++)
    {
      Exp[i] = Typ[i] = 0;
    }
    while ((TestNbr[0] & 1) == 0)
    { /* N even */
      if (Exp[NbrFactors] == 0)
      {
        PD[NbrFactors] = BigInt2;
      }
      Exp[NbrFactors]++;
      DivBigNbrByLong(TestNbr, 2, TestNbr, NumberLength);
    }
    if (Exp[NbrFactors] != 0)
    {
      NbrFactors++;
    }
    while (RemDivBigNbrByLong(TestNbr, 3, NumberLength) == 0)
    {
      if (Type == 1)
      {
        checkExpParity ^= true;
      }
      if (Exp[NbrFactors] == 0)
      {
        PD[NbrFactors] = BigInt3;
      }
      Exp[NbrFactors]++;
      DivBigNbrByLong(TestNbr, 3, TestNbr, NumberLength);
    }
    if (checkExpParity)
    {
      return -1; /* Discard it */
    }
    if (Exp[NbrFactors] != 0)
    {
      NbrFactors++;
    }
    Div = 5;
    TestComp = (long)TestNbr[0] + ((long)TestNbr[1] << 31);
    if (TestComp < 0)
    {
      TestComp = 10000 * DosALa31;
    }
    else
    {
      for (i = 2; i < NumberLength; i++)
      {
        if (TestNbr[i] != 0)
        {
          TestComp = 10000 * DosALa31;
          break;
        }
      }
    }
    while (Div < 131072)
    {
      if (Div % 3 != 0)
      {
        while (RemDivBigNbrByLong(TestNbr, Div, NumberLength) == 0)
        {
          if (Type == 1 && Div % 4 == 3)
          {
            checkExpParity ^= true;
          }
          if (Exp[NbrFactors] == 0)
          {
            PD[NbrFactors] = BigInteger.valueOf(Div);
          }
          Exp[NbrFactors]++;
          DivBigNbrByLong(TestNbr, Div, TestNbr, NumberLength);
          TestComp = (long)TestNbr[0] + ((long)TestNbr[1] << 31);
          if (TestComp < 0)
          {
            TestComp = 10000 * DosALa31;
          }
          else
          {
            for (i = 2; i < NumberLength; i++)
            {
              if (TestNbr[i] != 0)
              {
                TestComp = 10000 * DosALa31;
                break;
              }
            }
          } /* end while */
        }
        if (checkExpParity)
        {
          return -1; /* Discard it */
        }
        if (Exp[NbrFactors] != 0)
        {
          NbrFactors++;
        }
      }
      Div += 2;
      if (TestComp < Div * Div && TestComp != 1)
      {
        if (Type == 1 && TestComp % 4 == 3)
        {
          return -1; /* Discard it */
        }
        if (Exp[NbrFactors] != 0)
        {
          NbrFactors++;
        }
        PD[NbrFactors] = BigInteger.valueOf(TestComp);
        Exp[NbrFactors] = 1;
        TestComp = 1;
        NbrFactors++;
        break;
      }
    } /* end while */
    return TestComp;
  }

  int PowerCheck(final int i)
  {
    long New;
    final int maxExpon = (PD[i].bitLength() - 1) / 17;
    int h, j;
    long modulus;
    int intLog2N;
    double log2N;
    BigInteger root, rootN1, rootN, dif, nextroot;
    final int prime2310x1[] =
      { 2311, 4621, 9241, 11551, 18481, 25411, 32341, 34651, 43891, 50821 };
    // Primes of the form 2310x+1.
    boolean expon2 = true, expon3 = true, expon5 = true;
    boolean expon7 = true, expon11 = true;
    for (h = 0; h < prime2310x1.length; h++)
    {
      final long testprime = prime2310x1[h];
      final long mod = PD[i].mod(BigInteger.valueOf(testprime)).intValue();
      if (expon2 && modPow(mod, testprime / 2, testprime) > 1)
        expon2 = false;
      if (expon3 && modPow(mod, testprime / 3, testprime) > 1)
        expon3 = false;
      if (expon5 && modPow(mod, testprime / 5, testprime) > 1)
        expon5 = false;
      if (expon7 && modPow(mod, testprime / 7, testprime) > 1)
        expon7 = false;
      if (expon11 && modPow(mod, testprime / 11, testprime) > 1)
        expon11 = false;
    }
    boolean ProcessExpon[] = new boolean[maxExpon + 1];
    boolean primes[] = new boolean[2 * maxExpon + 3];
    for (h = 2; h <= maxExpon; h++)
    {
      ProcessExpon[h] = true;
    }
    for (h = 2; h < primes.length; h++)
    {
      primes[h] = true;
    }
    for (h = 2; h * h < primes.length; h++)
    { // Generation of primes
      for (j = h * h; j < primes.length; j += h)
      { // using Eratosthenes sieve
        primes[j] = false;
      }
    }
    for (h = 13; h < primes.length; h++)
    {
      if (primes[h])
      {
        int processed = 0;
        for (j = 2 * h + 1; j < primes.length; j += 2 * h)
        {
          if (primes[j])
          {
            modulus = PD[i].mod(BigInteger.valueOf(j)).longValue();
            if (modPow(modulus, j / h, j) > 1)
            {
              for (j = h; j <= maxExpon; j += h)
              {
                ProcessExpon[j] = false;
              }
              break;
            }
          }
          if (++processed > 10)
            break;
        }
      }
    }
    for (int Exponent = maxExpon; Exponent >= 2; Exponent--)
    {
      if (Exponent % 2 == 0 && !expon2)
        continue; // Not a square
      if (Exponent % 3 == 0 && !expon3)
        continue; // Not a cube
      if (Exponent % 5 == 0 && !expon5)
        continue; // Not a fifth power
      if (Exponent % 7 == 0 && !expon7)
        continue; // Not a 7th power
      if (Exponent % 11 == 0 && !expon11)
        continue; // Not an 11th power
      if (!ProcessExpon[Exponent])
        continue;
      New = System.currentTimeMillis();
      if (OldTimeElapsed >= 0
        && OldTimeElapsed / 1000 != (OldTimeElapsed + New - Old) / 1000)
      {
        OldTimeElapsed += New - Old;
        Old = New;
        labelStatus.setText(
          "Time elapsed: "
            + GetDHMS(OldTimeElapsed / 1000)
            + "    Power exponent: "
            + Exponent);
//        Thread.yield();
        if (TerminateThread)
        {
          throw new ArithmeticException();
        }
      }
      intLog2N = PD[i].bitLength() - 1;
      log2N =
        intLog2N
          + Math.log(PD[i].shiftRight(intLog2N - 32).add(BigInt1).doubleValue())
            / Math.log(2)
          - 32;
      log2N /= Exponent;
      if (log2N < 32)
      {
        root = BigInteger.valueOf((long) Math.exp(log2N * Math.log(2)));
      }
      else
      {
        intLog2N = (int) Math.floor(log2N) - 32;
        root =
          BigInteger
            .valueOf((long) Math.exp((log2N - intLog2N) * Math.log(2)) + 10)
            .shiftLeft(intLog2N);
      }
      for (;;)
      {
        rootN1 = root.pow(Exponent - 1);
        rootN = root.multiply(rootN1);
        dif = PD[i].subtract(rootN);
        if (dif.signum() == 0)
        { // Perfect power
          PD[i] = root;
          Exp[i] *= Exponent;
          return 1;
        }
        nextroot =
          dif
            .add(BigInt1)
            .divide(BigInteger.valueOf(Exponent).multiply(rootN1))
            .add(root)
            .subtract(BigInt1);
        if (root.compareTo(nextroot) <= 0)
          break; // Not a perfect power
        root = nextroot;
      }
    }
    return 0;
  }

  // Perform Lehman algorithm
  static BigInteger Lehman(final BigInteger nbr, final int k)
  {
    final long bitsSqr[] = { 0x0000000000000003L, // 3
      0x0000000000000013L, // 5
      0x0000000000000017L, // 7
      0x000000000000023BL, // 11
      0x000000000000161BL, // 13
      0x000000000001A317L, // 17
      0x0000000000030AF3L, // 19
      0x000000000005335FL, // 23
      0x0000000013D122F3L, // 29
      0x00000000121D47B7L, // 31
      0x000000165E211E9BL, // 37
      0x000001B382B50737L, // 41
      0x0000035883A3EE53L, // 43
      0x000004351B2753DFL, // 47
      0x0012DD703303AED3L, // 53
      0x022B62183E7B92BBL, // 59
      0x1713E6940A59F23BL, // 61
    };
    final int primes[] =
      { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 };
    int nbrs[] = new int[17];
    int diffs[] = new int[17];
    int i, j, m;
    int intLog2N;
    double log2N;
    BigInteger root, rootN, dif, nextroot;
    BigInteger bM, a, c, r, sqr, val;
    if (!nbr.testBit(0))
    { // nbr Even
      r = BigInt0;
      m = 1;
      bM = BigInt1;
    }
    else
    {
      if (k % 2 == 0)
      { // k Even
        r = BigInt1;
        m = 2;
        bM = BigInt2;
      }
      else
      { // k Odd
        r = BigInteger.valueOf(k).add(nbr).and(BigInt3);
        m = 4;
        bM = BigInteger.valueOf(4);
      }
    }
    sqr = nbr.multiply(BigInteger.valueOf(k)).shiftLeft(2);
    intLog2N = sqr.bitLength() - 1;
    log2N =
      intLog2N
        + Math.log(sqr.shiftRight(intLog2N - 32).add(BigInt1).doubleValue())
          / Math.log(2)
        - 32;
    log2N /= 2;
    if (log2N < 32)
    {
      root = BigInteger.valueOf((long) Math.exp(log2N * Math.log(2)));
    }
    else
    {
      intLog2N = (int) Math.floor(log2N) - 32;
      root =
        BigInteger
          .valueOf((long) Math.exp((log2N - intLog2N) * Math.log(2)) + 10)
          .shiftLeft(intLog2N);
    }
    for (;;)
    {
      rootN = root.multiply(root);
      dif = sqr.subtract(rootN);
      if (dif.signum() == 0)
      { // Perfect power
        break;
      }
      nextroot =
        dif.add(BigInt1).divide(BigInt2.multiply(root)).add(root).subtract(
          BigInt1);
      if (root.compareTo(nextroot) <= 0)
        break; // Not a perfect power
      root = nextroot;
    }
    a = root;
    while (!a.mod(bM).equals(r) ||
           a.multiply(a).compareTo(sqr)<0)
    {
      a = a.add(BigInt1);
    }
    c = a.multiply(a).subtract(sqr);
    for (i = 0; i < 17; i++)
    {
      final BigInteger prime = BigInteger.valueOf(primes[i]);
      nbrs[i] = c.mod(prime).intValue();
      diffs[i] = bM.multiply(a.shiftLeft(1).add(bM)).mod(prime).intValue();
    }
    for (j = 0; j < 10000; j++)
    {
      for (i = 0; i < 17; i++)
      {
        if ((bitsSqr[i] & (1L << nbrs[i])) == 0)
        { // Not a perfect square
          break;
        }
      }
      if (i == 17)
      { // Test for perfect square
        val = a.add(BigInteger.valueOf(m * j));
        c = val.multiply(val).subtract(sqr);
        intLog2N = c.bitLength() - 1;
        log2N =
          intLog2N
            + Math.log(c.shiftRight(intLog2N - 32).add(BigInt1).doubleValue())
              / Math.log(2)
            - 32;
        log2N /= 2;
        if (log2N < 32)
        {
          root = BigInteger.valueOf((long) Math.exp(log2N * Math.log(2)));
        }
        else
        {
          intLog2N = (int) Math.floor(log2N) - 32;
          root =
            BigInteger
              .valueOf((long) Math.exp((log2N - intLog2N) * Math.log(2)) + 10)
              .shiftLeft(intLog2N);
        }
        for (;;)
        {
          rootN = root.multiply(root);
          dif = c.subtract(rootN);
          if (dif.signum() == 0)
          { // Perfect power -> factor found
            root = nbr.gcd(val.add(root));
            if (root.compareTo(BigInteger.valueOf(10000)) > 0)
            {
              return root; // Return non-trivial found
            }
          }
          nextroot =
            dif.add(BigInt1).divide(BigInt2.multiply(root)).add(root).subtract(
              BigInt1);
          if (root.compareTo(nextroot) <= 0)
            break; // Not a perfect power
          root = nextroot;
        }
      }
      for (i = 0; i < 17; i++)
      {
        nbrs[i] = (nbrs[i] + diffs[i]) % primes[i];
        diffs[i] = (diffs[i] + 2 * m * m) % primes[i];
      }
    }
    return BigInt1; // Factor not found
  }

  final void PowerPM1Check()
  {
    if (onlyFactoring)
    {
      boolean plus1 = false;
      boolean minus1 = false;
      int Exponent = 0;
      int i, j;
      int modulus;
      long i2;
      final int mod9 = NumberToFactor.mod(BigInteger.valueOf(9L)).intValue();
      final int maxExpon = NumberToFactor.bitLength();
      final double logar =
        (maxExpon - 32) * Math.log(2)
          + Math.log(NumberToFactor.shiftRight(maxExpon - 32).longValue());
      boolean ProcessExpon[] = new boolean[maxExpon + 1];
      boolean primes[] = new boolean[2 * maxExpon + 3];
      for (i = 2; i <= maxExpon; i++)
      {
        ProcessExpon[i] = true;
      }
      for (i = 2; i < primes.length; i++)
      {
        primes[i] = true;
      }
      for (i = 2; i * i < primes.length; i++)
      { // Generation of primes
        for (j = i * i; j < primes.length; j += i)
        {
          primes[j] = false;
        }
      }
      // If the number +/- 1 is multiple of a prime but not a multiple
      // of its square then the number +/- 1 cannot be a perfect power.
      for (i = 2; i < primes.length; i++)
      {
        if (primes[i])
        {
          i2 = (long) i * (long) i;
          modulus = NumberToFactor.mod(BigInteger.valueOf(i)).intValue();
          if (modulus == 1
            && NumberToFactor.mod(BigInteger.valueOf(i2)).longValue() != 1L)
          {
            plus1 = true; // NumberFactor cannot be a power + 1
          }
          if (modulus == i - 1
            && NumberToFactor.mod(BigInteger.valueOf(i2)).longValue() != i2 - 1L)
          {
            minus1 = true; // NumberFactor cannot be a power - 1
          }
          if (minus1 && plus1)
            return;
          if (!ProcessExpon[i / 2])
          {
            continue;
          }
          if (modulus > (plus1 ? 1 : 2) && modulus < (minus1 ? i - 1 : i - 2))
          {
            for (j = i / 2; j <= maxExpon; j += i / 2)
            {
              ProcessExpon[j] = false;
            }
          }
          else
          {
            if (modulus == i - 2)
            {
              for (j = i - 1; j <= maxExpon; j += i - 1)
              {
                ProcessExpon[j] = false;
              }
            }
          }
        }
      }
      for (j = 2; j < 100; j++)
      {
        double u = logar / Math.log(j) + .000005;
        Exponent = (int) Math.floor(u);
        if (u - Exponent > .00001)
          continue;
        if (Exponent % 3 == 0 && mod9 > 2 && mod9 < 7)
          continue;
        if (!ProcessExpon[Exponent])
          continue;
        if (ProcessExponent(Exponent))
          return;
      }
      for (; Exponent >= 2; Exponent--)
      {
        if (Exponent % 3 == 0 && mod9 > 2 && mod9 < 7)
          continue;
        if (!ProcessExpon[Exponent])
          continue;
        if (ProcessExponent(Exponent))
          return;
      }
    }
  }

  private boolean ProcessExponent(int Exponent)
  {
    BigInteger NFp1, NFm1, root, rootN1, rootN, rootbak;
    BigInteger nextroot, dif;
    int intLog2N;
    double log2N;
    long New = System.currentTimeMillis();
    if (OldTimeElapsed >= 0
      && OldTimeElapsed / 1000 != (OldTimeElapsed + New - Old) / 1000)
    {
      OldTimeElapsed += New - Old;
      Old = New;
      labelStatus.setText(
        "Time elapsed: "
          + GetDHMS(OldTimeElapsed / 1000)
          + "    Power +/- 1 exponent: "
          + Exponent);
//      Thread.yield();
      if (TerminateThread)
      {
        throw new ArithmeticException();
      }
    }
    NFp1 = NumberToFactor.add(BigInt1);
    NFm1 = NumberToFactor.subtract(BigInt1);
    intLog2N = NFp1.bitLength() - 1;
    log2N =
      intLog2N
        + Math.log(NFp1.shiftRight(intLog2N - 32).add(BigInt1).doubleValue())
          / Math.log(2)
        - 32;
    log2N /= Exponent;
    if (log2N < 32)
    {
      root = BigInteger.valueOf((long) Math.exp(log2N * Math.log(2)));
    }
    else
    {
      intLog2N = (int) Math.floor(log2N) - 32;
      root =
        BigInteger
          .valueOf((long) Math.exp((log2N - intLog2N) * Math.log(2)) + 10)
          .shiftLeft(intLog2N);
    }
    rootbak = root;
    for (;;)
    {
      rootN1 = root.pow(Exponent - 1);
      rootN = root.multiply(rootN1);
      dif = NFp1.subtract(rootN);
      if (dif.signum() == 0)
      { // Perfect power
        Cunningham(root, Exponent, BigInt1.negate(), PD[NbrFactors - 1]);
        return true;
      }
      nextroot =
        dif
          .add(BigInt1)
          .divide(BigInteger.valueOf(Exponent).multiply(rootN1))
          .add(root)
          .subtract(BigInt1);
      if (root.compareTo(nextroot) <= 0)
        break; // Not a perfect power
      root = nextroot;
    }
    root = rootbak;
    for (;;)
    {
      rootN1 = root.pow(Exponent - 1);
      rootN = root.multiply(rootN1);
      dif = NFm1.subtract(rootN);
      if (dif.signum() == 0)
      { // Perfect power
        Cunningham(root, Exponent, BigInt1, PD[NbrFactors - 1]);
        return true;
      }
      nextroot =
        dif
          .add(BigInt1)
          .divide(BigInteger.valueOf(Exponent).multiply(rootN1))
          .add(root)
          .subtract(BigInt1);
      if (root.compareTo(nextroot) <= 0)
        break; // Not a perfect power
      root = nextroot;
    }
    return false;
  }

  final void LucasCheck()
  {
    int i, j;
    if (onlyFactoring)
    {
      if (NumberToFactor.bitLength() > 32)
      {
        int maxExpon = NumberToFactor.bitLength();
        double logar =
          (maxExpon - 32) * Math.log(2)
            + Math.log(NumberToFactor.shiftRight(maxExpon - 32).longValue());
        logar = logar / 0.481211825059603; // index of L
        if (logar + .000005 - Math.floor(logar + .000005) > .00001)
          return;
      }
      BigInteger LucasPrev = BigInteger.valueOf(-1);
      BigInteger LucasAct = BigInt2;
      BigInteger LucasNext;
      i = 0;
      for (;;)
      {
        j = LucasAct.compareTo(NumberToFactor);
        if (j == 0)
        {
          FactorLucas(i, PD[NbrFactors - 1]);
          return;
        }
        if (j > 0)
        {
          return;
        }
        LucasNext = LucasPrev.add(LucasAct);
        LucasPrev = LucasAct;
        LucasAct = LucasNext;
        i++;
      }
    }
  }

  final void FibonacciCheck()
  {
    int i, j;
    if (onlyFactoring)
    {
      if (NumberToFactor.bitLength() > 32)
      {
        int maxExpon = NumberToFactor.bitLength();
        double logar =
          (maxExpon - 32) * Math.log(2)
            + Math.log(NumberToFactor.shiftRight(maxExpon - 32).longValue());
        logar = (logar + 0.80471895621705) / 0.481211825059603; // index of F.
        if (logar + .000005 - Math.floor(logar + .000005) > .00001)
          return;
      }
      BigInteger FibonPrev = BigInt1;
      BigInteger FibonAct = BigInt0;
      BigInteger FibonNext;
      i = 0;
      for (;;)
      {
        j = FibonAct.compareTo(NumberToFactor);
        if (j == 0)
        {
          FactorFibonacci(i, PD[NbrFactors - 1]);
          return;
        }
        if (j > 0)
        {
          return;
        }
        FibonNext = FibonPrev.add(FibonAct);
        FibonPrev = FibonAct;
        FibonAct = FibonNext;
        i++;
      }
    }
  }

  // Prime checking routine
  // Return codes: 0 = Number is prime.
  //               1 = Number is composite.
  int AprtCle(BigInteger N)
  {
    int i, j, G, H, I, J, K, P, Q, T, U, W, X;
    int IV, InvX, LEVELnow, NP, PK, PL, PM, SW, VK, TestedQs, TestingQs;
    int QQ, T1, T3, U1, U3, V1, V3;
    int LengthN, LengthS;
    long Mask;
    double dS;
    String primalityString = "";

    lowerTextArea.setText("Starting Prime Check routine.");
    NumberLength = BigNbrToBigInt(N, TestNbr);
    GetYieldFrequency();
    GetMontgomeryParms();
    if (!Computing3Squares)
    {
      textAreaContents = "";
      StringToLabel = "Testing primality of ";
      insertBigNbr(N);
      addStringToLabel("(" + N.toString().length() + " digits)");
      primalityString =
        textAreaContents + StringToLabel + "\nAPRT-CLE progress: ";
    }
    j = PK = PL = PM = 0;
    for (I = 0; I < NumberLength; I++)
    {
      biS[I] = 0;
      for (J = 0; J < PWmax; J++)
      {
        aiJX[J][I] = 0;
      }
    }
    GetPrimes2Test : for (i = 0; i < LEVELmax; i++)
    {
      biS[0] = 2;
      for (I = 1; I < NumberLength; I++)
      {
        biS[I] = 0;
      }
      for (j = 0; j < aiNQ[i]; j++)
      {
        Q = aiQ[j];
        U = aiT[i] * Q;
        do
        {
          U /= Q;
          MultBigNbrByLong(biS, Q, biS, NumberLength);
        }
        while (U % Q == 0);

        // Exit loop if S^2 > N.

        if (CompareSquare(biS, TestNbr) > 0)
        {
          break GetPrimes2Test;
        }
      } /* End for j */
    } /* End for i */
    if (i == LEVELmax)
    { /* too big */
      return ProbabilisticPrimeTest(N);
    }
    LEVELnow = i;
    TestingQs = j;
    T = aiT[LEVELnow];
    NP = aiNP[LEVELnow];

    MainStart : for (;;)
    {
      for (i = 0; i < NP; i++)
      {
        P = aiP[i];
        SW = TestedQs = 0;
        Q = W = (int) BigNbrModLong(TestNbr, P * P);
        for (J = P - 2; J > 0; J--)
        {
          W = (W * Q) % (P * P);
        }
        if (P > 2 && W != 1)
        {
          SW = 1;
        }
        for (;;)
        {
          for (j = TestedQs; j <= TestingQs; j++)
          {
            Q = aiQ[j] - 1;
            G = aiG[j];
            K = 0;
            while (Q % P == 0)
            {
              K++;
              Q /= P;
            }
            Q = aiQ[j];
            if (K == 0)
            {
              continue;
            }
            if (!Computing3Squares)
            {
              lowerTextArea.setText(
                primalityString
                  + "P = "
                  + P
                  + ",  Q = "
                  + Q
                  + "  ("
                  + (i * (TestingQs + 1) + j) * 100 / (NP * (TestingQs + 1))
                  + "%)");
            }
            PM = 1;
            for (I = 1; I < K; I++)
            {
              PM = PM * P;
            }
            PL = (P - 1) * PM;
            PK = P * PM;
            J = 1;
            for (I = 1; I < Q; I++)
            {
              J = J * G % Q;
              aiIndx[J] = I;
            }
            J = 1;
            for (I = 1; I <= Q - 2; I++)
            {
              J = J * G % Q;
              aiF[I] = aiIndx[(Q + 1 - J) % Q];
            }
            for (I = 0; I < PK; I++)
            {
              for (J = 0; J < NumberLength; J++)
              {
                aiJ0[I][J] = aiJ1[I][J] = 0;
              }
            }
            if (P > 2)
            {
              JacobiSum(1, 1, P, PK, PL, PM, Q);
            }
            else
            {
              if (K != 1)
              {
                JacobiSum(1, 1, P, PK, PL, PM, Q);
                for (I = 0; I < PK; I++)
                {
                  for (J = 0; J < NumberLength; J++)
                  {
                    aiJW[I][J] = 0;
                  }
                }
                if (K != 2)
                {
                  for (I = 0; I < PM; I++)
                  {
                    for (J = 0; J < NumberLength; J++)
                    {
                      aiJW[I][J] = aiJ0[I][J];
                    }
                  }
                  JacobiSum(2, 1, P, PK, PL, PM, Q);
                  for (I = 0; I < PM; I++)
                  {
                    for (J = 0; J < NumberLength; J++)
                    {
                      aiJS[I][J] = aiJ0[I][J];
                    }
                  }
                  JS_JW(PK, PL, PM, P);
                  NormalizeJS(PK, PL, PM, P);
                  for (I = 0; I < PM; I++)
                  {
                    for (J = 0; J < NumberLength; J++)
                    {
                      aiJ1[I][J] = aiJS[I][J];
                    }
                  }
                  JacobiSum(3 << (K - 3), 1 << (K - 3), P, PK, PL, PM, Q);
                  for (J = 0; J < NumberLength; J++)
                  {
                    for (I = 0; I < PK; I++)
                    {
                      aiJW[I][J] = 0;
                    }
                    for (I = 0; I < PM; I++)
                    {
                      aiJS[I][J] = aiJ0[I][J];
                    }
                  }
                  JS_2(PK, PL, PM, P);
                  NormalizeJS(PK, PL, PM, P);
                  for (I = 0; I < PM; I++)
                  {
                    for (J = 0; J < NumberLength; J++)
                    {
                      aiJ2[I][J] = aiJS[I][J];
                    }
                  }
                }
              }
            }
            for (J = 0; J < NumberLength; J++)
            {
              aiJ00[0][J] = aiJ01[0][J] = MontgomeryMultR1[J];
              for (I = 1; I < PK; I++)
              {
                aiJ00[I][J] = aiJ01[I][J] = 0;
              }
            }
            VK = (int) BigNbrModLong(TestNbr, PK);
            for (I = 1; I < PK; I++)
            {
              if (I % P != 0)
              {
                U1 = 1;
                U3 = I;
                V1 = 0;
                V3 = PK;
                while (V3 != 0)
                {
                  QQ = U3 / V3;
                  T1 = U1 - V1 * QQ;
                  T3 = U3 - V3 * QQ;
                  U1 = V1;
                  U3 = V3;
                  V1 = T1;
                  V3 = T3;
                }
                aiInv[I] = (U1 + PK) % PK;
              }
              else
              {
                aiInv[I] = 0;
              }
            }
            if (P != 2)
            {
              for (IV = 0; IV <= 1; IV++)
              {
                for (X = 1; X < PK; X++)
                {
                  for (I = 0; I < PK; I++)
                  {
                    for (J = 0; J < NumberLength; J++)
                    {
                      aiJS[I][J] = aiJ0[I][J];
                    }
                  }
                  if (X % P == 0)
                  {
                    continue;
                  }
                  if (IV == 0)
                  {
                    LongToBigNbr(X, biExp, NumberLength);
                  }
                  else
                  {
                    LongToBigNbr(VK * X / PK, biExp, NumberLength);
                    if (VK * X / PK == 0)
                    {
                      continue;
                    }
                  }
                  JS_E(PK, PL, PM, P);
                  for (I = 0; I < PK; I++)
                  {
                    for (J = 0; J < NumberLength; J++)
                    {
                      aiJW[I][J] = 0;
                    }
                  }
                  InvX = aiInv[X];
                  for (I = 0; I < PK; I++)
                  {
                    J = I * InvX % PK;
                    AddBigNbrModN(aiJW[J], aiJS[I], aiJW[J], TestNbr,
                                  NumberLength);
                  }
                  NormalizeJW(PK, PL, PM, P);
                  if (IV == 0)
                  {
                    for (I = 0; I < PK; I++)
                    {
                      for (J = 0; J < NumberLength; J++)
                      {
                        aiJS[I][J] = aiJ00[I][J];
                      }
                    }
                  }
                  else
                  {
                    for (I = 0; I < PK; I++)
                    {
                      for (J = 0; J < NumberLength; J++)
                      {
                        aiJS[I][J] = aiJ01[I][J];
                      }
                    }
                  }
                  JS_JW(PK, PL, PM, P);
                  if (IV == 0)
                  {
                    for (I = 0; I < PK; I++)
                    {
                      for (J = 0; J < NumberLength; J++)
                      {
                        aiJ00[I][J] = aiJS[I][J];
                      }
                    }
                  }
                  else
                  {
                    for (I = 0; I < PK; I++)
                    {
                      for (J = 0; J < NumberLength; J++)
                      {
                        aiJ01[I][J] = aiJS[I][J];
                      }
                    }
                  }
                } /* end for X */
              } /* end for IV */
            }
            else
            {
              if (K == 1)
              {
                MultBigNbrByLongModN(MontgomeryMultR1, Q, aiJ00[0],
                                     TestNbr, NumberLength);
                for (J = 0; J < NumberLength; J++)
                {
                  aiJ01[0][J] = MontgomeryMultR1[J];
                }
              }
              else
              {
                if (K == 2)
                {
                  if (VK == 1)
                  {
                    for (J = 0; J < NumberLength; J++)
                    {
                      aiJ01[0][J] = MontgomeryMultR1[J];
                    }
                  }
                  for (J = 0; J < NumberLength; J++)
                  {
                    aiJS[0][J] = aiJ0[0][J];
                    aiJS[1][J] = aiJ0[1][J];
                  }
                  JS_2(PK, PL, PM, P);
                  if (VK == 3)
                  {
                    for (J = 0; J < NumberLength; J++)
                    {
                      aiJ01[0][J] = aiJS[0][J];
                      aiJ01[1][J] = aiJS[1][J];
                    }
                  }
                  MultBigNbrByLongModN(aiJS[0], Q, aiJ00[0], TestNbr,
                                       NumberLength);
                  MultBigNbrByLongModN(aiJS[1], Q, aiJ00[1], TestNbr,
                                       NumberLength);
                }
                else
                {
                  for (IV = 0; IV <= 1; IV++)
                  {
                    for (X = 1; X < PK; X += 2)
                    {
                      for (I = 0; I <= PM; I++)
                      {
                        for (J = 0; J < NumberLength; J++)
                        {
                          aiJS[I][J] = aiJ1[I][J];
                        }
                      }
                      if (X % 8 == 5 || X % 8 == 7)
                      {
                        continue;
                      }
                      if (IV == 0)
                      {
                        LongToBigNbr(X, biExp, NumberLength);
                      }
                      else
                      {
                        LongToBigNbr(VK * X / PK, biExp, NumberLength);
                        if (VK * X / PK == 0)
                        {
                          continue;
                        }
                      }
                      JS_E(PK, PL, PM, P);
                      for (I = 0; I < PK; I++)
                      {
                        for (J = 0; J < NumberLength; J++)
                        {
                          aiJW[I][J] = 0;
                        }
                      }
                      InvX = aiInv[X];
                      for (I = 0; I < PK; I++)
                      {
                        J = I * InvX % PK;
                        AddBigNbrModN(aiJW[J], aiJS[I], aiJW[J], TestNbr,
                                      NumberLength);
                      }
                      NormalizeJW(PK, PL, PM, P);
                      if (IV == 0)
                      {
                        for (I = 0; I < PK; I++)
                        {
                          for (J = 0; J < NumberLength; J++)
                          {
                            aiJS[I][J] = aiJ00[I][J];
                          }
                        }
                      }
                      else
                      {
                        for (I = 0; I < PK; I++)
                        {
                          for (J = 0; J < NumberLength; J++)
                          {
                            aiJS[I][J] = aiJ01[I][J];
                          }
                        }
                      }
                      NormalizeJS(PK, PL, PM, P);
                      JS_JW(PK, PL, PM, P);
                      if (IV == 0)
                      {
                        for (I = 0; I < PK; I++)
                        {
                          for (J = 0; J < NumberLength; J++)
                          {
                            aiJ00[I][J] = aiJS[I][J];
                          }
                        }
                      }
                      else
                      {
                        for (I = 0; I < PK; I++)
                        {
                          for (J = 0; J < NumberLength; J++)
                          {
                            aiJ01[I][J] = aiJS[I][J];
                          }
                        }
                      }
                    } /* end for X */
                    if (IV == 0 || VK % 8 == 1 || VK % 8 == 3)
                    {
                      continue;
                    }
                    for (I = 0; I < PM; I++)
                    {
                      for (J = 0; J < NumberLength; J++)
                      {
                        aiJW[I][J] = aiJ2[I][J];
                        aiJS[I][J] = aiJ01[I][J];
                      }
                    }
                    for (; I < PK; I++)
                    {
                      for (J = 0; J < NumberLength; J++)
                      {
                        aiJW[I][J] = aiJS[I][J] = 0;
                      }
                    }
                    JS_JW(PK, PL, PM, P);
                    for (I = 0; I < PM; I++)
                    {
                      for (J = 0; J < NumberLength; J++)
                      {
                        aiJ01[I][J] = aiJS[I][J];
                      }
                    }
                  } /* end for IV */
                }
              }
            }
            for (I = 0; I < PL; I++)
            {
              for (J = 0; J < NumberLength; J++)
              {
                aiJS[I][J] = aiJ00[I][J];
              }
            }
            for (; I < PK; I++)
            {
              for (J = 0; J < NumberLength; J++)
              {
                aiJS[I][J] = 0;
              }
            }
            DivBigNbrByLong(TestNbr, PK, biExp, NumberLength);
            JS_E(PK, PL, PM, P);
            for (I = 0; I < PK; I++)
            {
              for (J = 0; J < NumberLength; J++)
              {
                aiJW[I][J] = 0;
              }
            }
            for (I = 0; I < PL; I++)
            {
              for (J = 0; J < PL; J++)
              {
                MontgomeryMult(aiJS[I], aiJ01[J], biTmp);
                AddBigNbrModN(biTmp, aiJW[(I + J) % PK], aiJW[(I + J) % PK], TestNbr,
                              NumberLength);
              }
            }
            NormalizeJW(PK, PL, PM, P);
            MatchingRoot : do
            {
              H = -1;
              W = 0;
              for (I = 0; I < PL; I++)
              {
                if (!BigNbrIsZero(aiJW[I]))
                {
                  if (H == -1
                    && BigNbrAreEqual(aiJW[I], MontgomeryMultR1))
                  {
                    H = I;
                  }
                  else
                  {
                    H = -2;
                    AddBigNbrModN(aiJW[I], MontgomeryMultR1, biTmp, TestNbr,
                                  NumberLength);
                    if (BigNbrIsZero(biTmp))
                    {
                      W++;
                    }
                  }
                }
              }
              if (H >= 0)
              {
                break MatchingRoot;
              }
              if (W != P - 1)
              {
                return 1; /* Not prime */
              }
              for (I = 0; I < PM; I++)
              {
                AddBigNbrModN(aiJW[I], MontgomeryMultR1, biTmp, TestNbr, NumberLength);
                if (BigNbrIsZero(biTmp))
                {
                  break;
                }
              }
              if (I == PM)
              {
                return 1; /* Not prime */
              }
              for (J = 1; J <= P - 2; J++)
              {
                AddBigNbrModN(aiJW[I + J * PM], MontgomeryMultR1, biTmp, TestNbr,
                              NumberLength);
                if (!BigNbrIsZero(biTmp))
                {
                  return 1; /* Not prime */
                }
              }
              H = I + PL;
            }
            while (false);
            if (SW == 1 || H % P == 0)
            {
              continue;
            }
            if (P != 2)
            {
              SW = 1;
              continue;
            }
            if (K == 1)
            {
              if ((TestNbr[0] & 3) == 1)
              {
                SW = 1;
              }
              continue;
            }

            // if (Q^((N-1)/2) mod N != N-1), N is not prime.

            MultBigNbrByLongModN(MontgomeryMultR1, Q, biTmp, TestNbr,
                                 NumberLength);
            for (I = 0; I < NumberLength; I++)
            {
              biR[I] = biTmp[I];
            }
            I = NumberLength - 1;
            Mask = 0x40000000;
            while ((TestNbr[I] & Mask) == 0)
            {
              Mask /= 2;
              if (Mask == 0)
              {
                I--;
                Mask = 0x40000000L;
              }
            }
            do
            {
              Mask /= 2;
              if (Mask == 0)
              {
                I--;
                Mask = 0x40000000L;
              }
              MontgomeryMult(biR, biR, biT);
              for (J = 0; J < NumberLength; J++)
              {
                biR[J] = biT[J];
              }
              if ((TestNbr[I] & Mask) != 0)
              {
                MontgomeryMult(biR, biTmp, biT);
                for (J = 0; J < NumberLength; J++)
                {
                  biR[J] = biT[J];
                }
              }
            }
            while (I > 0 || Mask > 2);
            AddBigNbrModN(biR, MontgomeryMultR1, biTmp, TestNbr, NumberLength);
            if (!BigNbrIsZero(biTmp))
            {
              return 1; /* Not prime */
            }
            SW = 1;
          } /* end for j */
          if (SW == 0)
          {
            TestedQs = TestingQs + 1;
            if (TestingQs < aiNQ[LEVELnow] - 1)
            {
              TestingQs++;
              Q = aiQ[TestingQs];
              U = T * Q;
              do
              {
                MultBigNbrByLong(biS, Q, biS, NumberLength);
                U /= Q;
              }
              while (U % Q == 0);
              continue; /* Retry */
            }
            LEVELnow++;
            if (LEVELnow == LEVELmax)
            {
              return ProbabilisticPrimeTest(N); /* Cannot tell */
            }
            T = aiT[LEVELnow];
            NP = aiNP[LEVELnow];
            biS[0] = 2;
            for (J = 1; J < NumberLength; J++)
            {
              biS[J] = 0;
            }
            for (J = 0; J <= aiNQ[LEVELnow]; J++)
            {
              Q = aiQ[J];
              U = T * Q;
              do
              {
                MultBigNbrByLong(biS, Q, biS, NumberLength);
                U /= Q;
              }
              while (U % Q == 0);
              if (CompareSquare(biS, TestNbr) > 0)
              {
                TestingQs = J;
                continue MainStart; /* Retry from the beginning */
              }
            } /* end for J */
            return ProbabilisticPrimeTest(N); /* Program error */
          } /* end if */
          break;
        }             /* end for (;;) */
      } /* end for i */

      // Final Test

      LengthN = NumberLength;
      for (I = 0; I < NumberLength; I++)
      {
        biN[I] = TestNbr[I];
        TestNbr[I] = biS[I];
        biR[I] = 0;
      }
      for (;;)
      {
        if (TestNbr[NumberLength - 1] != 0)
        {
          break;
        }
        NumberLength--;
      }
      dN = (double) TestNbr[NumberLength - 1];
      if (NumberLength > 1)
      {
        dN += (double) TestNbr[NumberLength - 2] / dDosALa31;
      }
      if (NumberLength > 2)
      {
        dN += (double) TestNbr[NumberLength - 3] / dDosALa62;
      }
      LengthS = NumberLength;
      dS = dN;
      MontgomeryMultR1[0] = 1;
      for (I = 1; I < NumberLength; I++)
      {
        MontgomeryMultR1[I] = 0;
      }

      biR[0] = 1;
      BigNbrModN(biN, LengthN, biT); /* Compute N mod S */
      for (J = 1; J <= T; J++)
      {
        MultBigNbrModN(biR, biT, biTmp, TestNbr, NumberLength);
        for (i = NumberLength - 1; i > 0; i--)
        {
          if (biTmp[i] != 0)
          {
            break;
          }
        }
        if (i == 0 && biTmp[0] != 1)
        {
          return 0; /* Number is prime */
        }
        for (;;)
        {
          if (biTmp[NumberLength - 1] != 0)
          {
            break;
          }
          NumberLength--;
        }
        for (I = 0; I < NumberLength; I++)
        {
          TestNbr[I] = biTmp[I];
        }
        dN = (double) TestNbr[NumberLength - 1];
        if (NumberLength > 1)
        {
          dN += (double) TestNbr[NumberLength - 2] / dDosALa31;
        }
        if (NumberLength > 2)
        {
          dN += (double) TestNbr[NumberLength - 3] / dDosALa62;
        }
        for (i = NumberLength - 1; i > 0; i--)
        {
          if (TestNbr[i] != biTmp[i])
          {
            break;
          }
        }
        if (TestNbr[i] > biTmp[i])
        {
          BigNbrModN(biN, LengthN, biTmp); /* Compute N mod R */
          if (BigNbrIsZero(biTmp))
          { /* If N is multiple of R.. */
            return 1; /* Number is composite */
          }
        }
        dN = dS;
        NumberLength = LengthS;
        for (I = 0; I < NumberLength; I++)
        {
          biR[I] = TestNbr[I];
          TestNbr[I] = biS[I];
        }
      } /* End for J */
      return 0; /* Number is prime */
    }
  }

  // Prime checking routine
  // Return codes: 0 = Number is prime.
  //               1 = Number is composite.
  int ProbabilisticPrimeTest(BigInteger N)
  {
    long Base, Q;
    int baseNbr, nbrBases, exp, index, j, k;
    long mask;

    lowerTextArea.setText("Starting Rabin probabilistic prime check routine.");
    NumberLength = BigNbrToBigInt(N, TestNbr);
    GetYieldFrequency();
    exp = N.subtract(BigInt1).getLowestSetBit();
    GetMontgomeryParms();
    Base = 1;
    nbrBases = N.bitLength() / 2;
    for (baseNbr = 0; baseNbr < nbrBases; baseNbr++)
    {
      if (Base < 3)
      {
        Base++;
      }
      else
      {
        calculate_new_prime4 : for (;;)
        {
          Base += 2;
          for (Q = 3; Q * Q <= Base; Q += 2)
          { /* Check if Base is prime */
            if (Base % Q == 0)
            {
              continue calculate_new_prime4; /* Composite */
            }
          }
          break; /* Prime found */
        }
      } /* end if */
      lowerTextArea.setText(
        "Rabin probabilistic prime check routine\n\nBase used: "
          + Base
          + " ("
          + baseNbr * 100 / nbrBases
          + "%)");
      System.arraycopy(MontgomeryMultR1, 0, biN, 0, NumberLength);
      index = NumberLength - 1;
      mask = 0x40000000L;
      for (k = NumberLength * 31; k > exp; k--)
      {
        MontgomeryMult(biN, biN, biT);
        if ((TestNbr[index] & mask) != 0)
        {
          MultBigNbrByLongModN(biT, Base, biT, TestNbr, NumberLength);
        }
        System.arraycopy(biT, 0, biN, 0, NumberLength);
        mask >>= 1;
        if (mask == 0)
        {
          index--;
          mask = 0x40000000L;
        }
      }
      for (j = 0; j < NumberLength; j++)
      {
        if (biN[j] != MontgomeryMultR1[j])
        {
          break;
        }
      }
      if (j == NumberLength)
      {
        continue;
      } /* Probable prime, go to next base */
      for (k = 0; k < exp; k++)
      {
        AddBigNbr(biN, MontgomeryMultR1, biT, NumberLength);
          /* If temp is congruent to -1, the number is prob. prime */
        for (j = 0; j < NumberLength; j++)
        {
          if (biT[j] != TestNbr[j])
          {
            break;
          }
        }
        if (j == NumberLength)
        {
          break;
        } /* Probable prime, go to next base */
        MontgomeryMult(biN, biN, biT);
          /* Check whether square equals 1 */
        for (j = 0; j < NumberLength; j++)
        {
          if (biT[j] != MontgomeryMultR1[j])
          {
            break;
          }
        }
        if (j == NumberLength)
        { // Check whether the number can be factored
          // by computing gcd(temp-1, N).
          SubtractBigNbrModN(biN, MontgomeryMultR1, biU, TestNbr, NumberLength);
          GcdBigNbr(biU, TestNbr, biV, NumberLength);
          for (j = 0; j < NumberLength; j++)
          {
            if (biV[j] != BigNbr1[j])
            {
              break;
            }
          }
          if (j < NumberLength)
          {  // biV is a proper factor.
            probableFactor = BigIntToBigNbr(biV, NumberLength);
            return 1; /* Composite number */
          }          
        }
        System.arraycopy(biT, 0, biN, 0, NumberLength);
      }
      if (k == exp)
      {
        return 1; /* Composite number */
      }
    }
    return 0; /* Indicate probable prime */
  }

  public String InitSIQSStrings(BigInteger NbrToFactor, int SieveLimit)
  {
    textAreaContents = "";
    StringToLabel = "Factoring ";
    insertBigNbr(NbrToFactor);
    addStringToLabel("(" + NbrToFactor.toString().length() + " digits)");
    String SIQSInfoText =
      textAreaContents
        + StringToLabel
        + "\n\nSIQS parameters: "
        + nbrPrimes
        + " primes, sieve limit: "
        + SieveLimit;
    lowerTextArea.setText(
      SIQSInfoText + "\nSearching for Knuth-Schroeppel multiplier...");
    return SIQSInfoText;
  }

  public void showMatrixSize(String SIQSInfoText, int rows, int cols)
  {
    lowerTextArea.setText(
             SIQSInfoText +"\nSolving "+rows+"x"+cols+
             " congruence matrix using Block Lanczos algorithm");
  }

  public String getMultAndFactorBase(int multiplier, long FactorBase)
  {
    return "\nMultiplier: " + multiplier + ", factor base: " + FactorBase;
  }

  public void ShowSIQSInfo(long time, int congruencesFound,
                           int matrixBLength, long t)
  {
    if (time > 5 && congruencesFound > 10)
    {
      long u = time * (matrixBLength - congruencesFound) /
                 congruencesFound / 1000;
      labelStatus.setText(
            GetDHMS(t)
            + "     Congruences found: "
            + congruencesFound
            + " ("
            + (int) ((float) (congruencesFound * 100) / (float) matrixBLength)
            + "%)   End sieve in "
            + GetDHMS(u/2));
    }
    else
    {
      labelStatus.setText(
            GetDHMS(t)
            + "     Congruences found: "
            + congruencesFound
            + " ("
            + (int) ((float) (congruencesFound * 100) / (float) matrixBLength)
            + "%)");
    }
  }

  public void saveSIQSStatistics(long polynomialsSieved, long trialDivisions,
                          long smoothsFound, long totalPartials,
                          long partialsFound, long ValuesSieved)
  {
    this.polynomialsSieved = polynomialsSieved;
    this.trialDivisions = trialDivisions;
    this.smoothsFound = smoothsFound;
    this.totalPartials = totalPartials;
    this.partialsFound = partialsFound;
    this.ValuesSieved = ValuesSieved;
  }

  public boolean getTerminateThread()
  {
    return TerminateThread;
  }

  public void writeLowerPane(String data)
  {
    lowerTextArea.setText(data);
  }

  public void setOld(long old)
  {
    Old = old;
  }

  public long getOld()
  {
    return Old;
  }

  public void setOldTimeElapsed(long oldTimeElapsed)
  {
    OldTimeElapsed = oldTimeElapsed;
  }

  public long getOldTimeElapsed()
  {
    return OldTimeElapsed;
  }
  // Compare Nbr1^2 vs. Nbr2
  int CompareSquare(int Nbr1[], int Nbr2[])
  {
    int I, k;

    for (I = NumberLength - 1; I > 0; I--)
    {
      if (Nbr1[I] != 0)
      {
        break;
      }
    }
    k = NumberLength / 2;
    if (NumberLength % 2 == 0)
    {
      if (I >= k)
      {
        return 1;
      } // Nbr1^2 > Nbr2
      if (I < k - 1 || biS[k - 1] < 65536)
      {
        return -1;
      } // Nbr1^2 < Nbr2
    }
    else
    {
      if (I < k)
      {
        return -1;
      } // Nbr1^2 < Nbr2
      if (I > k || biS[k] >= 65536)
      {
        return 1;
      } // Nbr1^2 > Nbr2
    }
    MultBigNbr(biS, biS, biTmp, NumberLength);
    SubtractBigNbr(biTmp, TestNbr, biTmp, NumberLength);
    if (BigNbrIsZero(biTmp))
    {
      return 0;
    } // Nbr1^2 == Nbr2
    if (biTmp[NumberLength - 1] >= 0)
    {
      return 1;
    } // Nbr1^2 > Nbr2
    return -1; // Nbr1^2 < Nbr2
  }

  // Perform JS <- JS ^ E

  void JS_E(int PK, int PL, int PM, int P)
  {
    int I, J, K;
    long Mask;

    for (I = NumberLength - 1; I > 0; I--)
    {
      if (biExp[I] != 0)
      {
        break;
      }
    }
    if (I == 0 && biExp[0] == 1)
    {
      return;
    } // Return if E == 1
    for (K = 0; K < PL; K++)
    {
      for (J = 0; J < NumberLength; J++)
      {
        aiJW[K][J] = aiJS[K][J];
      }
    }
    Mask = 0x40000000L;
    for (;;)
    {
      if ((biExp[I] & Mask) != 0)
      {
        break;
      }
      Mask /= 2;
    }
    do
    {
      JS_2(PK, PL, PM, P);
      Mask /= 2;
      if (Mask == 0)
      {
        Mask = 0x40000000L;
        I--;
      }
      if ((biExp[I] & Mask) != 0)
      {
        JS_JW(PK, PL, PM, P);
      }
    }
    while (I > 0 || Mask != 1);
  }

  // Perform JS <- JS * JW

  void JS_JW(int PK, int PL, int PM, int P)
  {
    int I, J, K;
    for (I = 0; I < PL; I++)
    {
      for (J = 0; J < PL; J++)
      {
        K = (I + J) % PK;
        MontgomeryMult(aiJS[I], aiJW[J], biTmp);
        AddBigNbrModN(aiJX[K], biTmp, aiJX[K], TestNbr, NumberLength);
      }
    }
    for (I = 0; I < PK; I++)
    {
      for (J = 0; J < NumberLength; J++)
      {
        aiJS[I][J] = aiJX[I][J];
        aiJX[I][J] = 0;
      }
    }
    NormalizeJS(PK, PL, PM, P);
  }

  // Perform JS <- JS ^ 2

  void JS_2(int PK, int PL, int PM, int P)
  {
    int I, J, K;
    for (I = 0; I < PL; I++)
    {
      K = 2 * I % PK;
      MontgomeryMult(aiJS[I], aiJS[I], biTmp);
      AddBigNbrModN(aiJX[K], biTmp, aiJX[K], TestNbr, NumberLength);
      AddBigNbrModN(aiJS[I], aiJS[I], biT, TestNbr, NumberLength);
      for (J = I + 1; J < PL; J++)
      {
        K = (I + J) % PK;
        MontgomeryMult(biT, aiJS[J], biTmp);
        AddBigNbrModN(aiJX[K], biTmp, aiJX[K], TestNbr, NumberLength);
      }
    }
    for (I = 0; I < PK; I++)
    {
      for (J = 0; J < NumberLength; J++)
      {
        aiJS[I][J] = aiJX[I][J];
        aiJX[I][J] = 0;
      }
    }
    NormalizeJS(PK, PL, PM, P);
  }

  // Normalize coefficient of JS
  void NormalizeJS(int PK, int PL, int PM, int P)
  {
    int I, J;
    for (I = PL; I < PK; I++)
    {
      if (!BigNbrIsZero(aiJS[I]))
      {
        for (J = 0; J < NumberLength; J++)
        {
          biT[J] = aiJS[I][J];
        }
        for (J = 1; J < P; J++)
        {
          SubtractBigNbrModN(aiJS[I - J * PM], biT, aiJS[I - J * PM],
                             TestNbr, NumberLength);
        }
        for (J = 0; J < NumberLength; J++)
        {
          aiJS[I][J] = 0;
        }
      }
    }
  }

  // Normalize coefficient of JW
  void NormalizeJW(int PK, int PL, int PM, int P)
  {
    int I, J;
    for (I = PL; I < PK; I++)
    {
      if (!BigNbrIsZero(aiJW[I]))
      {
        for (J = 0; J < NumberLength; J++)
        {
          biT[J] = aiJW[I][J];
        }
        for (J = 1; J < P; J++)
        {
          SubtractBigNbrModN(aiJW[I - J * PM], biT, aiJW[I - J * PM],
                             TestNbr, NumberLength);
        }
        for (J = 0; J < NumberLength; J++)
        {
          aiJW[I][J] = 0;
        }
      }
    }
  }

  void JacobiSum(int A, int B, int P, int PK, int PL, int PM, int Q)
  {
    int I, J, K;

    for (I = 0; I < PL; I++)
    {
      for (J = 0; J < NumberLength; J++)
      {
        aiJ0[I][J] = 0;
      }
    }
    for (I = 1; I <= Q - 2; I++)
    {
      J = (A * I + B * aiF[I]) % PK;
      if (J < PL)
      {
        AddBigNbrModN(aiJ0[J], MontgomeryMultR1, aiJ0[J], TestNbr,
                      NumberLength);
      }
      else
      {
        for (K = 1; K < P; K++)
        {
          SubtractBigNbrModN(
            aiJ0[J - K * PM],
            MontgomeryMultR1,
            aiJ0[J - K * PM],
            TestNbr, NumberLength);
        }
      }
    }
  }

  BigInteger fnECM(BigInteger N, int FactorIndex)
  {
    int I, J, Pass, Qaux;
    long L1, L2, LS, P, Q, IP, Paux = 1;
    int[] A0 = new int[NLen];
    int[] A02 = new int[NLen];
    int[] A03 = new int[NLen];
    int[] AA = new int[NLen];
    int[] DX = new int[NLen];
    int[] DZ = new int[NLen];
    int[] GD = new int[NLen];
    int[] M = new int[NLen];
    int[] TX = new int[NLen];
    fieldTX = TX;
    int[] TZ = new int[NLen];
    fieldTZ = TZ;
    int[] UX = new int[NLen];
    fieldUX = UX;
    int[] UZ = new int[NLen];
    fieldUZ = UZ;
    int[] W1 = new int[NLen];
    int[] W2 = new int[NLen];
    int[] W3 = new int[NLen];
    int[] W4 = new int[NLen];
    int[] WX = new int[NLen];
    int[] WZ = new int[NLen];
    int[] X = new int[NLen];
    int[] Z = new int[NLen];
    int[] Aux1 = new int[NLen];
    fieldAux1 = Aux1;
    int[] Aux2 = new int[NLen];
    fieldAux2 = Aux2;
    int[] Aux3 = new int[NLen];
    fieldAux3 = Aux3;
    int[] Aux4 = new int[NLen];
    fieldAux4 = Aux4;
    int[] Xaux = new int[NLen];
    int[] Zaux = new int[NLen];
    int[][] root = new int[480][NLen];
    byte[] sieve = new byte[23100];
    byte[] sieve2310 = new byte[2310];
    int[] sieveidx = new int[480];
    String UpperLine;
    String LowerLine;
    String primalityString;
    int Prob, i, j, u;
    BigInteger NN;

    fieldAA = AA;
    NumberLength = BigNbrToBigInt(N, TestNbr);
    GetYieldFrequency();
    GetMontgomeryParms();
    for (I = 0; I < NumberLength; I++)
    {
      M[I] = DX[I] = DZ[I] = W3[I] = W4[I] = GD[I] = 0;
    }
    textAreaContents = "";
    StringToLabel = "Factoring ";
    insertBigNbr(N);
    addStringToLabel("(" + N.toString().length() + " digits)");
    EC--;
    SmallPrime[0] = 2;
    P = 3;
    indexM = 1;
    for (indexM = 1; indexM < SmallPrime.length; indexM++)
    {
      SmallPrime[indexM] = (int) P; /* Store prime */
      calculate_new_prime1 : for (;;)
      {
        P += 2;
        for (Q = 3; Q * Q <= P; Q += 2)
        { /* Check if P is prime */
          if (P % Q == 0)
          {
            continue calculate_new_prime1; /* Composite */
          }
        }
        break; /* Prime found */
      }
    }
    foundByLehman = false;
    do
    {
      new_curve : for (;;)
      {
        if (NextEC > 0)
        {
          EC = NextEC;
          NextEC = -1;
          if (EC >= TYP_SIQS)
          {
            return BigInt1;
          }
        }
        else
        {
          EC++;
          NN = Lehman(NumberToFactor, EC);
          if (!NN.equals(BigInt1))
          {                // Factor found.
            foundByLehman = true; 
            return NN;
          }
          L1 = N.toString().length();   // Get number of digits.
          if (L1 > 30 && L1 <= 90 &&    // If between 30 and 90 digits...
              (digitsInGroup & 0x400) == 0)
          {                             // Switch to SIQS checkbox is set.
            int limit = limits[((int)L1 - 31) / 5];
            if (EC % 50000000 >= limit && !forcedECM)
            {                           // Switch to SIQS.
              EC += TYP_SIQS;
              return BigInt1;
            }
          }
        }
        Typ[FactorIndex] = EC;
        L1 = 2000;
        L2 = 200000;
        LS = 45;
        Paux = EC;
        nbrPrimes = 303; /* Number of primes less than 2000 */
        if (EC > 25)
        {
          if (EC < 326)
          {
            L1 = 50000;
            L2 = 5000000;
            LS = 224;
            Paux = EC - 24;
            nbrPrimes = 5133; /* Number of primes less than 50000 */
          }
          else
          {
            if (EC < 2000)
            {
              L1 = 1000000;
              L2 = 100000000;
              LS = 1001;
              Paux = EC - 299;
              nbrPrimes = 78498; /* Number of primes less than 1000000 */
            }
            else
            {
              L1 = 11000000;
              L2 = 1100000000;
              LS = 3316;
              Paux = EC - 1900;
              nbrPrimes = 726517; /* Number of primes less than 11000000 */
            }
          }
        }
        primalityString =
          textAreaContents
            + StringToLabel
            + "\nLimit (B1="
            + L1
            + "; B2="
            + L2
            + ")    Curve ";
        UpperLine = "Digits in factor:   ";
        LowerLine = "Probability:        ";
        for (I = 0; I < 6; I++)
        {
          UpperLine += "    >= " + (I * 5 + 15);
          Prob =
            (int) Math.round(
              100
                * (1 - Math.exp(- ((double) L1 * (double) Paux) / ProbArray[I])));
          if (Prob == 100)
          {
            LowerLine += "    100% ";
          }
          else
          {
            if (Prob >= 10)
            {
              LowerLine += "     " + Prob + "% ";
            }
            else
            {
              LowerLine += "      " + Prob + "% ";
            }
          }
        } /* end for */
        lowerTextArea.setText(
          primalityString + EC + "\n" + UpperLine + "\n" + LowerLine);
        LongToBigNbr(2 * (EC + 1), Aux1, NumberLength);
        LongToBigNbr(3 * (EC + 1) * (EC + 1) - 1, Aux2, NumberLength);
        ModInvBigNbr(Aux2, Aux2, TestNbr, NumberLength);
        MultBigNbrModN(Aux1, Aux2, Aux3, TestNbr, NumberLength);
        MultBigNbrModN(Aux3, MontgomeryMultR1, A0, TestNbr, NumberLength);
        MontgomeryMult(A0, A0, A02);
        MontgomeryMult(A02, A0, A03);
        SubtractBigNbrModN(A03, A0, Aux1, TestNbr, NumberLength);
        MultBigNbrByLongModN(A02, 9, Aux2, TestNbr, NumberLength);
        SubtractBigNbrModN(Aux2, MontgomeryMultR1, Aux2, TestNbr, NumberLength);
        MontgomeryMult(Aux1, Aux2, Aux3);
        if (BigNbrIsZero(Aux3))
        {
          continue;
        }
        MultBigNbrByLongModN(A0, 4, Z, TestNbr, NumberLength);
        MultBigNbrByLongModN(A02, 6, Aux1, TestNbr, NumberLength);
        SubtractBigNbrModN(MontgomeryMultR1, Aux1, Aux1, TestNbr, NumberLength);
        MontgomeryMult(A02, A02, Aux2);
        MultBigNbrByLongModN(Aux2, 3, Aux2, TestNbr, NumberLength);
        SubtractBigNbrModN(Aux1, Aux2, Aux1, TestNbr, NumberLength);
        MultBigNbrByLongModN(A03, 4, Aux2, TestNbr, NumberLength);
        ModInvBigNbr(Aux2, Aux2, TestNbr, NumberLength);
        MontgomeryMult(Aux2, MontgomeryMultAfterInv, Aux3);
        MontgomeryMult(Aux1, Aux3, A0);
        AddBigNbrModN(A0, MontgomeryMultR2, Aux1, TestNbr, NumberLength);
        LongToBigNbr(4, Aux2, NumberLength);
        ModInvBigNbr(Aux2, Aux3, TestNbr, NumberLength);
        MultBigNbrModN(Aux3, MontgomeryMultR1, Aux2, TestNbr, NumberLength);
        MontgomeryMult(Aux1, Aux2, AA);
        MultBigNbrByLongModN(A02, 3, Aux1, TestNbr, NumberLength);
        AddBigNbrModN(Aux1, MontgomeryMultR1, X, TestNbr, NumberLength);
        /**************/
        /* First step */
        /**************/
        System.arraycopy(X, 0, Xaux, 0, NumberLength);
        System.arraycopy(Z, 0, Zaux, 0, NumberLength);
        System.arraycopy(MontgomeryMultR1, 0, GcdAccumulated, 0, NumberLength);
        for (Pass = 0; Pass < 2; Pass++)
        {
          /* For powers of 2 */
          indexPrimes = 0;
          StepECM = 1;
          for (I = 1; I <= L1; I <<= 1)
          {
            duplicate(X, Z, X, Z);
          }
          for (I = 3; I <= L1; I *= 3)
          {
            duplicate(W1, W2, X, Z);
            add3(X, Z, X, Z, W1, W2, X, Z);
          }

          if (Pass == 0)
          {
            MontgomeryMult(GcdAccumulated, Z, Aux1);
            System.arraycopy(Aux1, 0, GcdAccumulated, 0, NumberLength);
          }
          else
          {
            GcdBigNbr(Z, TestNbr, GD, NumberLength);
            if (!BigNbrAreEqual(GD, BigNbr1))
            {
              break new_curve;
            }
          }

          /* for powers of odd primes */

          indexM = 1;
          do
          {
            indexPrimes++;
            P = SmallPrime[indexM];
            for (IP = P; IP <= L1; IP *= P)
            {
              prac((int) P, X, Z, W1, W2, W3, W4);
            }
            indexM++;
            if (Pass == 0)
            {
              MontgomeryMult(GcdAccumulated, Z, Aux1);
              System.arraycopy(Aux1, 0, GcdAccumulated, 0, NumberLength);
            }
            else
            {
              GcdBigNbr(Z, TestNbr, GD, NumberLength);
              if (!BigNbrAreEqual(GD, BigNbr1))
              {
                break new_curve;
              }
            }
          }
          while (SmallPrime[indexM - 1] <= LS);
          P += 2;

          /* Initialize sieve2310[n]: 1 if gcd(P+2n,2310) > 1, 0 otherwise */
          u = (int) P;
          for (i = 0; i < 2310; i++)
          {
            sieve2310[i] =
              (u % 3 == 0
                || u % 5 == 0
                || u % 7 == 0
                || u % 11 == 0 ? (byte) 1 : (byte) 0);
            u += 2;
          }
          do
          {
            /* Generate sieve */
            GenerateSieve((int) P, sieve, sieve2310, SmallPrime);

            /* Walk through sieve */

            for (i = 0; i < 23100; i++)
            {
              if (sieve[i] != 0)
                continue; /* Do not process composites */
              if (P + 2 * i > L1)
                break;
              indexPrimes++;
              prac((int) (P + 2 * i), X, Z, W1, W2, W3, W4);
              if (Pass == 0)
              {
                MontgomeryMult(GcdAccumulated, Z, Aux1);
                System.arraycopy(Aux1, 0, GcdAccumulated, 0, NumberLength);
              }
              else
              {
                GcdBigNbr(Z, TestNbr, GD, NumberLength);
                if (!BigNbrAreEqual(GD, BigNbr1))
                {
                  break new_curve;
                }
              }
            }
            P += 46200;
          }
          while (P < L1);
          if (Pass == 0)
          {
            if (BigNbrIsZero(GcdAccumulated))
            { // If GcdAccumulated is
              System.arraycopy(Xaux, 0, X, 0, NumberLength);
              System.arraycopy(Zaux, 0, Z, 0, NumberLength);
              continue; // multiple of TestNbr, continue.
            }
            GcdBigNbr(GcdAccumulated, TestNbr, GD, NumberLength);
            if (!BigNbrAreEqual(GD, BigNbr1))
            {
              break new_curve;
            }
            break;
          }
        } /* end for Pass */

        /******************************************************/
        /* Second step (using improved standard continuation) */
        /******************************************************/
        StepECM = 2;
        j = 0;
        for (u = 1; u < 2310; u += 2)
        {
          if (u % 3 == 0 || u % 5 == 0 || u % 7 == 0 || u % 11 == 0)
          {
            sieve2310[u / 2] = (byte) 1;
          }
          else
          {
            sieve2310[(sieveidx[j++] = u / 2)] = (byte) 0;
          }
        }
        System.arraycopy(sieve2310, 0, sieve2310, 1155, 1155);
        System.arraycopy(X, 0, Xaux, 0, NumberLength); // (X:Z) -> Q (output
        System.arraycopy(Z, 0, Zaux, 0, NumberLength); //         from step 1)
        for (Pass = 0; Pass < 2; Pass++)
        {
          System.arraycopy(
            MontgomeryMultR1,
            0,
            GcdAccumulated,
            0,
            NumberLength);
          System.arraycopy(X, 0, UX, 0, NumberLength);
          System.arraycopy(Z, 0, UZ, 0, NumberLength); // (UX:UZ) -> Q
          ModInvBigNbr(Z, Aux2, TestNbr, NumberLength);
          MontgomeryMult(Aux2, MontgomeryMultAfterInv, Aux1);
          MontgomeryMult(Aux1, X, root[0]); // root[0] <- X/Z (Q)
          J = 0;
          AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, W1);
          SubtractBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, W2);
          MontgomeryMult(W1, W2, TX);
          SubtractBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, AA, Aux2);
          AddBigNbrModN(Aux2, W2, Aux3, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux3, TZ); // (TX:TZ) -> 2Q
          SubtractBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
          AddBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux2, W1);
          AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
          SubtractBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux2, W2);
          AddBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, Aux2);
          MontgomeryMult(Aux2, UZ, X);
          SubtractBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, Aux2);
          MontgomeryMult(Aux2, UX, Z); // (X:Z) -> 3Q
          for (I = 5; I < 2310; I += 2)
          {
            System.arraycopy(X, 0, WX, 0, NumberLength);
            System.arraycopy(Z, 0, WZ, 0, NumberLength);
            SubtractBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
            AddBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
            MontgomeryMult(Aux1, Aux2, W1);
            AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
            SubtractBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
            MontgomeryMult(Aux1, Aux2, W2);
            AddBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
            MontgomeryMult(Aux1, Aux1, Aux2);
            MontgomeryMult(Aux2, UZ, X);
            SubtractBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
            MontgomeryMult(Aux1, Aux1, Aux2);
            MontgomeryMult(Aux2, UX, Z); // (X:Z) -> 5Q, 7Q, ...
            if (Pass == 0)
            {
              MontgomeryMult(GcdAccumulated, Aux1, Aux2);
              System.arraycopy(Aux2, 0, GcdAccumulated, 0, NumberLength);
            }
            else
            {
              GcdBigNbr(Aux1, TestNbr, GD, NumberLength);
              if (!BigNbrAreEqual(GD, BigNbr1))
              {
                break new_curve;
              }
            }
            if (I == 1155)
            {
              System.arraycopy(X, 0, DX, 0, NumberLength);
              System.arraycopy(Z, 0, DZ, 0, NumberLength); // (DX:DZ) -> 1155Q
            }
            if (I % 3 != 0 && I % 5 != 0 && I % 7 != 0 && I % 11 != 0)
            {
              J++;
              ModInvBigNbr(Z, Aux2, TestNbr, NumberLength);
              MontgomeryMult(Aux2, MontgomeryMultAfterInv, Aux1);
              MontgomeryMult(Aux1, X, root[J]); // root[J] <- X/Z
            }
            System.arraycopy(WX, 0, UX, 0, NumberLength); // (UX:UZ) <-
            System.arraycopy(WZ, 0, UZ, 0, NumberLength); // Previous (X:Z)
          } /* end for I */
          AddBigNbrModN(DX, DZ, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, W1);
          SubtractBigNbrModN(DX, DZ, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, W2);
          MontgomeryMult(W1, W2, X);
          SubtractBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, AA, Aux2);
          AddBigNbrModN(Aux2, W2, Aux3, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux3, Z);
          System.arraycopy(X, 0, UX, 0, NumberLength);
          System.arraycopy(Z, 0, UZ, 0, NumberLength); // (UX:UZ) -> 2310Q
          AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, W1);
          SubtractBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, W2);
          MontgomeryMult(W1, W2, TX);
          SubtractBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, AA, Aux2);
          AddBigNbrModN(Aux2, W2, Aux3, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux3, TZ); // (TX:TZ) -> 2*2310Q
          SubtractBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
          AddBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux2, W1);
          AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
          SubtractBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux2, W2);
          AddBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, Aux2);
          MontgomeryMult(Aux2, UZ, X);
          SubtractBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
          MontgomeryMult(Aux1, Aux1, Aux2);
          MontgomeryMult(Aux2, UX, Z); // (X:Z) -> 3*2310Q
          Qaux = (int) (L1 / 4620);
          maxIndexM = (int) (L2 / 4620);
          for (indexM = 0; indexM <= maxIndexM; indexM++)
          {
            if (indexM >= Qaux)
            { // If inside step 2 range... 
              if (indexM == 0)
              {
                ModInvBigNbr(UZ, Aux2, TestNbr, NumberLength);
                MontgomeryMult(Aux2, MontgomeryMultAfterInv, Aux3);
                MontgomeryMult(UX, Aux3, Aux1); // Aux1 <- X/Z (2310Q)
              }
              else
              {
                ModInvBigNbr(Z, Aux2, TestNbr, NumberLength);
                MontgomeryMult(Aux2, MontgomeryMultAfterInv, Aux3);
                MontgomeryMult(X, Aux3, Aux1); // Aux1 <- X/Z (3,5,*
              } //              2310Q)

              /* Generate sieve */
              if (indexM % 10 == 0 || indexM == Qaux)
              {
                GenerateSieve(
                  indexM / 10 * 46200 + 1,
                  sieve,
                  sieve2310,
                  SmallPrime);
              }
              /* Walk through sieve */
              J = 1155 + (indexM % 10) * 2310;
              for (i = 0; i < 480; i++)
              {
                j = sieveidx[i]; // 0 < J < 1155
                if (sieve[J + j] != 0 && sieve[J - 1 - j] != 0)
                {
                  continue; // Do not process if both are composite numbers.
                }
                SubtractBigNbrModN(Aux1, root[i], M, TestNbr, NumberLength);
                MontgomeryMult(GcdAccumulated, M, Aux2);
                System.arraycopy(Aux2, 0, GcdAccumulated, 0, NumberLength);
              }
              if (Pass != 0)
              {
                GcdBigNbr(GcdAccumulated, TestNbr, GD, NumberLength);
                if (!BigNbrAreEqual(GD, BigNbr1))
                {
                  break new_curve;
                }
              }
            }
            if (indexM != 0)
            { // Update (X:Z)
              System.arraycopy(X, 0, WX, 0, NumberLength);
              System.arraycopy(Z, 0, WZ, 0, NumberLength);
              SubtractBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
              AddBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
              MontgomeryMult(Aux1, Aux2, W1);
              AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
              SubtractBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
              MontgomeryMult(Aux1, Aux2, W2);
              AddBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
              MontgomeryMult(Aux1, Aux1, Aux2);
              MontgomeryMult(Aux2, UZ, X);
              SubtractBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
              MontgomeryMult(Aux1, Aux1, Aux2);
              MontgomeryMult(Aux2, UX, Z);
              System.arraycopy(WX, 0, UX, 0, NumberLength);
              System.arraycopy(WZ, 0, UZ, 0, NumberLength);
            }
          } // end for Q
          if (Pass == 0)
          {
            if (BigNbrIsZero(GcdAccumulated))
            { // If GcdAccumulated is
              System.arraycopy(Xaux, 0, X, 0, NumberLength);
              System.arraycopy(Zaux, 0, Z, 0, NumberLength);
              continue; // multiple of TestNbr, continue.
            }
            GcdBigNbr(GcdAccumulated, TestNbr, GD, NumberLength);
            if (BigNbrAreEqual(GD, TestNbr))
            {
              break;
            }
            if (!BigNbrAreEqual(GD, BigNbr1))
            {
              break new_curve;
            }
            break;
          }
        } /* end for Pass */
        performLehman = true;
      }       /* End curve calculation */
    }
    while (BigNbrAreEqual(GD, TestNbr));
    lowerTextArea.setText("");
    StepECM = 0; /* do not show pass number on screen */
    return BigIntToBigNbr(GD, NumberLength);
  }

  /**********************/
  /* Auxiliary routines */
  /**********************/

  static void GenerateSieve(
    int initial,
    byte[] sieve,
    byte[] sieve2310,
    int[] SmallPrime)
  {
    int i, j, Q;
    for (i = 0; i < 23100; i += 2310)
    {
      System.arraycopy(sieve2310, 0, sieve, i, 2310);
    }
    j = 5;
    Q = 13; /* Point to prime 13 */
    do
    {
      if (initial > Q * Q)
      {
        for (i = (int) (((long) initial * ((Q - 1) / 2)) % Q);
          i < 23100;
          i += Q)
        {
          sieve[i] = 1; /* Composite */
        }
      }
      else
      {
        i = Q * Q - initial;
        if (i < 46200)
        {
          for (i = i / 2; i < 23100; i += Q)
          {
            sieve[i] = 1; /* Composite */
          }
        }
        else
        {
          break;
        }
      }
      Q = SmallPrime[++j];
    }
    while (Q < 5000);
  }

  private void GetYieldFrequency()
  {
    yieldFreq = 1000000 / (NumberLength * NumberLength);
    if (yieldFreq > 100000)
      yieldFreq = yieldFreq / 100000 * 100000;
    else if (yieldFreq > 10000)
      yieldFreq = yieldFreq / 10000 * 10000;
    else if (yieldFreq > 1000)
      yieldFreq = yieldFreq / 1000 * 1000;
    else if (yieldFreq > 100)
      yieldFreq = yieldFreq / 100 * 100;
  }

  static int BigNbrToBigInt(BigInteger N, int TestNbr[])
  {
    byte[] Result;
    long[] Temp;
    int i, j, mask;
    long p;
    int NumberLength;

    Result = N.toByteArray();
    NumberLength = (Result.length * 8 + 30)/31;
    Temp = new long[NumberLength+1];
    j = 0;
    mask = 1;
    p = 0;
    for (i = Result.length - 1; i >= 0; i--)
    {
      p += mask * (long) (Result[i] >= 0 ? Result[i] : Result[i] + 256);
      mask <<= 8;
      if (mask == 0)
      {                        // Overflow
        Temp[j++] = p;
        mask = 1;
        p = 0;
      }
    }
    Temp[j] = p;
    Convert32To31Bits(Temp, TestNbr, NumberLength);
    if (TestNbr[NumberLength - 1] > Mi)
    {
      TestNbr[NumberLength] = 0;
      NumberLength++;
    }
    TestNbr[NumberLength] = 0;
    return NumberLength;
  }

  /*********************************************************/
  /* Start of code "borrowed" from Paul Zimmermann's ECM4C */
  /*********************************************************/
  final static int ADD = 6; /* number of multiplications in an addition */
  final static int DUP = 5; /* number of multiplications in a duplicate */

  /* returns the number of modular multiplications */
  static int lucas_cost(int n, double v)
  {
    int c, d, e, r;

    d = n;
    r = (int) ((double) d / v + 0.5);
    if (r >= n)
      return (ADD * n);
    d = n - r;
    e = 2 * r - n;
    c = DUP + ADD; /* initial duplicate and final addition */
    while (d != e)
    {
      if (d < e)
      {
        r = d;
        d = e;
        e = r;
      }
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
      { /* condition 1 */
        r = (2 * d - e) / 3;
        e = (2 * e - d) / 3;
        d = r;
        c += 3 * ADD; /* 3 additions */
      }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
      { /* condition 2 */
        d = (d - e) / 2;
        c += ADD + DUP; /* one addition, one duplicate */
      }
      else if (d <= (4 * e))
      { /* condition 3 */
        d -= e;
        c += ADD; /* one addition */
      }
      else if ((d + e) % 2 == 0)
      { /* condition 4 */
        d = (d - e) / 2;
        c += ADD + DUP; /* one addition, one duplicate */
      }
      else if (d % 2 == 0)
      { /* condition 5 */
        d /= 2;
        c += ADD + DUP; /* one addition, one duplicate */
      }
      else if (d % 3 == 0)
      { /* condition 6 */
        d = d / 3 - e;
        c += 3 * ADD + DUP; /* three additions, one duplicate */
      }
      else if ((d + e) % 3 == 0)
      { /* condition 7 */
        d = (d - 2 * e) / 3;
        c += 3 * ADD + DUP; /* three additions, one duplicate */
      }
      else if ((d - e) % 3 == 0)
      { /* condition 8 */
        d = (d - e) / 3;
        c += 3 * ADD + DUP; /* three additions, one duplicate */
      }
      else if (e % 2 == 0)
      { /* condition 9 */
        e /= 2;
        c += ADD + DUP; /* one addition, one duplicate */
      }
    }
    return (c);
  }

  /* computes nP from P=(x:z) and puts the result in (x:z). Assumes n>2. */
  void prac(
    int n,
    int[] x,
    int[] z,
    int[] xT,
    int[] zT,
    int[] xT2,
    int[] zT2)
  {
    int d, e, r, i;
    int[] t;
    int[] xA = x, zA = z;
    int[] xB = fieldAux1, zB = fieldAux2;
    int[] xC = fieldAux3, zC = fieldAux4;
    double v[] =
      {
        1.61803398875,
        1.72360679775,
        1.618347119656,
        1.617914406529,
        1.612429949509,
        1.632839806089,
        1.620181980807,
        1.580178728295,
        1.617214616534,
        1.38196601125 };

    /* chooses the best value of v */
    r = lucas_cost(n, v[0]);
    i = 0;
    for (d = 1; d < 10; d++)
    {
      e = lucas_cost(n, v[d]);
      if (e < r)
      {
        r = e;
        i = d;
      }
    }
    d = n;
    r = (int) ((double) d / v[i] + 0.5);
    /* first iteration always begins by Condition 3, then a swap */
    d = n - r;
    e = 2 * r - n;
    System.arraycopy(xA, 0, xB, 0, NumberLength); // B = A
    System.arraycopy(zA, 0, zB, 0, NumberLength);
    System.arraycopy(xA, 0, xC, 0, NumberLength); // C = A
    System.arraycopy(zA, 0, zC, 0, NumberLength);
    duplicate(xA, zA, xA, zA); /* A=2*A */
    while (d != e)
    {
      if (d < e)
      {
        r = d;
        d = e;
        e = r;
        t = xA;
        xA = xB;
        xB = t;
        t = zA;
        zA = zB;
        zB = t;
      }
      /* do the first line of Table 4 whose condition qualifies */
      if (4 * d <= 5 * e && ((d + e) % 3) == 0)
      { /* condition 1 */
        r = (2 * d - e) / 3;
        e = (2 * e - d) / 3;
        d = r;
        add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T = f(A,B,C) */
        add3(xT2, zT2, xT, zT, xA, zA, xB, zB); /* T2 = f(T,A,B) */
        add3(xB, zB, xB, zB, xT, zT, xA, zA); /* B = f(B,T,A) */
        t = xA;
        xA = xT2;
        xT2 = t;
        t = zA;
        zA = zT2;
        zT2 = t; /* swap A and T2 */
      }
      else if (4 * d <= 5 * e && (d - e) % 6 == 0)
      { /* condition 2 */
        d = (d - e) / 2;
        add3(xB, zB, xA, zA, xB, zB, xC, zC); /* B = f(A,B,C) */
        duplicate(xA, zA, xA, zA); /* A = 2*A */
      }
      else if (d <= (4 * e))
      { /* condition 3 */
        d -= e;
        add3(xT, zT, xB, zB, xA, zA, xC, zC); /* T = f(B,A,C) */
        t = xB;
        xB = xT;
        xT = xC;
        xC = t;
        t = zB;
        zB = zT;
        zT = zC;
        zC = t; /* circular permutation (B,T,C) */
      }
      else if ((d + e) % 2 == 0)
      { /* condition 4 */
        d = (d - e) / 2;
        add3(xB, zB, xB, zB, xA, zA, xC, zC); /* B = f(B,A,C) */
        duplicate(xA, zA, xA, zA); /* A = 2*A */
      }
      else if (d % 2 == 0)
      { /* condition 5 */
        d /= 2;
        add3(xC, zC, xC, zC, xA, zA, xB, zB); /* C = f(C,A,B) */
        duplicate(xA, zA, xA, zA); /* A = 2*A */
      }
      else if (d % 3 == 0)
      { /* condition 6 */
        d = d / 3 - e;
        duplicate(xT, zT, xA, zA); /* T1 = 2*A */
        add3(xT2, zT2, xA, zA, xB, zB, xC, zC); /* T2 = f(A,B,C) */
        add3(xA, zA, xT, zT, xA, zA, xA, zA); /* A = f(T1,A,A) */
        add3(xT, zT, xT, zT, xT2, zT2, xC, zC); /* T1 = f(T1,T2,C) */
        t = xC;
        xC = xB;
        xB = xT;
        xT = t;
        t = zC;
        zC = zB;
        zB = zT;
        zT = t; /* circular permutation (C,B,T) */
      }
      else if ((d + e) % 3 == 0)
      { /* condition 7 */
        d = (d - 2 * e) / 3;
        add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T1 = f(A,B,C) */
        add3(xB, zB, xT, zT, xA, zA, xB, zB); /* B = f(T1,A,B) */
        duplicate(xT, zT, xA, zA);
        add3(xA, zA, xA, zA, xT, zT, xA, zA); /* A = 3*A */
      }
      else if ((d - e) % 3 == 0)
      { /* condition 8 */
        d = (d - e) / 3;
        add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T1 = f(A,B,C) */
        add3(xC, zC, xC, zC, xA, zA, xB, zB); /* C = f(A,C,B) */
        t = xB;
        xB = xT;
        xT = t;
        t = zB;
        zB = zT;
        zT = t; /* swap B and T */
        duplicate(xT, zT, xA, zA);
        add3(xA, zA, xA, zA, xT, zT, xA, zA); /* A = 3*A */
      }
      else if (e % 2 == 0)
      { /* condition 9 */
        e /= 2;
        add3(xC, zC, xC, zC, xB, zB, xA, zA); /* C = f(C,B,A) */
        duplicate(xB, zB, xB, zB); /* B = 2*B */
      }
    }
    add3(x, z, xA, zA, xB, zB, xC, zC);
  }

  /* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
       using 5/6 mul, 6 add/sub and 6 mod. One assumes that Q-R=P or R-Q=P where P=(x:z).
       Uses the following global variables:
       - n : number to factor
       - x, z : coordinates of P
       - u, v, w : auxiliary variables
  Modifies: x3, z3, u, v, w.
  (x3,z3) may be identical to (x2,z2) and to (x,z)
  */
  void add3(
    int[] x3,
    int[] z3,
    int[] x2,
    int[] z2,
    int[] x1,
    int[] z1,
    int[] x,
    int[] z)
  {
    int[] t = fieldTX;
    int[] u = fieldTZ;
    int[] v = fieldUX;
    int[] w = fieldUZ;
    SubtractBigNbrModN(x2, z2, v, TestNbr, NumberLength); // v = x2-z2
    AddBigNbrModN(x1, z1, w, TestNbr, NumberLength);      // w = x1+z1
    MontgomeryMult(v, w, u);       // u = (x2-z2)*(x1+z1)
    AddBigNbrModN(x2, z2, w, TestNbr, NumberLength);      // w = x2+z2
    SubtractBigNbrModN(x1, z1, t, TestNbr, NumberLength); // t = x1-z1
    MontgomeryMult(t, w, v);       // v = (x2+z2)*(x1-z1)
    AddBigNbrModN(u, v, t, TestNbr, NumberLength);        // t = 2*(x1*x2-z1*z2)
    MontgomeryMult(t, t, w);       // w = 4*(x1*x2-z1*z2)^2
    SubtractBigNbrModN(u, v, t, TestNbr, NumberLength);   // t = 2*(x2*z1-x1*z2)
    MontgomeryMult(t, t, v);       // v = 4*(x2*z1-x1*z2)^2
    if (BigNbrAreEqual(x, x3))
    {
      System.arraycopy(x, 0, u, 0, NumberLength);
      System.arraycopy(w, 0, t, 0, NumberLength);
      MontgomeryMult(z, t, w);
      MontgomeryMult(v, u, z3);
      System.arraycopy(w, 0, x3, 0, NumberLength);
    }
    else
    {
      MontgomeryMult(w, z, x3); // x3 = 4*z*(x1*x2-z1*z2)^2
      MontgomeryMult(x, v, z3); // z3 = 4*x*(x2*z1-x1*z2)^2
    }
  }

  /* computes 2P=(x2:z2) from P=(x1:z1), with 5 mul, 4 add/sub, 5 mod.
       Uses the following global variables:
       - n : number to factor
       - b : (a+2)/4 mod n
       - u, v, w : auxiliary variables
  Modifies: x2, z2, u, v, w
  */
  void duplicate(int[] x2, int[] z2, int[] x1, int[] z1)
  {
    int[] u = fieldUZ;
    int[] v = fieldTX;
    int[] w = fieldTZ;
    AddBigNbrModN(x1, z1, w, TestNbr, NumberLength);      // w = x1+z1
    MontgomeryMult(w, w, u);       // u = (x1+z1)^2
    SubtractBigNbrModN(x1, z1, w, TestNbr, NumberLength); // w = x1-z1
    MontgomeryMult(w, w, v);       // v = (x1-z1)^2
    MontgomeryMult(u, v, x2);      // x2 = u*v = (x1^2 - z1^2)^2
    SubtractBigNbrModN(u, v, w, TestNbr, NumberLength);   // w = u-v = 4*x1*z1
    MontgomeryMult(fieldAA, w, u);
    AddBigNbrModN(u, v, u, TestNbr, NumberLength);        // u = (v+b*w)
    MontgomeryMult(w, u, z2);      // z2 = (w*u)
  }
  /* End of code "borrowed" from Paul Zimmermann's ECM4C */

  static BigInteger BigIntToBigNbr(int[] GD, int NumberLength)
  {
    byte[] Result;
    long[] Temp;
    int i, NL;
    long digit;

    Temp = new long[NumberLength];
    Convert31To32Bits(GD, Temp, NumberLength);
    NL = NumberLength * 4;
    Result = new byte[NL];
    for (i = 0; i < NumberLength; i++)
    {
      digit = Temp[i];
      Result[NL - 1 - 4 * i] = (byte) (digit & 0xFF);
      Result[NL - 2 - 4 * i] = (byte) (digit / 0x100 & 0xFF);
      Result[NL - 3 - 4 * i] = (byte) (digit / 0x10000 & 0xFF);
      Result[NL - 4 - 4 * i] = (byte) (digit / 0x1000000 & 0xFF);
    }
    return (new BigInteger(Result));
  }

  static void LongToBigNbr(long Nbr, int Out[], int NumberLength)
  {
    int i;

    Out[0] = (int)(Nbr & 0x7FFFFFFF);
    Out[1] = (int)((Nbr >> 31) & 0x7FFFFFFF);
    for (i = 2; i < NumberLength; i++)
    {
      Out[i] = (Nbr < 0 ? 0x7FFFFFFF : 0);
    }
  }

  boolean BigNbrIsZero(int Nbr[])
  {
    for (int i = 0; i < NumberLength; i++)
    {
      if (Nbr[i] != 0)
      {
        return false;
      }
    }
    return true;
  }

  boolean BigNbrAreEqual(int Nbr1[], int Nbr2[])
  {
    for (int i = 0; i < NumberLength; i++)
    {
      if (Nbr1[i] != Nbr2[i])
      {
        return false;
      }
    }
    return true;
  }

  static void ChSignBigNbr(int Nbr[], int NumberLength)
  {
    int carry = 0;
    for (int i = 0; i < NumberLength; i++)
    {
      carry = (carry >> 31) - Nbr[i];
      Nbr[i] = carry & 0x7FFFFFFF;
    }
  }

  static void AddBigNbr(int Nbr1[], int Nbr2[], int Sum[],
                        int NumberLength)
  {
    long carry = 0;
    for (int i = 0; i < NumberLength; i++)
    {
      carry = (carry >> 31) + (long)Nbr1[i] + (long)Nbr2[i];
      Sum[i] = (int)(carry & 0x7FFFFFFFL);
    }
  }

  static void SubtractBigNbr(int Nbr1[], int Nbr2[], int Diff[],
                             int NumberLength)
  {
    long carry = 0;
    for (int i = 0; i < NumberLength; i++)
    {
      carry = (carry >> 31) + (long)Nbr1[i] - (long)Nbr2[i];
      Diff[i] = (int)(carry & 0x7FFFFFFFL);
    }
  }

  void AddBigNbr32(long Nbr1[], long Nbr2[], long Sum[], int NumberLength)
  {
    long carry = 0;
    for (int i = 0; i < NumberLength; i++)
    {
      carry = (carry >> 32) + Nbr1[i] + Nbr2[i];
      Sum[i] = carry & 0xFFFFFFFFL;
    }
  }

  void SubtractBigNbr32(long Nbr1[], long Nbr2[], long Diff[], int NumberLength)
  {
    long carry = 0;
    for (int i = 0; i < NumberLength; i++)
    {
      carry = (carry >> 32) + Nbr1[i] - Nbr2[i];
      Diff[i] = carry & 0xFFFFFFFFL;
    }
  }

  static void MultBigNbr(int Nbr1[], int Nbr2[], int Prod[], int NumberLength)
  {
    long MaxUInt = 0x7FFFFFFFL;
    long carry, Pr;
    int i, j;
    carry = Pr = 0;
    for (i = 0; i < NumberLength; i++)
    {
      Pr = carry & MaxUInt;
      carry >>>= 31;
      for (j = 0; j <= i; j++)
      {
        Pr += (long)Nbr1[j] * (long)Nbr2[i - j];
        carry += (Pr >>> 31);
        Pr &= MaxUInt;
      }
      Prod[i] = (int)Pr;
    }
  }

  static void MultBigNbrByLong(int Nbr1[], long Nbr2, int Prod[],
                               int NumberLength)
  {
    long MaxUInt = 0x7FFFFFFFL;
    long Pr;
    int i;
    Pr = 0;
    for (i = 0; i < NumberLength; i++)
    {
      Pr = (Pr >> 31) + Nbr2 * (long)Nbr1[i];
      Prod[i] = (int)(Pr & MaxUInt);
    }
  }

  long BigNbrModLong(int Nbr1[], long Nbr2)
  {
    int i;
    long Rem = 0;

    for (i = NumberLength - 1; i >= 0; i--)
    {
      Rem = ((Rem << 31) + (long)Nbr1[i]) % Nbr2;
    }
    return Rem;
  }

  static void AddBigNbrModN(int Nbr1[], int Nbr2[], int Sum[],
                            int TestNbr[], int NumberLength)
  {
    long MaxUInt = 0x7FFFFFFFL;
    long carry = 0;
    int i;

    for (i = 0; i < NumberLength; i++)
    {
      carry = (carry >> 31) + (long)Nbr1[i] + (long)Nbr2[i] - (long)TestNbr[i];
      Sum[i] = (int)(carry & MaxUInt);
    }
    if (carry < 0)
    {
      carry = 0;
      for (i = 0; i < NumberLength; i++)
      {
        carry = (carry >> 31) + (long)Sum[i] + (long)TestNbr[i];
        Sum[i] = (int)(carry & MaxUInt);
      }
    }
  }

  static void SubtractBigNbrModN(int Nbr1[], int Nbr2[], int Diff[],
                                 int TestNbr[], int NumberLength)
  {
    long MaxUInt = 0x7FFFFFFFL;
    long carry = 0;
    int i;

    for (i = 0; i < NumberLength; i++)
    {
      carry = (carry >> 31) + (long)Nbr1[i] - (long)Nbr2[i];
      Diff[i] = (int)(carry & MaxUInt);
    }
    if (carry < 0)
    {
      carry = 0;
      for (i = 0; i < NumberLength; i++)
      {
        carry = (carry >> 31) + (long)Diff[i] + (long)TestNbr[i];
        Diff[i] = (int)(carry & MaxUInt);
      }
    }
  }

  void MontgomeryMult(int Nbr1[], int Nbr2[], int Prod[])
  {
    long New;
    int NumberLength = this.NumberLength;
    String EcmInfo;

    if (TerminateThread)
    {
      throw new ArithmeticException();
    }
    if (lModularMult >= 0)
    {
      lModularMult++;
      if ((lModularMult % yieldFreq) == 0)
      {
        Thread.yield();
        New = System.currentTimeMillis();
        if (OldTimeElapsed >= 0
          && OldTimeElapsed / 1000 != (OldTimeElapsed + New - Old) / 1000)
        {
          OldTimeElapsed += New - Old;
          Old = New;
          switch (StepECM)
          {
            case 1 :
              EcmInfo = "Step 1: " + (indexPrimes * 100 / nbrPrimes) + "%";
              break;
            case 2 :
              EcmInfo =
                "Step 2: "
                  + (maxIndexM == 0 ? 0 : (indexM * 100 / maxIndexM))
                  + "%";
              break;
            default :
              EcmInfo = "";
          }
          labelStatus.setText(
            "Time elapsed: "
              + GetDHMS(OldTimeElapsed / 1000)
              + "    mod mult: "
              + (lModularMult >= 0 ? String.valueOf(lModularMult) : "I don't know")
              + "   "
              + EcmInfo);
        }
      }
    }
    switch (NumberLength)
    {
      case 2 :
        MontgomeryMult2(Nbr1, Nbr2, Prod);
        break;
      case 3 :
        MontgomeryMult3(Nbr1, Nbr2, Prod);
        break;
      case 4 :
        MontgomeryMult4(Nbr1, Nbr2, Prod);
        break;
      case 5 :
        MontgomeryMult5(Nbr1, Nbr2, Prod);
        break;
      case 6 :
        MontgomeryMult6(Nbr1, Nbr2, Prod);
        break;
      case 7 :
        MontgomeryMult7(Nbr1, Nbr2, Prod);
        break;
      case 8 :
        MontgomeryMult8(Nbr1, Nbr2, Prod);
        break;
      case 9 :
        MontgomeryMult9(Nbr1, Nbr2, Prod);
        break;
      case 10 :
        MontgomeryMult10(Nbr1, Nbr2, Prod);
        break;
      case 11 :
        MontgomeryMult11(Nbr1, Nbr2, Prod);
        break;
      default:
        if (KARATSUBA_ENABLED == false || NumberLength < KARATSUBA_CUTOFF)
        {
          LargeMontgomeryMult(Nbr1, Nbr2, Prod);
        }
        else
        {
          KaratsubaMontgomeryMult(Nbr1, Nbr2, Prod);
        }
        break;
    }
  }

  void MontgomeryMult2(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1;
    Prod0 = Prod1 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    for (int i=0; i<2; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = Pr >>> 31;
    }
    if (Prod1 > TestNbr1 || Prod1 == TestNbr1 && (Prod0 >= TestNbr0))
    {
      Prod0 = (Pr = Prod0 - TestNbr0) & 0x7FFFFFFFL;
      Prod1 = ((Pr >> 31) + Prod1 - TestNbr1) & 0x7FFFFFFFL;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
  }

  void MontgomeryMult3(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1, Prod2;
    Prod0 = Prod1 = Prod2 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long TestNbr2 = TestNbr[2];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    long Nbr2_2 = Nbr2[2];
    for (int i=0; i<3; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr2 + Nbr * Nbr2_2 + Prod2) & 0x7FFFFFFFL;
      Prod2 = Pr >>> 31;
    }
    if (Prod2 > TestNbr2
       || Prod2 == TestNbr2
       && (Prod1 > TestNbr1 || Prod1 == TestNbr1 && (Prod0 >= TestNbr0)))
    {
      Prod0 = (Pr = Prod0 - TestNbr0) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >> 31) + Prod1 - TestNbr1) & 0x7FFFFFFFL;
      Prod2 = ((Pr >> 31) + Prod2 - TestNbr2) & 0x7FFFFFFFL;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
  }

  void MontgomeryMult4(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1, Prod2, Prod3;
    Prod0 = Prod1 = Prod2 = Prod3 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long TestNbr2 = TestNbr[2];
    long TestNbr3 = TestNbr[3];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    long Nbr2_2 = Nbr2[2];
    long Nbr2_3 = Nbr2[3];
    for (int i=0; i<4; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr2 + Nbr * Nbr2_2 + Prod2) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr3 + Nbr * Nbr2_3 + Prod3) & 0x7FFFFFFFL;
      Prod3 = Pr >>> 31;
    }
    if (Prod3 > TestNbr3
      || Prod3 == TestNbr3
      && (Prod2 > TestNbr2
      || Prod2 == TestNbr2
      && (Prod1 > TestNbr1 || Prod1 == TestNbr1 && (Prod0 >= TestNbr0))))
    {
      Prod0 = (Pr = Prod0 - TestNbr0) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >> 31) + Prod1 - TestNbr1) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >> 31) + Prod2 - TestNbr2) & 0x7FFFFFFFL;
      Prod3 = ((Pr >> 31) + Prod3 - TestNbr3) & 0x7FFFFFFFL;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
    Prod[3] = (int)Prod3;
  }

  void MontgomeryMult5(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1, Prod2, Prod3, Prod4;
    Prod0 = Prod1 = Prod2 = Prod3 = Prod4 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long TestNbr2 = TestNbr[2];
    long TestNbr3 = TestNbr[3];
    long TestNbr4 = TestNbr[4];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    long Nbr2_2 = Nbr2[2];
    long Nbr2_3 = Nbr2[3];
    long Nbr2_4 = Nbr2[4];
    for (int i=0; i<5; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr2 + Nbr * Nbr2_2 + Prod2) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr3 + Nbr * Nbr2_3 + Prod3) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr4 + Nbr * Nbr2_4 + Prod4) & 0x7FFFFFFFL;
      Prod4 = Pr >>> 31;
    }
    if (Prod4 > TestNbr4
       || Prod4 == TestNbr4
       && (Prod3 > TestNbr3
       || Prod3 == TestNbr3
       && (Prod2 > TestNbr2
       || Prod2 == TestNbr2
       && (Prod1 > TestNbr1 || Prod1 == TestNbr1 && (Prod0 >= TestNbr0)))))
    {
      Prod0 = (Pr = Prod0 - TestNbr0) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >> 31) + Prod1 - TestNbr1) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >> 31) + Prod2 - TestNbr2) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >> 31) + Prod3 - TestNbr3) & 0x7FFFFFFFL;
      Prod4 = ((Pr >> 31) + Prod4 - TestNbr4) & 0x7FFFFFFFL;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
    Prod[3] = (int)Prod3;
    Prod[4] = (int)Prod4;
  }

  void MontgomeryMult6(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1, Prod2, Prod3, Prod4, Prod5;
    Prod0 = Prod1 = Prod2 = Prod3 = Prod4 = Prod5 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long TestNbr2 = TestNbr[2];
    long TestNbr3 = TestNbr[3];
    long TestNbr4 = TestNbr[4];
    long TestNbr5 = TestNbr[5];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    long Nbr2_2 = Nbr2[2];
    long Nbr2_3 = Nbr2[3];
    long Nbr2_4 = Nbr2[4];
    long Nbr2_5 = Nbr2[5];
    for (int i=0; i<6; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr2 + Nbr * Nbr2_2 + Prod2) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr3 + Nbr * Nbr2_3 + Prod3) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr4 + Nbr * Nbr2_4 + Prod4) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr5 + Nbr * Nbr2_5 + Prod5) & 0x7FFFFFFFL;
      Prod5 = Pr >>> 31;
    }
    if (Prod5 > TestNbr5
       || Prod5 == TestNbr5
       && (Prod4 > TestNbr4
       || Prod4 == TestNbr4
       && (Prod3 > TestNbr3
       || Prod3 == TestNbr3
       && (Prod2 > TestNbr2
       || Prod2 == TestNbr2
       && (Prod1 > TestNbr1
       || Prod1 == TestNbr1
       && (Prod0 >= TestNbr0))))))
    {
      Prod0 = (Pr = Prod0 - TestNbr0) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >> 31) + Prod1 - TestNbr1) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >> 31) + Prod2 - TestNbr2) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >> 31) + Prod3 - TestNbr3) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >> 31) + Prod4 - TestNbr4) & 0x7FFFFFFFL;
      Prod5 = ((Pr >> 31) + Prod5 - TestNbr5) & 0x7FFFFFFFL;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
    Prod[3] = (int)Prod3;
    Prod[4] = (int)Prod4;
    Prod[5] = (int)Prod5;
  }

  void MontgomeryMult7(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr;
    int MontDig;
    int Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6;
    Prod0 = Prod1 = Prod2 = Prod3 = Prod4 = Prod5 = Prod6 = 0;
    int TestNbr0 = (int)TestNbr[0];
    int TestNbr1 = (int)TestNbr[1];
    int TestNbr2 = (int)TestNbr[2];
    int TestNbr3 = (int)TestNbr[3];
    int TestNbr4 = (int)TestNbr[4];
    int TestNbr5 = (int)TestNbr[5];
    int TestNbr6 = (int)TestNbr[6];
    int Nbr2_0 = (int)Nbr2[0];
    int Nbr2_1 = (int)Nbr2[1];
    int Nbr2_2 = (int)Nbr2[2];
    int Nbr2_3 = (int)Nbr2[3];
    int Nbr2_4 = (int)Nbr2[4];
    int Nbr2_5 = (int)Nbr2[5];
    int Nbr2_6 = (int)Nbr2[6];
    int Sum;
    for (int i=0; i<7; i++)
    {
      Pr = (Nbr = Nbr1[i]) * (long)Nbr2_0 + Prod0;
      MontDig = ((int) Pr * (int)MontgomeryMultN) & 0x7FFFFFFF;
      Prod0 = (int)(Pr = ((MontDig * (long)TestNbr0 + Pr) >>> 31) +
              MontDig * (long)TestNbr1 + Nbr * (long)Nbr2_1 + Prod1) & 0x7FFFFFFF;
      Prod1 = (int)(Pr = (Pr >>> 31) +
              MontDig * (long)TestNbr2 + Nbr * (long)Nbr2_2 + Prod2) & 0x7FFFFFFF;
      Prod2 = (int)(Pr = (Pr >>> 31) +
              MontDig * (long)TestNbr3 + Nbr * (long)Nbr2_3 + Prod3) & 0x7FFFFFFF;
      Prod3 = (int)(Pr = (Pr >>> 31) +
              MontDig * (long)TestNbr4 + Nbr * (long)Nbr2_4 + Prod4) & 0x7FFFFFFF;
      Prod4 = (int)(Pr = (Pr >>> 31) +
              MontDig * (long)TestNbr5 + Nbr * (long)Nbr2_5 + Prod5) & 0x7FFFFFFF;
      Prod5 = (int)(Pr = (Pr >>> 31) +
              MontDig * (long)TestNbr6 + Nbr * (long)Nbr2_6 + Prod6) & 0x7FFFFFFF;
      Prod6 = (int)(Pr >>> 31);
    }
    if (Prod6 > TestNbr6
       || Prod6 == TestNbr6
       && (Prod5 > TestNbr5
       || Prod5 == TestNbr5
       && (Prod4 > TestNbr4
       || Prod4 == TestNbr4
       && (Prod3 > TestNbr3
       || Prod3 == TestNbr3
       && (Prod2 > TestNbr2
       || Prod2 == TestNbr2
       && (Prod1 > TestNbr1
       || Prod1 == TestNbr1
       && (Prod0 >= TestNbr0)))))))
    {
      Prod0 = (Sum = Prod0 - TestNbr0) & 0x7FFFFFFF;
      Prod1 = (Sum = (Sum >> 31) + (Prod1 - TestNbr1)) & 0x7FFFFFFF;
      Prod2 = (Sum = (Sum >> 31) + (Prod2 - TestNbr2)) & 0x7FFFFFFF;
      Prod3 = (Sum = (Sum >> 31) + (Prod3 - TestNbr3)) & 0x7FFFFFFF;
      Prod4 = (Sum = (Sum >> 31) + (Prod4 - TestNbr4)) & 0x7FFFFFFF;
      Prod5 = (Sum = (Sum >> 31) + (Prod5 - TestNbr5)) & 0x7FFFFFFF;
      Prod6 = ((Sum >> 31) + (Prod6 - TestNbr6)) & 0x7FFFFFFF;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
    Prod[3] = (int)Prod3;
    Prod[4] = (int)Prod4;
    Prod[5] = (int)Prod5;
    Prod[6] = (int)Prod6;
  }

  void MontgomeryMult8(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6, Prod7;
    Prod0 = Prod1 = Prod2 = Prod3 = Prod4 = Prod5 = Prod6 = Prod7 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long TestNbr2 = TestNbr[2];
    long TestNbr3 = TestNbr[3];
    long TestNbr4 = TestNbr[4];
    long TestNbr5 = TestNbr[5];
    long TestNbr6 = TestNbr[6];
    long TestNbr7 = TestNbr[7];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    long Nbr2_2 = Nbr2[2];
    long Nbr2_3 = Nbr2[3];
    long Nbr2_4 = Nbr2[4];
    long Nbr2_5 = Nbr2[5];
    long Nbr2_6 = Nbr2[6];
    long Nbr2_7 = Nbr2[7];
    for (int i=0; i<8; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr2 + Nbr * Nbr2_2 + Prod2) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr3 + Nbr * Nbr2_3 + Prod3) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr4 + Nbr * Nbr2_4 + Prod4) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr5 + Nbr * Nbr2_5 + Prod5) & 0x7FFFFFFFL;
      Prod5 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr6 + Nbr * Nbr2_6 + Prod6) & 0x7FFFFFFFL;
      Prod6 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr7 + Nbr * Nbr2_7 + Prod7) & 0x7FFFFFFFL;
      Prod7 = Pr >>> 31;
    }
    if (Prod7 > TestNbr7
       || Prod7 == TestNbr7
       && (Prod6 > TestNbr6
       || Prod6 == TestNbr6
       && (Prod5 > TestNbr5
       || Prod5 == TestNbr5
       && (Prod4 > TestNbr4
       || Prod4 == TestNbr4
       && (Prod3 > TestNbr3
       || Prod3 == TestNbr3
       && (Prod2 > TestNbr2
       || Prod2 == TestNbr2
       && (Prod1 > TestNbr1
       || Prod1 == TestNbr1
       && (Prod0 >= TestNbr0))))))))
    {
      Prod[0] = (int)((Pr = Prod0 - TestNbr0) & 0x7FFFFFFFL);
      Prod[1] = (int)((Pr = (Pr >> 31) + Prod1 - TestNbr1) & 0x7FFFFFFFL);
      Prod[2] = (int)((Pr = (Pr >> 31) + Prod2 - TestNbr2) & 0x7FFFFFFFL);
      Prod[3] = (int)((Pr = (Pr >> 31) + Prod3 - TestNbr3) & 0x7FFFFFFFL);
      Prod[4] = (int)((Pr = (Pr >> 31) + Prod4 - TestNbr4) & 0x7FFFFFFFL);
      Prod[5] = (int)((Pr = (Pr >> 31) + Prod5 - TestNbr5) & 0x7FFFFFFFL);
      Prod[6] = (int)((Pr = (Pr >> 31) + Prod6 - TestNbr6) & 0x7FFFFFFFL);
      Prod[7] = (int)(((Pr >> 31) + Prod7 - TestNbr7) & 0x7FFFFFFFL);
      return;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
    Prod[3] = (int)Prod3;
    Prod[4] = (int)Prod4;
    Prod[5] = (int)Prod5;
    Prod[6] = (int)Prod6;
    Prod[7] = (int)Prod7;
  }

  void MontgomeryMult9(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6, Prod7, Prod8;
    Prod0 =
          Prod1 = Prod2 = Prod3 = Prod4 = Prod5 = Prod6 = Prod7 = Prod8 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long TestNbr2 = TestNbr[2];
    long TestNbr3 = TestNbr[3];
    long TestNbr4 = TestNbr[4];
    long TestNbr5 = TestNbr[5];
    long TestNbr6 = TestNbr[6];
    long TestNbr7 = TestNbr[7];
    long TestNbr8 = TestNbr[8];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    long Nbr2_2 = Nbr2[2];
    long Nbr2_3 = Nbr2[3];
    long Nbr2_4 = Nbr2[4];
    long Nbr2_5 = Nbr2[5];
    long Nbr2_6 = Nbr2[6];
    long Nbr2_7 = Nbr2[7];
    long Nbr2_8 = Nbr2[8];
    for (int i=0; i<9; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr2 + Nbr * Nbr2_2 + Prod2) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr3 + Nbr * Nbr2_3 + Prod3) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr4 + Nbr * Nbr2_4 + Prod4) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr5 + Nbr * Nbr2_5 + Prod5) & 0x7FFFFFFFL;
      Prod5 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr6 + Nbr * Nbr2_6 + Prod6) & 0x7FFFFFFFL;
      Prod6 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr7 + Nbr * Nbr2_7 + Prod7) & 0x7FFFFFFFL;
      Prod7 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr8 + Nbr * Nbr2_8 + Prod8) & 0x7FFFFFFFL;
      Prod8 = Pr >>> 31;
    }
    if (Prod8 > TestNbr8
       || Prod8 == TestNbr8
       && (Prod7 > TestNbr7
       || Prod7 == TestNbr7
       && (Prod6 > TestNbr6
       || Prod6 == TestNbr6
       && (Prod5 > TestNbr5
       || Prod5 == TestNbr5
       && (Prod4 > TestNbr4
       || Prod4 == TestNbr4
       && (Prod3 > TestNbr3
       || Prod3 == TestNbr3
       && (Prod2 > TestNbr2
       || Prod2 == TestNbr2
       && (Prod1 > TestNbr1
       || Prod1 == TestNbr1
       && (Prod0 >= TestNbr0)))))))))
    {
      Prod0 = (Pr = Prod0 - TestNbr0) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >> 31) + Prod1 - TestNbr1) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >> 31) + Prod2 - TestNbr2) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >> 31) + Prod3 - TestNbr3) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >> 31) + Prod4 - TestNbr4) & 0x7FFFFFFFL;
      Prod5 = (Pr = (Pr >> 31) + Prod5 - TestNbr5) & 0x7FFFFFFFL;
      Prod6 = (Pr = (Pr >> 31) + Prod6 - TestNbr6) & 0x7FFFFFFFL;
      Prod7 = (Pr = (Pr >> 31) + Prod7 - TestNbr7) & 0x7FFFFFFFL;
      Prod8 = ((Pr >> 31) + Prod8 - TestNbr8) & 0x7FFFFFFFL;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
    Prod[3] = (int)Prod3;
    Prod[4] = (int)Prod4;
    Prod[5] = (int)Prod5;
    Prod[6] = (int)Prod6;
    Prod[7] = (int)Prod7;
    Prod[8] = (int)Prod8;
  }

  void MontgomeryMult10(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6, Prod7, Prod8,
         Prod9;
    Prod0 = Prod1 =
         Prod2 = Prod3 = Prod4 = Prod5 = Prod6 = Prod7 = Prod8 = Prod9 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long TestNbr2 = TestNbr[2];
    long TestNbr3 = TestNbr[3];
    long TestNbr4 = TestNbr[4];
    long TestNbr5 = TestNbr[5];
    long TestNbr6 = TestNbr[6];
    long TestNbr7 = TestNbr[7];
    long TestNbr8 = TestNbr[8];
    long TestNbr9 = TestNbr[9];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    long Nbr2_2 = Nbr2[2];
    long Nbr2_3 = Nbr2[3];
    long Nbr2_4 = Nbr2[4];
    long Nbr2_5 = Nbr2[5];
    long Nbr2_6 = Nbr2[6];
    long Nbr2_7 = Nbr2[7];
    long Nbr2_8 = Nbr2[8];
    long Nbr2_9 = Nbr2[9];
    for (int i=0; i<10; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr2 + Nbr * Nbr2_2 + Prod2) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr3 + Nbr * Nbr2_3 + Prod3) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr4 + Nbr * Nbr2_4 + Prod4) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr5 + Nbr * Nbr2_5 + Prod5) & 0x7FFFFFFFL;
      Prod5 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr6 + Nbr * Nbr2_6 + Prod6) & 0x7FFFFFFFL;
      Prod6 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr7 + Nbr * Nbr2_7 + Prod7) & 0x7FFFFFFFL;
      Prod7 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr8 + Nbr * Nbr2_8 + Prod8) & 0x7FFFFFFFL;
      Prod8 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr9 + Nbr * Nbr2_9 + Prod9) & 0x7FFFFFFFL;
      Prod9 = Pr >>> 31;
    }
    if (Prod9 > TestNbr9
       || Prod9 == TestNbr9
       && (Prod8 > TestNbr8
       || Prod8 == TestNbr8
       && (Prod7 > TestNbr7
       || Prod7 == TestNbr7
       && (Prod6 > TestNbr6
       || Prod6 == TestNbr6
       && (Prod5 > TestNbr5
       || Prod5 == TestNbr5
       && (Prod4 > TestNbr4
       || Prod4 == TestNbr4
       && (Prod3 > TestNbr3
       || Prod3 == TestNbr3
       && (Prod2 > TestNbr2
       || Prod2 == TestNbr2
       && (Prod1 > TestNbr1
       || Prod1 == TestNbr1
       && (Prod0 >= TestNbr0))))))))))
    {
      Prod0 = (Pr = Prod0 - TestNbr0) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >> 31) + Prod1 - TestNbr1) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >> 31) + Prod2 - TestNbr2) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >> 31) + Prod3 - TestNbr3) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >> 31) + Prod4 - TestNbr4) & 0x7FFFFFFFL;
      Prod5 = (Pr = (Pr >> 31) + Prod5 - TestNbr5) & 0x7FFFFFFFL;
      Prod6 = (Pr = (Pr >> 31) + Prod6 - TestNbr6) & 0x7FFFFFFFL;
      Prod7 = (Pr = (Pr >> 31) + Prod7 - TestNbr7) & 0x7FFFFFFFL;
      Prod8 = (Pr = (Pr >> 31) + Prod8 - TestNbr8) & 0x7FFFFFFFL;
      Prod9 = ((Pr >> 31) + Prod9 - TestNbr9) & 0x7FFFFFFFL;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
    Prod[3] = (int)Prod3;
    Prod[4] = (int)Prod4;
    Prod[5] = (int)Prod5;
    Prod[6] = (int)Prod6;
    Prod[7] = (int)Prod7;
    Prod[8] = (int)Prod8;
    Prod[9] = (int)Prod9;
  }

  void MontgomeryMult11(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6, Prod7, Prod8,
         Prod9, Prod10;
    Prod0 = Prod1 = Prod2 = Prod3 =
            Prod4 = Prod5 = Prod6 = Prod7 = Prod8 = Prod9 = Prod10 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long TestNbr2 = TestNbr[2];
    long TestNbr3 = TestNbr[3];
    long TestNbr4 = TestNbr[4];
    long TestNbr5 = TestNbr[5];
    long TestNbr6 = TestNbr[6];
    long TestNbr7 = TestNbr[7];
    long TestNbr8 = TestNbr[8];
    long TestNbr9 = TestNbr[9];
    long TestNbr10 = TestNbr[10];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    long Nbr2_2 = Nbr2[2];
    long Nbr2_3 = Nbr2[3];
    long Nbr2_4 = Nbr2[4];
    long Nbr2_5 = Nbr2[5];
    long Nbr2_6 = Nbr2[6];
    long Nbr2_7 = Nbr2[7];
    long Nbr2_8 = Nbr2[8];
    long Nbr2_9 = Nbr2[9];
    long Nbr2_10 = Nbr2[10];
    for (int i=0; i<11; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr2 + Nbr * Nbr2_2 + Prod2) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr3 + Nbr * Nbr2_3 + Prod3) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr4 + Nbr * Nbr2_4 + Prod4) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr5 + Nbr * Nbr2_5 + Prod5) & 0x7FFFFFFFL;
      Prod5 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr6 + Nbr * Nbr2_6 + Prod6) & 0x7FFFFFFFL;
      Prod6 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr7 + Nbr * Nbr2_7 + Prod7) & 0x7FFFFFFFL;
      Prod7 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr8 + Nbr * Nbr2_8 + Prod8) & 0x7FFFFFFFL;
      Prod8 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr9 + Nbr * Nbr2_9 + Prod9) & 0x7FFFFFFFL;
      Prod9 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr10 + Nbr * Nbr2_10 + Prod10) & 0x7FFFFFFFL;
      Prod10 = Pr >>> 31;
    }
    if (Prod10 > TestNbr10
       || Prod10 == TestNbr10
       && (Prod9 > TestNbr9
       || Prod9 == TestNbr9
       && (Prod8 > TestNbr8
       || Prod8 == TestNbr8
       && (Prod7 > TestNbr7
       || Prod7 == TestNbr7
       && (Prod6 > TestNbr6
       || Prod6 == TestNbr6
       && (Prod5 > TestNbr5
       || Prod5 == TestNbr5
       && (Prod4 > TestNbr4
       || Prod4 == TestNbr4
       && (Prod3 > TestNbr3
       || Prod3 == TestNbr3
       && (Prod2 > TestNbr2
       || Prod2 == TestNbr2
       && (Prod1 > TestNbr1
       || Prod1 == TestNbr1
       && (Prod0 >= TestNbr0)))))))))))
    {
      Prod0 = (Pr = Prod0 - TestNbr0) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >> 31) + Prod1 - TestNbr1) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >> 31) + Prod2 - TestNbr2) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >> 31) + Prod3 - TestNbr3) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >> 31) + Prod4 - TestNbr4) & 0x7FFFFFFFL;
      Prod5 = (Pr = (Pr >> 31) + Prod5 - TestNbr5) & 0x7FFFFFFFL;
      Prod6 = (Pr = (Pr >> 31) + Prod6 - TestNbr6) & 0x7FFFFFFFL;
      Prod7 = (Pr = (Pr >> 31) + Prod7 - TestNbr7) & 0x7FFFFFFFL;
      Prod8 = (Pr = (Pr >> 31) + Prod8 - TestNbr8) & 0x7FFFFFFFL;
      Prod9 = (Pr = (Pr >> 31) + Prod9 - TestNbr9) & 0x7FFFFFFFL;
      Prod10 = ((Pr >> 31) + Prod10 - TestNbr10) & 0x7FFFFFFFL;
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
    Prod[3] = (int)Prod3;
    Prod[4] = (int)Prod4;
    Prod[5] = (int)Prod5;
    Prod[6] = (int)Prod6;
    Prod[7] = (int)Prod7;
    Prod[8] = (int)Prod8;
    Prod[9] = (int)Prod9;
    Prod[10] = (int)Prod10;
  }  

  void LargeMontgomeryMult(int Nbr1[], int Nbr2[], int Prod[])
  {
    long Pr, Nbr, MontDig;
    long Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6, Prod7, Prod8,
         Prod9, Prod10;
    Prod0 = Prod1 = Prod2 = Prod3 =
            Prod4 = Prod5 = Prod6 = Prod7 = Prod8 = Prod9 = Prod10 = 0;
    long TestNbr0 = TestNbr[0];
    long TestNbr1 = TestNbr[1];
    long TestNbr2 = TestNbr[2];
    long TestNbr3 = TestNbr[3];
    long TestNbr4 = TestNbr[4];
    long TestNbr5 = TestNbr[5];
    long TestNbr6 = TestNbr[6];
    long TestNbr7 = TestNbr[7];
    long TestNbr8 = TestNbr[8];
    long TestNbr9 = TestNbr[9];
    long TestNbr10 = TestNbr[10];
    long Nbr2_0 = Nbr2[0];
    long Nbr2_1 = Nbr2[1];
    long Nbr2_2 = Nbr2[2];
    long Nbr2_3 = Nbr2[3];
    long Nbr2_4 = Nbr2[4];
    long Nbr2_5 = Nbr2[5];
    long Nbr2_6 = Nbr2[6];
    long Nbr2_7 = Nbr2[7];
    long Nbr2_8 = Nbr2[8];
    long Nbr2_9 = Nbr2[9];
    long Nbr2_10 = Nbr2[10];
    int j;
    for (j = 11; j < NumberLength; j++)
    {
      Prod[j] = 0;
    }
    for (int i=0; i<NumberLength; i++)
    {
      Pr = (Nbr = Nbr1[i]) * Nbr2_0 + Prod0;
      MontDig = ((int) Pr * MontgomeryMultN) & 0x7FFFFFFFL;
      Prod0 = (Pr = ((MontDig * TestNbr0 + Pr) >>> 31) +
              MontDig * TestNbr1 + Nbr * Nbr2_1 + Prod1) & 0x7FFFFFFFL;
      Prod1 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr2 + Nbr * Nbr2_2 + Prod2) & 0x7FFFFFFFL;
      Prod2 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr3 + Nbr * Nbr2_3 + Prod3) & 0x7FFFFFFFL;
      Prod3 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr4 + Nbr * Nbr2_4 + Prod4) & 0x7FFFFFFFL;
      Prod4 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr5 + Nbr * Nbr2_5 + Prod5) & 0x7FFFFFFFL;
      Prod5 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr6 + Nbr * Nbr2_6 + Prod6) & 0x7FFFFFFFL;
      Prod6 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr7 + Nbr * Nbr2_7 + Prod7) & 0x7FFFFFFFL;
      Prod7 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr8 + Nbr * Nbr2_8 + Prod8) & 0x7FFFFFFFL;
      Prod8 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr9 + Nbr * Nbr2_9 + Prod9) & 0x7FFFFFFFL;
      Prod9 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr10 + Nbr * Nbr2_10 + Prod10) & 0x7FFFFFFFL;
      Prod10 = (Pr = (Pr >>> 31) +
              MontDig * TestNbr[11] + Nbr * Nbr2[11] + Prod[11]) & 0x7FFFFFFFL;
      for (j = 12; j < NumberLength; j++)
      {
        Prod[j-1] = (int)((Pr = (Pr >>> 31) +
                 MontDig * TestNbr[j] + Nbr * Nbr2[j] + Prod[j]) & 0x7FFFFFFFL);
      }
      Prod[j - 1] = (int)(Pr >>> 31);
    }
    Prod[0] = (int)Prod0;
    Prod[1] = (int)Prod1;
    Prod[2] = (int)Prod2;
    Prod[3] = (int)Prod3;
    Prod[4] = (int)Prod4;
    Prod[5] = (int)Prod5;
    Prod[6] = (int)Prod6;
    Prod[7] = (int)Prod7;
    Prod[8] = (int)Prod8;
    Prod[9] = (int)Prod9;
    Prod[10] = (int)Prod10;
    for (j = NumberLength - 1; j >= 0; j--)
    {
      if (Prod[j] != TestNbr[j])
      {
        break;
      }
    } 
    if (j < 0 || j >= 0 && Prod[j] >= TestNbr[j])
    {
      Pr = 0;
      for (j = 0; j < NumberLength; j++)
      {
        Prod[j] = (int)((Pr = (Pr >> 31) +
                             (long)Prod[j] - (long)TestNbr[j]) & 0x7FFFFFFFL);
      }
    }
  }

          // Algorithm REDC: x is the product of both numbers
          //  m = (x mod R) * N' mod R
          //  t = (x + m*N) / R
          //  if t < N
          //     return t
          //  else return t - N
  void KaratsubaMontgomeryMult(int Nbr1[], int Nbr2[], int Prod[])
  {
    int carry = 0;
    int i;
    for (i=0; i<karatLength; i++)
    {
      arrNbr[i] = Nbr1[i];
      arrNbr[karatLength+i] = (int)Nbr2[i];
    }
    KaratsubaMultiply(0, karatLength, 2*karatLength);        // Get x.
    System.arraycopy(arrNbr, 0,                              // Source
                     arrNbrM, 0,                             // Destination
                     2*karatLength);                         // Length
    System.arraycopy(montgKaratsubaArr, 0,                   // Source
                     arrNbr, 0,                              // Destination
                     karatLength);                           // Length
    KaratsubaMultiply(0, karatLength, 2*karatLength);        // Get m.
    System.arraycopy(arrNbr, 0,                              // Source
                     arrNbr, 0,                              // Destination
                     karatLength);                           // Length
    System.arraycopy(KaratsubaTestNbr, 0,                    // Source
                     arrNbr, 0,                              // Destination
                     karatLength);                           // Length
    KaratsubaMultiply(0, karatLength, 2*karatLength);        // Get m*N.
    for (i=0; i<2*karatLength; i++)
    {
      carry += arrNbr[i] + arrNbrM[i];
      if (carry < 0)
      {                         // Bit 31 is set.
        arrNbr[i] = carry - (-0x80000000);
        carry = 1;
      }
      else
      {
        arrNbr[i] = carry;
        carry = 0;
      }
    }                                                        // t is got.
    if (carry==0)
    {
      for (i=2*karatLength-1; i>=karatLength; i--)
      {
        if (arrNbr[i] > KaratsubaTestNbr[i])
        {
          carry = 1;
          break;
        }
        if (arrNbr[i] < KaratsubaTestNbr[i])
        {
          break;
        }
      }
    }
    if (carry != 0)
    {
      carry = 0;
      for (i=2*karatLength-1; i>=karatLength; i--)
      {
        carry += arrNbr[i] - KaratsubaTestNbr[i];
        if (carry < 0)
        {                         // Bit 31 is set.
          arrNbr[i] = carry - (-0x80000000);
        }
        else
        {
          arrNbr[i] = carry;
        }
      }
    }
    for (i=karatLength; i<2*karatLength; i++)
    {
      Prod[i-karatLength] = arrNbr[i];
    }
  }

  void GetMontgomeryParms()
  {
    int NumberLength = this.NumberLength;
    int N, x, i, j, k, div;
    int length;

    dN = (double) TestNbr[NumberLength - 1];
    if (NumberLength > 1)
    {
      dN += (double) TestNbr[NumberLength - 2] / dDosALa31;
    }
    if (NumberLength > 2)
    {
      dN += (double) TestNbr[NumberLength - 3] / dDosALa62;
    }

    x = N = (int) TestNbr[0]; // 2 least significant bits of inverse correct.
    x = x * (2 - N * x); // 4 least significant bits of inverse correct.
    x = x * (2 - N * x); // 8 least significant bits of inverse correct.
    x = x * (2 - N * x); // 16 least significant bits of inverse correct.
    x = x * (2 - N * x); // 32 least significant bits of inverse correct.
    MontgomeryMultN = (-x) & 0x7FFFFFFF;
    if (KARATSUBA_ENABLED && NumberLength >= KARATSUBA_CUTOFF)
    {
      length = NumberLength;
      div = 1;
      while (length > KARATSUBA_CUTOFF)
      {
        div *= 2;
        length = (length+1)/2;
      }
      karatLength = length * div;
      for (i=NumberLength; i<karatLength; i++)
      {
        TestNbr[i]=0;
      }
      NumberLength = karatLength;
      for (i=0; i<NumberLength; i++)
      {
        KaratsubaTestNbr[i] = (int)TestNbr[i];
      }
      montgKaratsubaArr[0] = (int)MontgomeryMultN;
      j=1;
      for (i=2; i<=length; i<<=1)
      {
          // Compute x <- x * (2 + N * x)
        System.arraycopy(KaratsubaTestNbr, 0, arrNbr, 0, j);
        System.arraycopy(montgKaratsubaArr, 0, arrNbr, j, j);
        NormalMultiply(0, j);
        if ((arrNbrAux[0] += 2) < 0)
        {
          arrNbrAux[0] &= 0x7FFFFFFF;
          for (k=1; k<i; k++)
          {
            if (++arrNbrAux[k] >= 0)
            {
              break;
            }
            arrNbrAux[k] = 0;
          }
        }
        System.arraycopy(arrNbrAux, 0, arrNbr, 0, i);
        System.arraycopy(montgKaratsubaArr, 0, arrNbr, i, i);
        NormalMultiply(0, i);
        System.arraycopy(arrNbrAux, 0, montgKaratsubaArr, 0, i);
        j<<=1;
      }
    }
    j = NumberLength;
    MontgomeryMultR1[j] = 1;
    do
    {
      MontgomeryMultR1[--j] = 0;
    }
    while (j > 0);
    AdjustModN(MontgomeryMultR1, TestNbr, NumberLength);
    MultBigNbrModN(MontgomeryMultR1, MontgomeryMultR1, MontgomeryMultR2,
                   TestNbr, NumberLength);
    MontgomeryMult(MontgomeryMultR2, MontgomeryMultR2, MontgomeryMultAfterInv);
    AddBigNbrModN(MontgomeryMultR1, MontgomeryMultR1, MontgomeryMultR2,
                  TestNbr, NumberLength);
  }

  private int KaratsubaAdd(int idxAddend1, int idxAddend2,
                           int idxSum, int length)
  {
    int carry = 0;
    for (int i=0; i<length; i++)
    {
      carry += arrNbr[idxAddend1 + i] + arrNbr[idxAddend2 + i];
      if (carry < 0)
      {                         // Bit 31 is set.
        arrNbr[idxSum + i] = carry - (-0x80000000);
        carry = 1;
      }
      else
      {
        arrNbr[idxSum + i] = carry;
        carry = 0;
      }
    }
    return carry;
  }

    // The return value is the sign: true: negative.
    // In result the absolute value of the difference is computed.
  private boolean absSubtract(int idxMinuend, int idxSubtrahend,
                              int idxResult, int length)
  {
    boolean sign = false;
    int carry = 0;
    int [] arrNbr = this.arrNbr;
    int i;
    for (i=length-1; i>=0; i--)
    {
      if (arrNbr[idxMinuend + i] != arrNbr[idxSubtrahend + i])
      {
        break;
      }
    }
    if (i>=0 && arrNbr[idxMinuend + i] < arrNbr[idxSubtrahend + i])
    {
      sign = true;
      i = idxMinuend;
      idxMinuend = idxSubtrahend;
      idxSubtrahend = i;
    }
    for (i=0; i<length; i++)
    {
      carry = arrNbr[idxMinuend + i] - arrNbr[idxSubtrahend + i] - carry;
      if (carry < 0)
      {
        arrNbr[idxResult + i] = carry + (-0x80000000);
        carry = 1;
      }
      else
      {
        arrNbr[idxResult + i] = carry;
        carry = 0;
      }
    }
    return sign;
  }

  private void NormalMultiply(int idxFactor1, int length)
  {
    long sum, carry, prod;
    int destIndex, j, offs, tmp;
    int idxFactor2 = idxFactor1 + length;
    sum = 0;
    int lastDestIndex = 2*length-1;
    for (destIndex=0; destIndex<lastDestIndex; destIndex++)
    {
      carry = 0;
      offs = idxFactor2+destIndex;
      if (destIndex<length)
      {
        for (j=destIndex; j>0; j-=2)
        {
          prod = (long)arrNbr[idxFactor1+j]*(long)arrNbr[offs-j] +
                 (long)arrNbr[idxFactor1+j-1]*(long)arrNbr[offs-j+1];
          sum += prod & 0x7FFFFFFFL;
          carry += prod >>> 31;
        }
        if (j==0)
        {
          prod = (long)arrNbr[idxFactor1+j]*(long)arrNbr[offs-j];
          sum += prod & 0x7FFFFFFFL;
          carry += prod >>> 31;
        }
      }
      else
      {
        for (j=length-1; j>destIndex-length; j-=2)
        {
          prod = (long)arrNbr[idxFactor1+j]*(long)arrNbr[offs-j] +
                 (long)arrNbr[idxFactor1+j-1]*(long)arrNbr[offs-j+1];
          sum += prod & 0x7FFFFFFFL;
          carry += prod >>> 31;
        }
        if (j==destIndex-length)
        {
          prod = (long)arrNbr[idxFactor1+j]*(long)arrNbr[offs-j];
          sum += prod & 0x7FFFFFFFL;
          carry += prod >>> 31;
        }
      }
      tmp = (int)(sum>>>31);
      arrNbrAux[destIndex] = (int)sum & 0x7FFFFFFF;
      sum = carry + tmp;
    }
    arrNbrAux[destIndex] = (int)sum;
    System.arraycopy(arrNbrAux, 0, arrNbr, idxFactor1, 2*length);
    return;
  }

  private void KaratsubaMultiply(int idxFactor1, int length, int endIndex)
  {
    int [] arrNbr = this.arrNbr;
    int idxFactor2 = idxFactor1 + length;
    int i, offs, carry, carry2, tmp;
    int middle;
    boolean sign;
    if (length <= KARATSUBA_CUTOFF)
    {
      // Check if one of the factors is equal to zero.
      for (i=length-1; i>=0; i--)
      {
        if (arrNbr[idxFactor1+i] != 0)
        {
          break;
        }
      }
      if (i>=0)
      {
        for (i=length-1; i>=0; i--)
        {
          if (arrNbr[idxFactor2+i] != 0)
          {
            break;
          }
        }
      }
      if (i<0)
      {              // One of the factors is equal to zero.
        for (i=length-1; i>=0; i--)
        {
          arrNbr[idxFactor1+i] = arrNbr[idxFactor2+i] = 0;
        }
        return;
      }
      NormalMultiply(idxFactor1, length);
    }
        // Length > KARATSUBA_CUTOFF: Use Karatsuba multiplication.
        // At this moment the order is: High1, Low1, High2, Low2.
        // Exchange low part of first factor with high part of 2nd factor.
    int halfLength = length/2;
    for (i=idxFactor1+halfLength; i<idxFactor2; i++)
    {
      tmp = arrNbr[i];
      arrNbr[i] = arrNbr[i+halfLength];
      arrNbr[i+halfLength] = tmp;
    }
        // At this moment the order is: High1, High2, Low1, Low2.
        // Get absolute values of (High1-Low1) and (Low2-High2) and the signs.
    sign = absSubtract(idxFactor1, idxFactor2, endIndex, halfLength);
    sign = absSubtract(idxFactor2+halfLength, idxFactor1+halfLength,
                        endIndex+halfLength, halfLength) != sign;
    middle = endIndex;
    endIndex += length;
                                      // Multiply both low parts.
    KaratsubaMultiply(idxFactor1, halfLength, endIndex);
                                      // Multiply both high parts.
    KaratsubaMultiply(idxFactor2, halfLength, endIndex);
    KaratsubaMultiply(middle, halfLength, endIndex);
    if (sign)
    {            // (High1-Low1) * (Low2-High2) is negative.
      if (absSubtract(idxFactor1, middle, middle, length))
      {          // Result is still negative.
        absSubtract(idxFactor2, middle, middle, length);
        carry2 = 0;
      }
      else
      {
        carry2 = KaratsubaAdd(idxFactor2, middle, middle, length);
      }
    }
    else
    {            // (High1-Low1) * (Low2-High2) is non-negative.
      carry2 = KaratsubaAdd(idxFactor1, middle, middle, length);
      carry2 += KaratsubaAdd(idxFactor2, middle, middle, length);
    }
    carry = 0;
    offs = idxFactor1+halfLength;
    for (i=0; i<length; i++)
    {
      carry += arrNbr[offs+i] + arrNbr[middle+i];
      if (carry < 0)
      {
        arrNbr[offs+i] = carry - (-0x80000000);
        carry = 1;
      }
      else
      {
        arrNbr[offs+i] = carry;
        carry = 0;
      }
    }
    arrNbr[idxFactor1+halfLength+i] += carry + carry2;
    if (arrNbr[idxFactor1+halfLength+i] < 0)
    {
      arrNbr[idxFactor1+halfLength+i] -= (-0x80000000);
      for (i=halfLength+1; i<length; i++)
      {
        if (++arrNbr[idxFactor1+i] >= 0)
        {
          break;
        }
        arrNbr[idxFactor1+i] = 0;
      }
    }
  }

  void BigNbrModN(int Nbr[], int Length, int Mod[])
  {
    int i, j;
    for (i = 0; i < NumberLength; i++)
    {
      Mod[i] = Nbr[i + Length - NumberLength];
    }
    Mod[i] = 0;
    AdjustModN(Mod, TestNbr, NumberLength);
    for (i = Length - NumberLength - 1; i >= 0; i--)
    {
      for (j = NumberLength; j > 0; j--)
      {
        Mod[j] = Mod[j - 1];
      }
      Mod[0] = Nbr[i];
      AdjustModN(Mod, TestNbr, NumberLength);
    }
  }

  static void MultBigNbrModN(int Nbr1[], int Nbr2[], int Prod[],
                             int TestNbr[], int NumberLength)
  {
    long MaxUInt = 0x7FFFFFFFL;
    int i, j;
    long Pr, Nbr;

    if (NumberLength >= 2 &&
        TestNbr[NumberLength-1]==0 && TestNbr[NumberLength-2]<0x40000000)
    {
      NumberLength--;
    }
    i = NumberLength;
    do
    {
      Prod[--i] = 0;
    }
    while (i > 0);
    i = NumberLength;
    do
    {
      Nbr = Nbr1[--i];
      j = NumberLength;
      do
      {
        Prod[j] = Prod[j - 1];
        j--;
      }
      while (j > 0);
      Prod[0] = 0;
      Pr = 0;
      for (j = 0; j < NumberLength; j++)
      {
        Pr = (Pr >>> 31) + Nbr * Nbr2[j] + Prod[j];
        Prod[j] = (int)(Pr & MaxUInt);
      }
      Prod[j] += (Pr >>> 31);
      AdjustModN(Prod, TestNbr, NumberLength);
    }
    while (i > 0);
  }

  static void MultBigNbrByLongModN(int Nbr1[], long Nbr2, int Prod[],
                                   int TestNbr[], int NumberLength)
  {
    long MaxUInt = 0x7FFFFFFFL;
    long Pr;
    int j;

    if (NumberLength>=2 &&
    		TestNbr[NumberLength-1]==0 && TestNbr[NumberLength-2]<0x40000000)
    {
      NumberLength--;
    }
    Pr = 0;
    for (j = 0; j < NumberLength; j++)
    {
      Pr = (Pr >>> 31) + Nbr2 * Nbr1[j];
      Prod[j] = (int)(Pr & MaxUInt);
    }
    Prod[j] = (int)(Pr >>> 31);
    AdjustModN(Prod, TestNbr, NumberLength);
  }

  static void AdjustModN(int Nbr[], int TestNbr[], int NumberLength)
  {
    long MaxUInt = 0x7FFFFFFFL;
    long TrialQuotient;
    long carry;
    int i;
    double dAux, dN;

    dN = (double) TestNbr[NumberLength - 1];
    if (NumberLength > 1)
    {
      dN += (double) TestNbr[NumberLength - 2] / dDosALa31;
    }
    if (NumberLength > 2)
    {
      dN += (double) TestNbr[NumberLength - 3] / dDosALa62;
    }
    dAux =
      (double) Nbr[NumberLength] * dDosALa31 + (double) Nbr[NumberLength - 1];
    if (NumberLength > 1)
    {
      dAux += (double) Nbr[NumberLength - 2] / dDosALa31;
    }
    TrialQuotient = (long) (dAux / dN) + 3;
    if (TrialQuotient >= DosALa32)
    {
      carry = 0;
      for (i = 0; i < NumberLength; i++)
      {
        carry = Nbr[i + 1] - (TrialQuotient >>> 31) * TestNbr[i] - carry;
        Nbr[i + 1] = (int)(carry & MaxUInt);
        carry = (MaxUInt - carry) >>> 31;
      }
      TrialQuotient &= MaxUInt;
    }
    carry = 0;
    for (i = 0; i < NumberLength; i++)
    {
      carry = Nbr[i] - TrialQuotient * TestNbr[i] - carry;
      Nbr[i] = (int)(carry & MaxUInt);
      carry = (MaxUInt - carry) >>> 31;
    }
    Nbr[NumberLength] -= (int)carry;
    while ((Nbr[NumberLength] & MaxUInt) != 0)
    {
      carry = 0;
      for (i = 0; i < NumberLength; i++)
      {
        carry += (long)Nbr[i] + (long)TestNbr[i];
        Nbr[i] = (int)(carry & MaxUInt);
        carry >>= 31;
      }
      Nbr[NumberLength] += (int)carry;
    }
  }

  static void DivBigNbrByLong(int Dividend[], long Divisor, int Quotient[],
                              int NumberLength)
  {
    int i;
    boolean ChSignDivisor = false;
    long Divid, Rem = 0;

    if (Divisor < 0)
    {                            // If divisor is negative...
      ChSignDivisor = true;      // Indicate to change sign at the end and
      Divisor = -Divisor;        // convert divisor to positive.
    }
    if (Dividend[i = NumberLength - 1] >= 0x40000000)
    {                            // If dividend is negative...
      Rem = Divisor - 1;
    }
    for ( ; i >= 0; i--)
    {
      Divid = Dividend[i] + (Rem << 31);
      Rem = Divid - (Quotient[i] = (int)(Divid / Divisor))*Divisor;
    }
    if (ChSignDivisor)
    {                            // Change sign if divisor is negative.
                                 // Convert divisor to positive.
      ChSignBigNbr(Quotient, NumberLength);
    }
  }

  static long RemDivBigNbrByLong(int Dividend[], long Divisor,
                                 int NumberLength)
  {
    int i;
    long Rem = 0;
    long Mod2_31;
    int divis = (int)(Divisor < 0?-Divisor:Divisor);
    if (Divisor < 0)
    {                            // If divisor is negative...
      Divisor = -Divisor;        // Convert divisor to positive.
    }
    Mod2_31 = ((-2147483648)-divis)%divis;  // 2^31 mod divis.
    if (Dividend[i = NumberLength - 1] >= 0x40000000)
    {                            // If dividend is negative...
      Rem = Divisor - 1;
    }
    for ( ; i >= 0; i--)
    {
      Rem = Rem * Mod2_31 + Dividend[i];
      do
      {
        Rem = (Rem >> 31)*Mod2_31+(Rem & 0x7FFFFFFF);
      } while (Rem > 0x1FFFFFFFFL);
    }
    return Rem % divis;
  }

  // Gcd calculation:
  // Step 1: Set k<-0, and then repeatedly set k<-k+1, u<-u/2, v<-v/2
  //         zero or more times until u and v are not both even.
  // Step 2: If u is odd, set t<-(-v) and go to step 4. Otherwise set t<-u.
  // Step 3: Set t<-t/2
  // Step 4: If t is even, go back to step 3.
  // Step 5: If t>0, set u<-t, otherwise set v<-(-t).
  // Step 6: Set t<-u-v. If t!=0, go back to step 3.
  // Step 7: The GCD is u*2^k.
  public void GcdBigNbr(int Nbr1[], int Nbr2[], int Gcd[], int NumberLength)
  {
    int i, k;

    System.arraycopy(Nbr1, 0, CalcAuxGcdU, 0, NumberLength);
    System.arraycopy(Nbr2, 0, CalcAuxGcdV, 0, NumberLength);
    for (i = 0; i < NumberLength; i++)
    {
      if (CalcAuxGcdU[i] != 0)
      {
        break;
      }
    }
    if (i == NumberLength)
    {
      System.arraycopy(CalcAuxGcdV, 0, Gcd, 0, NumberLength);
      return;
    }
    for (i = 0; i < NumberLength; i++)
    {
      if (CalcAuxGcdV[i] != 0)
      {
        break;
      }
    }
    if (i == NumberLength)
    {
      System.arraycopy(CalcAuxGcdU, 0, Gcd, 0, NumberLength);
      return;
    }
    if (CalcAuxGcdU[NumberLength - 1] >= 0x40000000L)
    {
      ChSignBigNbr(CalcAuxGcdU, NumberLength);
    }
    if (CalcAuxGcdV[NumberLength - 1] >= 0x40000000L)
    {
      ChSignBigNbr(CalcAuxGcdV, NumberLength);
    }
    k = 0;
    while ((CalcAuxGcdU[0] & 1) == 0 && (CalcAuxGcdV[0] & 1) == 0)
    { // Step 1
      k++;
      DivBigNbrByLong(CalcAuxGcdU, 2, CalcAuxGcdU, NumberLength);
      DivBigNbrByLong(CalcAuxGcdV, 2, CalcAuxGcdV, NumberLength);
    }
    if ((CalcAuxGcdU[0] & 1) == 1)
    { // Step 2
      System.arraycopy(CalcAuxGcdV, 0, CalcAuxGcdT, 0, NumberLength);
      ChSignBigNbr(CalcAuxGcdT, NumberLength);
    }
    else
    {
      System.arraycopy(CalcAuxGcdU, 0, CalcAuxGcdT, 0, NumberLength);
    }
    do
    {
      while ((CalcAuxGcdT[0] & 1) == 0)
      { // Step 4
        DivBigNbrByLong(CalcAuxGcdT, 2, CalcAuxGcdT, NumberLength); // Step 3
      }
      if (CalcAuxGcdT[NumberLength - 1] < 0x40000000L)
      { // Step 5
        System.arraycopy(CalcAuxGcdT, 0, CalcAuxGcdU, 0, NumberLength);
      }
      else
      {
        System.arraycopy(CalcAuxGcdT, 0, CalcAuxGcdV, 0, NumberLength);
        ChSignBigNbr(CalcAuxGcdV, NumberLength);
      }                                                // Step 6
      SubtractBigNbr(CalcAuxGcdU, CalcAuxGcdV, CalcAuxGcdT, NumberLength);
      for (i = 0; i < NumberLength; i++)
      {
        if (CalcAuxGcdT[i] != 0)
        {
          break;
        }
      }
    }
    while (i != NumberLength);
    System.arraycopy(CalcAuxGcdU, 0, Gcd, 0, NumberLength); // Step 7
    while (k > 0)
    {
      AddBigNbr(Gcd, Gcd, Gcd, NumberLength);
      k--;
    }
  }

  String BigNbrToString(int Nbr[])
  {
    long Rem;
    int i;
    boolean ChSign = false;
    String nbrOutput = "";

    if (Nbr[NumberLength - 1] >= 0x40000000)
    {
      ChSignBigNbr(Nbr, NumberLength);
      ChSign = true;
    }
    System.arraycopy(Nbr, 0, CalcBigNbr, 0, NumberLength);
    do
    {
      Rem = 0;
      for (i = NumberLength - 1; i >= 0; i--)
      {
        CalcBigNbr[i] += Rem << 31;
        Rem = CalcBigNbr[i] % Mi;
        CalcBigNbr[i] /= Mi;
      }
      nbrOutput = String.valueOf(Rem + Mi).substring(1) + nbrOutput;
      for (i = 0; i < NumberLength; i++)
      {
        if (CalcBigNbr[i] != 0)
        {
          break;
        }
      }
    }
    while (i < NumberLength);
    while (nbrOutput.charAt(0) == '0' && nbrOutput.length() > 1)
    {
      nbrOutput = nbrOutput.substring(1);
    }
    if (ChSign)
    {
      ChSignBigNbr(Nbr, NumberLength);
      nbrOutput = "-" + nbrOutput;
    }
    return nbrOutput;
  }

  static void Convert31To32Bits(int[] nbr31bits, long[] nbr32bits,
                                int NumberLength)
  {
    int i, j, k;
    i = 0;
    for (j = -1; j < NumberLength; j++)
    {
      k = i%31;
      if (k == 0)
      {
        j++;
      }
      if (j == NumberLength)
      {
        break;
      }
      if (j == NumberLength-1)
      {
        nbr32bits[i] = nbr31bits[j] >> k;
      }
      else
      {
        nbr32bits[i] = ((nbr31bits[j] >> k) |
                        (nbr31bits[j+1] << (31-k))) & 0xFFFFFFFFL;
      }
      i++;
    }
    for (; i<NumberLength; i++)
    {
      nbr32bits[i] = 0;
    }
  }

  static void Convert32To31Bits(long [] nbr32bits, int [] nbr31bits,
                                int NumberLength)
  {
    int i, j, k;
    j = 0;
    nbr32bits[NumberLength] = 0;
    for (i = 0; i < NumberLength; i++)
    {
      k = i & 0x0000001F;
      if (k == 0)
      {
        nbr31bits[i] = (int)(nbr32bits[j] & 0x7FFFFFFF);
      }
      else
      {
        nbr31bits[i] = (int)(((nbr32bits[j] >> (32-k)) |
                              (nbr32bits[j+1] << k)) & 0x7FFFFFFF);
        j++;
      }
    }
  }

  /***********************************************************************/
  /* NAME: ModInvBigNbr                                                  */
  /*                                                                     */
  /* PURPOSE: Find the inverse multiplicative modulo v.                  */
  /*                                                                     */
  /* The algorithm terminates with u1 = u^(-1) MOD v.                    */
  /***********************************************************************/
  public void ModInvBigNbr(int[] a, int[] inv, int[] b, int NumberLength)
  {
    int i;
    int Dif, E;
    int st1, st2;
    long Yaa, Yab; // 2^E * A'     = Yaa A + Yab B
    long Yba, Ybb; // 2^E * B'     = Yba A + Ybb B
    long Ygb0; // 2^E * Mu'    = Yaa Mu + Yab Gamma + Ymb0 B0
    long Ymb0; // 2^E * Gamma' = Yba Mu + Ybb Gamma + Ygb0 B0
    int Iaa, Iab, Iba, Ibb;
    long Tmp1, Tmp2, Tmp3, Tmp4, Tmp5;
    int B0l;
    int invB0l;
    int Al, Bl, T1, Gl, Ml, P;
    long carry1, carry2, carry3, carry4;
    int Yaah, Yabh, Ybah, Ybbh;
    int Ymb0h, Ygb0h;
    long Pr1, Pr2, Pr3, Pr4, Pr5, Pr6, Pr7;
    long[] B = this.CalcAuxModInvBB;
    long[] CalcAuxModInvA = this.CalcAuxModInvA;
    long[] CalcAuxModInvB = this.CalcAuxModInvB;
    long[] CalcAuxModInvMu = this.CalcAuxModInvMu;
    long[] CalcAuxModInvGamma = this.CalcAuxModInvGamma;

    if (NumberLength >= 2 &&
        b[NumberLength-1] == 0 && b[NumberLength-2] < 0x40000000)
    {
      NumberLength--;
    }
    Convert31To32Bits(a, CalcAuxModInvA, NumberLength);
    Convert31To32Bits(b, CalcAuxModInvB, NumberLength);
    System.arraycopy(CalcAuxModInvB, 0, B, 0, NumberLength);
    B0l = (int)B[0];
    invB0l = B0l; // 2 least significant bits of inverse correct.
    invB0l = invB0l * (2 - B0l * invB0l); // 4 LSB of inverse correct.
    invB0l = invB0l * (2 - B0l * invB0l); // 8 LSB of inverse correct.
    invB0l = invB0l * (2 - B0l * invB0l); // 16 LSB of inverse correct.
    invB0l = invB0l * (2 - B0l * invB0l); // 32 LSB of inverse correct.
    for (i = NumberLength - 1; i >= 0; i--)
    {
      CalcAuxModInvGamma[i] = 0;
      CalcAuxModInvMu[i] = 0;
    }
    CalcAuxModInvMu[0] = 1;
    Dif = 0;
    outer_loop : for (;;)
    {
      Iaa = Ibb = 1;
      Iab = Iba = 0;
      Al = (int) CalcAuxModInvA[0];
      Bl = (int) CalcAuxModInvB[0];
      E = 0;
      P = 1;
      if (Bl == 0)
      {
        for (i = NumberLength - 1; i >= 0; i--)
        {
          if (CalcAuxModInvB[i] != 0)
            break;
        }
        if (i < 0)
          break; // Go out of loop if CalcAuxModInvB = 0
      }
      for (;;)
      {
        T1 = 0;
        while ((Bl & 1) == 0)
        {
          if (E == 31)
          {
            Yaa = Iaa;
            Yab = Iab;
            Yba = Iba;
            Ybb = Ibb;
            Gl = (int) CalcAuxModInvGamma[0];
            Ml = (int) CalcAuxModInvMu[0];
            Dif++;
            T1++;
            Yaa <<= T1;
            Yab <<= T1;
            Ymb0 = (- (int) Yaa * Ml - (int) Yab * Gl) * invB0l;
            Ygb0 = (-Iba * Ml - Ibb * Gl) * invB0l;
            carry1 = carry2 = carry3 = carry4 = 0;
            Yaah = (int) (Yaa >> 32);
            Yabh = (int) (Yab >> 32);
            Ybah = (int) (Yba >> 32);
            Ybbh = (int) (Ybb >> 32);
            Ymb0h = (int) (Ymb0 >> 32);
            Ygb0h = (int) (Ygb0 >> 32);
            Yaa &= 0xFFFFFFFFL;
            Yab &= 0xFFFFFFFFL;
            Yba &= 0xFFFFFFFFL;
            Ybb &= 0xFFFFFFFFL;
            Ymb0 &= 0xFFFFFFFFL;
            Ygb0 &= 0xFFFFFFFFL;

            st1 = Yaah * 6 + Yabh * 2 + Ymb0h;
            st2 = Ybah * 6 + Ybbh * 2 + Ygb0h;
            for (i = 0; i < NumberLength; i++)
            {
              Pr1 = Yaa * (Tmp1 = CalcAuxModInvMu[i]);
              Pr2 = Yab * (Tmp2 = CalcAuxModInvGamma[i]);
              Pr3 = Ymb0 * (Tmp3 = B[i]);
              Pr4 =
                (Pr1 & 0xFFFFFFFFL)
                  + (Pr2 & 0xFFFFFFFFL)
                  + (Pr3 & 0xFFFFFFFFL)
                  + carry3;
              Pr5 = Yaa * (Tmp4 = CalcAuxModInvA[i]);
              Pr6 = Yab * (Tmp5 = CalcAuxModInvB[i]);
              Pr7 = (Pr5 & 0xFFFFFFFFL) + (Pr6 & 0xFFFFFFFFL) + carry1;
              switch (st1)
              {
                case -9 :
                  carry3 = -Tmp1 - Tmp2 - Tmp3;
                  carry1 = -Tmp4 - Tmp5;
                  break;
                case -8 :
                  carry3 = -Tmp1 - Tmp2;
                  carry1 = -Tmp4 - Tmp5;
                  break;
                case -7 :
                  carry3 = -Tmp1 - Tmp3;
                  carry1 = -Tmp4;
                  break;
                case -6 :
                  carry3 = -Tmp1;
                  carry1 = -Tmp4;
                  break;
                case -5 :
                  carry3 = -Tmp1 + Tmp2 - Tmp3;
                  carry1 = -Tmp4 + Tmp5;
                  break;
                case -4 :
                  carry3 = -Tmp1 + Tmp2;
                  carry1 = -Tmp4 + Tmp5;
                  break;
                case -3 :
                  carry3 = -Tmp2 - Tmp3;
                  carry1 = -Tmp5;
                  break;
                case -2 :
                  carry3 = -Tmp2;
                  carry1 = -Tmp5;
                  break;
                case -1 :
                  carry3 = -Tmp3;
                  carry1 = 0;
                  break;
                case 0 :
                  carry3 = 0;
                  carry1 = 0;
                  break;
                case 1 :
                  carry3 = Tmp2 - Tmp3;
                  carry1 = Tmp5;
                  break;
                case 2 :
                  carry3 = Tmp2;
                  carry1 = Tmp5;
                  break;
                case 3 :
                  carry3 = Tmp1 - Tmp2 - Tmp3;
                  carry1 = Tmp4 - Tmp5;
                  break;
                case 4 :
                  carry3 = Tmp1 - Tmp2;
                  carry1 = Tmp4 - Tmp5;
                  break;
                case 5 :
                  carry3 = Tmp1 - Tmp3;
                  carry1 = Tmp4;
                  break;
                case 6 :
                  carry3 = Tmp1;
                  carry1 = Tmp4;
                  break;
                case 7 :
                  carry3 = Tmp1 + Tmp2 - Tmp3;
                  carry1 = Tmp4 + Tmp5;
                  break;
                case 8 :
                  carry3 = Tmp1 + Tmp2;
                  carry1 = Tmp4 + Tmp5;
                  break;
              }
              carry3 += (Pr1 >>> 32) + (Pr2 >>> 32) + (Pr3 >>> 32) + (Pr4 >> 32);
              carry1 += (Pr5 >>> 32) + (Pr6 >>> 32) + (Pr7 >> 32);
              if (i > 0)
              {
                CalcAuxModInvMu[i - 1] = Pr4 & 0xFFFFFFFFL;
                CalcAuxModInvA[i - 1] = Pr7 & 0xFFFFFFFFL;
              }
              Pr1 = Yba * Tmp1;
              Pr2 = Ybb * Tmp2;
              Pr3 = Ygb0 * Tmp3;
              Pr4 =
                (Pr1 & 0xFFFFFFFFL)
                  + (Pr2 & 0xFFFFFFFFL)
                  + (Pr3 & 0xFFFFFFFFL)
                  + carry4;
              Pr5 = Yba * Tmp4;
              Pr6 = Ybb * Tmp5;
              Pr7 = (Pr5 & 0xFFFFFFFFL) + (Pr6 & 0xFFFFFFFFL) + carry2;
              switch (st2)
              {
                case -9 :
                  carry4 = -Tmp1 - Tmp2 - Tmp3;
                  carry2 = -Tmp4 - Tmp5;
                  break;
                case -8 :
                  carry4 = -Tmp1 - Tmp2;
                  carry2 = -Tmp4 - Tmp5;
                  break;
                case -7 :
                  carry4 = -Tmp1 - Tmp3;
                  carry2 = -Tmp4;
                  break;
                case -6 :
                  carry4 = -Tmp1;
                  carry2 = -Tmp4;
                  break;
                case -5 :
                  carry4 = -Tmp1 + Tmp2 - Tmp3;
                  carry2 = -Tmp4 + Tmp5;
                  break;
                case -4 :
                  carry4 = -Tmp1 + Tmp2;
                  carry2 = -Tmp4 + Tmp5;
                  break;
                case -3 :
                  carry4 = -Tmp2 - Tmp3;
                  carry2 = -Tmp5;
                  break;
                case -2 :
                  carry4 = -Tmp2;
                  carry2 = -Tmp5;
                  break;
                case -1 :
                  carry4 = -Tmp3;
                  carry2 = 0;
                  break;
                case 0 :
                  carry4 = 0;
                  carry2 = 0;
                  break;
                case 1 :
                  carry4 = Tmp2 - Tmp3;
                  carry2 = Tmp5;
                  break;
                case 2 :
                  carry4 = Tmp2;
                  carry2 = Tmp5;
                  break;
                case 3 :
                  carry4 = Tmp1 - Tmp2 - Tmp3;
                  carry2 = Tmp4 - Tmp5;
                  break;
                case 4 :
                  carry4 = Tmp1 - Tmp2;
                  carry2 = Tmp4 - Tmp5;
                  break;
                case 5 :
                  carry4 = Tmp1 - Tmp3;
                  carry2 = Tmp4;
                  break;
                case 6 :
                  carry4 = Tmp1;
                  carry2 = Tmp4;
                  break;
                case 7 :
                  carry4 = Tmp1 + Tmp2 - Tmp3;
                  carry2 = Tmp4 + Tmp5;
                  break;
                case 8 :
                  carry4 = Tmp1 + Tmp2;
                  carry2 = Tmp4 + Tmp5;
                  break;
              }
              carry4 += (Pr1 >>> 32) + (Pr2 >>> 32) + (Pr3 >>> 32) + (Pr4 >> 32);
              carry2 += (Pr5 >>> 32) + (Pr6 >>> 32) + (Pr7 >> 32);
              if (i > 0)
              {
                CalcAuxModInvGamma[i - 1] = Pr4 & 0xFFFFFFFFL;
                CalcAuxModInvB[i - 1] = Pr7 & 0xFFFFFFFFL;
              }
            }

            if ((int) CalcAuxModInvA[i - 1] < 0)
            {
              carry1 -= Yaa;
              carry2 -= Yba;
            }
            if ((int) CalcAuxModInvB[i - 1] < 0)
            {
              carry1 -= Yab;
              carry2 -= Ybb;
            }
            if ((int) CalcAuxModInvMu[i - 1] < 0)
            {
              carry3 -= Yaa;
              carry4 -= Yba;
            }
            if ((int) CalcAuxModInvGamma[i - 1] < 0)
            {
              carry3 -= Yab;
              carry4 -= Ybb;
            }
            CalcAuxModInvA[i - 1] = carry1 & 0xFFFFFFFFL;
            CalcAuxModInvB[i - 1] = carry2 & 0xFFFFFFFFL;
            CalcAuxModInvMu[i - 1] = carry3 & 0xFFFFFFFFL;
            CalcAuxModInvGamma[i - 1] = carry4 & 0xFFFFFFFFL;
            continue outer_loop;
          }
          Bl >>= 1;
          Dif++;
          E++;
          P *= 2;
          T1++;
        } /* end while */
        Iaa <<= T1;
        Iab <<= T1;
        if (Dif >= 0)
        {
          Dif = -Dif;
          if (((Al + Bl) & 3) == 0)
          {
            T1 = Iba;
            Iba += Iaa;
            Iaa = T1;
            T1 = Ibb;
            Ibb += Iab;
            Iab = T1;
            T1 = Bl;
            Bl += Al;
            Al = T1;
          }
          else
          {
            T1 = Iba;
            Iba -= Iaa;
            Iaa = T1;
            T1 = Ibb;
            Ibb -= Iab;
            Iab = T1;
            T1 = Bl;
            Bl -= Al;
            Al = T1;
          }
        }
        else
        {
          if (((Al + Bl) & 3) == 0)
          {
            Iba += Iaa;
            Ibb += Iab;
            Bl += Al;
          }
          else
          {
            Iba -= Iaa;
            Ibb -= Iab;
            Bl -= Al;
          }
        }
        Dif--;
      }
    }
    if (CalcAuxModInvA[0] != 1)
    {
      SubtractBigNbr32(B, CalcAuxModInvMu, CalcAuxModInvMu, NumberLength);
    }
    if ((int) CalcAuxModInvMu[i = NumberLength - 1] < 0)
    {
      AddBigNbr32(B, CalcAuxModInvMu, CalcAuxModInvMu, NumberLength);
    }
    for (; i >= 0; i--)
    {
      if (B[i] != CalcAuxModInvMu[i])
        break;
    }
    if (i < 0 || B[i] < CalcAuxModInvMu[i])
    { // If B < Mu
      SubtractBigNbr32(CalcAuxModInvMu, B, CalcAuxModInvMu, NumberLength); // Mu <- Mu - B
    }
    Convert32To31Bits(CalcAuxModInvMu, inv, NumberLength);
  }

  void FactorFibonacci(int Index, BigInteger BigOriginal)
  {
    int Index2, k;
    BigInteger Nro1;
    if (onlyFactoring)
    {
      NroFact = 1;
      Factores[0] = BigOriginal;
      Index2 = Index;
      while (Index2 % 2 == 0)
      {
        Index2 /= 2;
      }
      k = 1; /* Factor F(Index2) */
      while (k * k <= Index2)
      {
        if (Index2 % k == 0)
        {
          Nro1 = Fibonacci(k).gcd(BigOriginal);
          InsertFactor(Nro1);
          InsertFactor(BigOriginal.divide(Nro1));
          Nro1 = Fibonacci(Index / k).gcd(BigOriginal);
          InsertFactor(Nro1);
          InsertFactor(BigOriginal.divide(Nro1));
        }
        k += 2;
      }
      Index2 = Index;
      while (Index2 % 2 == 0)
      {
        Index2 /= 2;
        InsertLucasFactor(Index2, BigOriginal);
      }
      SortFactors();
    }
  }

  void FactorLucas(int Index, BigInteger BigOriginal)
  {
    NroFact = 1;
    if (onlyFactoring)
    {
      Factores[0] = BigOriginal;
      InsertLucasFactor(Index, BigOriginal);
      SortFactors();
    }
  }

  void InsertLucasFactor(int Index, BigInteger BigOriginal)
  {
    int k;
    if (onlyFactoring)
    {
      BigInteger Fibo = BigInt0;
      BigInteger BigInt5 = BigInteger.valueOf(5);
      BigInteger Nro1;
      k = 1; /* Factor L(Index) */
      while (k * k <= Index)
      {
        if (Index % k == 0)
        {
          Nro1 = Lucas(k).gcd(BigOriginal);
          InsertFactor(Nro1);
          if (k % 5 == 0)
          {
            Fibo = Fibonacci(k);
            InsertFactor(
              BigInt5.multiply(Fibo).subtract(BigInt5).multiply(Fibo).add(
                BigInt1));
            InsertFactor(
              BigInt5.multiply(Fibo).add(BigInt5).multiply(Fibo).add(BigInt1));
          }
          else
          {
            InsertFactor(Nro1);
            InsertFactor(BigOriginal.divide(Nro1));
          }
          Nro1 = Lucas(Index / k).gcd(BigOriginal);
          InsertFactor(Nro1);
          if ((Index / k) % 5 == 0)
          {
            Fibo = Fibonacci(Index / k);
            InsertFactor(
              BigInt5.multiply(Fibo).subtract(BigInt5).multiply(Fibo).add(
                BigInt1));
            InsertFactor(
              BigInt5.multiply(Fibo).add(BigInt5).multiply(Fibo).add(BigInt1));
          }
          else
          {
            InsertFactor(Nro1);
            InsertFactor(BigOriginal.divide(Nro1));
          }
        }
        k++;
      }
    }
  }

  BigInteger Fibonacci(int Index)
  {
    int i;
    if (onlyFactoring)
    {
      BigInteger FibonPrev = BigInt1;
      BigInteger FibonAct = BigInt0;
      BigInteger FibonNext;
      for (i = 1; i <= Index; i++)
      {
        FibonNext = FibonPrev.add(FibonAct);
        FibonPrev = FibonAct;
        FibonAct = FibonNext;
      }
      return FibonAct;
    }
    else
    {
      return null;
    }
  }

  BigInteger Lucas(int Index)
  {
    int i;
    if (onlyFactoring)
    {
      BigInteger LucasPrev = BigInteger.valueOf(-1);
      BigInteger LucasAct = BigInt2;
      BigInteger LucasNext;
      for (i = 1; i <= Index; i++)
      {
        LucasNext = LucasPrev.add(LucasAct);
        LucasPrev = LucasAct;
        LucasAct = LucasNext;
      }
      return LucasAct;
    }
    else
    {
      return null;
    }
  }

  void Cunningham(
    BigInteger BigBase,
    int Expon,
    BigInteger BigIncre,
    BigInteger BigOriginal)
  {
    int Expon2, k;
    int Incre;
    BigInteger Nro1;

    if (onlyFactoring)
    {
      Incre = BigIncre.intValue();
      NroFact = 1;
      Factores[0] = BigOriginal;
      Expon2 = Expon;
      Nro1 = BigBase.pow(Expon2).add(BigIncre);
      InsertFactor(Nro1);
      if ((digitsInGroup & 0x1000) == 0 && Nro1.bitLength() > 60)
      {
        try
        { // Get known primitive factors.
          lowerTextArea.setText(
            "Requesting known primitive factors from Web server.");
          URL script =
            new URL(
              getDocumentBase(),
              "factors.pl?base="
                + BigBase.intValue()
                + "&expon="
                + Expon
                + "&type="
                + (Incre == 1 ? "p" : "m"));
          BufferedInputStream buffer =
            new BufferedInputStream(script.openStream());
          BufferedReader in = new BufferedReader(new InputStreamReader(buffer));
          String factorsAscii = in.readLine();
          in.close();
          if (factorsAscii.length() > 0)
          { // Factors found in server.
            int indexFactors = 0;
            int newIndexFactors = 0;
            do
            { // Loop through factors.
              newIndexFactors = factorsAscii.indexOf('*', indexFactors);
              if (newIndexFactors > 0)
              {
                InsertFactor(
                  new BigInteger(
                    factorsAscii.substring(indexFactors, newIndexFactors)));
              }
              else
              {
                InsertFactor(
                  new BigInteger(factorsAscii.substring(indexFactors)));
              }
              indexFactors = newIndexFactors + 1;
            }
            while (indexFactors > 0);
          }
        }
        catch (Exception e)
        {
        }
      }
      while (Expon2 % 2 == 0 && Incre == -1)
      {
        Expon2 /= 2;
        InsertFactor(BigBase.pow(Expon2).add(BigInt1));
        InsertAurifFactors(BigBase, Expon2, 1);
      }
      k = 1;
      while (k * k <= Expon)
      {
        if (Expon % k == 0)
        {
          if (k % 2 != 0)
          { /* Only for odd exponent */
            Nro1 = BigBase.pow(Expon / k).add(BigIncre).gcd(BigOriginal);
            InsertFactor(Nro1);
            InsertFactor(BigOriginal.divide(Nro1));
            InsertAurifFactors(BigBase, Expon / k, Incre);
          }
          if ((Expon / k) % 2 != 0)
          { /* Only for odd exponent */
            Nro1 = BigBase.pow(k).add(BigIncre).gcd(BigOriginal);
            InsertFactor(Nro1);
            InsertFactor(BigOriginal.divide(Nro1));
            InsertAurifFactors(BigBase, k, Incre);
          }
        }
        k++;
      }
      SortFactors();
    }
  }

  void InsertAurifFactors(BigInteger BigBase, int Expon, int Incre)
  {
    int t1, t2, t3, N, N1, q, L, j, k, Base;
    if (BigBase.compareTo(BigInteger.valueOf(386)) <= 0)
    {
      Base = BigBase.intValue();
      if (Expon % 2 == 0 && Incre == -1)
      {
        do
        {
          Expon /= 2;
        }
        while (Expon % 2 == 0);
        Incre = Base % 4 - 2;
      }
      if (Expon % Base == 0
        && Expon / Base % 2 != 0
        && ((Base % 4 != 1 && Incre == 1) || (Base % 4 == 1 && Incre == -1)))
      {
        N = Base;
        if (N % 4 == 1)
        {
          N1 = N;
        }
        else
        {
          N1 = 2 * N;
        }
        DegreeAurif = Totient(N1) / 2;
        for (k = 1; k <= DegreeAurif; k += 2)
        {
          AurifQ[k] = JacobiSymbol(N, k);
        }
        for (k = 2; k <= DegreeAurif; k += 2)
        {
          t1 = k; // Calculate t2 = gcd(k, N1)
          t2 = N1;
          while (t1 != 0)
          {
            t3 = t2 % t1;
            t2 = t1;
            t1 = t3;
          }
          AurifQ[k] = Moebius(N1 / t2) * Totient(t2) * Cos((N - 1) * k);
        }
        Gamma[0] = Delta[0] = 1;
        for (k = 1; k <= DegreeAurif / 2; k++)
        {
          Gamma[k] = Delta[k] = 0;
          for (j = 0; j < k; j++)
          {
            Gamma[k] =
              Gamma[k]
                + N * AurifQ[2 * k
                - 2 * j
                - 1] * Delta[j]
                - AurifQ[2 * k
                - 2 * j] * Gamma[j];
            Delta[k] =
              Delta[k]
                + AurifQ[2 * k
                + 1
                - 2 * j] * Gamma[j]
                - AurifQ[2 * k
                - 2 * j] * Delta[j];
          }
          Gamma[k] /= 2 * k;
          Delta[k] = (Delta[k] + Gamma[k]) / (2 * k + 1);
        }
        for (k = DegreeAurif / 2 + 1; k <= DegreeAurif; k++)
        {
          Gamma[k] = Gamma[DegreeAurif - k];
        }
        for (k = (DegreeAurif + 1) / 2; k < DegreeAurif; k++)
        {
          Delta[k] = Delta[DegreeAurif - k - 1];
        }
        q = Expon / Base;
        L = 1;
        while (L * L <= q)
        {
          if (q % L == 0)
          {
            GetAurifeuilleFactor(L, BigBase);
            if (q != L * L)
            {
              GetAurifeuilleFactor(q / L, BigBase);
            }
          }
          L += 2;
        }
      }
    }
  }
  // Sort the factors
  void SortFactors()
  {
    int j, k;
    BigInteger Nro1;
    for (k = 0; k < NroFact - 1; k++)
    {
      for (j = k + 1; j < NroFact; j++)
      {
        if (Factores[k].compareTo(Factores[j]) > 0)
        {
          Nro1 = Factores[k];
          Factores[k] = Factores[j];
          Factores[j] = Nro1;
        }
      }
    }
    for (k = 0; k < NroFact; k++)
    {
      PD[NbrFactors + k - 1] = Factores[k];
      Exp[NbrFactors + k - 1] = 1;
      Typ[NbrFactors + k - 1] = -1; /* Unknown */
    }
    NbrFactors += k - 1;
  }

  int JacobiSymbol(int M, int Q)
  {
    int k, t1, t2, t3, jacobi;
    if (onlyFactoring)
    {

      // Calculate gcd(M,Q)

      t1 = M;
      t2 = Q;
      while (t1 != 0)
      {
        t3 = t2 % t1;
        t2 = t1;
        t1 = t3;
      }
      if (t2 > 1)
      {
        return 0;
      }
      jacobi = 1;
      while (Q % 2 == 0)
      {
        Q /= 2;
      }
      if (Q % 3 == 0)
      {
        do
        {
          jacobi = (jacobi * M) % 3;
          Q /= 3;
        }
        while (Q % 3 == 0);
        jacobi = (jacobi + 1) % 3 - 1;
      }

      k = 5;
      while (k * k <= Q)
      {
        if (k % 3 != 0)
        {
          while (Q % k == 0)
          {
            Q /= k;
            jacobi = (jacobi + k) % k;
            for (t1 = (k - 1) / 2; t1 > 0; t1--)
            {
              jacobi = jacobi * M % k;
            }
            jacobi = (jacobi + 1) % k - 1;
          }
        }
        k += 2;
      }
      if (Q > 1)
      {
        jacobi = (jacobi + Q) % Q;
        for (t1 = (Q - 1) / 2; t1 > 0; t1--)
        {
          jacobi = jacobi * M % Q;
        }
        jacobi = (jacobi + 1) % Q - 1;
      }
      return jacobi;
    }
    else
    {
      return 0;
    }
  }

  static int Cos(int N)
  {
    switch (N % 8)
    {
      case 0 :
        return 1;
      case 4 :
        return -1;
    }
    return 0;
  }

  int Totient(int N)
  {
    int totient, q, k;

    if (onlyFactoring)
    {
      totient = q = N;
      if (q % 2 == 0)
      {
        totient /= 2;
        do
        {
          q /= 2;
        }
        while (q % 2 == 0);
      }
      if (q % 3 == 0)
      {
        totient = totient * 2 / 3;
        do
        {
          q /= 3;
        }
        while (q % 3 == 0);
      }
      k = 5;
      while (k * k <= q)
      {
        if (k % 3 != 0 && q % k == 0)
        {
          totient = totient * (k - 1) / k;
          do
          {
            q /= k;
          }
          while (q % k == 0);
        }
        k += 2;
      }
      if (q > 1)
      {
        totient = totient * (q - 1) / q;
      }
      return totient;
    }
    else
    {
      return 0;
    }
  }

  int Moebius(int N)
  {
    int moebius, q, k;

    if (onlyFactoring)
    {
      moebius = 1;
      q = N;
      if (q % 2 == 0)
      {
        moebius = -moebius;
        q /= 2;
        if (q % 2 == 0)
        {
          return 0;
        }
      }
      if (q % 3 == 0)
      {
        moebius = -moebius;
        q /= 3;
        if (q % 3 == 0)
        {
          return 0;
        }
      }
      k = 5;
      while (k * k <= q)
      {
        if (k % 3 != 0)
        {
          while (q % k == 0)
          {
            moebius = -moebius;
            q /= k;
            if (q % k == 0)
            {
              return 0;
            }
          }
        }
        k += 2;
      }
      if (q > 1)
      {
        moebius = -moebius;
      }
      return moebius;
    }
    else
    {
      return 0;
    }
  }

  void InsertFactor(BigInteger N)
  {
    int g;

    for (g = NroFact - 1; g >= 0; g--)
    {
      Factores[NroFact] = Factores[g].gcd(N);
      if (!Factores[NroFact].equals(BigInt1) &&
          !Factores[NroFact].equals(Factores[g]))
      {
        Factores[g] = Factores[g].divide(Factores[NroFact]);
        NroFact++;
      }
    }
  }

  void GetAurifeuilleFactor(int L, BigInteger BigBase)
  {
    BigInteger X, Csal, Dsal, Nro1;
    int k;

    if (onlyFactoring)
    {
      X = BigBase.pow(L);
      Csal = Dsal = BigInt1;
      for (k = 1; k < DegreeAurif; k++)
      {
        Csal = Csal.multiply(X).add(BigInteger.valueOf(Gamma[k]));
        Dsal = Dsal.multiply(X).add(BigInteger.valueOf(Delta[k]));
      }
      Csal = Csal.multiply(X).add(BigInteger.valueOf(Gamma[k]));
      Nro1 = Dsal.multiply(BigBase.pow((L + 1) / 2));
      InsertFactor(Csal.add(Nro1));
      InsertFactor(Csal.subtract(Nro1));
    }
  }

  boolean ComputeFourSquares(BigInteger PD[], int Exp[])
  {
    if (onlyFactoring)
    {
      int indexPrimes;
      BigInteger p, q, K, Mult1, Mult2, Mult3, Mult4;
      BigInteger Tmp, Tmp1, Tmp2, Tmp3, M1, M2, M3, M4;

      Quad1 = BigInt1; /* 1 = 1^2 + 0^2 + 0^2 + 0^2 */
      Quad2 = BigInt0;
      Quad3 = BigInt0;
      Quad4 = BigInt0;
      for (indexPrimes = NbrFactors - 1; indexPrimes >= 0; indexPrimes--)
      {
        if (Exp[indexPrimes] % 2 == 0)
        {
          continue;
        }
        p = PD[indexPrimes];
        q = p.subtract(BigInt1); /* q = p-1 */
        if (p.equals(BigInt2))
        {
          Mult1 = BigInt1; /* 2 = 1^2 + 1^2 + 0^2 + 0^2 */
          Mult2 = BigInt1;
          Mult3 = BigInt0;
          Mult4 = BigInt0;
        }
        else
        { /* Prime not 2 */
          if (!p.testBit(1))
          { /* if p = 1 (mod 4) */
            K = BigInt1;
            do
            { // Compute Mult1 = sqrt(-1) mod p
              K = K.add(BigInt1);
              Mult1 = K.modPow(q.shiftRight(2), p);
            }
            while (Mult1.equals(BigInt1) || Mult1.equals(q));
            if (!Mult1.multiply(Mult1).mod(p).equals(q))
            {
              return false; /* The number is not prime */
            }
            Mult2 = BigInt1;
            for (;;)
            {
              K = Mult1.multiply(Mult1).add(Mult2.multiply(Mult2)).divide(p);
              if (K.equals(BigInt1))
              {
                Mult3 = BigInt0;
                Mult4 = BigInt0;
                break;
              }
              if (p.mod(K).signum() == 0)
              {
                return false; /* The number is not prime */
              }
              M1 = Mult1.mod(K);
              M2 = Mult2.mod(K);
              if (M1.compareTo(K.shiftRight(1)) > 0)
              {
                M1 = M1.subtract(K);
              }
              if (M2.compareTo(K.shiftRight(1)) > 0)
              {
                M2 = M2.subtract(K);
              }
              Tmp = Mult1.multiply(M1).add(Mult2.multiply(M2)).divide(K);
              Mult2 = Mult1.multiply(M2).subtract(Mult2.multiply(M1)).divide(K);
              Mult1 = Tmp;
            } /* end while */
          } /* end p = 1 (mod 4) */
          else
          { /* if p = 3 (mod 4) */
            // Compute Mult1 and Mult2 so Mult1^2 + Mult2^2 = -1 (mod p)
            Mult1 = BigInt0;
            do
            {
              Mult1 = Mult1.add(BigInt1);
            }
            while (BigInt1
              .negate()
              .subtract(Mult1.multiply(Mult1))
              .modPow(q.shiftRight(1), p)
              .compareTo(BigInt1)
              > 0);
            Mult2 =
              BigInt1.negate().subtract(Mult1.multiply(Mult1)).modPow(
                p.add(BigInt1).shiftRight(2),
                p);
            Mult3 = BigInt1;
            Mult4 = BigInt0;
            for (;;)
            {
              K =
                Mult1
                  .multiply(Mult1)
                  .add(Mult2.multiply(Mult2))
                  .add(Mult3.multiply(Mult3))
                  .add(Mult4.multiply(Mult4))
                  .divide(p);
              if (K.equals(BigInt1))
              {
                break;
              }
              if (!K.testBit(0))
              { // If K is even ...
                if (Mult1.add(Mult2).testBit(0))
                {
                  if (!Mult1.add(Mult3).testBit(0))
                  {
                    Tmp = Mult2;
                    Mult2 = Mult3;
                    Mult3 = Tmp;
                  }
                  else
                  {
                    Tmp = Mult2;
                    Mult2 = Mult4;
                    Mult4 = Tmp;
                  }
                } // At this moment Mult1+Mult2 = even, Mult3+Mult4 = even
                Tmp1 = Mult1.add(Mult2).shiftRight(1);
                Tmp2 = Mult1.subtract(Mult2).shiftRight(1);
                Tmp3 = Mult3.add(Mult4).shiftRight(1);
                Mult4 = Mult3.subtract(Mult4).shiftRight(1);
                Mult3 = Tmp3;
                Mult2 = Tmp2;
                Mult1 = Tmp1;
                continue;
              } /* end if k is even */
              M1 = Mult1.mod(K);
              M2 = Mult2.mod(K);
              M3 = Mult3.mod(K);
              M4 = Mult4.mod(K);
              if (M1.compareTo(K.shiftRight(1)) > 0)
              {
                M1 = M1.subtract(K);
              }
              if (M2.compareTo(K.shiftRight(1)) > 0)
              {
                M2 = M2.subtract(K);
              }
              if (M3.compareTo(K.shiftRight(1)) > 0)
              {
                M3 = M3.subtract(K);
              }
              if (M4.compareTo(K.shiftRight(1)) > 0)
              {
                M4 = M4.subtract(K);
              }
              Tmp1 =
                Mult1
                  .multiply(M1)
                  .add(Mult2.multiply(M2))
                  .add(Mult3.multiply(M3))
                  .add(Mult4.multiply(M4))
                  .divide(K);
              Tmp2 =
                Mult1
                  .multiply(M2)
                  .subtract(Mult2.multiply(M1))
                  .add(Mult3.multiply(M4))
                  .subtract(Mult4.multiply(M3))
                  .divide(K);
              Tmp3 =
                Mult1
                  .multiply(M3)
                  .subtract(Mult3.multiply(M1))
                  .subtract(Mult2.multiply(M4))
                  .add(Mult4.multiply(M2))
                  .divide(K);
              Mult4 =
                Mult1
                  .multiply(M4)
                  .subtract(Mult4.multiply(M1))
                  .add(Mult2.multiply(M3))
                  .subtract(Mult3.multiply(M2))
                  .divide(K);
              Mult3 = Tmp3;
              Mult2 = Tmp2;
              Mult1 = Tmp1;
            } /* end while */
          } /* end if p = 3 (mod 4) */
        } /* end prime not 2 */
        Tmp1 =
          Mult1.multiply(Quad1).add(Mult2.multiply(Quad2)).add(
            Mult3.multiply(Quad3)).add(
            Mult4.multiply(Quad4));
        Tmp2 =
          Mult1
            .multiply(Quad2)
            .subtract(Mult2.multiply(Quad1))
            .add(Mult3.multiply(Quad4))
            .subtract(Mult4.multiply(Quad3));
        Tmp3 =
          Mult1
            .multiply(Quad3)
            .subtract(Mult3.multiply(Quad1))
            .subtract(Mult2.multiply(Quad4))
            .add(Mult4.multiply(Quad2));
        Quad4 =
          Mult1
            .multiply(Quad4)
            .subtract(Mult4.multiply(Quad1))
            .add(Mult2.multiply(Quad3))
            .subtract(Mult3.multiply(Quad2));
        Quad3 = Tmp3;
        Quad2 = Tmp2;
        Quad1 = Tmp1;
      } /* end for indexPrimes */
      for (indexPrimes = 0; indexPrimes < NbrFactors; indexPrimes++)
      {
        p = PD[indexPrimes].pow(Exp[indexPrimes] / 2);
        Quad1 = Quad1.multiply(p);
        Quad2 = Quad2.multiply(p);
        Quad3 = Quad3.multiply(p);
        Quad4 = Quad4.multiply(p);
      }
      Quad1 = Quad1.abs();
      Quad2 = Quad2.abs();
      Quad3 = Quad3.abs();
      Quad4 = Quad4.abs();
      // Sort squares
      if (Quad1.compareTo(Quad2) < 0)
      {
        Tmp = Quad1;
        Quad1 = Quad2;
        Quad2 = Tmp;
      }
      if (Quad1.compareTo(Quad3) < 0)
      {
        Tmp = Quad1;
        Quad1 = Quad3;
        Quad3 = Tmp;
      }
      if (Quad1.compareTo(Quad4) < 0)
      {
        Tmp = Quad1;
        Quad1 = Quad4;
        Quad4 = Tmp;
      }
      if (Quad2.compareTo(Quad3) < 0)
      {
        Tmp = Quad2;
        Quad2 = Quad3;
        Quad3 = Tmp;
      }
      if (Quad2.compareTo(Quad4) < 0)
      {
        Tmp = Quad2;
        Quad2 = Quad4;
        Quad4 = Tmp;
      }
      if (Quad3.compareTo(Quad4) < 0)
      {
        Tmp = Quad3;
        Quad3 = Quad4;
        Quad4 = Tmp;
      }
    }
    return true;
  }

  private long modPow(long NbrMod, long Expon, long currentPrime)
  {
    long Power = 1;
    long Square = NbrMod;
    while (Expon != 0)
    {
      if ((Expon & 1) == 1)
      {
        Power = (Power * Square) % currentPrime;
      }
      Square = (Square * Square) % currentPrime;
      Expon /= 2;
    }
    return Power;
  }

  public String StartFactorExprBatch(String inputStr, int type)
  {
    batchFinished = false;
    this.inputStr = inputStr;
    batchPrime = (type == 1);
    calcThread = new Thread(this);
    calcThread.start();
    return "";
  }

  void BatchThread()
  {
    outputStr = new StringBuffer();
    String expr, firstN, nextN, endExpression, ExpressionToCompute;
    int StrLen = inputStr.length();
    int index = 0;
    int counter;
    int newIndex, newIndex2, newIndex3;
    BigInteger N;
    BigInteger ExpressionResult[] = new BigInteger[1];
    textNumber.setEditable(false);
    labelTop.setText("Number to factor:");
    while (index < StrLen)
    {
      newIndex = inputStr.indexOf('\n', index);
      if (newIndex < 0)
      {                // No line feed found, it means last line.
        newIndex = StrLen;
      }
      expr = inputStr.substring(index, newIndex).trim();
      index = newIndex + 1;
      if (expr.length() == 0)
      {                      // Replicate blank line.
        outputStr.append('\n');
      }
      else if (expr.charAt(0) == '#')
      {                      // Replicate comment.
        outputStr.append(expr);
        outputStr.append('\n');
      }
      else if (expr.charAt(0) == 'x')
      {                      // Loop command.
        newIndex = expr.indexOf(';');
        if (newIndex < 0)
        {                    // Semicolon not found.
          outputStr.append(expr);
          outputStr.append(": three semicolons expected but none found.\n");
          continue;
        }
        if (expr.length() < 3)
        {
          outputStr.append(expr);
          outputStr.append(": initial expression too short\n");
          continue;
        }
        firstN = expr.substring(2, newIndex);
        if (expr.charAt(1) != '=')
        {
          outputStr.append(expr);
          outputStr.append(": Syntax error on first expression.\n");
          continue;
        }
        if (!evaluateExpression(1, firstN, ExpressionResult))
        {
          continue;
        }
        N = ExpressionResult[0];
        newIndex2 = expr.indexOf(';', newIndex+1);
        if (newIndex2 < 0)
        {                    // Semicolon not found.
          outputStr.append(expr);
          outputStr.append(": three semicolons expected but only one found\n");
          continue;
        }
        nextN = expr.substring(newIndex+1, newIndex2);
        if (nextN.length() < 3 ||
            nextN.charAt(0) != 'x' || nextN.charAt(1) != '=')
        {
          outputStr.append(expr);
          outputStr.append(": Syntax error on next expression.\n");
          continue;
        }
        nextN = nextN.substring(2);
        newIndex3 = expr.indexOf(';', newIndex2+1);
        if (newIndex3 < 0)
        {                    // Semicolon not found.
          outputStr.append(expr);
          outputStr.append(": three semicolons expected but only two found\n");
          continue;
        }
        endExpression = expr.substring(newIndex2+1, newIndex3);
        ExpressionToCompute = expr.substring(newIndex3+1);
        counter = 1;
        for (;;)
        {
          if (!computeExpression(endExpression, N, counter,
                                 ExpressionResult))
          {
            break;
          }
          if (ExpressionResult[0].signum() > 0)
          {                  // Go out if the end expression is > 0.
            break;
          }
          if (!computeExpression(ExpressionToCompute, N, counter,
                                 ExpressionResult))
          {
            break;
          }
          factorExpression(ExpressionResult);
          if (!computeExpression(nextN, N, counter,
                                 ExpressionResult))
          {
            break;
          }
          N = ExpressionResult[0];
          counter++;
        }
      }
      else
      {
        if (!evaluateExpression(0, expr, ExpressionResult))
        {
          continue;
        }
        factorExpression(ExpressionResult);
      }
    }
    batchFinished = true;
    onlyFactoring = true;
    textNumber.setEditable(true);
    textNumber.setText("");
    upperTextArea.setText("");
    lowerTextArea.setText("");
    labelStatus.setText("");
    labelTop.setText
     ("Type number or numerical expression to factor here and press Return:");
  }

  private boolean computeExpression(String expr, BigInteger N, int counter,
                                    BigInteger [] ExpressionResult)
  {
    StringBuffer varsReplaced = new StringBuffer();
    for (int i=0; i<expr.length(); i++)
    {
      if (expr.charAt(i) == 'x')
      {
        varsReplaced.append(N.toString());
      }
      else if (expr.charAt(i) == 'i')
      {
        varsReplaced.append(counter);
      }
      else
      {
        varsReplaced.append(expr.charAt(i));
      }
    }
    return evaluateExpression(1, varsReplaced.toString(), ExpressionResult);
  }

  private boolean evaluateExpression(int type, String expr,
                                     BigInteger [] ExpressionResult)
  {
    int ExpressionRC;
    try
    {
      ExpressionRC =
        expression.ComputeExpression(
              expr,
              type,
              ExpressionResult);
    }
    catch (OutOfMemoryError e)
    {
      outputStr.append(expr);
      outputStr.append(": Out of memory\n");
      return false;
    }
    catch (ArithmeticException e)
    {
      return false;
    }
    if (ExpressionRC > 1)
    {
      outputStr.append(expr);
      outputStr.append(": "+expressionText[-2 - ExpressionRC]);
      outputStr.append('\n');
      return false;
    }
    return true;
  }

  private void factorExpression(BigInteger [] ExpressionResult)
  {
    int i;
    NumberToFactor = ExpressionResult[0];
    BigInteger AbsNumberToFactor = NumberToFactor.abs();
    if (batchPrime)
    {
      if (AbsNumberToFactor.compareTo(BigInt3) <= 0)
      {
        NbrFactors = 0;     // Indicate it is prime (2 or 3).
      }
      else if (BigInt3.modPow(AbsNumberToFactor.subtract(BigInt1),
                    AbsNumberToFactor).equals(BigInt1))
      {                // Pseudoprime
        if (AbsNumberToFactor.bitLength() < 34)
        {                    // Small number
          long modulus = AbsNumberToFactor.longValue();
          if (modulus % 2 == 0)
          {
            NbrFactors = 1;     // Indicate it is composite.
          }
          else
          {
            NbrFactors = 0;     // Indicate prime in advance.
            for (long Div = 3; Div*Div <= modulus; Div += 2)
            {
              if (modulus % Div == 0)
              {
                NbrFactors = 1;     // Indicate it is composite.
                break;
              }
            }
          }
        }
        else
        {                      // Large number.
          textNumber.setText(NumberToFactor.toString());
          startNewFactorization(true); // Request complete factorization
        }
      }
      else
      {                        // Pseudoprime test indicate composite.
        NbrFactors = 1;        // Indicate number composite.
      }
    }
    else
    {
      if (AbsNumberToFactor.compareTo(BigInt1) <= 0)
      {
      }
      else if (AbsNumberToFactor.bitLength() < 34)
      {                    // Small number
        long modulus = AbsNumberToFactor.longValue();
        int Expon = 0;
        i = 0;
        while (modulus % 2 == 0)
        {
          Expon++;
          modulus /= 2;
        }
        if (Expon > 0)
        {
          PD[0] = BigInt2;
          Exp[0] = Expon;
          i++;
        }
        long Div = 3;
        while (Div*Div <= modulus)
        {
          Expon = 0;
          while (modulus % Div == 0)
          {
            Expon++;
            modulus /= Div;
          }
          if (Expon > 0)
          {
            PD[i] = BigInteger.valueOf(Div);
            Exp[i] = Expon;
            i++;
          }
          Div += 2;
        }
        if (modulus > 1)
        {
          PD[i] = BigInteger.valueOf(modulus);
          Exp[i] = 1;
          i++;
        }
        NbrFactors = i;
      }
      else
      {
        textNumber.setText(NumberToFactor.toString());
        startNewFactorization(true);     // Request complete factorization
        if (numberIsNegative)
        {
          NumberToFactor = NumberToFactor.negate();
        }
      }
    }
    outputStr.append(NumberToFactor.toString());
    if (batchPrime)
    {               // Test primality
      if (AbsNumberToFactor.signum() == 0)
      {
        outputStr.append(" is a zero\n");
      }
      if (AbsNumberToFactor.compareTo(BigInt1) == 0)
      {
        outputStr.append(" is a unit\n");
      }
      else if (NbrFactors == 0)
      {
        outputStr.append(" is prime\n");
      }
      else
      {
        outputStr.append(" is composite\n");
      }
    }
    else
    {               // Perform complete factorization
      if (AbsNumberToFactor.compareTo(BigInt1) <= 0)
      {
        outputStr.append(" = "+NumberToFactor.toString()+"\n");
      }
      else
      {
        if (NumberToFactor.signum() < 0)
        {
          outputStr.append(" = -1 * ");
        }
        else
        {
          outputStr.append(" = ");
        }
        for (i=0; i<NbrFactors; i++)
        {
          if (i > 0)
          {
            outputStr.append(" * ");             
          }
          outputStr.append(PD[i].toString());
          if (Exp[i] > 1)
          {
            outputStr.append('^');
            outputStr.append(Exp[i]);
          }
        }
        outputStr.append('\n');
      }
    }
  }

  public String resultBatch()
  {
    if (batchFinished)
    {
      return outputStr.toString();
    }
    return "";
  }
  class Command implements ActionListener
  {
    int id;
    ecm appletEcm;

    public Command(int id, ecm appletEcm)
    {
      this.id = id;
      this.appletEcm = appletEcm;
    }

    public void actionPerformed(ActionEvent e)
    {
      int nbrEntered;
      switch (id)
      {
        case ACTION_TEXT_NBR :
          if (calcThread != null && calcThread.isAlive())
          {
            new AlertContinue(appletEcm); // Pop up alert
          }
          else
          {
            startNewFactorization(false);    // Start new factorization.
          }
          break;
        case ACTION_TEXT_CURVE :
        case ACTION_BTN_CURVE :
          try
          {
            nbrEntered = Integer.parseInt(textCurve.getText());
          }
          catch (Exception exc)
          {
            nbrEntered = -1;
          }
          if (nbrEntered >= 0 && nbrEntered < 100000000)
          {
            NextEC = nbrEntered;
            if (calcThread != null)
            {
              TerminateThread = true;
              try
              {
                // Wait until the factorization thread dies
                calcThread.join();
              }
              catch (InterruptedException ie)
              {
              }
            }
            if (NextEC == 0)
            {
              NextEC = TYP_SIQS + EC;
              forcedECM = false;
            }
            else
            {
              forcedECM = true;
            }
            calcThread = new Thread(appletEcm); // Start factorization thread
            calcThread.start();
          }
          break;
        case ACTION_TEXT_FACTOR :
        case ACTION_BTN_FACTOR :
          try
          {
            InputFactor = new BigInteger(textFactor.getText().trim());
          }
          catch (Exception exc)
          {
            InputFactor = BigInt0;
          }

          if (InputFactor.compareTo(BigInt1) > 0)
          {
            if (calcThread != null)
            {
              TerminateThread = true;
              try
              {
                // Wait until the factorization thread dies
                calcThread.join();
              }
              catch (InterruptedException ie)
              {
              }
            }
            InsertNewFactor(InputFactor);
            NextEC = EC;
            calcThread = new Thread(appletEcm); // Start factorization thread
            calcThread.start();
          }
          break;
      } // end switch
    } // end method
  } // end inner class

} /* end applet */

class AlertContinue extends Frame
{
  static final long serialVersionUID = 30L;

  public AlertContinue(ecm applet)
  {
    super("Warning!");
    Label label1, label2;
    Button buttonYes, buttonNo;
    ActionListener al = new AlertActionListener(applet);  
    setSize(400, 210);
    setVisible(true);
    setLayout(null);
    label1 = new Label("A factorization is still in progress!", Label.CENTER);
    label2 =
      new Label("Do you want to stop it and start a new one?", Label.CENTER);
    label1.setFont(new Font("Helvetica", Font.PLAIN, 12));
    label2.setFont(new Font("Helvetica", Font.PLAIN, 12));
    buttonYes = new Button("Yes");
    buttonNo = new Button("No");
    label1.setBounds(new Rectangle(25, 40, 350, 20));
    label2.setBounds(new Rectangle(25, 70, 350, 20));
    buttonYes.setBounds(new Rectangle(100, 130, 40, 30));
    buttonYes.setActionCommand("Yes");
    buttonYes.addActionListener(al);
    buttonNo.setBounds(new Rectangle(260, 130, 40, 30));
    buttonNo.setActionCommand("No");
    buttonNo.addActionListener(al);
    add(label1);
    add(label2);
    add(buttonYes);
    add(buttonNo);
    buttonNo.requestFocus();
    validate();
  }
}

class AlertActionListener implements ActionListener
{
  private ecm appletEcm;
  public AlertActionListener(ecm myApplet)
  {
    appletEcm = myApplet;
  }

  public void actionPerformed(ActionEvent e)
  {
    String s = e.getActionCommand();
    if (s.equals("Yes"))
    {
      appletEcm.startNewFactorization(false);  // Start new factorization.
    }
    ((Frame)((Component)e.getSource()).getParent()).dispose();
  }

} /* end class */

class PanelWindowListener implements WindowListener
{
  public void windowClosing(WindowEvent arg0)
  {
    System.exit(0);
  }

  public void windowOpened(WindowEvent arg0)
  {
  }

  public void windowClosed(WindowEvent arg0)
  {
  }

  public void windowIconified(WindowEvent arg0)
  {
  }

  public void windowDeiconified(WindowEvent arg0)
  {
  }

  public void windowActivated(WindowEvent arg0)
  {
  }

  public void windowDeactivated(WindowEvent arg0)
  {
  }
}
