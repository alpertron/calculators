import java.awt.*;
import java.math.*;

interface factorApplet
{
  String InitSIQSStrings(BigInteger NbrToFactor, int SieveLimit);
  String getMultAndFactorBase(int multiplier, long FactorBase);
  void showMatrixSize(String SIQSInfoText, int rows, int cols);
  void ShowSIQSInfo(long time, int congruencesFound, int matrixBLength, long t);
  void ModInvBigNbr(int[] a, int[] inv, int[] b, int NumberLength);
  void GcdBigNbr(int Nbr1[], int Nbr2[], int Gcd[], int NumberLength);
  void saveSIQSStatistics(long polynomialsSieved, long trialDivisions,
                          long smoothsFound, long totalPartials,
                          long partialsFound, long ValuesSieved);
  boolean getTerminateThread();
  void writeLowerPane(String data);
  void setOld(long old);
  long getOld();
  void setOldTimeElapsed(long old);
  long getOldTimeElapsed();
}
