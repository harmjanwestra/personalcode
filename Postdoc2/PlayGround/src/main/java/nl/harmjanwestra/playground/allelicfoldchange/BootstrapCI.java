package nl.harmjanwestra.playground.allelicfoldchange;

import org.apache.commons.math3.distribution.NormalDistribution;

public class BootstrapCI {

	/* Program to accompany the article
   "The bias-corrected and accelerated (BCa) bootstrap interval"
   by Rick Wicklin, published 12JUL2017 on The DO Loop blog,
   http://blogs.sas.com/content/iml/2017/07/12/bca-bootstrap-interval.html

   Chernick and LaBudde (p. 78) state,
   "Unfortunately, in small to moderate samples for asymmetric or heavy-tailed
   distributions, the percentile method is not very good and so modifications
   are required to improve it."

   (p. 85) "The BCa method incorporates a parameter that Efron named the
   'acceleration constant.'  It is based on the third moment and corrects for skewness."

   (p. 86) "The method has a second-order accuracy property."

   Taken from: https://blogs.sas.com/content/iml/files/2017/07/bootBCa.txt
   */

	/*****************************************************/


	NormalDistribution normalDistribution = new NormalDistribution();

	public static void main(String[] args) {
		NormalDistribution normalDistribution = new NormalDistribution();
		System.out.println(normalDistribution.inverseCumulativeProbability(.975));
	}

	// expect samples on rows
	public double[][] JackSampMat(double[][] x) {
		int n = x.length;
		double[][] sample = new double[n - 1][];
		for (int i = 0; i < n - 1; i++) {
			sample[i] = x[i];
		}
		return sample;
	}


	/* compute bias-correction factor from the proportion of bootstrap estimates
   that are less than the observed estimate
*/
	public double bootBC(double[][] bootEst, double[][] est) {
		int B = bootEst.length * bootEst[0].length;
		double propLess = ctless(bootEst, est) / B;
		double z0 = normalDistribution.inverseCumulativeProbability(propLess);
		return z0;
	}

	private double ctless(double[][] bootEst, double[][] est) {
		int ctr = 0;
		for (int r = 0; r < bootEst.length; r++) {
			for (int c = 0; c < bootEst.length; c++) {
				if (bootEst[r][c] < est[r][c]) {
					ctr++;
				}
			}
		}
		return ctr;
	}


	/* compute acceleration factor, which is related to the skewness of bootstrap estimates.
   Use jackknife replicates to estimate.
*/
	public void bootAccel(double[][] x) {
		double[][] M = JackSampMat(x);
		double[] jStat = EvalStat(M); // eval returns double[]
//		double jackEst = mean(jStat);
//		double num = sum((jackEst - jStat)##3); // ## means ^3
//		double den = sum((jackEst - jStat)##2);
//		double ahat = num / (6 * den##(3 / 2));        /* ahat based on jackknife ==> not random */
//		return ahat;

	}

	private double[] mean(double[][] x) {
		double[] mean = new double[x.length];
		for (int r = 0; r < x.length; r++) {
			for (int c = 0; c < x[r].length; c++) {
				mean[r] += x[r][c];
			}
			mean[r] /= x[r].length;
		}
		return mean;
	}

	private double[][] transpose(double[] x) {
		double[][] t = new double[1][x.length];
		for (int d = 0; d < x.length; d++) {
			t[0][d] = x[d];
		}
		return t;
	}


	// let's assume our data is not multinomial for now
	public double[] EvalStat(double[][] M) {
		double[] estimates = new double[M[0].length];
		return estimates;
	}
}
