package nl.harmjanwestra.playground.regressionslides;


import org.apache.commons.math3.distribution.NormalDistribution;

public class Main {

	public static void main(String[] args) {


		int n = 1000;

		double intercept = 0;
		double interceptPresent = 0;
		boolean noisePresent = true;
		double[] y = null;
		double[] x = null;
		double beta = 0.3;
		boolean grouped = true;
		boolean groupedOld = false;

//		int figure = 1;
//
//		if(figure == 1){
//			System.out.println("Val\tX\tY");
//			double vary = 1;
//			NormalDistribution dy = new NormalDistribution(0, vary);
//			y = dy.sample(n);
//			x = new double[n];
//			for (int i = 0; i < n; i++) {
//				x[i] = y[i] * 3;
//				System.out.println(i + "\t" + x[i] + "\t" + y[i]);
//			}
//		}


		if (!grouped) {
			double varx = 1;
			double vary = 1;
			NormalDistribution dx = new NormalDistribution(0, varx);
			NormalDistribution dnoise = new NormalDistribution(0, 1);
			// NormalDistribution dy = new NormalDistribution(0, vary);

			x = dx.sample(n);
//			y = dy.sample(n);
			y = new double[n];
			double[] noise = dnoise.sample(n);
			double meanY = 2;
			double sumy = 0;
			for (int i = 0; i < n; i++) {
				y[i] = meanY + (x[i] * beta);
				if (noisePresent) {
					y[i] += noise[i];
				}
//			x[i] -= intercept;
				sumy += y[i];
			}

			sumy /= n;


			System.out.println("Val\tX\tY");
			for (int i = 0; i < n; i++) {
//			x[i] = (y[i] / beta) + (noisePresent * noise[i]) + (intercept * interceptPresent);
				System.out.println(i + "\t" + x[i] + "\t" + y[i]);
			}
			System.out.println();
			System.out.println(sumy);
		} else if (groupedOld) {
			int[] groups = new int[n];
			int ngroups = 5;
			int samplesPerGroup = n / ngroups;
			int sctr = 0;
			int gctr = 0;
			for (int i = 0; i < groups.length; i++) {
				if (sctr == samplesPerGroup) {
					gctr++;
					sctr = 0;
				}
				groups[i] = gctr;
				sctr++;
//				System.out.println(i + "\t" + groups[i]);
			}
			x = new double[n];
			y = new double[n];
			NormalDistribution dnoise = new NormalDistribution(0, 0.5);
			double[] noise = dnoise.sample(n);

			boolean groupDependentX = false;
			boolean groupDependentVarianceX = false;
			boolean groupDependentVarianceY = false;

			boolean betaPerGroup = true;

			System.out.println("Val\tX\tY\tGroup");
			for (int i = 0; i < n; i++) {
				NormalDistribution dx = null;
				if (groupDependentX) {
					if (groupDependentVarianceX) {
						dx = new NormalDistribution(-groups[i] * 2, 1 + (groups[i]));
					} else {
						dx = new NormalDistribution(-groups[i] * 2, 1);
					}
				} else {
					dx = new NormalDistribution(0, 1);
				}

				NormalDistribution dy = new NormalDistribution(groups[i] * 2, 0.0001);
				double ddx = dx.sample();
				double ddy = dy.sample();

				x[i] = ddx;
				y[i] = ddy;

				if (betaPerGroup) {
					double groupBeta = 0 + groups[i] * 0.2;
					y[i] += (x[i] * groupBeta);
				} else {
					y[i] += (x[i] * beta);
				}

				if (noisePresent) {
					if (groupDependentVarianceY) {
						NormalDistribution dnoiseGroup = new NormalDistribution(0, 1 + (groups[i]));
						double noiseval = dnoiseGroup.sample();
						y[i] += noiseval;
					} else {
						y[i] += noise[i];
					}
				}
				System.out.println(i + "\t" + x[i] + "\t" + y[i] + "\t" + groups[i]);
			}

//			String header = "Val";
//			for (int i = 0; i < ngroups; i++) {
//				header += "\tgroup" + i;
//			}
//			System.out.println(header);


		} else {
			int[] groups = new int[n];
			int ngroups = 5;
			int samplesPerGroup = n / ngroups;
			int sctr = 0;
			int gctr = 0;
			for (int i = 0; i < groups.length; i++) {
				if (sctr == samplesPerGroup) {
					gctr++;
					sctr = 0;
				}
				groups[i] = gctr;
				sctr++;
//				System.out.println(i + "\t" + groups[i]);
			}
			x = new double[n];
			y = new double[n];

			System.out.println("Val\tX\tY\tGroup");

			boolean groupDependentMeanX = false;
			boolean groupDependentErrorY = true;
			boolean groupDependentSlope = true;

			double groupScalarMeanX = -4;
			double groupScalarMeanY = 4;
			double groupScalarVarY = 1;
			double groupScalarSlopeScalar = 2;
			double groupIndependentSlope = 0.3;

			for (int i = 0; i < n; i++) {
				NormalDistribution dx = null;
				if (!groupDependentMeanX) {
					dx = new NormalDistribution(0, 1);
				} else {
					dx = new NormalDistribution(groups[i] * groupScalarMeanX, 1);
				}
				NormalDistribution dnoise = null;
				if (groupDependentErrorY) {
					dnoise = new NormalDistribution(0, 1 + (groupScalarVarY * groups[i]));
				} else {
					dnoise = new NormalDistribution(0, 1);
				}

				x[i] = dx.sample();
				double meanY = groups[i] * groupScalarMeanY;
				double noise = dnoise.sample();
				double slope = groupIndependentSlope;
				if(groupDependentSlope){
					slope = groups[i] * groupIndependentSlope;
				}
				y[i] = meanY + x[i] * slope + noise;

				System.out.println(i + "\t" + x[i] + "\t" + y[i] + "\t" + groups[i]);
			}
		}


	}


}

