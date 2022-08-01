package nl.harmjanwestra.playground;

import JSci.maths.statistics.FDistribution;

public class Fdist {

	public static void main(String[] args) {

		double rss1 = 316.78546485815684;
		double rss2 = 311.63860948506397;
		int df1 = 9;
		int df2 = 10;
		int n = 2638;
		int dfd = n - df2;
		int dfn = df2 - df1;
		double f_value = ((rss1 - rss2) / dfn) / (rss2 / dfd);


		FDistribution f = new FDistribution(dfn, dfd);
		double p = f.cumulative(f_value);

		System.out.println(p);
		System.out.println(1 - p);


		double sumOfSquaresModelA = rss1;
		int degreesOfFreedomA = n-df1;
		int degreesOfFreedomB = n-df2;
		double sumOfSquaresModelB = rss2;

		double meanSquareError = sumOfSquaresModelA / degreesOfFreedomA;
		int degreesOfFreedomDifference = Math.abs(degreesOfFreedomB - degreesOfFreedomA);
		double meanSquareErrorDiff = Math.abs((sumOfSquaresModelA - sumOfSquaresModelB) / (degreesOfFreedomDifference));
		double Fval = meanSquareErrorDiff / meanSquareError;
		FDistribution Fdist = new FDistribution(degreesOfFreedomDifference, degreesOfFreedomB);
		double pval = 1 - Fdist.cumulative(Fval);
		System.out.println(pval);
	}
}
