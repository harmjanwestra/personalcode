import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Heterogeneity;
import umcg.genetica.math.stats.ZScores;

public class WeightedCorrelation2 {

	public static void main(String[] args) {


		double[] x = new double[100];
		double[] y = new double[100];


		for (int i = 0; i < x.length; i++) {
			x[i] = Math.random() + Math.random() + Math.random() + Math.random();
			y[i] = Math.random() + Math.random() + Math.random() + Math.random() + x[i] * 0.3;
		}

		double meanx = Descriptives.mean(x);
		double meany = Descriptives.mean(y);
		double varx = Descriptives.variance(x);
		double vary = Descriptives.variance(y);


		Correlation c = new Correlation();
		double r = Correlation.correlate(x, y);
		for (int nrMissing = 0; nrMissing < 50; nrMissing += 10) {
			double meanxAct = meanx;
			double meanyAct = meany;
			double propmissing = (double) nrMissing / x.length;
			double varActX = Descriptives.variance(x, meanxAct) / (1 - propmissing);
			double varActY = Descriptives.variance(y, meanyAct) / (1 - propmissing);
			double r2 = Correlation.correlate(x, y, meanxAct, meanyAct, varActX, varActY);
			double r3 = WeightedCorrelation2.correlate(x, y, meanxAct, meanyAct, varx, vary, x.length - nrMissing);
			double r4 = WeightedCorrelation2.correlate(x, y, meanxAct, meanyAct, varx, vary, x.length + nrMissing);
			System.out.println(r + "\t" + r2 + "\t" + r3 + "\t" + r4);
		}


//		Heterogeneity j

	}

	public static double mean(double[] x) {
		return Descriptives.mean(x);
	}

	public static double correlate(double[] a1, double[] a2, double mean1, double mean2, double var1, double var2, double n) {
		double denom = Math.sqrt(var1 * var2);
		if (denom == 0.0) {
			return var1 == 0.0 && var2 == 0.0 ? 1.0 : 0.0;
		} else if (a1.length != a2.length) {
			throw new IllegalArgumentException("Arrays must have the same length : " + a1.length + ", " + a2.length);
		} else {
			double ans = 0.0;

			for (int i = 0; i < a1.length; ++i) {
				ans += (a1[i] - mean1) * (a2[i] - mean2);
			}

			return ans / (double) (n - 1) / denom;
		}
	}
}
