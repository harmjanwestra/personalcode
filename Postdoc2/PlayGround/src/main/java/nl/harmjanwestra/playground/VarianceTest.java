package nl.harmjanwestra.playground;

import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Heterogeneity;

public class VarianceTest {
	public static void main(String[] args) {
		double[] data = new double[1000];
		for (int i = 0; i < data.length; i++) {
			data[i] = Math.random();
			double v = Math.random();
			if (v < 0.5) {
				data[i] *= -1;
			}
			System.out.println(data[i]);
		}
		System.out.println();
		System.out.println(Descriptives.variance(data));

	}
}
