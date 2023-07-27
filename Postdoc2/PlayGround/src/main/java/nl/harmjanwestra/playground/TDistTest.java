package nl.harmjanwestra.playground;

import cern.jet.random.tdouble.StudentT;
import cern.jet.random.tdouble.engine.DRand;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import cern.jet.stat.tdouble.Probability;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;

public class TDistTest {

	public static void main(String[] args) {
		DoubleRandomEngine randomEngine = new DRand();
		int nrSamples = 100;
		StudentT tDistColt = new StudentT((double) (nrSamples - 2), randomEngine);

		for (double spearman = -0.9; spearman < 1; spearman += 0.1) {
			double t = spearman / Math.sqrt((1.0 - spearman * spearman) / (double) (nrSamples - 2));
			double p = 1;
			double porig = 1;
			double z = 0;
			if (t < 0) {
				p = tDistColt.cdf(t);
				porig = tDistColt.cdf(t);
			} else {
				System.out.println("flip: " + spearman + "\t" + t);
				p = tDistColt.cdf(-t);
				porig = tDistColt.cdf(t);
			}
			double ppdf = tDistColt.pdf(t);
			if (p < 2.0E-323) {
				p = 2.0E-323;
			}
			z = Probability.normalInverse(p);
			double pz = ZScores.zToP(z);
			System.out.println(spearman + "\tt: " + t + "\tp: " + p + "\tporig: " + porig + "\tz: " + z + "\tppdf: " + ppdf + "\tpz" + pz);
		}
	}

}
