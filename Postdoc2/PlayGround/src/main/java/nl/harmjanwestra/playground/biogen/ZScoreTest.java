package nl.harmjanwestra.playground.biogen;

import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;

public class ZScoreTest {


    public static void main(String[] args) {

        Correlation.correlationToZScore(8000);
        double z = Correlation.convertCorrelationToZScore(8000, 0.05);
        double p = ZScores.zToP(z);
        System.out.println(z + "\t" + p);
        System.out.println(z + "\t" + (p/2));

    }
}
