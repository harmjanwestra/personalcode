package nl.harmjanwestra.playground;

import cern.jet.stat.tdouble.Probability;
import umcg.genetica.math.stats.ZScores;

public class Logit {

    public static void main(String[] args) {
        double x = 30;
        double df = 5;
        double p = Probability.chiSquareComplemented(df, x);
        System.out.println(p);


    }
}
