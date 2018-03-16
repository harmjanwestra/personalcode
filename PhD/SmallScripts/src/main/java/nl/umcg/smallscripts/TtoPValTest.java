/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import cern.jet.random.StudentT;
import umcg.genetica.math.stats.TTest;

/**
 *
 * @author harmjan
 */
public class TtoPValTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        int df = 14000;

        for (double t = 0; t < 10000; t++) {
            StudentT tDistColt = new StudentT(df, (new cern.jet.random.engine.DRand()));

            double tTestPValue1 = tDistColt.cdf(t);

            if (tTestPValue1 > 0.5) {
                tTestPValue1 = 1 - tTestPValue1;
            }
            tTestPValue1 *= 2;

            System.out.println(t + "\t" + tTestPValue1);
        }

    }
}
