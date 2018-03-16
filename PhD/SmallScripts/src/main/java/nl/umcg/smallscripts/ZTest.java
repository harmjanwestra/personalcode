/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author harmjan
 */
public class ZTest {
    /*
     * double weightedZ = 0;
     int totalSample = 0;
     for (int d = 0; d < datasetZ.length; d++) {
     if (datasetZ[d] != null) {
     weightedZ += Math.sqrt(datasetWeights[d]) * datasetZ[d];
     totalSample += datasetWeights[d];
     }
     }

     //        System.out.println("WZ: " + weightedZ);
     //        System.out.println("TotalSample: " + totalSample);

     double hetSum = 0;
     int hetDf = 0;
     for (int d = 0; d < datasetZ.length; d++) {
     if (datasetZ[d] != null) {
     double expectedZ = Math.sqrt(datasetWeights[d]) * weightedZ / totalSample;

     hetSum += (datasetZ[d] - expectedZ) * (datasetZ[d] - expectedZ);
     //                System.out.println(d + "\tz:\t" + datasetZ[d] + "\tez:\t" + expectedZ + "\tdiff:\t" + (datasetZ[d] - expectedZ) + "\tHetSum:\t" + hetSum);
     hetDf++;
     }
     }
     */

    public static void main(String[] args) {
//        int nrSamples = 2000;
//        double obsZ = 23.9718d;
//        int obsZInt = (int) Math.round(obsZ * 10000d + 1000000d);
//
//        cern.jet.random.StudentT tDistColt = new cern.jet.random.StudentT(nrSamples - 2, (new cern.jet.random.engine.DRand()));
//        double[] zScoreToCorrelation = new double[2000001];
//        for (int corrInt = 0; corrInt < 2000001; corrInt++) {
//            double correlation = (double) (corrInt - 1000000) / 1000000d;
//            double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (nrSamples - 2)));
//            double pValue = 0;
//            double zScore = 0;
//            if (t < 0) {
//                pValue = tDistColt.cdf(t);
//                if (pValue < 2.0E-323) {
//                    pValue = 2.0E-323;
//                }
//                zScore = cern.jet.stat.Probability.normalInverse(pValue);
//            } else {
//                pValue = tDistColt.cdf(-t); //Take two sided P-Value
//                if (pValue < 2.0E-323) {
//                    pValue = 2.0E-323;
//                }
//                zScore = -cern.jet.stat.Probability.normalInverse(pValue);
//            }
//            int zScoreInt = (int) Math.round(zScore * 10000d + 1000000d);
//            if (zScoreInt < 0) {
//                zScoreInt = 0;
//            }
//            if (zScoreInt > 2000000) {
//                zScoreInt = 2000000;
//            }
//            zScoreToCorrelation[zScoreInt] = correlation;
//        }
//
//        double correlation = zScoreToCorrelation[obsZInt];
//
//        nrSamples = 5311;
//        tDistColt = new cern.jet.random.StudentT(nrSamples - 2, (new cern.jet.random.engine.DRand()));
//        double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (nrSamples - 2)));
//        double pValue = 0;
//        double zScore = 0;
//        if (t < 0) {
//            pValue = tDistColt.cdf(t);
//            if (pValue < 2.0E-323) {
//                pValue = 2.0E-323;
//            }
//            zScore = cern.jet.stat.Probability.normalInverse(pValue);
//        } else {
//            pValue = tDistColt.cdf(-t); //Take two sided P-Value
//            if (pValue < 2.0E-323) {
//                pValue = 2.0E-323;
//            }
//            zScore = -cern.jet.stat.Probability.normalInverse(pValue);
//        }
//        System.out.println(zScore);

        
//        ZScores.convertZScoreToCorrelation(5000, 3.5);

    }
}
