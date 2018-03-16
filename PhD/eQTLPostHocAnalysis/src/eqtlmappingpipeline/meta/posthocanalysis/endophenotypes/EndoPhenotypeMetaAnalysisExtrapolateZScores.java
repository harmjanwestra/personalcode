/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.endophenotypes;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class EndoPhenotypeMetaAnalysisExtrapolateZScores {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            EndoPhenotypeMetaAnalysisExtrapolateZScores p = new EndoPhenotypeMetaAnalysisExtrapolateZScores();
            String metafile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Meta/MetaMerged-40PCsTransEQTLs.txt";
            String outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Meta/MetaMerged-40PCsTransEQTLs-Extrapolated.txt";
            int newNrSamples = 5311;
            p.run(metafile, outfile, newNrSamples);
//            p.test();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String fileIn, String fileOut, int nrNewSamples) throws IOException {
        TextFile tf = new TextFile(fileIn, TextFile.R);
        TextFile tfo = new TextFile(fileOut, TextFile.W);
        String[] elems = tf.readLineElems(TextFile.tab); // skip header
        String newHeader = Strings.concat(elems, Strings.tab, 0, 6);
        for (int i = 6; i < elems.length; i += 6) {
            newHeader += "\t" + elems[i + 3];
            newHeader += "\t" + elems[i + 4];
        }
        tfo.writeln(newHeader);

        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String row = Strings.concat(elems, Strings.tab, 0, 6);
//            System.out.println(row);
            for (int i = 6; i < elems.length; i += 6) {
                double origP = Double.parseDouble(elems[i + 3]);
                double z = Double.parseDouble(elems[i + 4]);
                int nrSamples = Integer.parseInt(elems[i + 5]);

                // extrapolate
                System.out.println(origP + "\t" + z);
                double correlation3 = ZScores.zScoreToCorrelation(z, nrSamples);
                double newZ = ZScores.extrapolateZScore(nrSamples, nrNewSamples, z);
                double newP = ZScores.zToP(newZ);

                System.out.println("p " + origP + "\tz " + z + "\tn " + nrSamples + "\tnewz " + newZ + "\tnewp " + newP + "\tc " + correlation3);
                row += "\t" + newP + "\t" + newZ;

                
            }
            tfo.writeln(row);
//            System.exit(0);
            elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();
        tfo.close();
    }

    public void test() throws IOException {

        int nrSamples = 10000;
        cern.jet.random.StudentT tDistColt = new cern.jet.random.StudentT(nrSamples - 2, (new cern.jet.random.engine.DRand()));
        JSci.maths.statistics.TDistribution tDist = new JSci.maths.statistics.TDistribution(nrSamples - 2);

        System.out.println("CorrelationInitial\tCorrespondingZ\tCorrespondingP\tTValueBasedOnZ\tCorrelationBasedOnTValue\tMyCorr\tnewZ\tnewP");
        for (double correlation = -0.6; correlation <= 0.6; correlation += 0.01d) {

            //For a given correlation calculate the Z-Score:
            double t = correlation / (Math.sqrt((1 - correlation * correlation) / (double) (nrSamples - 2)));
            double pValue = 0;
            double obsZ = 0;
            if (t < 0) {
                pValue = tDistColt.cdf(t);
                if (pValue < 2.0E-323) {
                    pValue = 2.0E-323;
                }
                obsZ = cern.jet.stat.Probability.normalInverse(pValue);
            } else {
                pValue = tDistColt.cdf(-t); //Take two sided P-Value
                if (pValue < 2.0E-323) {
                    pValue = 2.0E-323;
                }
                obsZ = -cern.jet.stat.Probability.normalInverse(pValue);
            }

            //For a given Z-Score calculate the correlation:
            double obsPValue = 1;
            double obsT = 0;
            if (obsZ < 0) {
                obsPValue = cern.jet.stat.Probability.normal(obsZ);
                if (obsPValue == 0) {
                    obsPValue = Double.MIN_VALUE;
                }
                obsT = -tDist.inverse(obsPValue);
            } else {
                obsPValue = cern.jet.stat.Probability.normal(-obsZ);
                if (obsPValue == 0) {
                    obsPValue = Double.MIN_VALUE;
                }
                obsT = tDist.inverse(obsPValue);
            }
            double correlation2 = Math.sqrt(obsT * obsT / (obsT * obsT + nrSamples));
            if (obsT > 0) {
                correlation2 = -correlation2;
            }
            
            double correlation3 = ZScores.zScoreToCorrelation(obsZ, nrSamples);
            double newZ = ZScores.extrapolateZScore(nrSamples, nrSamples, obsZ);

            //Compare the previous correlation and the newly calculated correlation, are they the same?
            System.out.println(correlation + "\t" + obsZ + "\t" + (obsPValue*2) + "\t" + obsT + "\t" + correlation2 + "\t" + correlation3 + "\t" + newZ + "\t" + ZScores.zToP(newZ));


        }

    }
}
