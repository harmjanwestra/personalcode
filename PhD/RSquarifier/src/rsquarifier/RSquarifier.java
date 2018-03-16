/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package rsquarifier;

import java.io.IOException;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author harmjan
 */
public class RSquarifier {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            String infile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLProbesFDR0.05-ProbeLevel.txt.gz";
            String outfile = "/Volumes/iSnackHD/Data/Projects/2013-03-RandomEffectMetaAnalysis/";
            Gpio.createDir(outfile);
            outfile += "2013-03-25-RSquaresForMetaAnalysis.txt";

            TextFile tf = new TextFile(infile, TextFile.R);
            TextFile tfout = new TextFile(outfile, TextFile.W);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            tfout.writeln("Dataset Samplesizes: 891,963,1240,229,762,509,611");
            tfout.writeln("Pval\tSNP\tProbe\tMetaZ\tMetaR2\tEGCUT Z\tEGCUT R2\tSHIP_TREND Z\tSHIP_TREND R2\tGroningen-HT12 Z\tGroningen-HT12 R2\tGroningen-H8v2 Z\tGroningen-H8v2 R2\tRotterdam Z\tRotterdam R2\tDILGOM Z\tDILGOM R2\tINCHIANTI Z\tINCHIANTI R2");
            int numWithAll9Ds = 0;
            while (elems != null) {

                String pval = elems[0];
                String snp = elems[1];
                String probe = elems[4];

                String metazstr = elems[10];
                String datasets = elems[11];
                String datasetszstr = elems[12];
                String dastasetsize = elems[13];
                String hugo = elems[15];

                String[] zelems = datasetszstr.split(",");
                String[] selems = dastasetsize.split(",");
                int ctr = 0;
                double[] zs = new double[zelems.length - 2];
                int[] samplesizes = new int[zelems.length - 2];

                boolean allsamplespresent = true;
                int totalSamples = 0;
                int additionalsamples = 0;
                for (int i = 0; i < zelems.length; i++) {
//                    System.out.println(z);
                    String zstr = zelems[i];

                    if (i < zelems.length - 2) {
                        if (zstr.equals("-")) {
                            ctr++;
                            allsamplespresent = false;
                        } else {
                            zs[i] = Double.parseDouble(zstr);
                            samplesizes[i] = Integer.parseInt(selems[i]);
                            totalSamples += samplesizes[i];
                        }
                    } else {
                        if (zstr.equals("-")) {
                        } else {
                            additionalsamples += Integer.parseInt(selems[i]);
                        }
                    }

                }

                Correlation.correlationToZScore(totalSamples);
                if (allsamplespresent) {

                    // perform r2 conversion
                    double metaz = Double.parseDouble(metazstr);
                    double metar = ZScores.zScoreToCorrelation(metaz, totalSamples + additionalsamples);
                    metar *= metar;
                    tfout.append(pval);
                    tfout.append("\t");
                    tfout.append(snp);
                    tfout.append("\t");
                    tfout.append(probe);
                    tfout.append("\t");
                    tfout.append(metazstr);
                    tfout.append("\t");
                    tfout.append("" + metar);
                    for (int i = 0; i < zs.length; i++) {
                        
                        double r = ZScores.zScoreToCorrelation(zs[i], samplesizes[i]);
                        double r2 = r * r;

                        double diff = Math.abs(r2-metar);
                        double ratio = metar/r2;
                        tfout.append("\t" + zs[i]);
                        tfout.append("\t" + r2);
                        tfout.append("\t" + diff);

                    }

                    tfout.append("\n");
                    numWithAll9Ds++;
                }


                elems = tf.readLineElems(TextFile.tab);
            }


            tfout.close();
            tf.close();

            System.out.println(numWithAll9Ds);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
