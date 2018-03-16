/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author harmjan
 */
public class AnandAndiappanDataRewrite {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String infile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-11-19-AnandNeutrophils/2013-07-11-HT12v3SNPProbeCombos_replicated_r2_effects_only.txt";
        String outfile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-11-19-AnandNeutrophils/2013-07-11-HT12v3SNPProbeCombos-eQTLFile.txt";
        String expfile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-11-19-AnandNeutrophils/MeanExpression.txt";

        try {

            TextFile tf = new TextFile(infile, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            TextFile tfOut = new TextFile(outfile, TextFile.W);
            TextFile tfOut2 = new TextFile(expfile, TextFile.W);
            tfOut.writeln("PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR");
            while (elems != null) {
// original	rsid	array_address_id	probe_id	symbol	refseq	p_regression	geno1	geno2	geno3	mean1	mean2	mean3	mean1_log2	mean2_log2	mean3_log2	count1	count2	count3	comments
                String rs = elems[0];
                String probe = elems[2];
                String symbol = elems[4];
                double p = Double.parseDouble(elems[6]);
                String geno1 = elems[7].substring(0);
                String geno2 = elems[9].substring(0);

                int sum = 0;
                double meansum = 0;
                for (int q = 0; q < 3; q++) {
                    sum += Integer.parseInt(elems[16 + q]);
                    meansum += Double.parseDouble(elems[10] + q);
                }
                meansum /= 3;

                double multiplier = 1.0d / Math.log10(2.0d);
                meansum = (Math.log10(meansum) * multiplier);
                double z = Math.abs(ZScores.pToZ(p / 2d));

                double corr = ZScores.zScoreToCorrelation(z, sum);

                String output = p + "\t" + rs + "\t0\t0\t" + probe + "\t0\t0\tcis\t" + geno1 + "/" + geno2 + "\t" + geno1 + "\t" + z + "\tNeutrophils\t" + z + "\t" + sum + "\t" + meansum + "\t0\t" + symbol + "\t" + corr + "\t0\t0\t0\t0";
                tfOut.writeln(output);

                tfOut2.writeln(probe + "\t" + meansum + "\t" + 0);

                elems = tf.readLineElems(TextFile.tab);
            }
            tfOut.close();
            tfOut2.close();
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
