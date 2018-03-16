/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package rsquarifier;

import java.io.IOException;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ReperformMetaAnalysis {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        /*
         * 0	PValue
         1	SNPName
         2	SNPChr
         3	SNPChrPos
         4	ProbeName
         5	ProbeChr
         6	ProbeCenterChrPos
         7	CisTrans
         8	SNPType
         9	AlleleAssessed
         10	OverallZScore
         11	DatasetsWhereSNPProbePairIsAvailableAndPassesQC
         12	DatasetsZScores
         13	DatasetsNrSamples
         14	IncludedDatasetsMeanProbeExpression
         15	IncludedDatasetsProbeExpressionVariance
         16	HGNCName
         17	IncludedDatasetsCorrelationCoefficient
         18	Meta-Beta
         19	(SE)
         20	Beta
         21	(SE)
         22	FoldChange
         */
        try {
            double[] weights = new double[]{1.149324968, 0.903980852, 0.742652579, 1.404928207, 1.00382799, 1.048156609, 0.747128795, 1, 1};
            for (int perm = 0; perm < 11; perm++) {
                Gpio.createDir("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCReWeightedZ-Batch2/");
                String infile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz";
                String outfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCReWeightedZ-Batch2/eQTLs.txt.gz";
                String outfile2 = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCReWeightedZ-Batch2/eQTLs-Diff.txt.gz";

                if (perm > 0) {
                    infile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/PermutedEQTLsPermutationRound" + perm + ".txt.gz";
                    outfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCReWeightedZ-Batch2/PermutedEQTLsPermutationRound" + perm + ".txt.gz";
                    outfile2 = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCReWeightedZ-Batch2/PermutedEQTLsPermutationRound" + perm + "-Diff.txt.gz";
                }

                System.out.println(infile);
                TextFile tf = new TextFile(infile, TextFile.R);
                TextFile tfout = new TextFile(outfile, TextFile.W);
                TextFile tfout2 = new TextFile(outfile2, TextFile.W);

                String header = tf.readLine();
                tfout.writeln(header);
                tfout2.writeln(header);
                int ctr = 0;
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {
                    String datasetZ = elems[12];
                    String datasetS = elems[13];

                    String[] zelems = datasetZ.split(",");
                    String[] sElems = datasetS.split(",");

                    double[] z = new double[zelems.length];
                    int[] samples = new int[zelems.length];
                    for (int i = 0; i < zelems.length; i++) {
                        if (zelems[i].equals("-")) {
                            samples[i] = 0;
                            z[i] = Double.NaN;
                        } else {
                            samples[i] = Integer.parseInt(sElems[i]);
                            z[i] = Double.parseDouble(zelems[i]);
                        }
                    }

                    
                    double metaZ = ZScores.getWeightedZ(z, samples, weights);
                    double metaP = Descriptives.convertZscoreToPvalue(metaZ);
                    
                    String origP = elems[0];
                    String origZ = elems[10];
                    elems[0] = "" + metaP;
                    elems[10] = "" + metaZ;

                    String output = Strings.concat(elems, Strings.tab);
                    if (ctr < 1000) {
                        elems[8] = origP;
                        elems[9] = origZ;
                        tfout2.writeln(Strings.concat(elems, Strings.tab));
                    }

                    if (ctr % 100000 == 0) {
                        System.out.println(ctr + " lines parsed..");
                    }
                    tfout.writeln(output);
                    ctr++;
                    elems = tf.readLineElems(TextFile.tab);
                }
                tfout2.close();
                tfout.close();
                tf.close();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
