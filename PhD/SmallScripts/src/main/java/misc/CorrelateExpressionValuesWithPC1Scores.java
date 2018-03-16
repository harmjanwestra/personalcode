/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CorrelateExpressionValuesWithPC1Scores {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            String[] pc1files = new String[]{
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/MarjoleinVisit3/2013-10-10-ExpressionDataForInteractionTerms/Covariate.txt-DirectionCorrected.txt",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/2013-10-10-BloodHT12DataForInteractionTerms/Covariate.txt-DirectionCorrected.txt",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/2013-10-10-EGCUTExpressionDataForInteractionTerms/Covariate.txt-DirectionCorrected.txt"};
            String[] datasets = new String[]{
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/MarjoleinVisit3/2013-10-10-ExpressionDataForInteractionTerms/ExpressionData.QuantileNormalized.Log2Transformed-HT12v3.txt.gz",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/2013-10-10-BloodHT12DataForInteractionTerms/ExpressionData.QuantileNormalized.Log2Transformed.txt.gz",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/2013-10-10-EGCUTExpressionDataForInteractionTerms/ExpressionData.QuantileNormalized.Log2Transformed.txt.gz"
            };

            String outfilename = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/CorrelationsWithPC1-MetaZ.txt";

            HashSet<String> probes = new HashSet<String>();

            int[] sampleSizes = new int[pc1files.length];
            ArrayList<HashMap<String, Double>> zScoresPerDataset = new ArrayList<HashMap<String, Double>>();
            for (int d = 0; d < datasets.length; d++) {
                HashMap<String, Double> zscoreperprobe = new HashMap<String, Double>();
                DoubleMatrixDataset<String, String> expressionMatrix = new DoubleMatrixDataset<String, String>(datasets[d]);
                DoubleMatrixDataset<String, String> pcScoreFile = new DoubleMatrixDataset<String, String>(pc1files[d]);
                pcScoreFile.transposeDataset();

                double[] vals = new double[pcScoreFile.nrCols];
                for (int col = 0; col < pcScoreFile.nrCols; col++) {
                    String sample = pcScoreFile.colObjects.get(col);
                    Integer id = expressionMatrix.hashCols.get(sample);
                    vals[id] = pcScoreFile.rawData[0][col];
                }

                for (int row = 0; row < expressionMatrix.nrRows; row++) {
                    String probe = expressionMatrix.rowObjects.get(row);
                    double correlation = JSci.maths.ArrayMath.correlation(expressionMatrix.rawData[row], vals);
                    Correlation.correlationToZScore(vals.length);
                    double z = Correlation.convertCorrelationToZScore(vals.length, correlation);
                    zscoreperprobe.put(probe, z);
                    probes.add(probe);
                }

                zScoresPerDataset.add(zscoreperprobe);
                sampleSizes[d] = vals.length;

            }

            System.out.println(probes.size() +" probes detected");
            
            TextFile out = new TextFile(outfilename, TextFile.W);
            out.writeln("Probe\tMetaZ\tIndividualZScores");
            for (String probe : probes) {
                double[] zScores = new double[pc1files.length];
                int[] samples = new int[pc1files.length];
                for (int d = 0; d < zScores.length; d++) {
                    HashMap<String, Double> zs = zScoresPerDataset.get(d);
                    Double z = zs.get(probe);
                    if (z == null) {
                        zScores[d] = Double.NaN;
                        samples[d] = 0;
                    } else {
                        zScores[d] = z;
                        samples[d] = sampleSizes[d];
                    }
                }
                double metaZ = ZScores.getWeightedZ(zScores, sampleSizes);
                out.writeln(probe + "\t" + metaZ +"\t"+Strings.concat(zScores, Strings.comma));
            }
            out.close();


        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
