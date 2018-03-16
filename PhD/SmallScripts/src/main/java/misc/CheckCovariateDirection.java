/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class CheckCovariateDirection {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String[] covariateFiles = new String[]{
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/MarjoleinVisit3/2013-10-10-ExpressionDataForInteractionTerms/Covariate.txt",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/2013-10-10-BloodHT12DataForInteractionTerms/Covariate.txt",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/2013-10-10-EGCUTExpressionDataForInteractionTerms/Covariate.txt"};
            String[] expressionDatasets = new String[]{
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/MarjoleinVisit3/2013-10-10-ExpressionDataForInteractionTerms/ExpressionData.QuantileNormalized.Log2Transformed.txt.gz",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/2013-10-10-BloodHT12DataForInteractionTerms/ExpressionData.QuantileNormalized.Log2Transformed.txt.gz",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/2013-10-10-EGCUTExpressionDataForInteractionTerms/ExpressionData.QuantileNormalized.Log2Transformed.txt.gz"
            };

            for (int i = 0; i < covariateFiles.length; i++) {
                DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(expressionDatasets[i]);
                String s = covariateFiles[i];
                TextFile tf = new TextFile(s, TextFile.R);
                int nrPos = 0;
                int nrNeg = 0;
                String header = tf.readLine();
                String[] elems = tf.readLineElems(TextFile.tab);
                double[] pcScores = new double[ds.nrCols];
                while (elems != null) {
                    String sample = elems[0];
                    Double d = Double.parseDouble(elems[1]);
                    Integer id = ds.hashCols.get(sample);
                    pcScores[id] = d;
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();

                int nrProbesCorrelatingPositively = 0;
                for (int row = 0; row < ds.nrRows; row++) {
                    double[] data = ds.rawData[row];
                    double corr = JSci.maths.ArrayMath.correlation(pcScores, data);
                    if (corr >= 0) {
                        nrProbesCorrelatingPositively++;
                    } else {
                        nrProbesCorrelatingPositively--;
                    }
                }

                TextFile out = new TextFile(s + "-DirectionCorrected.txt", TextFile.W);

                for (int j = 0; j < pcScores.length; j++) {
                    if (nrProbesCorrelatingPositively < 0) {
                        out.writeln(ds.colObjects.get(j) + "\t" + (-pcScores[j]));
                    } else {
                        out.writeln(ds.colObjects.get(j) + "\t" + (pcScores[j]));
                    }
                }
                out.close();

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
