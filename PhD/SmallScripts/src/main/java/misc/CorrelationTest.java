/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;

/**
 *
 * @author harmjan
 */
public class CorrelationTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String fileName = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodH8v2/BloodH8v2OriginalExpressionDataCorrectedFor4GWASPCs/ExpressionData.txt.gz";

        try {
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(fileName);
            ds.transposeDataset();
            double[][] correlationMatrix1 = new double[ds.nrRows][ds.nrRows];
            double[][] correlationMatrix2 = new double[ds.nrRows][ds.nrRows];

            Correlation ourcorrelations = new Correlation();

            for (int row = 0; row < ds.nrRows; row++) {
                double[] data1 = ds.rawData[row];
                for (int row2 = 0; row2 < ds.nrRows; row2++) {
                    double[] data2 = ds.rawData[row2];
                    double jscicorr = JSci.maths.ArrayMath.correlation(data1, data2);
                    double ourcorr = Correlation.correlate(data1, data2);
//                    if (jscicorr != ourcorr) {
                    System.out.println(jscicorr + "\t" + ourcorr + "\t" + Math.abs(ourcorr - jscicorr));
//                    }
                }
            }



        } catch (IOException e) {
            e.printStackTrace();

        }
    }
}
