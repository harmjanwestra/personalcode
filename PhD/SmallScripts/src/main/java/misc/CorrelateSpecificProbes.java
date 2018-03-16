/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class CorrelateSpecificProbes {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String matrix1 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/2013-10-10-BloodHT12DataForInteractionTerms/ExpressionData.QuantileNormalized.Log2Transformed.txt.gz";
        String matrix2 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/2013-10-10-BloodHT12DataForInteractionTerms/ExpressionData.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors.txt.gz";
        HashSet<String> probes = null;
        HashSet<String> probes2 = new HashSet<String>();
        probes2.add("Comp1");
        try {
            CorrelateSpecificProbes c = new CorrelateSpecificProbes();
            c.run(matrix1, matrix2, probes, probes2);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String matrix1, String matrix2, HashSet<String> probes, HashSet<String> probes2) throws IOException {

        DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>(matrix1, probes);
        DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(matrix2, null, probes2);
        ds2.transposeDataset();

        Correlation correlation = new Correlation();

        for (int row1 = 0; row1 < ds1.nrRows; row1++) {

            for (int row2 = 0; row2 < ds2.nrRows; row2++) {
                ArrayList<Double> values1 = new ArrayList<Double>();
                ArrayList<Double> values2 = new ArrayList<Double>();
                for (int col1 = 0; col1 < ds1.nrCols; col1++) {
                    String ind = ds1.colObjects.get(col1);
                    Integer indid2 = ds2.hashCols.get(ind);
                    if (indid2 != null) {
                        values1.add(ds1.rawData[row1][col1]);
                        values2.add(ds2.rawData[row2][indid2]);
                    }
                }
                double[] vals1 = Primitives.toPrimitiveArr(values1.toArray(new Double[0]));
                double[] vals2 = Primitives.toPrimitiveArr(values2.toArray(new Double[0]));
                double r = Correlation.correlate(vals1, vals2);
                System.out.println(ds1.rowObjects.get(row1) + "\t" + ds2.rowObjects.get(row2) + "\t" + r);
            }
        }

    }
}
