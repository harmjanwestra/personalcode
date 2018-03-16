/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class MatrixGetSubSetOfColumns {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String file1 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-Groningen/Normalization/ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.CovariatesRemoved.txt.gz";
        String file2 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.40PCAsOverSamplesRemoved.txt.gz";
        String query = "2000128";
        String query2 = "130717";

        try {

            DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>(file1);



            DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(file2);
            ds2.recalculateHashMaps();
            Integer rowId11 = ds1.hashRows.get(query);
            Integer rowId12 = ds1.hashRows.get(query2);
            Integer rowId21 = ds2.hashRows.get(query);
            Integer rowId22 = ds2.hashRows.get(query2);



            for (int i = 0; i < ds2.nrCols; i++) {
                String sample = ds2.colObjects.get(i);
                Integer otherDsSampleId = ds2.hashCols.get(sample);
                if (otherDsSampleId != null) {
                    System.out.println(sample + "\t" + ds1.rawData[rowId11][i] + "\t" + ds1.rawData[rowId12][i] + "\t" + ds2.rawData[rowId21][otherDsSampleId] + "\t" + ds2.rawData[rowId22][otherDsSampleId]);

                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}