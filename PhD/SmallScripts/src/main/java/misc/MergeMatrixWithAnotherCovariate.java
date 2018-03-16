/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.Gpio;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class MergeMatrixWithAnotherCovariate {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String matrix1 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
            String matrix2 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/MetaAnalysis-PC1Only/MetaAnalysisZScoreMatrix.txt";
            String matrixOut = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/MetaAnalysis-PC1OnlyMerged/MetaAnalysisZScoreMatrix.txt";
            Gpio.createDir("/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/MetaAnalysis-PC1OnlyMerged/");

            DoubleMatrixDataset<String, String> metamatrix1 = new DoubleMatrixDataset<String, String>(matrix1);
            DoubleMatrixDataset<String, String> metamatrix2 = new DoubleMatrixDataset<String, String>(matrix2);

            for (int col = 0; col < metamatrix1.nrCols; col++) {
                String colObj = metamatrix1.colObjects.get(col);
                if (metamatrix2.hashCols.get(colObj) == null) {
                    System.err.println("Missing eQTL: " + colObj);
                }
            }
//            
            double[][] fullnewmatrix = new double[metamatrix1.nrRows + metamatrix2.nrRows][metamatrix1.nrCols];
            System.out.println("size: " + (metamatrix1.nrRows + metamatrix2.nrRows) + "x" + metamatrix1.nrCols);
            for (int row = 0; row < metamatrix1.nrRows; row++) {
                fullnewmatrix[row] = metamatrix1.rawData[row];
            }

            int rowCtr = metamatrix1.nrRows;
            for (int row = 0; row < metamatrix2.nrRows; row++) {
                for (int col = 0; col < metamatrix2.nrCols; col++) {
                    String colname = metamatrix2.colObjects.get(col);
                    Integer id = metamatrix1.hashCols.get(colname);
                    if (id != null) {
                        fullnewmatrix[rowCtr][id] = metamatrix2.rawData[row][col];
                    }
                }
                metamatrix1.rowObjects.add("PC1DirectionCorrected" + metamatrix2.rowObjects.get(row));
                rowCtr++;
            }

            DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>();
            dsout.rawData = fullnewmatrix;
            dsout.colObjects = metamatrix1.colObjects;
            dsout.rowObjects = metamatrix1.rowObjects;
            dsout.recalculateHashMaps();
            dsout.save(matrixOut);




        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
