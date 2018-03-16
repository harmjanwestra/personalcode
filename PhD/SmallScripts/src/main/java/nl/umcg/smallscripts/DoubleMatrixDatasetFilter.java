/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.smallscripts;

import java.io.IOException;
import java.util.Set;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class DoubleMatrixDatasetFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String filterRows = null;
        String filterCols = "/Volumes/iSnackHD/Data/Projects/JosephPowell/Data2/EGCUT/EGCUTExpSamples.txt";
        String matrixIn = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/RawExpressionData/ExpressionData.txt.gz";
        String matrixOut = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/2013-12-13-CellTypeProxy/ExpressionData.txt.gz";

        DoubleMatrixDatasetFilter filter = new DoubleMatrixDatasetFilter();
        try {
            filter.filter(filterRows, filterCols, matrixIn, matrixOut);
        } catch (IOException e) {

        }
    }

    private void filter(String filterRows, String filterCols, String matrixIn, String matrixOut) throws IOException {
        Set<String> rows = null;
        Set<String> cols = null;

        if (filterCols != null) {
            TextFile tf = new TextFile(filterCols, TextFile.R);
            cols = tf.readAsSet(0, TextFile.tab);
            tf.close();
        }

        if (filterRows != null) {
            TextFile tf = new TextFile(filterRows, TextFile.R);
            rows = tf.readAsSet(0, TextFile.tab);
            tf.close();
        }

        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(matrixIn, rows, cols);
        ds.save(matrixOut);
        
    }

}
