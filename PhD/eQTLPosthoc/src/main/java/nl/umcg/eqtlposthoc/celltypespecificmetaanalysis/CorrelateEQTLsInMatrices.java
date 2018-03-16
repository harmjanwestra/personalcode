/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class CorrelateEQTLsInMatrices {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String f1 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/Groningen/CellTypeSpecificityMatrix.binary";
            String f2 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/EGCUT/CellTypeSpecificityMatrix.binary";
            String fo = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/MetaAnalysisWithoutEGCUT/GroningenVsEGCUTComparisonOnCovariates.txt";
            CorrelateEQTLsInMatrices c = new CorrelateEQTLsInMatrices();
            c.correlate(f1, false, f2, false, fo);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void correlate(String ds1, boolean transpose1, String ds2, boolean transpose2, String outputfile) throws IOException {
        DoubleMatrixDataset<String, String> matrix1 = new DoubleMatrixDataset<String, String>(ds1);
        DoubleMatrixDataset<String, String> matrix2 = new DoubleMatrixDataset<String, String>(ds2);
        
        HashSet<String> idsToAvoid = new HashSet<String>();
        idsToAvoid.add("CellTypeZScore");
        idsToAvoid.add("MainEffectZScore");
        idsToAvoid.add("CellTypeInteractionZScore");
        idsToAvoid.add("CellTypeSNPZScore");

        if (transpose1) {
            matrix1.transposeDataset();
        }

        if (transpose2) {
            matrix2.transposeDataset();
        }

        // get the row intersect..
        ArrayList<String> rowNames = new ArrayList<String>();
        for (int row = 0; row < matrix1.nrRows; row++) {
            String rowName = matrix1.rowObjects.get(row);
            if (!idsToAvoid.contains(rowName)) {
                Integer id = matrix2.hashRows.get(rowName);
                if (id != null) {
                    rowNames.add(rowName);
                }
            }

        }

        TextFile tfOut = new TextFile(outputfile, TextFile.W);
        for (int r = 0; r < rowNames.size(); r++) {
            String rowName = rowNames.get(r);
            Integer row = matrix1.hashRows.get(rowName);
            Integer id = matrix2.hashRows.get(rowName);
            if (id != null) {
                // now get the overlapping columns
                ArrayList<Double> valsX = new ArrayList<Double>();
                ArrayList<Double> valsY = new ArrayList<Double>();

                for (int col = 0; col < matrix1.nrCols; col++) {
                    String colName = matrix1.colObjects.get(col);
                    if (!idsToAvoid.contains(colName)) {
                        Integer id2 = matrix2.hashCols.get(colName);
                        if (id2 != null) {
                            valsX.add(matrix1.rawData[row][col]);
                            valsY.add(matrix2.rawData[id][id2]);
                        }
                    }

                }

                double[] xarr = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                double[] yarr = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));

                double corr = JSci.maths.ArrayMath.correlation(xarr, yarr);
                tfOut.writeln(rowName+"\t"+corr+"\t"+xarr.length);
            }
        }

        tfOut.close();

//        DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<String, String>();
//        output.rawData = correlationmatrix;
//        output.colObjects = rowNames;
//        output.rowObjects = rowNames;
//        output.recalculateHashMaps();
//        output.save(outputfile);

    }
}
