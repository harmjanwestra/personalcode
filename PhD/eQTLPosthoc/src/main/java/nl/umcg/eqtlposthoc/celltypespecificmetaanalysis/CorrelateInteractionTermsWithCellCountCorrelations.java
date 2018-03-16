/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 * nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.CorrelateInteractionTermsWithCellCountCorrelations
 */
public class CorrelateInteractionTermsWithCellCountCorrelations {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        String metaMatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-31-MetaAnalysis5Datasets/MetaAnalysis/CellTypeSpecificityMatrix.binary";
        String cellCountMatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/PCAOnEGCUTEndophenotypeCorrelations/principalcomponents.txt";
        String output = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-31-MetaAnalysis5Datasets/CorrelationsWithCellCounts/MetaAnalysisZScoreMatrix-CorrelatedWithCellCountCorrelations.txt";

        CorrelateInteractionTermsWithCellCountCorrelations c = new CorrelateInteractionTermsWithCellCountCorrelations();
        try {
            c.run(metaMatrix, cellCountMatrix, output);
        } catch (IOException ex) {
            Logger.getLogger(CorrelateInteractionTermsWithCellCountCorrelations.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void run(String metaMatrixFile, String cellCountCorrelationMatrix, String outputfileName) throws IOException {
        DoubleMatrixDataset<String, String> metamatrix = new DoubleMatrixDataset<String, String>(metaMatrixFile); // eQTLs on columns, covariates on rows
        DoubleMatrixDataset<String, String> cellcountc = new DoubleMatrixDataset<String, String>(cellCountCorrelationMatrix); // endophenotypes on columns, covariates on rows

        Integer[] rtr = new Integer[cellcountc.nrRows];
        for (int r = 0; r < cellcountc.nrRows; r++) {

            String cov = cellcountc.rowObjects.get(r);
            Integer idInMetaMatrix = metamatrix.hashRows.get(cov);

            rtr[r] = idInMetaMatrix;


        }

        DoubleMatrixDataset<String, String> output = new DoubleMatrixDataset<String, String>();
        double[][] matrixOut = new double[metamatrix.nrCols][cellcountc.nrCols];

        Integer mainFx = metamatrix.hashRows.get("MainEffectZScore");
        double[] mainFxRow = metamatrix.rawData[mainFx];
        // correct for effect flips, per eQTL

        for (int col = 0; col < cellcountc.nrCols; col++) {
            for (int col2 = 0; col2 < metamatrix.nrCols; col2++) {
                ArrayList<Double> valX = new ArrayList<Double>();
                ArrayList<Double> valY = new ArrayList<Double>();
                for (int row = 0; row < cellcountc.nrRows; row++) {
                    if (rtr[row] != null) {
                        valX.add(cellcountc.rawData[row][col]);
                        if (mainFxRow[col2] < 0) {
                            valY.add(-1 * metamatrix.rawData[rtr[row]][col2]);
                        } else {
                            valY.add(metamatrix.rawData[rtr[row]][col2]);
                        }

                    }
                }
                double[] xArr = Primitives.toPrimitiveArr(valX.toArray(new Double[0]));
                double[] yArr = Primitives.toPrimitiveArr(valY.toArray(new Double[0]));
                double r = JSci.maths.ArrayMath.correlation(xArr, yArr);
                matrixOut[col2][col] = r;
            }
        }

        output.rawData = matrixOut;
        output.rowObjects = metamatrix.colObjects;
        output.colObjects = cellcountc.colObjects;

        output.recalculateHashMaps();
        output.save(outputfileName);
    }
}
