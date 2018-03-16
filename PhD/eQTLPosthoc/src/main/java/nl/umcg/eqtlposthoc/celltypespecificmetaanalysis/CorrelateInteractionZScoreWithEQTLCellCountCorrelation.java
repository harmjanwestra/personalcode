/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class CorrelateInteractionZScoreWithEQTLCellCountCorrelation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here
            String metamatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-24-12KEQTLs-WSHIP-TREND/MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
            String cellcountcorrelationmatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-24-12KEQTLs-WSHIP-TREND/MetaAnalysis/MetaAnalysisZScoreMatrix-CorrelatedWithCellCountCorrelations.txt";
            String output = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-25-PostHocAnalyses/MetaWSHIP-TREND/CorrelationsBetweenInteractionZScoresAndInteractionTermCellCountCorrrelations/";
            CorrelateInteractionZScoreWithEQTLCellCountCorrelation q = new CorrelateInteractionZScoreWithEQTLCellCountCorrelation();

            q.run(metamatrix, cellcountcorrelationmatrix, output);
        } catch (IOException ex) {
            Logger.getLogger(CorrelateInteractionZScoreWithEQTLCellCountCorrelation.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void run(String metaMatrixFileLoc, String cellcountcorrelationMatrixFileLoc, String output) throws IOException {
        DoubleMatrixDataset<String, String> metamatrix = new DoubleMatrixDataset<String, String>(metaMatrixFileLoc);
        DoubleMatrixDataset<String, String> ccmatrix = new DoubleMatrixDataset<String, String>(cellcountcorrelationMatrixFileLoc); // eQTLs on rows, cellcounts on columns

        SpearmansCorrelation spearman = new SpearmansCorrelation();

        Integer rowIdInteractionTerm = metamatrix.hashRows.get("CellTypeInteractionZScore");
        Integer mainEffectTerm = metamatrix.hashRows.get("MainEffectZScore");
        // MainEffectZScore
        for (int col = 0; col < ccmatrix.nrCols; col++) {
            String cellcount = ccmatrix.colObjects.get(col);
            ArrayList<Double> xArrl = new ArrayList<Double>();
            ArrayList<Double> yArrl = new ArrayList<Double>();
            for (int row = 0; row < ccmatrix.nrRows; row++) {
                String name = ccmatrix.rowObjects.get(row);
                Integer colIdInMetaMatrix = metamatrix.hashCols.get(name);
                if (colIdInMetaMatrix != null) {
                    xArrl.add(ccmatrix.rawData[row][col]);
                    double main = metamatrix.rawData[mainEffectTerm][colIdInMetaMatrix];
                    if (main < 0) {
                        yArrl.add(-metamatrix.rawData[rowIdInteractionTerm][colIdInMetaMatrix]);
                    } else {
                        yArrl.add(metamatrix.rawData[rowIdInteractionTerm][colIdInMetaMatrix]);
                    }

                }
            }
            double[] xarr = Primitives.toPrimitiveArr(xArrl.toArray(new Double[0]));
            double[] yarr = Primitives.toPrimitiveArr(yArrl.toArray(new Double[0]));
            double correlation = spearman.correlation(xarr, yarr);
            new ScatterPlot(500, 500, xarr, yarr, ScatterPlot.OUTPUTFORMAT.PNG, output + cellcount + "-" + correlation + ".png");

        }



    }
}
