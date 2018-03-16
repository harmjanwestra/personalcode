/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class MetaAnalysisMatrixFlipEffectDirections {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String dsLoc = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-24-12KEQTLs-WSHIP-TREND/MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
        String outputfile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-24-12KEQTLs-WSHIP-TREND/MetaAnalysis/MetaAnalysisZScoreMatrix-EffectsFlipped-WithProxyInteractionTerms.txt";
        try {
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(dsLoc);
            Integer rowIdMainFx = ds.hashRows.get("MainEffectZScore");


            DoubleMatrixDataset<String, String> dsOutput = new DoubleMatrixDataset<String, String>(ds.nrRows, ds.nrCols);


            int rowctr = 0;
            double[] mainFx = ds.rawData[rowIdMainFx];
            ArrayList<String> newRowObjects = new ArrayList<String>();
            for (int r = 0; r < ds.nrRows; r++) {
//                if (ds.rowObjects.get(r).equals("CellTypeInteractionZScore") || ds.rowObjects.get(r).equals("MainEffectZScore") || ds.rowObjects.get(r).equals("CellTypeZScore") || ds.rowObjects.get(r).equals("CellTypeSNPZScore")) {
//                    // do nothing
//                } else {
                    for (int c = 0; c < ds.nrCols; c++) {
                        if (mainFx[c] < 0 && !ds.rowObjects.get(r).equals("MainEffectZScore")) {
                            dsOutput.rawData[rowctr][c] = -ds.rawData[r][c];
                        } else {
                            dsOutput.rawData[rowctr][c] = ds.rawData[r][c];
                        }

                    }
                    newRowObjects.add(ds.rowObjects.get(r));
                    rowctr++;
//                }
            }

            dsOutput.rowObjects = newRowObjects;
            dsOutput.colObjects = ds.colObjects;
            dsOutput.recalculateHashMaps();
            dsOutput.save(outputfile);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
