/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ReplaceInteractionTermsWithProbeCellCountCorrelations {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String metaFile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/MetaAnalysis/CellCountCorrelationMatrixCorrelatedWithInteractionTerms.txt";
        String corrfile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/2013-07-15-EGCUTCellCountCorrelations/CellCountCorrelationMatrix.txtEndophenotypeVsExpressionCorrelationMatrix.txt";
        String output = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-24-12KEQTLs-WSHIP-TREND/MetaAnalysis/MetaAnalysisZScoreMatrix-ValuesReplacedBy.txt";
        
        ReplaceInteractionTermsWithProbeCellCountCorrelations c = new ReplaceInteractionTermsWithProbeCellCountCorrelations();
        try {
            c.run(metaFile, output, corrfile);
        } catch (IOException ex) {
            Logger.getLogger(ReplaceInteractionTermsWithProbeCellCountCorrelations.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void run(String metaFile, String output, String cellCountCorrelations) throws IOException {
        DoubleMatrixDataset<String, String> dsmeta = new DoubleMatrixDataset<String, String>(metaFile); // eQTLs on rows
        DoubleMatrixDataset<String, String> dscellcountcorrelation = new DoubleMatrixDataset<String, String>(cellCountCorrelations); // probes on rows
        
        DoubleMatrixDataset<String, String> outputds = new DoubleMatrixDataset<String, String>(dsmeta.nrRows, dscellcountcorrelation.nrCols);
        outputds.rowObjects = dsmeta.rowObjects;
        outputds.colObjects = dscellcountcorrelation.colObjects;
        for(int r=0;r<dsmeta.nrRows;r++){
            String eQTL = dsmeta.rowObjects.get(r);
            String probe = eQTL.split("-")[1];
            Integer rowId = dscellcountcorrelation.hashRows.get(probe);
            System.arraycopy(dscellcountcorrelation.rawData[rowId], 0, outputds.rawData[r], 0, dscellcountcorrelation.nrCols);
        }
        outputds.recalculateHashMaps();
        outputds.save(output);
        
        
    }
}
