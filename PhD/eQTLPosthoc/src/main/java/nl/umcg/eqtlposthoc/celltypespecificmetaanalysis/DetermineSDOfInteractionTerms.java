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
public class DetermineSDOfInteractionTerms {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here
            String matrixIn = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-31-MetaAnalysis5Datasets/MetaAnalysis/CellTypeSpecificityMatrix.binary";
            DetermineSDOfInteractionTerms d = new DetermineSDOfInteractionTerms();
            d.run(matrixIn);
        } catch (IOException ex) {
            Logger.getLogger(DetermineSDOfInteractionTerms.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.exit(0);
    }

    public void run(String matrix) throws IOException {
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(matrix);

        for (int row = 0; row < ds.nrRows; row++) {
            double sd = JSci.maths.ArrayMath.standardDeviation(ds.rawData[row]);
            double mean = JSci.maths.ArrayMath.mean(ds.rawData[row]);
            cern.colt.list.DoubleArrayList list = new cern.colt.list.DoubleArrayList(ds.rawData[row]);
            double skewness = cern.jet.stat.Descriptive.skew(list,mean,sd);
            double kurtosis = cern.jet.stat.Descriptive.kurtosis(list,mean,sd);
            
            System.out.println(ds.rowObjects.get(row) + "\t" + mean + "\t"+  sd + "\t" + skewness + "\t"  + kurtosis);
        }
    }
}
