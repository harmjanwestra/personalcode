/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class CheckForNaNs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String[] datasetFileDirs = new String[]{
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutputMT2/",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/CellTypeSpecificTesteQTLOutputMT/",
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/MarjoleinVisit3/CellTypeSpecificityTestOutputMT2/"
            };

            for (int d = 0; d < datasetFileDirs.length; d++) {
                System.out.println("");
                System.out.println("Loading dataset: " + d);
                String matrixIn = datasetFileDirs[d] + "CellTypeSpecificityMatrix.binary";

                DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(matrixIn);

                check(ds);
            }

        } catch (IOException e) {
            e.printStackTrace();

        }
    }

    public static void check(DoubleMatrixDataset<String, String> ds) {
        for (int r = 0; r < ds.nrRows; r++) {
            for (int c = 0; c < ds.nrCols; c++) {
                if (Double.isNaN(ds.rawData[r][c])) {
                    System.out.println(ds.rowObjects.get(r) + "\t" + ds.colObjects.get(c) + " is NaN!");
                }
                if (ds.rawData[r][c] == 0d) {
                    System.out.println(ds.rowObjects.get(r) + "\t" + ds.colObjects.get(c) + " is 0.0!");
                }
            }
        }
    }
}
