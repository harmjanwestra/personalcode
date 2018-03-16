/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class DoubleMatrixDatasetReplaceProbeNames {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String fileIn = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/ExpressionData-HT12v4/ExpressionData.txt.QuantileNormalized.Log2Transformed.txt.gz";
        String fileOut = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/ExpressionData-HT12v4/ExpressionData.txt.QuantileNormalized.Log2Transformed-HT12v4.txt";
        String pbt = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-18-ProbeAnnotationFile.txt";
        String source = "HT12v3.txt";
        String dest = "HT12v4.txt";

        try {
            DoubleMatrixDatasetReplaceProbeNames f = new DoubleMatrixDatasetReplaceProbeNames();
            f.run(fileIn, fileOut, pbt, source, dest);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String fileIn, String fileOut, String pbt, String source, String dest) throws IOException {
        ProbeTranslation pbtl = new ProbeTranslation();
        HashMap<String, String> probeToProbe = pbtl.getProbeTranslation(pbt, source, dest);

        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(fileIn);
        boolean[] includerow = new boolean[ds.nrRows];
        ArrayList<String> probes = new ArrayList<String>();
        int ctr = 0;
        for (int i = 0; i < ds.nrRows; i++) {
            String probe = ds.rowObjects.get(i);
            String otherprobe = probeToProbe.get(probe);

            if (probe.startsWith("HT12v4")) {
                includerow[i] = true;
                probes.add(probe.replace("HT12v4-", ""));
                ctr++;
            } else if (otherprobe != null && !otherprobe.equals("-")) {
                includerow[i] = true;
                probes.add(otherprobe);
                ctr++;
            }
        }

        if (ctr == 0) {
            System.err.println("ERROR!");
            System.exit(0);
        }

        double[][] matrix = new double[ctr][ds.nrCols];
        int ctr2 = 0;
        for (int i = 0; i < ds.nrRows; i++) {
            if (includerow[i]) {
                matrix[ctr2] = ds.rawData[i];
                ctr2++;
            }
        }

        ds.rawData = matrix;
        ds.rowObjects = probes;
        ds.recalculateHashMaps();
        ds.save(fileOut);
    }
}
