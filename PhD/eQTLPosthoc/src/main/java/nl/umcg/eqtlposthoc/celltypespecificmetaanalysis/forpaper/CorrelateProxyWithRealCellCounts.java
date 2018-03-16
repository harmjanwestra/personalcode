/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import edu.ufl.cise.colamd.tdouble.Dcolamd;
import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class CorrelateProxyWithRealCellCounts {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String proxyfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/2013-12-13-CellTypeProxy/CellTypeProxyFile.txt";
        String cellcountMatrix = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/EGCUTEndophenotypesValidSamples-NoSex.txt-asLoadedByNormalizer.txt";
        String query = "Neutrophils";
        try {
            TextFile tf = new TextFile(proxyfile, TextFile.R);

            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            HashMap<String, Double> proxyPerSample = new HashMap<String, Double>();
            while (elems != null) {
                String sample = elems[0];
                Double count = Double.parseDouble(elems[1]);
                proxyPerSample.put(sample, count);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(cellcountMatrix);

            Integer queryId = ds.hashRows.get(query);
            System.out.println("Sample\tProxy\tCellcount");
            if (queryId != null) {
                for (int col = 0; col < ds.nrCols; col++) {
                    String colName = ds.colObjects.get(col);
                    Double d = proxyPerSample.get(colName);
                    if (d != null) {
                        System.out.println(colName + "\t" + d + "\t" + ds.rawData[queryId][col]);
                    }
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
