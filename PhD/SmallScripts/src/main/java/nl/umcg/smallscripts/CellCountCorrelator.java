/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;

/**
 *
 * @author harmjan
 */
public class CellCountCorrelator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String probefile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/HT12v3SNPProbeCombos.txt";
            String dsloc = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTest/CellTypeSpecificProbeExpression.txt";
            String cellloc = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTest/CellTypeProxyFile.txt";
            String corrout = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTest/CorrelationOfAllCisWithProxyPhenotype3.txt";
            
            // load probes
            TextFile tf = new TextFile(probefile, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            HashSet<String> probes = new HashSet<String>();
            while (elems != null) {
                if (elems.length > 1) {
                    probes.add(elems[1]);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            // load cell counts
            HashMap<String, Double> sampleToCellCt = new HashMap<String, Double>();
            TextFile tf2 = new TextFile(cellloc, TextFile.R);
            String[] elems2 = tf2.readLineElems(TextFile.tab);
            while (elems2 != null) {
                if (elems2.length > 1) {
                    String sample = elems2[0];
                    try {
                        Double ct = Double.parseDouble(elems2[1]);
                        sampleToCellCt.put(sample, ct);
                    } catch (NumberFormatException e) {
                    }
                }
                elems2 = tf2.readLineElems(TextFile.tab);

            }
            tf2.close();


            // load exp data
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(dsloc);

            TextFile corrrrrr = new TextFile(corrout, TextFile.W);
            for (int r = 0; r < ds.nrRows; r++) {
                String probe = ds.rowObjects.get(r);
                if (probes.contains(probe)) {
                    ArrayList<Double> valsX = new ArrayList<Double>();
                    ArrayList<Double> valsY = new ArrayList<Double>();
                    for (int c = 0; c < ds.nrCols; c++) {
                        String sample = ds.colObjects.get(c);
                        Double cc = sampleToCellCt.get(sample);
                        if (cc != null) {
                            valsX.add(ds.rawData[r][c]);
                            valsY.add(cc);
                        }
                    }
                    double[] xArr = toPrimitiveArr(valsX.toArray(new Double[0]));
                    double[] yArr = toPrimitiveArr(valsY.toArray(new Double[0]));
                    double correlation = JSci.maths.ArrayMath.correlation(xArr, yArr);
                    corrrrrr.writeln(probe + "\t" + correlation);

                }
            }

            corrrrrr.close();

            //LD:D
            // tnx :D
        } catch (IOException e) {
        }
    }

    private static double[] toPrimitiveArr(Double[] toArray) {
        double[] arr = new double[toArray.length];
        for (int i = 0; i < toArray.length; i++) {
            arr[i] = toArray[i];
        }
        return arr;
    }
}
