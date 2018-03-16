/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class DannyTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String dannyvector = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/file.txt";
            String dannyOut = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/DannyFehrmann/";
            Gpio.createDir(dannyOut);
            String hjvector = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix-Meta.txt";
            String probetrans = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt";

            HashMap<String, String> wg6toht12 = new HashMap<String, String>();
            TextFile tf = new TextFile(probetrans, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                wg6toht12.put(elems[14], elems[5]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(dannyvector);

            HashMap<String, Double> probeToHJVector = new HashMap<String, Double>();

            TextFile tf2 = new TextFile(hjvector, TextFile.R);
            tf2.readLine();
            String[] elems2 = tf2.readLineElems(TextFile.tab);
            while (elems2 != null) {
                probeToHJVector.put(elems2[0], Double.parseDouble(elems2[1]));
                elems2 = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            for (int c = 0; c < ds.nrCols; c++) {
                ArrayList<Double> valsWG6 = new ArrayList<Double>();
                ArrayList<Double> valsHT12 = new ArrayList<Double>();
                int ctr = 0;
                TextFile tfout = new TextFile(dannyOut+ds.colObjects.get(c)+".txt", TextFile.W);
                tfout.writeln("WG6Probe\tHT12v3Probe\tWG6Val\tHT12v3Val");
                for (int r = 0; r < ds.nrRows; r++) {
                    String probe = ds.rowObjects.get(r);
                    String ht12 = wg6toht12.get(probe);
                    if (ht12 != null) {
                        Double val = probeToHJVector.get(ht12);
                        if (val != null) {
                            valsWG6.add(ds.rawData[r][c]);
                            valsHT12.add(val);
                        }
                        tfout.writeln(probe+"\t"+ht12+"\t"+ds.rawData[r][c]+"\t"+val);
                        ctr++;
                    }
                }
                tfout.close();

                System.out.println(ctr);

                double[] xArr = toPrimitiveArr(valsWG6.toArray(new Double[0]));
                double[] yArr = toPrimitiveArr(valsHT12.toArray(new Double[0]));
                
                


                double r = JSci.maths.ArrayMath.correlation(xArr, yArr);
                System.out.println(ds.colObjects.get(c) + "\t" + r);

            }



        } catch (IOException e) {
            e.printStackTrace();
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
