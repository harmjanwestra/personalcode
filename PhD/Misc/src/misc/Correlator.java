/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class Correlator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String d1 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix.txt";
            String d2 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix-Meta.txt";

            DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>(d1);
            ds1.transposeDataset();
            DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(d2);
            ds2.transposeDataset();



            for (int r = 0; r < ds1.nrRows; r++) {
                String snpprobe = ds1.rowObjects.get(r);
                Integer snpprobeds2 = ds2.hashRows.get(snpprobe);
                if (snpprobeds2 != null) {
                    ArrayList<Double> valx = new ArrayList<Double>();
                    ArrayList<Double> valy = new ArrayList<Double>();
                    for (int c = 0; c < ds1.nrCols-2; c++) {
                        String probe = ds1.colObjects.get(c);
                        Integer probe2 = ds2.hashCols.get(probe);

                        if (probe2 != null) {
                            valx.add(ds1.rawData[r][c]);
                            valy.add(ds2.rawData[snpprobeds2][probe2]);
                        }
                    }
                    double[] xArr = toPrimitiveArr(valx.toArray(new Double[0]));
                    double[] yArr = toPrimitiveArr(valy.toArray(new Double[0]));

                    double corr = JSci.maths.ArrayMath.correlation(xArr, yArr);
                    System.out.println(ds1.rowObjects.get(r) + "\t" + corr);

                }

            }




//            String d1 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix.txt";
//            String d2 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/RotterdamStudy/MarjoleinVisit3/CellTypeSpecificityTestOutput/CellTypeSpecificityMatrix-HT12v3.txt";
            boolean switchEffect2 = true;
//            String d3 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/CellTypeSpecificTesteQTLOutput2/CellTypeSpecificityMatrix.txt";
//            String dOut = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix-Meta.txt";
//
//
//            int n1 = 891;
//            int n2 = 755;
//            int n3 = 1220;
//
//
//            DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>(d1);
//            ds1.transposeDataset();
//            DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(d2);
//            ds2.transposeDataset();
//
//            DoubleMatrixDataset<String, String> ds3 = new DoubleMatrixDataset<String, String>(d3);
//            ds3.transposeDataset();
//
//            TextFile tfOut = new TextFile(dOut, TextFile.W);
//
//            for (int r = 0; r < ds1.nrRows; r++) {
//                double[] z = new double[3];
//                int[] ns = new int[3];
//
//
//
//                String snpprobe = ds1.rowObjects.get(r);
//                Integer snpprobeds2 = ds2.hashRows.get(snpprobe);
//                Integer snpprobeds3 = ds3.hashRows.get(snpprobe);
//
//                tfOut.writeln(snpprobe);
//                for (int c = 0; c < ds1.nrCols - 2; c++) {
//                    String probe = ds1.colObjects.get(c);
//                    if (snpprobeds2 != null) {
//                        Integer probe2 = ds2.hashCols.get(probe);
//                        if (probe2 != null) {
//                            z[1] = ds2.rawData[snpprobeds2][probe2];
//                            ns[1] = n2;
//                        } else {
//                            z[1] = Double.NaN;
//                            ns[1] = 0;
//                        }
//                    }
//
//                    if (snpprobeds3 != null) {
//                        Integer probe2 = ds3.hashCols.get(probe);
//                        if (probe2 != null) {
//                            z[2] = ds3.rawData[snpprobeds3][probe2];
//                            ns[2] = n3;
//                        } else {
//                            z[2] = Double.NaN;
//                            ns[2] = 0;
//                        }
//                    }
//
//                    // meta-analyze...
//                    double metaZ = ZScores.getWeightedZ(z, ns);
//                    tfOut.writeln(probe + "\t" + metaZ);
//                }
//            }
//
//            tfOut.close();






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
