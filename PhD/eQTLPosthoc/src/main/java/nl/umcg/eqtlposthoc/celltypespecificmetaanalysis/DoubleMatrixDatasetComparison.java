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
public class DoubleMatrixDatasetComparison {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here


        try {
            DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutputMT2/CellTypeSpecificityMatrix.binary");

            DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix.binary");


            boolean datamissing = false;
            for (int i = 0; i < ds1.nrRows; i++) {

                String rowName = ds1.rowObjects.get(i);
                Integer rowIdInOtherdataset = ds2.hashRows.get(rowName);
                if (rowIdInOtherdataset == null) {
                    System.out.println("Row: " + rowName + " not present in ds2");
                    datamissing = true;
                } else {
                    System.out.println(rowName + "\t" + rowIdInOtherdataset);
                }
            }

            for (int i = 0; i < ds2.nrRows; i++) {
                String rowName = ds2.rowObjects.get(i);
                Integer rowIdInOtherdataset = ds1.hashRows.get(rowName);
                if (rowIdInOtherdataset == null) {
                    System.out.println("Row: " + rowName + " not present in ds1");
                    datamissing = true;
                }
            }


            for (int i = 0; i < ds1.nrCols; i++) {
                String colName = ds1.colObjects.get(i);
                Integer colIdInOtherdataset = ds2.hashCols.get(colName);
                if (colIdInOtherdataset == null) {
                    System.out.println("Col: " + colName + " not present in ds2");
                    datamissing = true;
                }
            }

            for (int i = 0; i < ds2.nrCols; i++) {
                String colName = ds2.colObjects.get(i);
                Integer colIdInOtherdataset = ds1.hashCols.get(colName);
                if (colIdInOtherdataset == null) {
                    System.out.println("Col: " + colName + " not present in ds1");
                    datamissing = true;
                }
            }


            // now compare values
            if (!datamissing) {
                System.out.println("Data is not missing..");
                System.out.println(ds1.nrRows);
                System.out.println(ds1.nrCols);
                System.out.println(ds2.nrRows);
                System.out.println(ds2.nrCols);
                for (int i = 0; i < ds1.nrRows; i++) {
                    String rowName = ds1.rowObjects.get(i);
                    Integer rowIdInOtherdataset = ds2.hashRows.get(rowName);
                    if (rowIdInOtherdataset != null) {

                        for (int j = 0; j < ds1.nrCols; j++) {
                            String colName = ds1.colObjects.get(j);
                            Integer colIdInOtherDataset = ds2.hashCols.get(colName);

                            if(ds1.rawData[i][j] == 0d){
                                System.out.println("Value is 0.0 "+ds1.rowObjects.get(i)+"\t"+ds1.colObjects.get(j));
                            }
                            
                            double diff = Math.abs(ds1.rawData[i][j] - ds2.rawData[rowIdInOtherdataset][colIdInOtherDataset]);
                            if (diff > 0) {
                                System.out.println("Differences in values for\t" + rowName + "\t" + colName + "\t" + ds1.rawData[i][j] + "\t" + ds2.rawData[rowIdInOtherdataset][colIdInOtherDataset]);
                            }


                        }
                    } else {
                        System.out.println(rowName + " has null identifier in ds2");
                    }
                }
            } else {
                System.out.println("Data is missing...");
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
