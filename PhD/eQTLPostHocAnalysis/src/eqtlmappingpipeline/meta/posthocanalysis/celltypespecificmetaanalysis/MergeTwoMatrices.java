/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class MergeTwoMatrices {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {
            String file1 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/2013-07-15-EGCUTCellCountCorrelations/CellCountCorrelationMatrixCorrelatedWithInteractionTerms.txt";
            String file2 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/MetaAnalysisWithoutEGCUT/CellCountCorrelationMatrixCorrelatedWithInteractionTerms2.txt";
            String outfi = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/MetaAnalysisWithoutEGCUT/CellCountCorrelationMatrixCorrelatedWithInteractionTermsMergedWithEGCUT.txt";

            DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>(file1);
            DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(file2);

            MergeTwoMatrices m = new MergeTwoMatrices();

            DoubleMatrixDataset<String, String> ds3 = m.mergeUsingRowIds(ds1, "EGCUT_", ds2, "MetaWOEGCUT");
            ds3.save(outfi);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public DoubleMatrixDataset<String, String> mergeUsingRowIds(DoubleMatrixDataset<String, String> ds1, String prefix1, DoubleMatrixDataset<String, String> ds2, String prefix2) {
        // get overlapping rows
        HashSet<String> selectedRows = new HashSet<String>();
        for (int row = 0; row < ds1.nrRows; row++) {
            String name = ds1.rowObjects.get(row);
            Integer idInOther = ds2.hashRows.get(name);

            if (idInOther != null) {
                if (selectedRows.contains(name)) {
                    System.err.println("WARNING: " + name + " is a duplicate rowname in either dataset?");
                }
                selectedRows.add(name);
            }
        }
        
        if(selectedRows.isEmpty()){
            System.err.println("ERROR: no rows to merge!");
            return null;
        }

        DoubleMatrixDataset<String, String> ds3 = new DoubleMatrixDataset<String, String>();

        ArrayList<String> newRowNames = new ArrayList<String>();
        newRowNames.addAll(Arrays.asList(selectedRows.toArray(new String[0])));


        double[][] outputData = new double[selectedRows.size()][ds1.nrCols + ds2.nrCols];

        ds3.colObjects = new ArrayList<String>();
        for (int col = 0; col < ds1.colObjects.size(); col++) {
            ds3.colObjects.add(prefix1 + ds1.colObjects.get(col));
        }

        for (int col = 0; col < ds2.colObjects.size(); col++) {
            ds3.colObjects.add(prefix2 + ds2.colObjects.get(col));
        }

        ds3.rowObjects = newRowNames;

        for (int i = 0; i < selectedRows.size(); i++) {
            String name = newRowNames.get(i);
            Integer id1 = ds1.hashRows.get(name);
            Integer id2 = ds2.hashRows.get(name);
            System.arraycopy(ds1.rawData[id1], 0, outputData[i], 0, ds1.nrCols);
            System.arraycopy(ds2.rawData[id2], 0, outputData[i], ds1.nrCols, ds2.nrCols);
        }
        ds3.rawData = outputData;
        ds3.recalculateHashMaps();
        return ds3;
    }
}
