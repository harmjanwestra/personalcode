package nl.harmjanwestra.playground.mat;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;

public class MergeMatrices {


    public static void main(String[] args) {

        String mat1 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.CMM.TMM.filtered.txt.gz";
        String mat2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\ENA.geneCounts.2019-03-02.TMM.txt.gz";
        String colLimit = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\includeSample.CMC.TargetALS.GTEx.Braineac.ROSMAP.MSBB.MayoTCX.MayoCER.EUR.ENA.OutliersRemoved.txt";
        String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.txt.gz";

        try {

            MergeMatrices m = new MergeMatrices();
            m.mergeSharedRows(mat1, mat2, colLimit, null, out);

        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void mergeSharedColumns(String mat1, String mat2, String out) throws Exception {
        DoubleMatrixDataset<String, String> ds1 = DoubleMatrixDataset.loadDoubleData(mat1);
        DoubleMatrixDataset<String, String> ds2 = DoubleMatrixDataset.loadDoubleData(mat2);


        HashMap<String, Integer> colMap = new HashMap<String, Integer>();
        ArrayList<String> newCols = new ArrayList<String>();

        int ctr = 0;
        for (String col : ds1.getColObjects()) {
            if (ds2.containsCol(col)) {

                if (colMap.containsKey(col)) {
                    colMap.put(col, ctr);
                    newCols.add(col);
                    ctr++;
                }
            }
        }

        HashMap<String, Integer> rowMap = new HashMap<>();
        ArrayList<String> newRows = new ArrayList<String>();
        ctr = 0;
        for (String row : ds1.getRowObjects()) {
            if (!ds2.containsRow(row)) {
                newRows.add(row);
                rowMap.put(row, ctr);
                ctr++;
            }
        }

        System.out.println("Output will be " + rowMap.size() + " x " + colMap.size());


        DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<>();
        DoubleMatrix2D dmd = new DenseDoubleMatrix2D(rowMap.size(), colMap.size());


        merge(ds1, dmd, colMap, rowMap);
        merge(ds2, dmd, colMap, rowMap);
        dsout.setColObjects(newCols);
        dsout.setRowObjects(newRows);
        dsout.setMatrix(dmd);
        dsout.save(out);


    }

    public void mergeSharedRows(String mat1, String mat2, String limitCols, String limitRows, String out) throws Exception {
        DoubleMatrixDataset<String, String> ds1 = DoubleMatrixDataset.loadDoubleData(mat1);
        System.out.println(ds1.rows() + " x " + ds1.columns() + " - " + mat1);
        DoubleMatrixDataset<String, String> ds2 = DoubleMatrixDataset.loadDoubleData(mat2);
        System.out.println(ds2.rows() + " x " + ds2.columns() + " - " + mat2);

        HashSet<String> allowedCols = getSet(limitCols);
        HashSet<String> allowedRows = getSet(limitRows);

        HashMap<String, Integer> colMap = new HashMap<String, Integer>();
        ArrayList<String> newCols = new ArrayList<String>();

        int ctr = 0;
        for (String col : ds1.getColObjects()) {
            if (!colMap.containsKey(col)) {
                if (allowedCols == null || allowedCols.contains(col)) {
                    colMap.put(col, ctr);
                    newCols.add(col);
                    ctr++;
                }
            }
        }
        for (String col : ds2.getColObjects()) {
            if (!colMap.containsKey(col)) {
                if (allowedCols == null || allowedCols.contains(col)) {
                    colMap.put(col, ctr);
                    newCols.add(col);
                    ctr++;
                }
            }
        }


        HashMap<String, Integer> rowMap = new HashMap<>();
        ArrayList<String> newRows = new ArrayList<String>();
        ctr = 0;
        for (String row : ds1.getRowObjects()) {
            if (allowedRows == null || allowedRows.contains(row)) {
                if (ds2.containsRow(row)) {
                    newRows.add(row);
                    rowMap.put(row, ctr);
                    ctr++;
                }
            }
        }


        System.out.println("Output will be " + rowMap.size() + " x " + colMap.size());


        DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<>();
        DoubleMatrix2D dmd = new DenseDoubleMatrix2D(rowMap.size(), colMap.size());

        System.out.println("Merge ds 1");
        merge(ds1, dmd, colMap, rowMap);
        System.out.println("Merge ds 2");
        merge(ds2, dmd, colMap, rowMap);

        dsout.setColObjects(newCols);
        dsout.setRowObjects(newRows);
        dsout.setMatrix(dmd);
        dsout.save(out);
    }

    private HashSet<String> getSet(String limitCols) throws IOException {

        if (limitCols == null) {
            return null;
        } else {
            HashSet<String> set = new HashSet<>();
            TextFile tf = new TextFile(limitCols, TextFile.R);
            set.addAll(tf.readAsArrayList());
            tf.close();
            return set;

        }
    }

    private void merge(DoubleMatrixDataset<String, String> ds1,
                       DoubleMatrix2D dmd,
                       HashMap<String, Integer> colMap,
                       HashMap<String, Integer> rowMap) {
        int[] index = new int[ds1.columns()];
        for (int col = 0; col < ds1.columns(); col++) {
            String namecol = ds1.getColObjects().get(col);
            Integer id = colMap.get(namecol);
            if (id == null) {
                index[col] = -1;
            } else {
                index[col] = id;
            }
        }

        ProgressBar pb = new ProgressBar(ds1.rows());
        for (int row = 0; row < ds1.rows(); row++) {
            String name = ds1.getRowObjects().get(row);
            Integer newrowid = rowMap.get(name);
            if (newrowid != null) {
                for (int col = 0; col < ds1.columns(); col++) {
                    int newcolid = index[col];
                    if (newcolid != -1) {
                        dmd.setQuick(newrowid, newcolid, ds1.getElementQuick(row, col));
                    }
                }
            }
            pb.iterate();
        }
        pb.close();
    }
}
