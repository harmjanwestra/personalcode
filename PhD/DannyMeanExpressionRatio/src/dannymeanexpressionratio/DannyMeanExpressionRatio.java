/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package dannymeanexpressionratio;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class DannyMeanExpressionRatio {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            String file1 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/GEOAffyData/GSE6613vsGSE22886/expHGU133A_GSE22886.txt";
            String catFile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/GEOAffyData/GSE6613vsGSE22886/GSE22886.txt";
            String fileOut = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/GEOAffyData/GSE6613vsGSE22886/HJ/expHGU133A_GSE22886_mean.txt";
            boolean normalize = false;

            DoubleMatrixDataset<String, String> averageCellTypes = DannyMeanExpressionRatio.getAveragePerCat(file1, catFile, fileOut, normalize);

            file1 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/GEOAffyData/GSE6613vsGSE22886/expHGU133A_GSE6613_wholeBlood.txt";
            catFile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/GEOAffyData/GSE6613vsGSE22886/GSE6613.txt";
            fileOut = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/GEOAffyData/GSE6613vsGSE22886/HJ/expHGU133A_GSE6613_wholeBlood_mean.txt";
            DoubleMatrixDataset<String, String> averageWholeBlood = DannyMeanExpressionRatio.getAveragePerCat(file1, catFile, fileOut, normalize);


            double[][] ratios = new double[averageCellTypes.nrRows][averageWholeBlood.nrCols];

            
            
            for (int r = 0; r < averageCellTypes.nrRows; r++) {
                String probe = averageCellTypes.rowObjects.get(r);
                Integer rowIdWholeBlood = averageWholeBlood.hashRows.get(probe);
                if (rowIdWholeBlood != null) {
                    for (int c = 0; c < averageCellTypes.nrCols; c++) {
                        for (int c2 = 0; c2 < averageWholeBlood.nrCols; c2++) {
                            ratios[r][c2] = averageCellTypes.rawData[r][c] / averageWholeBlood.rawData[rowIdWholeBlood][c2];
                        }
                    }
                }
            }
            
            String fileOut2 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/GEOAffyData/GSE6613vsGSE22886/HJ/expHGU133A_GSE6613_GSE22886_patientTypeRatio.txt";
            
            DoubleMatrixDataset<String, String> out = new DoubleMatrixDataset<String, String>();
            out.rawData = ratios;
            out.rowObjects = averageCellTypes.rowObjects;
            out.colObjects = averageWholeBlood.colObjects;
            out.recalculateHashMaps();
            
            out.save(fileOut2);
            
            
            ratios = new double[averageCellTypes.nrRows][averageCellTypes.nrCols*averageWholeBlood.nrCols];

            String[] names = new String[averageCellTypes.nrCols*averageWholeBlood.nrCols];
            
            for (int r = 0; r < averageCellTypes.nrRows; r++) {
                String probe = averageCellTypes.rowObjects.get(r);
                Integer rowIdWholeBlood = averageWholeBlood.hashRows.get(probe);
                if (rowIdWholeBlood != null) {
                    for (int c = 0; c < averageCellTypes.nrCols; c++) {
                        for (int c2 = 0; c2 < averageWholeBlood.nrCols; c2++) {
                            ratios[r][c*c2] = averageCellTypes.rawData[r][c] / averageWholeBlood.rawData[rowIdWholeBlood][c2];
                            names[c*c2] = averageCellTypes.colObjects.get(c)+"-"+averageWholeBlood.colObjects.get(c2);
                        }
                    }
                }
            }
            
            fileOut2 = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/GEOAffyData/GSE6613vsGSE22886/HJ/expHGU133A_GSE6613_GSE22886_cellTypeRatio.txt";
            
            out = new DoubleMatrixDataset<String, String>();
            out.rawData = ratios;
            out.rowObjects = averageCellTypes.rowObjects;
            out.colObjects = Arrays.asList(names);
            out.recalculateHashMaps();
            
            out.save(fileOut2);
            




        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static DoubleMatrixDataset<String, String> getAveragePerCat(String expressionFile, String catFile, String fileOut, boolean normalize) throws IOException {
        HashMap<String, String> sampleToCat = new HashMap<String, String>();

        TextFile tf = new TextFile(catFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        HashSet<String> categories = new HashSet<String>();
        while (elems != null) {
            String sample = elems[0];
            String cat = elems[1];
            sampleToCat.put(sample, cat);
            categories.add(cat);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        ArrayList<String> categoryArr = new ArrayList<String>();
        HashMap<String, Integer> catToId = new HashMap<String, Integer>();
        int ctr = 0;
        for (String cat : categories) {
            catToId.put(cat, ctr);
            categoryArr.add(cat);
            ctr++;
        }

        // samples on cols, probes on rows
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(expressionFile, ",");
        double[][] averagesPerCat = new double[ds.rowObjects.size()][categories.size()];
        int[][] samplesPerCat = new int[ds.rowObjects.size()][categories.size()];


        // center and scale the columns if required.          
        if (normalize) {
        }

        for (int r = 0; r < ds.rawData.length; r++) {

            for (int c = 0; c < ds.rawData[r].length; c++) {

                String sample = ds.colObjects.get(c);
                String cat = sampleToCat.get(sample);
                Integer id = catToId.get(cat);

                if (id == null) {
                    System.err.println("No ID for sample: " + sample + "\t" + cat);
                } else {
                    averagesPerCat[r][id] += ds.rawData[r][c];
                    samplesPerCat[r][id]++;
                }
            }

        }

        for (int r = 0; r < averagesPerCat.length; r++) {
            for (int c = 0; c < averagesPerCat[r].length; c++) {
                averagesPerCat[r][c] /= samplesPerCat[r][c];
            }
        }


        DoubleMatrixDataset<String, String> averageOut = new DoubleMatrixDataset<String, String>();
        averageOut.rawData = averagesPerCat;
        averageOut.colObjects = categoryArr;
        averageOut.rowObjects = ds.rowObjects;
        averageOut.recalculateHashMaps();
        averageOut.save(fileOut);

        return averageOut;
    }
}
