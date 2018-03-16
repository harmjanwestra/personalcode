/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class JaccardIndex {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        JaccardIndex i = new JaccardIndex();
        try {
            String meta = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/MetaAnalysisZScoreMatrix.binary";
            String fdr = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-09-27-FDREstimates/MetaAnalysisZScoreMatrix.binary-FDR.binary";


            double fdrthreshold = 0.0;
            String outputdir = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-10-02-JaccardIndex/";
            Gpio.createDir(outputdir);
            String output = outputdir + "JaccardIndex.txt";
            i.run(meta, fdr, fdrthreshold, output);
            String output2 = output + "-Filtered.txt";
            i.filterJaccardIndexTable(output, output2);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String meta, String fdr, double fdrthreshold, String output) throws IOException {
        DoubleMatrixDataset<String, String> m = new DoubleMatrixDataset<String, String>(meta);
        DoubleMatrixDataset<String, String> f = new DoubleMatrixDataset<String, String>(meta);




        ArrayList<HashSet<Integer>> significantEffectsPerCovariate = new ArrayList<HashSet<Integer>>();



        HashSet<Integer> rowsWithMoreThan1000SignificantEffects = new HashSet<Integer>();
        for (int row = 0; row < m.nrRows; row++) {

            HashSet<Integer> significantCisEffects = new HashSet<Integer>();

            for (int e = 0; e < m.nrCols; e++) {
                if (f.rawData[row][e] <= fdrthreshold) {
                    significantCisEffects.add(e);
                }
            }
            significantEffectsPerCovariate.add(significantCisEffects);

            if (significantEffectsPerCovariate.size() >= 1000) {
                rowsWithMoreThan1000SignificantEffects.add(row);
            }
//            for (int row2 = 0; row2 < m.nrRows; row2++) {
//                HashSet<String> significantCisEffects2 = new HashSet<String>();
//
//                for (int e = 0; e < m.nrCols; e++) {
//                    if (f.rawData[row2][e] <= fdrthreshold) {
//                        significantCisEffects2.add(m.colObjects.get(e));
//                    }
//                }
//
//                // calculate jaccard
//            }

        }


        double[][] jaccardindex = new double[rowsWithMoreThan1000SignificantEffects.size()][rowsWithMoreThan1000SignificantEffects.size()];
        ProgressBar pb = new ProgressBar(m.nrRows);
        int row1Ctr = 0;
        ArrayList<String> rowObjs = new ArrayList<String>();
        for (int row = 0; row < m.nrRows; row++) {
            if (rowsWithMoreThan1000SignificantEffects.contains(row)) {
                HashSet<Integer> significantCisEffects = significantEffectsPerCovariate.get(row);
                int row2Ctr = 0;
                for (int row2 = 0; row2 < m.nrRows; row2++) {
                    if (rowsWithMoreThan1000SignificantEffects.contains(row2)) {
                        HashSet<Integer> significantCisEffects2 = significantEffectsPerCovariate.get(row2);
                        HashSet<Integer> union = new HashSet<Integer>();
                        union.addAll(significantCisEffects);
                        union.addAll(significantCisEffects2);

                        int intersect = 0;
                        for (Integer s : significantCisEffects) {
                            if (significantCisEffects2.contains(s)) {
                                intersect++;
                            }
                        }
                        double jaccard = (double) intersect / union.size();
                        jaccardindex[row1Ctr][row2Ctr] = jaccard;
                        row2Ctr++;
                    }
                }
                jaccardindex[row][row] = 1.0;
                row1Ctr++;
                rowObjs.add(m.rowObjects.get(row));
            }
            pb.set(row);
        }

        pb.close();

        DoubleMatrixDataset<String, String> outputfile = new DoubleMatrixDataset<String, String>();
        outputfile.rawData = jaccardindex;
        outputfile.colObjects = rowObjs;
        outputfile.rowObjects = rowObjs;
        outputfile.recalculateHashMaps();
        outputfile.save(output);


    }

    public void filterJaccardIndexTable(String in, String out) throws IOException {
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(in);
        HashSet<Integer> rowsToExclude = new HashSet<Integer>();

        int[] frequencyDistribution = new int[100];


        for (int row = 0; row < ds.nrRows; row++) {
            int nrNull = 0;


            for (int col = 0; col < ds.nrCols; col++) {
                double val = ds.rawData[row][col];
                if (!Double.isNaN(val)) {
                    int bin = (int) Math.round(val * frequencyDistribution.length);
                    if (bin == frequencyDistribution.length) {
                        bin--;
                    }
                    frequencyDistribution[bin]++;
                }
            }

            for (int col = 0; col < ds.nrCols; col++) {

                double val = ds.rawData[row][col];


                if (col != row) {
                    if (val <= 0.5 || Double.isNaN(val)) {
                        nrNull++;
                    }
                }
            }

            if (nrNull == ds.nrCols || nrNull == ds.nrCols - 1) {
                System.out.println("Excluding row: " + row + "\t" + ds.rowObjects.get(row));
                rowsToExclude.add(row);
            }
        }

        System.out.println("Rows remaining: " + (ds.nrRows - rowsToExclude.size()));

        for (int bin = 0; bin < frequencyDistribution.length; bin++) {
            System.out.println(((double) bin / frequencyDistribution.length) + "\t" + frequencyDistribution[bin]);
        }

    }
}
