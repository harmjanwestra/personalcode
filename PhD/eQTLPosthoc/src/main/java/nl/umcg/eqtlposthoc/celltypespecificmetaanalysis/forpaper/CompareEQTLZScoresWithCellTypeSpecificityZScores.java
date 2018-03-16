/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class CompareEQTLZScoresWithCellTypeSpecificityZScores {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String[] eqtlsFiles = new String[]{
            "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-09-04-CD4AndCD8Cells/CD4-eQTLs-HT12v3.txt.gz",
            "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-09-04-CD4AndCD8Cells/CD8-eQTLs-HT12v3.txt.gz",
            //            "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-06-BloodH8v2/eQTLs/eQTLs.txt.gz",
            // /Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-10-21-FairFaxBCells
            "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-09-04-StrangerLCL/eQTLs-medR.txt.gz",
            "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-10-21-FairFaxBCells/eQTLs.txt.gz",
            "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-10-21-FairFaxMonocytes/eQTLs.txt.gz",
            "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-11-19-AnandNeutrophils/2013-07-11-HT12v3SNPProbeCombos-eQTLFile.txt"// "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Unsorted/eQTLs-HT12v3.txt.gz"
        };

        String[] averageExpressionFiles = new String[]{
            "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-09-04-CD4AndCD8Cells/CD4-ExpressionMeanExpressionVariance-HT12v3.txt",
            "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-09-04-CD4AndCD8Cells/CD8-ExpressionMeanExpressionVariance-HT12v3.txt",
            //            "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-06-BloodH8v2/MeanExpresion.txt",

            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/MedianExpresionOverAllPopulations.txt",
            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/Fairfax/BCells/MeanExpresion.txt",
            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/Fairfax/Monocytes/MeanExpresion.txt",
            "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-11-19-AnandNeutrophils/MeanExpression.txt"//"/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/MeanExpresion.txt"
        };

        String[] datasetNames = new String[]{
            "CD4+ T-Cells",
            "CD8+ T-Cells",
            //            "Peripheral Blood (BloodH8v2)",

            "Lymphoblastoid cell lines",
            "B Cells",
            "Monocytes",
            "Neutrophils"// "Peripheral Blood"
        };

        boolean useCorrelationCoefficientAsEffectSize = true;

        String fdrFile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-09-27-FDREstimates/MetaAnalysisZScoreMatrix.binary-FDR.binary";
        String specificityfile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/MetaAnalysisZScoreMatrix.binary";
        String outputdir = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/Replication/2013-11-27-ReplicationInCellTypeSpecificDatasets/";
        String probeTranslationFile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";
        double pvaluethreshold = 1E-106;
        double fdrThreshold = 0.05;

        boolean onlyTopFx = false;
        try {
            CompareEQTLZScoresWithCellTypeSpecificityZScores q = new CompareEQTLZScoresWithCellTypeSpecificityZScores();
            q.run(eqtlsFiles, averageExpressionFiles, datasetNames, specificityfile, fdrFile, fdrThreshold, probeTranslationFile, pvaluethreshold, outputdir, onlyTopFx, useCorrelationCoefficientAsEffectSize);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String[] eqtlsFiles, String[] averageExpressionFiles, String[] datasetnames, String specificityFile, String fdrFile, double fdrThreshold, String probeTranslation, double pvalueThreshold, String outputdir, boolean onlyConsiderTopFx, boolean useCorrelationCoeff) throws IOException {
        outputdir = Gpio.formatAsDirectory(outputdir);
        Gpio.createDir(outputdir);

        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> probeToHUGO = pb.getProbeTranslation(probeTranslation, "HT12v3.txt", "Gene");
        DoubleMatrixDataset<String, String> dmd = new DoubleMatrixDataset<String, String>(specificityFile);
        DoubleMatrixDataset<String, String> dmdfdr = new DoubleMatrixDataset<String, String>(fdrFile);
        WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();

        ArrayList<HashMap<Pair<String, String>, Double>> eqtlsForAllDatasets = new ArrayList<HashMap<Pair<String, String>, Double>>();
        ArrayList<HashMap<String, Double>> eqtlsExpMeanForAllDatasets = new ArrayList<HashMap<String, Double>>();
        ArrayList<HashMap<String, Double>> eqtlsExpVariForAllDatasets = new ArrayList<HashMap<String, Double>>();
//        ArrayList<HashMap<Pair<String, String>, Double>> eqtlsExpVariForAllDatasets = new ArrayList<HashMap<Pair<String, String>, Double>>();

        for (int d = 0; d < datasetnames.length; d++) {
            String eqtlsFile = eqtlsFiles[d];
            System.out.println("Parsing file: " + eqtlsFile);
            TextFile tf = new TextFile(eqtlsFile, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            HashMap<Pair<String, String>, Double> eqtlsfordataset = new HashMap<Pair<String, String>, Double>();
            HashMap<String, Double> eqtlsfordatasetexpmean = new HashMap<String, Double>();
            HashMap<String, Double> eqtlsfordatasetexpvari = new HashMap<String, Double>();
//            HashMap<Pair<String, String>, Double> eqtlsfordatasetexpvari = new HashMap<Pair<String, String>, Double>();
            while (elems != null) {
                String snp = elems[1];
                String probe = elems[4];
                String z = elems[10];
                if (useCorrelationCoeff) {
                    z = elems[17];
                }
                double dz = Math.abs(Double.parseDouble(z));
                if (useCorrelationCoeff) {
                    // dz *= dz;
                }

                Pair<String, String> e = new Pair<String, String>(snp, probe);
                eqtlsfordataset.put(e, dz);

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile tf2 = new TextFile(averageExpressionFiles[d], TextFile.R);
            System.out.println("Parsing file: " + averageExpressionFiles[d]);
            elems = tf2.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems.length == 3) {
                    String probe = elems[0];
                    String val = elems[1];
                    String val2 = elems[2];
                    double em = Double.parseDouble(val);
                    double ev = Math.sqrt(Double.parseDouble(val2));
                    eqtlsfordatasetexpmean.put(probe, em);
                    eqtlsfordatasetexpvari.put(probe, ev);
                }
                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            System.out.println(eqtlsfordatasetexpmean.size() + " expression means loaded.");

            eqtlsForAllDatasets.add(eqtlsfordataset);
            eqtlsExpMeanForAllDatasets.add(eqtlsfordatasetexpmean);
            eqtlsExpVariForAllDatasets.add(eqtlsfordatasetexpvari);
        }

        TextFile outputfile = new TextFile(outputdir + "Comparison.txt", TextFile.W);

        String header = "Covariate\tGene";

        String yAxisName = "eQTL effect size (Z-Score)";
        for (String datasetname : datasetnames) {
            header += "\t" + datasetname + "-Specific-N";
            header += "\t" + datasetname + "-Generic-N";
            header += "\t" + datasetname + "-Aspecific-N";
            header += "\t" + datasetname + "-SpecificVsAspecific-P";
            header += "\t" + datasetname + "-SpecificVsAspecific-AUC";
            header += "\t" + datasetname + "-SpecificVsGeneric-P";
            header += "\t" + datasetname + "-SpecificVsGeneric-AUC";
            header += "\t" + datasetname + "-AspecificVsGeneric-P";
            header += "\t" + datasetname + "-AspecificVsGeneric-AUC";
        }
        outputfile.writeln(header);

        ProgressBar progressbar = new ProgressBar(dmd.nrRows);
        for (int covariate = 0; covariate < dmd.nrRows; covariate++) {
            String covariateName = dmd.rowObjects.get(covariate);

            if (covariateName.equals("CellTypeInteractionZScore")) {

                TextFile tfOut = new TextFile(outputdir + covariateName + "-Summary.txt", TextFile.W);
                String summaryHeader = "eQTL";
                for (String datasetname : datasetnames) {
                    summaryHeader += "\t" + datasetname + "-ZScore";
                    summaryHeader += "\t" + datasetname + "-AverageExpression";
                }
                tfOut.writeln(summaryHeader);

                boolean plotCovariate = false;

                double[][][] valuesPerDataset = new double[datasetnames.length][3][0];
                double[][][] valuesVariancePerDataset = new double[datasetnames.length][3][0];
                double[][][] valuesMeanExpPerDataset = new double[datasetnames.length][3][0];
                String[][] xLabels = new String[datasetnames.length][3];
                String[][] xLabelsMeanExp = new String[datasetnames.length][3];
                for (int d = 0; d < datasetnames.length; d++) {
                    xLabels[d][0] = "Lymphocyte specific";
                    xLabelsMeanExp[d][0] = "Lymphocyte specific";
                    xLabels[d][1] = "Generic";
                    xLabelsMeanExp[d][1] = "Generic";
                    xLabels[d][2] = "Neutrophil specific";
                    xLabelsMeanExp[d][2] = "Neutrophil specific";
                }

                String[] outputStrings = new String[dmd.nrCols];

                String outputStr = covariateName + "\t" + probeToHUGO.get(covariateName);

                for (int d = 0; d < datasetnames.length; d++) {
                    ArrayList<Double> cellTypeSpecificityZScores = new ArrayList<Double>();
                    ArrayList<Double> eQTLZScores = new ArrayList<Double>();

                    HashMap<Pair<String, String>, Double> eqtls = eqtlsForAllDatasets.get(d);
                    HashMap<String, Double> eqtlsMeanExp = eqtlsExpMeanForAllDatasets.get(d);
                    HashMap<String, Double> eqtlsVarianceExp = eqtlsExpVariForAllDatasets.get(d);

                    ArrayList<Double> neutrospecific = new ArrayList<Double>();
                    ArrayList<Double> neutroaspecific = new ArrayList<Double>();
                    ArrayList<Double> othereqtls = new ArrayList<Double>();

                    ArrayList<Double> neutrospecificmeanexp = new ArrayList<Double>();
                    ArrayList<Double> neutroaspecificmeanexp = new ArrayList<Double>();
                    ArrayList<Double> othereqtlsmeanexp = new ArrayList<Double>();

                    ArrayList<Double> neutrospecificexpvariance = new ArrayList<Double>();
                    ArrayList<Double> neutroaspecificexpvariance = new ArrayList<Double>();
                    ArrayList<Double> othereqtlsexpvariance = new ArrayList<Double>();

                    HashSet<String> visitedProbesA = new HashSet<String>();
                    HashSet<String> visitedProbesB = new HashSet<String>();
                    HashSet<String> visitedProbesC = new HashSet<String>();

                    for (int eqtlId = 0; eqtlId < dmd.nrCols; eqtlId++) {
                        String eqtl = dmd.colObjects.get(eqtlId);
                        if (outputStrings[eqtlId] == null) {
                            outputStrings[eqtlId] = eqtl;
                        }

                        String[] snpprobelems = eqtl.split("-");
                        String snp = snpprobelems[0];
                        String probe = snpprobelems[1];
                        double celltypeSpecificityZScore = dmd.rawData[covariate][eqtlId];
                        Pair<String, String> p = new Pair<String, String>(snp, probe);
                        Double eqtlZScore = eqtls.get(p);
                        Double eqtlMeanExp = eqtlsMeanExp.get(probe);
                        Double eqtlVariance = eqtlsVarianceExp.get(probe);

                        double fdr = dmdfdr.rawData[covariate][eqtlId];

                        if (eqtlZScore == null) {
                            outputStrings[eqtlId] += "\tNaN\tNaN";
                        } else {
                            outputStrings[eqtlId] += "\t" + eqtlZScore + "\t" + eqtlMeanExp;
                        }

                        if (eqtlZScore != null) {
                            if (eqtlMeanExp == null) {
                                System.err.println("ERROR: we have QTL, but no mean expression??? " + probe + "\tDataset: " + d + "\t" + datasetnames[d]);
                            } else {
                                if (celltypeSpecificityZScore >= 0 && fdr < fdrThreshold) {

                                    if (!onlyConsiderTopFx || !visitedProbesA.contains(probe)) {
                                        neutrospecific.add(eqtlZScore);
                                        neutrospecificmeanexp.add(eqtlMeanExp);
                                        neutrospecificexpvariance.add(eqtlVariance);
                                        visitedProbesA.add(probe);
                                    }
                                } else if (celltypeSpecificityZScore < 0 && fdr < fdrThreshold) {

                                    if (!onlyConsiderTopFx || !visitedProbesB.contains(probe)) {
                                        neutroaspecific.add(eqtlZScore);
                                        neutroaspecificmeanexp.add(eqtlMeanExp);
                                        neutroaspecificexpvariance.add(eqtlVariance);
                                        visitedProbesB.add(probe);
                                    }
                                } else {

                                    if (!onlyConsiderTopFx || !visitedProbesC.contains(probe)) {
                                        othereqtls.add(eqtlZScore);
                                        othereqtlsmeanexp.add(eqtlMeanExp);
                                        othereqtlsexpvariance.add(eqtlVariance);
                                        visitedProbesC.add(probe);
                                    }
                                }

                                cellTypeSpecificityZScores.add(celltypeSpecificityZScore);
                                eQTLZScores.add(eqtlZScore);
                            }
                        }
                    } // end iterate over cis-eQTLs

                    double[] specificEQTLZScores = Primitives.toPrimitiveArr(neutrospecific.toArray(new Double[0]));
                    double[] aspecificEQTLZScores = Primitives.toPrimitiveArr(neutroaspecific.toArray(new Double[0]));
                    double[] genericEQTLZScores = Primitives.toPrimitiveArr(othereqtls.toArray(new Double[0]));

                    double[] specificMeanExpression = Primitives.toPrimitiveArr(neutrospecificmeanexp.toArray(new Double[0]));
                    double[] aspecificMeanExpression = Primitives.toPrimitiveArr(neutroaspecificmeanexp.toArray(new Double[0]));
                    double[] genericMeanExpression = Primitives.toPrimitiveArr(othereqtlsmeanexp.toArray(new Double[0]));

                    double[] specificVariance = Primitives.toPrimitiveArr(neutrospecificexpvariance.toArray(new Double[0]));
                    double[] aspecificVariance = Primitives.toPrimitiveArr(neutroaspecificexpvariance.toArray(new Double[0]));
                    double[] genericVariance = Primitives.toPrimitiveArr(othereqtlsexpvariance.toArray(new Double[0]));

                    xLabels[d][0] += " (n=" + aspecificEQTLZScores.length + ")";
                    xLabels[d][1] += " (n=" + genericEQTLZScores.length + ")";
                    xLabels[d][2] += " (n=" + specificEQTLZScores.length + ")";

                    xLabelsMeanExp[d][0] += " (n=" + aspecificMeanExpression.length + ")";
                    xLabelsMeanExp[d][1] += " (n=" + genericMeanExpression.length + ")";
                    xLabelsMeanExp[d][2] += " (n=" + specificMeanExpression.length + ")";

                    valuesPerDataset[d][0] = aspecificEQTLZScores;
                    valuesPerDataset[d][1] = genericEQTLZScores;
                    valuesPerDataset[d][2] = specificEQTLZScores;

                    valuesMeanExpPerDataset[d][0] = aspecificMeanExpression;
                    valuesMeanExpPerDataset[d][1] = genericMeanExpression;
                    valuesMeanExpPerDataset[d][2] = specificMeanExpression;

                    valuesVariancePerDataset[d][0] = aspecificVariance;
                    valuesVariancePerDataset[d][1] = genericVariance;
                    valuesVariancePerDataset[d][2] = specificVariance;

                    double p1 = mwm.returnWilcoxonMannWhitneyPValue(specificEQTLZScores, aspecificEQTLZScores);
                    double auc1 = mwm.getAUC();
                    double p2 = mwm.returnWilcoxonMannWhitneyPValue(specificEQTLZScores, genericEQTLZScores);
                    double auc2 = mwm.getAUC();
                    double p3 = mwm.returnWilcoxonMannWhitneyPValue(aspecificEQTLZScores, genericEQTLZScores);
                    double auc3 = mwm.getAUC();

                    if (p1 < pvalueThreshold) {
//                    plotCovariate = true;
                    }

                    // Create scatterplot: x-axis == celltype specificity ZScore, y-axis == eQTL Zscore in tissue
                    if (plotCovariate || covariateName.equals("CellTypeInteractionZScore") //                            || covariateName.equals("7330546")
                            //                            || covariateName.equals("2650750")
                            //                            || covariateName.equals("4920202")
                            ) {
                        System.out.println("Plotting covariate: " + dmd.rowObjects.get(covariate));
                        double[] cellTypeSpecificityZScoresArray = Primitives.toPrimitiveArr(cellTypeSpecificityZScores.toArray(new Double[0]));
                        double[] eQTLZScoreArray = Primitives.toPrimitiveArr(eQTLZScores.toArray(new Double[0]));
                        String outputFileName = outputdir + covariateName + "-" + datasetnames[d] + "-" + probeToHUGO.get(covariateName);
                        if (!useCorrelationCoeff) {
                            new ScatterPlot(500, 500, eQTLZScoreArray, cellTypeSpecificityZScoresArray, ScatterPlot.OUTPUTFORMAT.PDF, outputFileName + "-Scatterplot.pdf");
                        }
                    }
                    outputStr += "\t" + specificEQTLZScores.length + "\t" + genericEQTLZScores.length + "\t" + aspecificEQTLZScores.length + "\t" + p1 + "\t" + auc1 + "\t" + p2 + "\t" + auc2 + "\t" + p3 + "\t" + auc3;
                }

                for (String s : outputStrings) {
                    tfOut.writeln(s);
                }
                tfOut.close();
                
                outputfile.writeln(outputStr);
                String outputFileName = outputdir + covariateName + "-" + probeToHUGO.get(covariateName);

                if (plotCovariate || covariateName.equals("CellTypeInteractionZScore") //                        || covariateName.equals("7330546")
                        //                        || covariateName.equals("2650750")
                        //                        || covariateName.equals("4920202")
                        ) {
                    ViolinBoxPlot vbp = new ViolinBoxPlot();
                    String plotYAxis = "eQTL effect size (Z-Score)";
                    if (useCorrelationCoeff) {
                        plotYAxis = "eQTL effect size (R)";
                    }
                    vbp.draw(valuesPerDataset, datasetnames, xLabels, plotYAxis, ViolinBoxPlot.Output.PDF, outputFileName + "-ViolinBoxPlot.pdf", false);
                    vbp.draw(valuesMeanExpPerDataset, datasetnames, xLabelsMeanExp, "Mean expression", ViolinBoxPlot.Output.PDF, outputFileName + "-ExpMean-ViolinBoxPlot.pdf", false);
                    vbp.draw(valuesVariancePerDataset, datasetnames, xLabelsMeanExp, "Expression standard deviation", ViolinBoxPlot.Output.PDF, outputFileName + "-ExpVariance-ViolinBoxPlot.pdf", false);
                }
            }

            progressbar.iterate();
        }
        progressbar.close();
        outputfile.close();
    }
}
