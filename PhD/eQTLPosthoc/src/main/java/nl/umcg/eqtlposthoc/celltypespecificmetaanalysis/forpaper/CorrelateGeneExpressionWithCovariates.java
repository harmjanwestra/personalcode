/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class CorrelateGeneExpressionWithCovariates {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String fdrMatrix = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-09-23-FDREstimates/MetaAnalysisZScoreMatrix.binary-FDR.binary";
        double fdrThreshold = 0.05;
        String metaAnalysisMatrix = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/MetaAnalysisZScoreMatrix.binary";
        String probeAnnotationFile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-18-ProbeAnnotationFile.txt";
        String platformString = "HT12v3.txt";
        String geneString = "Gene";
        String geneExpressionMatrix = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/SeboRNASeq/expression_table.genes.exonic_v69.0.3.rawCounts.log2.qn.txt.centeredscaled.txt";
        String outputdir = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-09-24-CorrelationsOfCovariatesWithGeneExpressionFDR0.05/";

        try {
            CorrelateGeneExpressionWithCovariates g = new CorrelateGeneExpressionWithCovariates();
            // g.runWithFDR(outputdir, metaAnalysisMatrix, probeAnnotationFile, platformString, geneString, geneExpressionMatrix, fdrMatrix, fdrThreshold);
            boolean onlyIncludePositiveZScores = true;
            g.runIterateOverCovariates(outputdir, metaAnalysisMatrix, probeAnnotationFile, platformString, geneString, geneExpressionMatrix, fdrMatrix, fdrThreshold, onlyIncludePositiveZScores);

            onlyIncludePositiveZScores = false;
            g.runIterateOverCovariates(outputdir, metaAnalysisMatrix, probeAnnotationFile, platformString, geneString, geneExpressionMatrix, fdrMatrix, fdrThreshold, onlyIncludePositiveZScores);

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void runIterateOverCovariates(String outputdir, String metaAnalysisMatrix, String probeAnnotationFile, String platformString, String geneString, String geneExpressionMatrix,
            String fdrMatrixFile, double fdrThreshold, boolean onlyincludepositiveZScores) throws IOException {
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaAnalysisMatrix);
        DoubleMatrixDataset<String, String> fdrMatrix = new DoubleMatrixDataset<String, String>(fdrMatrixFile);

        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> platformToHUGO = pbt.getProbeTranslation(probeAnnotationFile, platformString, geneString);

        DoubleMatrixDataset<String, String> expressionMatrix = new DoubleMatrixDataset<String, String>(geneExpressionMatrix);

        String outfilename = outputdir + "ResultsCorrelationOverCovariates.txt";
        if (onlyincludepositiveZScores) {
            outfilename = outputdir + "ResultsCorrelationOverCovariatesOnlyPositiveZScores.txt";
        }
        TextFile tfOut = new TextFile(outfilename, TextFile.W);

        String header = "Covariate\tHugo\tN\tNWithDuplicateGenes";
        for (String celltype : expressionMatrix.colObjects) {
            header += "\t" + celltype + "-SPEARMAN\t"
                    + celltype + "-SPEARMANWITHDUPLICATEGENES\t"
                    + celltype + "-ZScore\t"
                    + celltype + "-ZScoreWITHDUPLICATEGENES\t"
                    + celltype + "-PValue\t"
                    + celltype + "-PValueWITHDUPLICATEGENES";
        }


        tfOut.writeln(header);

        for (int covariate = 0; covariate < metaMatrix.nrRows; covariate++) {

            String covariateName = metaMatrix.rowObjects.get(covariate);
            String covariateGene = platformToHUGO.get(covariateName);

            Integer n1 = null;
            Integer n2 = null;

            String cellTypeOutput = "";

            SpearmansCorrelation c = new SpearmansCorrelation();
            // iterate over the celltypes
            for (int celltype = 0; celltype < expressionMatrix.nrCols; celltype++) {
                // iterate over the covariates..

                ArrayList<Double> valsX = new ArrayList<Double>();
                ArrayList<Double> valsY = new ArrayList<Double>();
                ArrayList<Double> valsXWDups = new ArrayList<Double>();
                ArrayList<Double> valsYWDups = new ArrayList<Double>();
                HashSet<String> visitedGenes = new HashSet<String>();
                for (int eqtlId = 0; eqtlId < metaMatrix.nrCols; eqtlId++) {
                    String eQTL = metaMatrix.colObjects.get(eqtlId);
                    String eQTLProbe = eQTL.split("-")[1];


                    if (fdrMatrix.rawData[covariate][eqtlId] < fdrThreshold && (!onlyincludepositiveZScores || (onlyincludepositiveZScores && metaMatrix.rawData[covariate][eqtlId] > 0))) {
                        String hugo = platformToHUGO.get(eQTLProbe);
                        if (hugo != null && !hugo.equals("-")) {
                            Integer rowIdInExpresionMatrix = expressionMatrix.hashRows.get(hugo);
                            if (rowIdInExpresionMatrix == null) {
                                // split the hugo.
                                String[] hugoElems = hugo.split(";");
                                for (String hugoElem : hugoElems) {
                                    if (rowIdInExpresionMatrix == null) {
                                        rowIdInExpresionMatrix = expressionMatrix.hashRows.get(hugoElem);
                                    }
                                }
                            }

                            if (rowIdInExpresionMatrix != null) {
                                if (!visitedGenes.contains(hugo)) {
                                    valsX.add(metaMatrix.rawData[covariate][eqtlId]);
                                    valsY.add(expressionMatrix.rawData[rowIdInExpresionMatrix][celltype]);
                                    visitedGenes.add(hugo);
                                }

                                valsXWDups.add(metaMatrix.rawData[covariate][eqtlId]);
                                valsYWDups.add(expressionMatrix.rawData[rowIdInExpresionMatrix][celltype]);
                            }
                        }
                    }
                }


                if (valsX.size() > 1) {
                    // correlate //
                    double[] x1 = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                    double[] y1 = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));
                    double[] x2 = Primitives.toPrimitiveArr(valsXWDups.toArray(new Double[0]));
                    double[] y2 = Primitives.toPrimitiveArr(valsYWDups.toArray(new Double[0]));

                    if (n1 == null) {
                        n1 = x1.length;
                    }

                    if (n2 == null) {
                        n2 = x2.length;
                    }

                    double sp1 = c.correlation(x1, y1);
                    double sp2 = c.correlation(x2, y2);
                    Correlation.correlationToZScore(x1.length);
                    Correlation.correlationToZScore(x2.length);
                    double z1 = Correlation.convertCorrelationToZScore(x1.length, sp1);
                    double z2 = Correlation.convertCorrelationToZScore(x2.length, sp2);
                    cellTypeOutput += "\t" + sp1 + "\t" + sp2 + "\t" + z1 + "\t" + z2 + "\t" + ZScores.zToP(z1) + "\t" + ZScores.zToP(z2);
                } else {
                    cellTypeOutput += "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + 1.0 + "\t" + 1.0;
                }
            }

            tfOut.writeln(covariateName + "\t" + covariateGene + "\t" + n1 + "\t" + n2 + cellTypeOutput);

        }

        tfOut.close();


    }

    public void runWithFDR(String outputdir, String metaAnalysisMatrix, String probeAnnotationFile, String platformString, String geneString, String geneExpressionMatrix,
            String fdrMatrixFile, double fdrThreshold) throws IOException {
        Gpio.createDir(outputdir);
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaAnalysisMatrix);
        DoubleMatrixDataset<String, String> fdrMatrix = new DoubleMatrixDataset<String, String>(fdrMatrixFile);

        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> platformToHUGO = pbt.getProbeTranslation(probeAnnotationFile, platformString, geneString);

        DoubleMatrixDataset<String, String> expressionMatrix = new DoubleMatrixDataset<String, String>(geneExpressionMatrix);


        TextFile tfOut = new TextFile(outputdir + "Results.txt", TextFile.W);
        String header = "eQTL\tHugo\tN\tNWithDuplicateGenes";
        for (String celltype : expressionMatrix.colObjects) {
            header += "\t" + celltype + "-SPEARMAN\t" + celltype + "-SPEARMANWITHDUPLICATEGENES\t" + celltype + "-ZScore\t" + celltype + "-ZScoreWITHDUPLICATEGENES";
        }

        tfOut.writeln(header);
        DecimalFormat df = new DecimalFormat("0.E00");

        for (int eqtlId = 0; eqtlId < metaMatrix.nrCols; eqtlId++) {
            String eQTL = metaMatrix.colObjects.get(eqtlId);
            String eQTLProbe = eQTL.split("-")[1];
            String eQTLProbeHugo = platformToHUGO.get(eQTLProbe);

            Integer n1 = null;
            Integer n2 = null;

            String cellTypeOutput = "";

            SpearmansCorrelation c = new SpearmansCorrelation();
            // iterate over the celltypes
            for (int celltype = 0; celltype < expressionMatrix.nrCols; celltype++) {
                // iterate over the covariates..

                ArrayList<Double> valsX = new ArrayList<Double>();
                ArrayList<Double> valsY = new ArrayList<Double>();
                ArrayList<Double> valsXWDups = new ArrayList<Double>();
                ArrayList<Double> valsYWDups = new ArrayList<Double>();
                HashSet<String> visitedGenes = new HashSet<String>();
                for (int covariate = 0; covariate < metaMatrix.nrRows; covariate++) {
                    if (fdrMatrix.rawData[covariate][eqtlId] < fdrThreshold) {
                        String hugo = platformToHUGO.get(metaMatrix.rowObjects.get(covariate));
                        if (hugo != null && !hugo.equals("-")) {
                            Integer rowIdInExpresionMatrix = expressionMatrix.hashRows.get(hugo);
                            if (rowIdInExpresionMatrix == null) {
                                // split the hugo.
                                String[] hugoElems = hugo.split(";");
                                for (String hugoElem : hugoElems) {
                                    if (rowIdInExpresionMatrix == null) {
                                        rowIdInExpresionMatrix = expressionMatrix.hashRows.get(hugoElem);
                                    }
                                }
                            }

                            if (rowIdInExpresionMatrix != null) {
                                if (!visitedGenes.contains(hugo)) {
                                    valsX.add(metaMatrix.rawData[covariate][eqtlId]);
                                    valsY.add(expressionMatrix.rawData[rowIdInExpresionMatrix][celltype]);
                                    visitedGenes.add(hugo);
                                }

                                valsXWDups.add(metaMatrix.rawData[covariate][eqtlId]);
                                valsYWDups.add(expressionMatrix.rawData[rowIdInExpresionMatrix][celltype]);
                            }
                        }
                    }
                }


                if (valsX.size() > 1) {
                    // correlate //
                    double[] x1 = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                    double[] y1 = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));
                    double[] x2 = Primitives.toPrimitiveArr(valsXWDups.toArray(new Double[0]));
                    double[] y2 = Primitives.toPrimitiveArr(valsYWDups.toArray(new Double[0]));

                    if (n1 == null) {
                        n1 = x1.length;
                    }

                    if (n2 == null) {
                        n2 = x2.length;
                    }

                    double sp1 = c.correlation(x1, y1);
                    double sp2 = c.correlation(x2, y2);
                    Correlation.correlationToZScore(x1.length);
                    Correlation.correlationToZScore(x2.length);
                    double z1 = Correlation.convertCorrelationToZScore(x1.length, sp1);
                    double z2 = Correlation.convertCorrelationToZScore(x2.length, sp2);
                    cellTypeOutput += "\t" + sp1 + "\t" + sp2 + "\t" + z1 + "\t" + z2;
                } else {
                    cellTypeOutput += "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + 0;
                }
            }

            tfOut.writeln(metaMatrix.colObjects.get(eqtlId) + "\t" + eQTLProbeHugo + "\t" + n1 + "\t" + n2 + cellTypeOutput);

        }

        tfOut.close();


    }

    public void run(String outputdir, String metaAnalysisMatrix, String probeAnnotationFile, String platformString, String geneString, String geneExpressionMatrix) throws IOException {
        Gpio.createDir(outputdir);
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaAnalysisMatrix);
        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> platformToHUGO = pbt.getProbeTranslation(probeAnnotationFile, platformString, geneString);

        DoubleMatrixDataset<String, String> expressionMatrix = new DoubleMatrixDataset<String, String>(geneExpressionMatrix);


        TextFile tfOut = new TextFile(outputdir + "Results.txt", TextFile.W);
        String header = "eQTL\tHugo\tN\tNWithDuplicateGenes";
        for (String celltype : expressionMatrix.colObjects) {
            header += "\t" + celltype + "-SPEARMAN\t" + celltype + "-SPEARMANWITHDUPLICATEGENES\t" + celltype + "-ZScore\t" + celltype + "-ZScoreWITHDUPLICATEGENES";
        }

        tfOut.writeln(header);
        DecimalFormat df = new DecimalFormat("0.E00");

        for (int e = 0; e < metaMatrix.nrCols; e++) {
            String eQTL = metaMatrix.colObjects.get(e);
            String eQTLProbe = eQTL.split("-")[1];
            String eQTLProbeHugo = platformToHUGO.get(eQTLProbe);

            Integer n1 = null;
            Integer n2 = null;

            String cellTypeOutput = "";

            SpearmansCorrelation c = new SpearmansCorrelation();
            // iterate over the celltypes
            for (int celltype = 0; celltype < expressionMatrix.nrCols; celltype++) {
                // iterate over the covariates..            
                ArrayList<Double> valsX = new ArrayList<Double>();
                ArrayList<Double> valsY = new ArrayList<Double>();
                ArrayList<Double> valsXWDups = new ArrayList<Double>();
                ArrayList<Double> valsYWDups = new ArrayList<Double>();
                HashSet<String> visitedGenes = new HashSet<String>();
                for (int covariate = 0; covariate < metaMatrix.nrRows; covariate++) {
                    String hugo = platformToHUGO.get(metaMatrix.rowObjects.get(covariate));
                    if (hugo != null && !hugo.equals("-")) {
                        Integer rowIdInExpresionMatrix = expressionMatrix.hashRows.get(hugo);
                        if (rowIdInExpresionMatrix == null) {
                            // split the hugo.
                            String[] hugoElems = hugo.split(";");
                            for (String hugoElem : hugoElems) {
                                if (rowIdInExpresionMatrix == null) {
                                    rowIdInExpresionMatrix = expressionMatrix.hashRows.get(hugoElem);
                                }
                            }
                        }

                        if (rowIdInExpresionMatrix != null) {
                            if (!visitedGenes.contains(hugo)) {
                                valsX.add(metaMatrix.rawData[covariate][e]);
                                valsY.add(expressionMatrix.rawData[rowIdInExpresionMatrix][celltype]);
                                visitedGenes.add(hugo);
                            }

                            valsXWDups.add(metaMatrix.rawData[covariate][e]);
                            valsYWDups.add(expressionMatrix.rawData[rowIdInExpresionMatrix][celltype]);

                        }

                    }
                }

                // correlate //

                double[] x1 = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                double[] y1 = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));
                double[] x2 = Primitives.toPrimitiveArr(valsXWDups.toArray(new Double[0]));
                double[] y2 = Primitives.toPrimitiveArr(valsYWDups.toArray(new Double[0]));

                if (n1 == null) {
                    n1 = x1.length;
                }

                if (n2 == null) {
                    n2 = x2.length;
                }

                double sp1 = c.correlation(x1, y1);
                double sp2 = c.correlation(x2, y2);
                Correlation.correlationToZScore(x1.length);
                Correlation.correlationToZScore(x2.length);
                double z1 = Correlation.convertCorrelationToZScore(x1.length, sp1);
                double z2 = Correlation.convertCorrelationToZScore(x2.length, sp2);
                cellTypeOutput += "\t" + sp1 + "\t" + sp2 + "\t" + z1 + "\t" + z2;
            }

            tfOut.writeln(metaMatrix.colObjects.get(e) + "\t" + eQTLProbeHugo + "\t" + n1 + "\t" + n2 + cellTypeOutput);

        }

        tfOut.close();


    }
}
