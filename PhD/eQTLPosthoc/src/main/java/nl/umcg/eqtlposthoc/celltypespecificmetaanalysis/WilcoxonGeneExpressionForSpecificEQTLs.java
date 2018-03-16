/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
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
public class WilcoxonGeneExpressionForSpecificEQTLs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        double zthreshold = 3;
        String metaAnalysisMatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
        String probeAnnotationFile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";
        String platformString = "HT12v3.txt";
        String geneString = "Gene";
        String geneExpressionMatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/SeboRNASeq/expression_table.genes.exonic_v69.0.3.rawCounts.log2.qn.txt.centeredscaled.txt";
        String outputdir = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/WilcoxonGeneExpressionForSpecificEQTLs/";
        WilcoxonGeneExpressionForSpecificEQTLs w = new WilcoxonGeneExpressionForSpecificEQTLs();
        try {
            w.run(metaAnalysisMatrix, probeAnnotationFile, platformString, geneString, geneExpressionMatrix, zthreshold, outputdir);
        } catch (IOException ex) {
            Logger.getLogger(WilcoxonGeneExpressionForSpecificEQTLs.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private void run(String metaAnalysisMatrix, String probeAnnotationFile, String platformString, String geneString, String geneExpressionMatrix, double zthreshold, String outputdir) throws IOException {

        Gpio.createDir(outputdir);
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaAnalysisMatrix);
        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> platformToHUGO = pbt.getProbeTranslation(probeAnnotationFile, platformString, geneString);

        DoubleMatrixDataset<String, String> expressionMatrix = new DoubleMatrixDataset<String, String>(geneExpressionMatrix);

        WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();
        TextFile tfOut = new TextFile(outputdir + "Results.txt", TextFile.W);


        String header = "Covariate\tHugo\tN(Z>=3)\tN(Z>=3Bg)\tN(Z<=-3)\tN(Z<=-3Bg)";
        for (String celltype : expressionMatrix.colObjects) {
            header += "\t" + celltype + "-WilcoxonP (Z>=3)\t" + celltype + "-WilcoxonAUC (Z>=3)\t" + celltype + "-WilcoxonP (Z<=-3)\t" + celltype + "-WilcoxonAUC (Z<=-3)";
        }

        tfOut.writeln(header);
        DecimalFormat df = new DecimalFormat("0.E00");

        // correlate the eQTLs with cell type gene expression

        for (int covariate = 0; covariate < metaMatrix.nrRows; covariate++) {
            String covariateName = metaMatrix.rowObjects.get(covariate);
            String covariateStr = covariateName + "\t" + platformToHUGO.get(covariateName);
            String outputStr = "";


            double[][][] vals = new double[expressionMatrix.nrCols][2][0];

            double[] pvals = new double[expressionMatrix.nrCols];

            boolean plotVBP = false;

            int nrpos = 0;
            int nrneg = 0;
            int nrposbg = 0;
            int nrnegbg = 0;
            for (int celltype = 0; celltype < expressionMatrix.nrCols; celltype++) {
                ArrayList<Double> valsPositive = new ArrayList<Double>();
                ArrayList<Double> valsPositiveBgSet = new ArrayList<Double>();
                ArrayList<Double> valsNegative = new ArrayList<Double>();
                ArrayList<Double> valsNegativeBgSet = new ArrayList<Double>();

                HashSet<String> visitedPositiveGenes = new HashSet<String>();
                HashSet<String> visitedPositiveBgGenes = new HashSet<String>();
                HashSet<String> visitedNegativeGenes = new HashSet<String>();
                HashSet<String> visitedNegativeBgGenes = new HashSet<String>();


                for (int eqtl = 0; eqtl < metaMatrix.nrCols; eqtl++) {
                    double z = metaMatrix.rawData[covariate][eqtl];
                    String eqtlProbe = metaMatrix.colObjects.get(eqtl).split("-")[1];
                    String hugo = platformToHUGO.get(eqtlProbe);
                    Integer expressionRowId = expressionMatrix.hashRows.get(hugo);

                    if (expressionRowId != null && !hugo.equals("-")) {
                        double expressionValue = expressionMatrix.rawData[expressionRowId][celltype];
                        if (z >= zthreshold && !visitedPositiveGenes.contains(hugo)) {
                            valsPositive.add(expressionValue);
                            visitedPositiveGenes.add(hugo);
                        }
                        if (z <= -zthreshold && !visitedNegativeGenes.contains(hugo)) {
                            valsNegative.add(expressionValue);
                            visitedNegativeGenes.add(hugo);
                        }
                    }
                }
                
                for (int eqtl = 0; eqtl < metaMatrix.nrCols; eqtl++) {
                    double z = metaMatrix.rawData[covariate][eqtl];
                    String eqtlProbe = metaMatrix.colObjects.get(eqtl).split("-")[1];
                    String hugo = platformToHUGO.get(eqtlProbe);
                    Integer expressionRowId = expressionMatrix.hashRows.get(hugo);

                    if (expressionRowId != null && !hugo.equals("-")) {
                        double expressionValue = expressionMatrix.rawData[expressionRowId][celltype];
                        if (!visitedPositiveBgGenes.contains(hugo) && z < zthreshold && !visitedPositiveGenes.contains(hugo)) {
                            valsPositiveBgSet.add(expressionValue);
                            visitedPositiveBgGenes.add(hugo);
                        }
                        if (!visitedNegativeBgGenes.contains(hugo) && z > -zthreshold && !visitedNegativeGenes.contains(hugo)) {
                            valsNegativeBgSet.add(expressionValue);
                            visitedNegativeBgGenes.add(hugo);
                        }
                    }
                }

                double[] xarrd = Primitives.toPrimitiveArr(valsPositive.toArray(new Double[0]));
                double[] xarrdbg = Primitives.toPrimitiveArr(valsPositiveBgSet.toArray(new Double[0]));
                double[] yarrd = Primitives.toPrimitiveArr(valsNegative.toArray(new Double[0]));
                double[] yarrdbg = Primitives.toPrimitiveArr(valsNegativeBgSet.toArray(new Double[0]));

                double p = mwm.returnWilcoxonMannWhitneyPValue(xarrd, xarrdbg);
                double auc = mwm.getAUC();
                double p2 = mwm.returnWilcoxonMannWhitneyPValue(yarrd, yarrdbg);
                double auc2 = mwm.getAUC();

                vals[celltype][0] = xarrd;
                vals[celltype][1] = xarrdbg;

                if (p < 1E-10) {
                    plotVBP = true;
                }

                nrpos = xarrd.length;
                nrposbg = xarrdbg.length;
                nrneg = yarrd.length;
                nrnegbg = yarrdbg.length;

                pvals[celltype] = p;
                outputStr += "\t" + p + "\t" + auc + "\t" + p2 + "\t" + auc2;

            }

            tfOut.writeln(covariateStr + "\t" + nrpos + "\t" + nrposbg + "\t" + nrneg + "\t" + nrnegbg + outputStr);

            if (plotVBP) {
                ViolinBoxPlot vbp = new ViolinBoxPlot();
                String[][] xLabels = new String[expressionMatrix.nrCols][2];
                for (int i = 0; i < expressionMatrix.nrCols; i++) {
                    xLabels[i][0] = "Z>=3";
                    xLabels[i][1] = "Z<=-3";
                }

                String[] plotnames = expressionMatrix.colObjects.toArray(new String[0]);
                double lowestpval = 1;
                for (int i = 0; i < plotnames.length; i++) {
                    if (pvals[i] < lowestpval) {
                        lowestpval = pvals[i];
                    }
                    plotnames[i] += " (" + df.format(pvals[i]) + ")";
                }


                int height = 500;
                int margin = 50;
                int plotwidth = 100;
                int betweenplotmargin = 25;
                int width = (expressionMatrix.nrCols * plotwidth) + (betweenplotmargin * expressionMatrix.nrCols - 1) + (2 * margin);
//                vbp.draw(vals, plotnames, xLabels, width, height, margin, betweenplotmargin, ViolinBoxPlot.Output.PDF, outputdir + df.format(lowestpval) + "-" + covariateStr.replace("\t", "-") + ".pdf");
            }
//            
        }
        tfOut.close();
    }
}
