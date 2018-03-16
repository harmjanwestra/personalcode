/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import umcg.genetica.gwas.Independifier;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class DetermineSignificanceUsingFDR {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String cellTypeSpecMatLoc = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/MetaAnalysis-PC1OnlyMerged/MetaAnalysisZScoreMatrix.txt";
        String probeTranslationFile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";
        String probeTranslationSource = "Probe";
        String probeTranslationDest = "HT12v3.txt";
        String geneColumnId = "Gene";
        String permutationDirectory = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Unsorted/";
        int nrPermutations = 10;
        double fdrTrheshold = 0.05;
        String outputfolder = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/MetaAnalysisFDREstimates-PC1OnlyMerged/";
        String referenceDs = "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/Merged/";
        int window = 250000;
        double r2 = 0.2;
        Double effectiveTests = 6134.249734365972;

        try {
            DetermineSignificanceUsingFDR blaat = new DetermineSignificanceUsingFDR();
            blaat.run(cellTypeSpecMatLoc, probeTranslationFile, probeTranslationSource, probeTranslationDest, geneColumnId, permutationDirectory, effectiveTests, nrPermutations, fdrTrheshold, outputfolder, referenceDs, window, r2);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String cellTypeSpecMatLoc, String probeTranslationFile, String probeTranslationSource, String probeTranslationDest, String geneColumnId,
            String permutationDirectory, Double effectiveTests, int nrPermutations, double fdrThreshold, String outputFolder, String referenceDataset, int window, double r2) throws IOException {

        Gpio.createDir(outputFolder);
        DoubleMatrixDataset<String, String> cellTypeSpecificityMatrix = new DoubleMatrixDataset<String, String>(cellTypeSpecMatLoc);

        // get all SNP probe combos
        HashSet<String> snpProbeCombos = new HashSet<String>();
        snpProbeCombos.addAll(cellTypeSpecificityMatrix.colObjects);

        // load probe translation file
        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> metaProbeToHT12v3 = pbt.getProbeTranslation(probeTranslationFile, probeTranslationSource, probeTranslationDest);
        HashMap<String, String> HT12v3ToGene = pbt.getProbeTranslation(probeTranslationFile, probeTranslationDest, geneColumnId);

        double nrOfEffectiveTests;
        if (effectiveTests == null) {
            // get permuted p-values from meta-analysis permutation files
            double[] lowestPValues = new double[nrPermutations];
            for (int permutation = 1; permutation < nrPermutations; permutation++) {

                String fileName = permutationDirectory + "PermutedEQTLsPermutationRound" + permutation + ".txt.gz";
                System.out.println("Parsing filename: " + fileName);
                TextFile tf = new TextFile(fileName, TextFile.R);
                tf.readLine();
                lowestPValues[permutation] = 1;
                String[] elems = tf.readLineElems(TextFile.tab);

                while (elems != null) {
                    String snp = elems[1];
                    String probe = metaProbeToHT12v3.get(elems[4]);

                    String snpprobe = snp + "-" + probe;
                    if (snpProbeCombos.contains(snpprobe)) {
                        double pval = Double.parseDouble(elems[0]);
                        if (pval < lowestPValues[permutation]) {
                            lowestPValues[permutation] = pval;
                        }
                    }
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
                System.out.println("Lowest pvalue: " + lowestPValues[permutation]);

            }
            double averagePValue = JSci.maths.ArrayMath.mean(lowestPValues);
            nrOfEffectiveTests = 0.5 / averagePValue;
            System.out.println("Average P: " + averagePValue + "\tNr. Effective Tests: " + nrOfEffectiveTests);
        } else {
            nrOfEffectiveTests = effectiveTests;
        }

        TriTyperGenotypeData referenceGenotypes = new TriTyperGenotypeData(referenceDataset);
        SNPLoader genotypeLoader = referenceGenotypes.createSNPLoader();
        Independifier independifier = new Independifier(referenceGenotypes, genotypeLoader);

        String[] snpsPerEQTL = new String[cellTypeSpecificityMatrix.colObjects.size()];
        String[] probesPerEQTL = new String[cellTypeSpecificityMatrix.colObjects.size()];
        for (int i = 0; i < cellTypeSpecificityMatrix.colObjects.size(); i++) {
            String[] eqtlelems = cellTypeSpecificityMatrix.colObjects.get(i).split("-");
            snpsPerEQTL[i] = eqtlelems[0];
            probesPerEQTL[i] = eqtlelems[1];
        }

        // now iterate over all rows (covariates)
        double[][] fdrmatrix = new double[cellTypeSpecificityMatrix.nrRows][cellTypeSpecificityMatrix.nrCols];
        for (int covariate = 0; covariate < cellTypeSpecificityMatrix.nrRows; covariate++) {
            double[] zScores = cellTypeSpecificityMatrix.rawData[covariate];
            double[] pvalues = new double[zScores.length];
            HashMap<Double, Integer> pvalueCtr = new HashMap<Double, Integer>();
            for (int eqtl = 0; eqtl < zScores.length; eqtl++) {
                double p = ZScores.zToP(zScores[eqtl]);

                Integer ctr = pvalueCtr.get(p);
                if (ctr == null) {
                    ctr = 0;
                }
                ctr++;
                pvalueCtr.put(p, ctr);
                pvalues[eqtl] = p;
            }

            // get the unique pvalues
            Set<Double> uniquePvalues = pvalueCtr.keySet();

            // sort the pvalues
            double[] pvaluecopy = Primitives.toPrimitiveArr(uniquePvalues.toArray(new Double[0]));
            Arrays.sort(pvaluecopy);
            HashMap<Double, Double> fdrPerPvalue = new HashMap<Double, Double>();
            int cumulativeEQTLCounter = 0;
            int nrSignificant = 0;
            for (int i = 0; i < pvaluecopy.length; i++) {
                double p = pvaluecopy[i];
                double expected = p * nrOfEffectiveTests;
                Integer nrOfEQTLsWithPvalue = pvalueCtr.get(p);
                cumulativeEQTLCounter += nrOfEQTLsWithPvalue;
                double fdr = expected / cumulativeEQTLCounter;
                if (fdr < fdrThreshold) {
                    nrSignificant++;
                }
                fdrPerPvalue.put(p, fdr);
            }

            // output significant results
            String covariateName = cellTypeSpecificityMatrix.rowObjects.get(covariate);
            String gene = HT12v3ToGene.get(covariateName);
            for (int eqtl = 0; eqtl < pvalues.length; eqtl++) {
                double p = pvalues[eqtl];
                double fdr = fdrPerPvalue.get(p);
                fdrmatrix[covariate][eqtl] = fdr;
            }

            if (nrSignificant > 0 && covariateName.equals("CellTypeInteractionZScore")) {
                TextFile outputFile = new TextFile(outputFolder + covariateName + "-" + gene + ".txt", TextFile.W);
                outputFile.writeln("EQTL\tGene\tZScore\tPvalue\tFDR");
                HashMap<String, Double> topFxA = new HashMap<String, Double>();
                HashMap<String, Double> topFxB = new HashMap<String, Double>();
                HashMap<String, Double> topFxC = new HashMap<String, Double>();
                for (int eqtl = 0; eqtl < pvalues.length; eqtl++) {
                    double p = pvalues[eqtl];
                    double fdr = fdrPerPvalue.get(p);
                    double z = zScores[eqtl];
                    String eQTLName = cellTypeSpecificityMatrix.colObjects.get(eqtl);
                    String probe = probesPerEQTL[eqtl];
                    String snp = snpsPerEQTL[eqtl];

                    if (fdr <= fdrThreshold) {
                        if (z < 0) {
                            Double z1 = topFxA.get(snp);
                            if (z1 == null) {
                                topFxA.put(snp, z);
                            } else {
                                if (z < z1) {
                                    topFxA.put(snp, z);
                                }
                            }
                        }
                        if (z > 0) {
                            Double z1 = topFxC.get(snp);
                            if (z1 == null) {
                                topFxC.put(snp, z);
                            } else {
                                if (z > z1) {
                                    topFxC.put(snp, z);
                                }
                            }
                        }
                    } else {
                        Double z1 = topFxB.get(snp);
                        if (z1 == null) {
                            topFxB.put(snp, z);
                        } else {
                            if (Math.abs(z) > Math.abs(z1)) {
                                topFxC.put(snp, z);
                            }
                        }
                    }
                    String probegene = HT12v3ToGene.get(probe);
                    outputFile.writeln(eQTLName + "\t" + probegene + "\t" + zScores[eqtl] + "\t" + pvalues[eqtl] + "\t" + fdr);
                }
                outputFile.close();

                // now we have three groups of SNPs. Prune them and output
                if (covariateName.equals("CellTypeInteractionZScore")) {
                    String[] topFxAPruned = independifier.independify(topFxA.keySet().toArray(new String[0]), r2, window);
                    String[] topFxBPruned = independifier.independify(topFxB.keySet().toArray(new String[0]), r2, window);
                    String[] topFxCPruned = independifier.independify(topFxC.keySet().toArray(new String[0]), r2, window);

                    String[] topFxASNPs = new String[topFxAPruned.length];
                    String[] topFxBSNPs = new String[topFxBPruned.length];
                    String[] topFxCSNPs = new String[topFxCPruned.length];

                    // take top effect per cluster
                    for (int i = 0; i < topFxAPruned.length; i++) {
                        String s = topFxAPruned[i];
                        String[] cluster = s.split(";");
                        if (cluster.length == 1) {
                            topFxASNPs[i] = s;
                        } else {
                            String topFx = null;
                            double maxZ = Double.MAX_VALUE;
                            for (String s2 : cluster) {
                                double z = topFxA.get(s2);
                                if (z < maxZ) {
                                    maxZ = z;
                                    topFx = s2;
                                }
                            }
                            topFxASNPs[i] = topFx;
                        }
                    }
                    for (int i = 0; i < topFxBPruned.length; i++) {
                        String s = topFxBPruned[i];
                        String[] cluster = s.split(";");
                        if (cluster.length == 1) {
                            topFxBSNPs[i] = s;
                        } else {
                            String topFx = null;
                            double maxZ = -Double.MAX_VALUE;
                            for (String s2 : cluster) {
                                double z = topFxB.get(s2);
                                if (z > maxZ) {
                                    maxZ = z;
                                    topFx = s2;
                                }
                            }
                            topFxBSNPs[i] = topFx;
                        }

                    }
                    for (int i = 0; i < topFxCPruned.length; i++) {
                        String s = topFxCPruned[i];
                        String[] cluster = s.split(";");
                        if (cluster.length == 1) {
                            topFxCSNPs[i] = s;
                        } else {
                            String topFx = null;
                            double maxZ = 0;
                            for (String s2 : cluster) {
                                double z = topFxC.get(s2);
                                if (Math.abs(z) > Math.abs(maxZ)) {
                                    maxZ = z;
                                    topFx = s2;
                                }
                            }
                            topFxCSNPs[i] = topFx;
                        }
                    }



                    TextFile outputFile2 = new TextFile(outputFolder + covariateName + "-" + gene + "-PrunedSNPsPerCategory.txt", TextFile.W);
                    outputFile2.writeln("SNP\tCategory\tKey\teQTLGenes");
                    for (int i = 0; i < topFxASNPs.length; i++) {
                        HashSet<String> eQTLGenes = new HashSet<String>();
                        HashSet<String> allSNPsInCategory = new HashSet<String>();
                        allSNPsInCategory.addAll(Arrays.asList(topFxASNPs));
                        for (int q = 0; q < snpsPerEQTL.length; q++) {
                            if (allSNPsInCategory.contains(snpsPerEQTL[i]) && pvalues[q] < fdrThreshold && zScores[i] < 0) {
                                String probe = probesPerEQTL[i];
                                String eqtlgene = HT12v3ToGene.get(probe);
                                if (!eqtlgene.equals("-")) {
                                    eQTLGenes.add(eqtlgene);
                                }
                            }
                        }
                        String eqtlgenes = Strings.concat(eQTLGenes.toArray(new String[0]), Strings.semicolon);
                        outputFile2.writeln(topFxASNPs[i] + "\tA\tZ<threshold\t" + eqtlgenes);
                    }

                    for (int i = 0; i < topFxBSNPs.length; i++) {
                        HashSet<String> eQTLGenes = new HashSet<String>();
                        HashSet<String> allSNPsInCategory = new HashSet<String>();
                        allSNPsInCategory.addAll(Arrays.asList(topFxBSNPs));
                        for (int q = 0; q < snpsPerEQTL.length; q++) {
                            if (allSNPsInCategory.contains(snpsPerEQTL[i]) && pvalues[q] > fdrThreshold) {
                                String probe = probesPerEQTL[i];
                                String eqtlgene = HT12v3ToGene.get(probe);
                                if (!eqtlgene.equals("-")) {
                                    eQTLGenes.add(eqtlgene);
                                }
                            }
                        }
                        String eqtlgenes = Strings.concat(eQTLGenes.toArray(new String[0]), Strings.semicolon);
                        outputFile2.writeln(topFxBSNPs[i] + "\tB\tZ<threshold\t" + eqtlgenes);
                    }

                    for (int i = 0; i < topFxCSNPs.length; i++) {
                        HashSet<String> eQTLGenes = new HashSet<String>();
                        HashSet<String> allSNPsInCategory = new HashSet<String>();
                        allSNPsInCategory.addAll(Arrays.asList(topFxCSNPs));
                        for (int q = 0; q < snpsPerEQTL.length; q++) {
                            if (allSNPsInCategory.contains(snpsPerEQTL[i]) && pvalues[q] < fdrThreshold && zScores[i] > 0) {
                                String probe = probesPerEQTL[i];
                                String eqtlgene = HT12v3ToGene.get(probe);
                                if (!eqtlgene.equals("-")) {
                                    eQTLGenes.add(eqtlgene);
                                }
                            }
                        }
                        String eqtlgenes = Strings.concat(eQTLGenes.toArray(new String[0]), Strings.semicolon);
                        outputFile2.writeln(topFxCSNPs[i] + "\tC\tZ>threshold\t" + eqtlgenes);
                    }
                    outputFile2.close();

                }


            }

        }

        DoubleMatrixDataset<String, String> outputMatrix = new DoubleMatrixDataset<String, String>();
        outputMatrix.rawData = fdrmatrix;
        outputMatrix.colObjects = cellTypeSpecificityMatrix.colObjects;
        outputMatrix.rowObjects = cellTypeSpecificityMatrix.rowObjects;
        outputMatrix.recalculateHashMaps();
        outputMatrix.save(outputFolder + "MetaAnalysisZScoreMatrix.binary-FDR.binary");

        genotypeLoader.close();

    }
}
