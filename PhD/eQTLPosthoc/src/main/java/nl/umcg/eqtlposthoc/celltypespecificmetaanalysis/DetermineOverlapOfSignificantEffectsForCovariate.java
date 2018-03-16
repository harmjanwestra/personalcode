/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class DetermineOverlapOfSignificantEffectsForCovariate {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String metaFile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/MetaAnalysis-PC1OnlyMerged/MetaAnalysisZScoreMatrix.txt";
            String fdrFile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/MetaAnalysisFDREstimates-PC1OnlyMerged/MetaAnalysisZScoreMatrix.binary-FDR.binary";
            String probeTranslationFile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";
            String output = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/2013-10-09-MED8CovariateOverlap-PC1OnlyMerged/";
            String query = "CellTypeInteractionZScore";
            DetermineOverlapOfSignificantEffectsForCovariate f = new DetermineOverlapOfSignificantEffectsForCovariate();
            f.run(metaFile, fdrFile, probeTranslationFile, output, query);
        } catch (IOException e) {
            e.printStackTrace();

        }
    }

    private void run(String metaFile, String fdrFile, String probeTranslationFile, String output, String query) throws IOException {
        Gpio.createDir(output);
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(metaFile);
        DoubleMatrixDataset<String, String> dsfdr = new DoubleMatrixDataset<String, String>(fdrFile);

        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> probeToGene = pbt.getProbeTranslation(probeTranslationFile, "HT12v3.txt", "Gene");

        ArrayList<Integer> queryProbes = new ArrayList<Integer>();
        for (int row = 0; row < ds.nrRows; row++) {
            String probeName = ds.rowObjects.get(row);
            String gene = probeToGene.get(probeName);
            if ((gene != null && gene.contains(query)) || probeName.contains(query)) {
                queryProbes.add(row);
            }
        }

        for (Integer row : queryProbes) {
            String probeName = ds.rowObjects.get(row);
            TextFile outputfile = new TextFile(output + probeName + "-" + query + ".txt", TextFile.W);

            outputfile.writeln("Covariate\tCovariateGene\tJaccardIndex\tTotalOverlapEQTLs\tOppositeEffectsEQTLs\tOverlapPositiveEQTLs\tOverlapNegativeEQTLs\tOverlapTotalGenes\tOverlapOppositeGenes\tOverlapPositiveEQTLGenes\tOverlapNegativeEQTLGenes");

            // determine significant effects
            HashSet<Integer> significantEffectPositive = new HashSet<Integer>();
            HashSet<Integer> significantEffectNegative = new HashSet<Integer>();
            for (int col = 0; col < ds.nrCols; col++) {
                double fdr = dsfdr.rawData[row][col];
                double z = ds.rawData[row][col];
                if (fdr < 0.05 && z >= 0) {
                    significantEffectPositive.add(col);
                } else if (fdr < 0.05 && z < 0) {
                    significantEffectNegative.add(col);
                }
            }

            // now iterate over all rows..
            for (int row2 = 0; row2 < ds.nrRows; row2++) {

                String probeName2 = ds.rowObjects.get(row2);
                outputfile.append(probeName2 + "\t" + probeToGene.get(probeName2));

                HashSet<Integer> significantEffectPositive2 = new HashSet<Integer>();
                HashSet<Integer> significantEffectNegative2 = new HashSet<Integer>();

                for (int col = 0; col < ds.nrCols; col++) {
                    double fdr = dsfdr.rawData[row2][col];
                    double z = ds.rawData[row2][col];
                    if (fdr < 0.05 && z >= 0) {
                        significantEffectPositive2.add(col);
                    } else if (fdr < 0.05 && z < 0) {
                        significantEffectNegative2.add(col);
                    }
                }

                HashSet<String> positiveOverlapEQTLs = new HashSet<String>();
                HashSet<String> negativeOverlapEQTLs = new HashSet<String>();
                HashSet<String> totalOverlapEQTLs = new HashSet<String>();

                HashSet<String> oppositeEffects = new HashSet<String>();

                HashSet<String> negativeOverlapGenes = new HashSet<String>();
                HashSet<String> positiveOverlapGenes = new HashSet<String>();
                HashSet<String> oppositeOverlapGenes = new HashSet<String>();
                HashSet<String> totalOverlapGenes = new HashSet<String>();

                for (Integer i : significantEffectPositive2) {
                    if (significantEffectPositive.contains(i)) {
                        String eqtl = ds.colObjects.get(i);
                        positiveOverlapEQTLs.add(eqtl);
                        String probe = eqtl.split("-")[1];
                        String gene = probeToGene.get(probe);
                        positiveOverlapGenes.add(gene);
                        totalOverlapEQTLs.add(eqtl);
                        totalOverlapGenes.add(gene);
                    } else if (significantEffectNegative.contains(i)) {
                        String eqtl = ds.colObjects.get(i);
                        String probe = eqtl.split("-")[1];
                        String gene = probeToGene.get(probe);
                        totalOverlapEQTLs.add(eqtl);
                        totalOverlapGenes.add(gene);
                        oppositeEffects.add(eqtl);
                        oppositeOverlapGenes.add(gene);
                    }
                }

                for (Integer i : significantEffectNegative2) {
                    if (significantEffectNegative.contains(i)) {

                        String eqtl = ds.colObjects.get(i);
                        negativeOverlapEQTLs.add(eqtl);
                        String probe = eqtl.split("-")[1];
                        String gene = probeToGene.get(probe);
                        negativeOverlapGenes.add(gene);
                        totalOverlapGenes.add(gene);
                        totalOverlapEQTLs.add(eqtl);
                    } else if (significantEffectPositive.contains(i)) {
                        String eqtl = ds.colObjects.get(i);
                        String probe = eqtl.split("-")[1];
                        String gene = probeToGene.get(probe);
                        totalOverlapEQTLs.add(eqtl);
                        totalOverlapGenes.add(gene);
                        oppositeEffects.add(eqtl);
                        oppositeOverlapGenes.add(gene);
                    }
                }

                HashSet<Integer> union = new HashSet<Integer>();
                union.addAll(significantEffectNegative);
                union.addAll(significantEffectNegative2);
                union.addAll(significantEffectPositive);
                union.addAll(significantEffectPositive2);
                double jaccard = (double) totalOverlapEQTLs.size() / union.size();

                if (row.equals(row2)) {
                    if (jaccard < 1d) {
                        System.err.println("ERROR! Jaccard < 1");
                        for (Integer i : significantEffectNegative) {
                            String eqtl = ds.colObjects.get(i);
                            if (!significantEffectNegative2.contains(i)) {
                                System.err.println("Missing negative: " + eqtl);
                            }
                            if (!totalOverlapEQTLs.contains(eqtl)) {
                                System.err.println("Missing in total (from neg) " + eqtl);
                            }
                        }

                        for (Integer i : significantEffectPositive) {
                            String eqtl = ds.colObjects.get(i);
                            if (!significantEffectPositive2.contains(i)) {
                                System.err.println("Missing positive: " + ds.colObjects.get(i));
                            }
                            if (!totalOverlapEQTLs.contains(eqtl)) {
                                System.err.println("Missing in total (from pos) " + eqtl);
                            }
                        }

                    }
                }

                // TotalOverlapEQTLs\tOverlapPositiveEQTLs\tOverlapNegativeEQTLs\tOverlapTotalGenes\tOverlapPositiveEQTLGenes\tOverlapNegativeEQTLGenes
                outputfile.append("\t" + jaccard + "\t" + totalOverlapEQTLs.size()
                        + "\t" + oppositeEffects.size()
                        + "\t" + positiveOverlapEQTLs.size() + "\t" + negativeOverlapEQTLs.size()
                        + "\t" + Strings.concat(totalOverlapGenes.toArray(new String[0]), Strings.comma)
                        + "\t" + Strings.concat(oppositeOverlapGenes.toArray(new String[0]), Strings.comma)
                        + "\t" + Strings.concat(positiveOverlapGenes.toArray(new String[0]), Strings.comma)
                        + "\t" + Strings.concat(negativeOverlapGenes.toArray(new String[0]), Strings.comma)
                        + "\n");
            }





            outputfile.close();
        }


    }
}
