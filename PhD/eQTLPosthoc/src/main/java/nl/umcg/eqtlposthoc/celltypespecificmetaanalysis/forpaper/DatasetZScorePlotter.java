/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.CellTypeSpecificMetaAnalsysisDatasetContainer;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;
import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class DatasetZScorePlotter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {
            DatasetZScorePlotter m = new DatasetZScorePlotter();
            boolean flipEffects = true;
            String output = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/2013-11-27-MetaAnalysisZScorePlots/";
            Gpio.createDir(output);
            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            String[] datasetFileDirs = new String[]{
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-12-GroningenHT12v3/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-23_SHIP-TREND/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-24_ROTTERDAM-STUDY/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-25-EGCUT/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-29-DILGOM/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-08-19-KORA-F4/",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-08-19-INCHIANTI/"
            };

            String[] datasetNames = new String[]{
                "GroningenHT12v3",
                "SHIP-TREND",
                "Rotterdam",
                "EGCUT",
                "DILGOM",
                "KORA-F4",
                "INCHIANTI"
            };

//            Double neutFXFDR = 2.605579924;
//            System.out.println("Nominal PValue: " + ZScores.zToP(neutFXFDR));
            Double neutFXFDR = null;

//            String[] datasetFileDirs = new String[]{
//                "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/EGCUTHT12v3-PC1Only/",
//                "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/GroningenBloodHT12v3-PC1Only/",
//                "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/RotterdamHT12v3-PC1Only/"
//            };
//
//            String[] datasetNames = new String[]{
//                "EGCUT",
//                "GroningenHT12v3",
//                "Rotterdam"
//            };
//            boolean[] platformIsHT12v3 = new boolean[]{true, true, false};
            boolean[] platformIsHT12v3 = new boolean[]{
                true,
                true,
                false,
                true,
                true,
                true,
                true
            };

            m.run(output, probetranslationfile, datasetFileDirs, datasetNames, platformIsHT12v3, false, flipEffects, neutFXFDR);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String output, String probetranslationfile, String[] datasetFileDirs, String[] datasetNames, boolean[] platformIsHT12v3, boolean allowNaNs, boolean flipeffects, Double neutrofxzscorethreshold) throws IOException {

        Gpio.createDir(output);

        String ht12v3String = "HumanHT-12_V3_0_R2_11283641_A.txt";
        String ht12v4String = "HumanHT-12_V4_0_R1_15002873_B.txt";

        String hugoStr = "HUGO";

        //       PROBE LOOKUP TABLES..
        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> ht12v4ToHT12v3 = new HashMap<String, String>();
        HashMap<String, String> ht12v3ToHugo = new HashMap<String, String>();
        ht12v4ToHT12v3 = pbt.getProbeTranslation(probetranslationfile, ht12v4String, ht12v3String);
        ht12v3ToHugo = pbt.getProbeTranslation(probetranslationfile, ht12v3String, hugoStr);

        HashSet<String> requiredRows = new HashSet<String>();
        requiredRows.add("MainEffectZScore");
        requiredRows.add("CellTypeInteractionZScore");

        // load the datasets..
        CellTypeSpecificMetaAnalsysisDatasetContainer container = new CellTypeSpecificMetaAnalsysisDatasetContainer(datasetFileDirs, platformIsHT12v3, ht12v4ToHT12v3, requiredRows);
        SNP[][] snps = container.getSnps();

        ArrayList<DoubleMatrixDataset<String, String>> datasets = container.getDatasets();
        HashSet<String> allSNPs = container.getAllSNPs();
        HashSet<String> allEQTL = container.getAllEQTL();

        // hash the SNPs for easy lookup..
        Integer[][] snpIndex = container.getSnpIndex(); // for each SNP in the meta analysis, this index points to a SNP in a dataset (in the jagged allSNPs array)
        HashMap<String, Integer> snpToId = container.getSnpToId();

        // hash the eQTLs as well..
        Integer[][] eQTLIndex = container.geteQTLIndex();
        HashMap<String, Integer> eQTLToId = container.geteQTLToId();
        String[] allEQTLArr = allEQTL.toArray(new String[0]);

        // hash the probes / covariates...
        HashSet<String> allCovariates = container.getAllCovariates();
        String[] allCovariateArr = container.getAllCovariateArr();
        Integer[][] covariateIndex = container.getCovariateIndex();
        HashMap<String, Integer> covariateToId = container.getCovariateToId();

        // done hashing
        // create the output
        ZScorePlot zsMainEffect = new ZScorePlot();
        ZScorePlot zsNeutEffect = new ZScorePlot();
        String[] datasetNames2 = new String[datasetNames.length + 1];
        System.arraycopy(datasetNames, 0, datasetNames2, 0, datasetNames.length);
        datasetNames2[datasetNames2.length - 1] = "Meta-Analysis";
        String tableOut = "";
        if (neutrofxzscorethreshold != null) {
            zsMainEffect.init(datasets.size() + 1, datasetNames2, true, output + "MainEffectZScoreComparison-FDRSignificant.pdf");
            zsNeutEffect.init(datasets.size() + 1, datasetNames2, true, output + "InteractionEffectZScoreComparison-FDRSignificant.pdf");
            tableOut = output + "ComparisonTable-FDRSignificant.txt";
        } else {
            zsMainEffect.init(datasets.size() + 1, datasetNames2, true, output + "MainEffectZScoreComparison.pdf");
            zsNeutEffect.init(datasets.size() + 1, datasetNames2, true, output + "InteractionEffectZScoreComparison.pdf");
            tableOut = output + "ComparisonTable.txt";
        }

        TextFile tableOutput = new TextFile(tableOut, TextFile.W);
        String header = "eQTL\tMetaMainZ\tMetaInteractionZ\tMetaMainR\tMetaInteractionR\tMetaAnalysisN";
        for (String ds : datasetNames) {
            header += "\t" + ds + "-MainZ\t" + ds + "-InteractionZ\t" + ds + "-N";
        }
        tableOutput.writeln(header);
        for (int e = 0; e < allEQTLArr.length; e++) {

            String eQTL = allEQTLArr[e];
            String snp = eQTL.split("-")[0];
            String probe = eQTL.split("-")[1];
            Integer snpId = snpToId.get(snp);

            String eQTLOutput = eQTL;

            if (snpId == null) {
                System.err.println("ERROR snpId " + snp + " is null for eQTL: " + eQTL);
            }

            String firstSNPAlleles = null;
            String firstSNPAssessedAllele = null;
            String firstSNPName = null;
            Boolean[] flipEffect = new Boolean[datasets.size()];

            String mafStr = "";
            String crStr = "";
            String nrCalledStr = "";
            String hwepStr = "";

            String dsNames = "";

            for (int d = 0; d < datasets.size(); d++) {
                Integer snpIdInDataset = snpIndex[d][snpId];
                if (snpIdInDataset != null) {
                    SNP snpObj = snps[d][snpIdInDataset];
                    String SNPAlleles = BaseAnnot.getAllelesDescription(snpObj.getAlleles());
                    String SNPAssessedAllele = BaseAnnot.toString(snpObj.getMinorAllele());
                    if (firstSNPAlleles == null) {

                        firstSNPAlleles = SNPAlleles;
                        firstSNPAssessedAllele = SNPAssessedAllele;
                        flipEffect[d] = false;
                        firstSNPName = snpObj.getName();
                    } else {
                        flipEffect[d] = BaseAnnot.flipalleles(firstSNPAlleles, firstSNPAssessedAllele, SNPAlleles, SNPAssessedAllele);

                        String secondSNPName = snpObj.getName();
                        if (!firstSNPName.equals(secondSNPName)) {
                            System.err.println("SNPs not linked together properly");
                        }
                    }

                    mafStr += snpObj.getMAF() + ";";
                    crStr += snpObj.getCR() + ";";
                    nrCalledStr += snpObj.getNrCalled() + ";";
                    hwepStr += snpObj.getHWEP() + ";";
                    dsNames += datasetNames[d] + ";";

                } else {
//                    System.out.println(snp + " not present in dataset: " + datasetNames[d]);
                    mafStr += "-;";
                    crStr += "-;";
                    nrCalledStr += "-;";
                    hwepStr += "-;";
                    dsNames += "-;";
                    flipEffect[d] = null;
                }
            }

            // first meta-analyze the main effect z-score
            Integer mainEffectId = covariateToId.get("MainEffectZScore");
            Integer neutEffectId = covariateToId.get("CellTypeInteractionZScore");

            Triple<double[], int[], Double> resultsMainEffect = metaAnalyze(mainEffectId, snpId, e, eQTL, flipEffect, snps, datasets, snpIndex, eQTLIndex, covariateIndex, allCovariateArr);
            Triple<double[], int[], Double> resultsNeutroEffect = metaAnalyze(neutEffectId, snpId, e, eQTL, flipEffect, snps, datasets, snpIndex, eQTLIndex, covariateIndex, allCovariateArr);

            int[] samples = resultsMainEffect.getMiddle();
            int[] samples2 = resultsNeutroEffect.getMiddle();
            int sum = 0;
            int sum2 = 0;
            for (int i = 0; i < samples.length; i++) {
                sum += samples[i];
                sum2 += samples2[i];
            }
            if (sum != sum2) {
                System.err.println("ERROR: differens sums of samples for eQTL: " + eQTL);
            }

            Double metaZMain = resultsMainEffect.getRight();
            Double metaZNeutro = resultsNeutroEffect.getRight();

            if (metaZMain < 0) {
                metaZNeutro *= -1;
            }

            double explainedVarianceMainEffect = ZScores.zScoreToCorrelation(metaZMain, sum);
            double explainedVarianceInteractionEffect = ZScores.zScoreToCorrelation(metaZNeutro, sum);
            eQTLOutput += "\t" + metaZMain + "\t" + metaZNeutro + "\t" + explainedVarianceMainEffect + "\t" + explainedVarianceInteractionEffect + "\t" + sum;
            
            for (int dataset1 = 0; dataset1 < resultsMainEffect.getLeft().length; dataset1++) {
                double mainEffect = resultsMainEffect.getLeft()[dataset1];
                double neutroEffect = resultsNeutroEffect.getLeft()[dataset1];

                int sampleSize = samples[dataset1];
                if (mainEffect < 0) {
                    neutroEffect *= -1;
                }
                eQTLOutput += "\t" + mainEffect + "\t" + neutroEffect + "\t" + sampleSize;

                if (!Double.isNaN(mainEffect) && !Double.isNaN(neutroEffect)) {

                    if (neutrofxzscorethreshold == null || Math.abs(metaZNeutro) > neutrofxzscorethreshold) {
                        for (int dataset2 = dataset1 + 1; dataset2 < resultsMainEffect.getLeft().length; dataset2++) {
                            double mainEffect2 = resultsMainEffect.getLeft()[dataset2];
                            double neutroEffect2 = resultsNeutroEffect.getLeft()[dataset2];
                            if (!Double.isNaN(mainEffect2) && !Double.isNaN(neutroEffect2)) {
                                if (mainEffect2 < 0) {
                                    neutroEffect2 *= -1;
                                }

                                zsMainEffect.draw(mainEffect, mainEffect2, dataset1, dataset2);
                                zsNeutEffect.draw(neutroEffect, neutroEffect2, dataset1, dataset2);
                            }
                        }

                        zsMainEffect.draw(mainEffect, metaZMain, dataset1, datasetNames2.length - 1);
                        zsNeutEffect.draw(neutroEffect, metaZNeutro, dataset1, datasetNames2.length - 1);
                    }
                }

            }
            tableOutput.writeln(eQTLOutput);
        }
        tableOutput.close();
        try {
            if (neutrofxzscorethreshold != null) {
                zsMainEffect.write(output + "MainEffectZScoreComparison-FDRSignificant.pdf");
                zsNeutEffect.write(output + "InteractionEffectZScoreComparison-FDRSignificant.pdf");
            } else {
                zsMainEffect.write(output + "MainEffectZScoreComparison.pdf");
                zsNeutEffect.write(output + "InteractionEffectZScoreComparison.pdf");
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private Triple<double[], int[], Double> metaAnalyze(Integer c, Integer snpId, Integer e, String eQTL, Boolean[] flipEffect, SNP[][] snps, ArrayList<DoubleMatrixDataset<String, String>> datasets, Integer[][] snpIndex, Integer[][] eQTLIndex, Integer[][] covariateIndex, String[] allCovariateArr) {

        double[] zScores = new double[datasets.size()];
        int totalSample = 0;
        int[] sampleSize = new int[datasets.size()];
        for (int d = 0; d < datasets.size(); d++) {
            Integer snpIdInDataset = snpIndex[d][snpId];
            Integer eQTLIdInDataset = eQTLIndex[d][e];
            Integer covariateIdInDataset = covariateIndex[d][c];

            if (snpIdInDataset != null && eQTLIdInDataset != null && covariateIdInDataset != null) {
                SNP snpObj = snps[d][snpIdInDataset];
                sampleSize[d] = snpObj.getNrCalled();
                double z = datasets.get(d).rawData[covariateIdInDataset][eQTLIdInDataset];
                if (sampleSize[d] > 0) {
                    if (Double.isNaN(z)) {
                        System.err.println("ERROR z is NaN " + eQTL + "\t" + allCovariateArr[c]);
                    }
                    if (flipEffect[d]) {
                        z *= -1;
                    }
                    zScores[d] = z;
                    totalSample += sampleSize[d];
                } else {
                    if (z != 0d) {
                        System.out.println("Call rate is 0, but eQTL has Z-score? " + eQTL + "\t" + allCovariateArr[c] + "\t" + z);
                    }
                    zScores[d] = Double.NaN;
                    sampleSize[d] = 0;
                }
            } else {
                zScores[d] = Double.NaN;
                sampleSize[d] = 0;
            }
        }

        double metaZ = ZScores.getWeightedZ(zScores, sampleSize);
        if (Double.isNaN(metaZ)) {
            for (int d = 0; d < datasets.size(); d++) {
                Integer snpIdInDataset = snpIndex[d][snpId];
                Integer eQTLIdInDataset = eQTLIndex[d][e];
                Integer covariateIdInDataset = covariateIndex[d][c];
            }
        } else {
            double P = ZScores.zToP(metaZ);
        }

        return new Triple<double[], int[], Double>(zScores, sampleSize, metaZ);

    }

}
