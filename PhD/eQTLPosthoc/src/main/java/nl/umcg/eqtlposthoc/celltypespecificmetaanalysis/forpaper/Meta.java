/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.CellTypeSpecificMetaAnalsysisDatasetContainer;
import nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.CheckForNaNs;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class Meta {

    /**
     * @param args the command line arguments
     */
    public static void mainMeta(String[] args) {
        // TODO code application logic here

        try {
            Meta m = new Meta();
            boolean flipEffects = true;
            String output = "/Volumes/iSnackHD/AeroFS/2013-11-20-MetaAnalysisForPlots/";
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
            HashSet<String> covariatesToTest = new HashSet<String>();
            covariatesToTest.add("CellTypeInteractionZScore");
            covariatesToTest.add("MainEffectZScore");
            boolean[] platformIsHT12v3 = new boolean[]{
                true,
                true,
                false,
                true,
                true,
                true,
                true
            };

            m.run(output, probetranslationfile, datasetFileDirs, datasetNames, platformIsHT12v3, false, flipEffects, covariatesToTest, false);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        // TODO code application logic here

        try {
            Meta m = new Meta();
            boolean flipEffects = true;
            String output = "/Volumes/iSnackHD/Data/Projects/LudeFranke/2014-02-12-TransInteractionTerms/Meta/";
            Gpio.createDir(output);
            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            String[] datasetFileDirs = new String[]{
                "/Volumes/iSnackHD/Data/Projects/LudeFranke/2014-02-12-TransInteractionTerms/EGCUT/interaction/",
                "/Volumes/iSnackHD/Data/Projects/LudeFranke/2014-02-12-TransInteractionTerms/Groningen/interaction/",
                "/Volumes/iSnackHD/Data/Projects/LudeFranke/2014-02-12-TransInteractionTerms/Rotterdam/interaction/"
            };

            String[] datasetNames = new String[]{
                "EGCUT",
                "GroningenHT12v3",
                "Rotterdam"
            };

//            boolean[] platformIsHT12v3 = new boolean[]{true, true, false};
            HashSet<String> covariatesToTest = null;
//            covariatesToTest.add("CellTypeInteractionZScore");
//            covariatesToTest.add("MainEffectZScore");
            boolean[] platformIsHT12v3 = new boolean[]{
                true,
                true,
                false,};

            m.run(output, probetranslationfile, datasetFileDirs, datasetNames, platformIsHT12v3, false, flipEffects, covariatesToTest, false);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String output, String probetranslationfile, String[] datasetFileDirs, String[] datasetNames, boolean[] platformIsHT12v3, boolean allowNaNs, boolean flipeffects, HashSet<String> covariatesToTest, boolean forcePresenceOfEQTL) throws IOException {

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

        // load the datasets..
        CellTypeSpecificMetaAnalsysisDatasetContainer container;
        if (covariatesToTest == null) {
            container = new CellTypeSpecificMetaAnalsysisDatasetContainer(datasetFileDirs, platformIsHT12v3, ht12v4ToHT12v3);
        } else {
            container = new CellTypeSpecificMetaAnalsysisDatasetContainer(datasetFileDirs, platformIsHT12v3, ht12v4ToHT12v3, covariatesToTest);
        }

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
        DoubleMatrixDataset<String, String> outputMatrix = new DoubleMatrixDataset<String, String>();
        DoubleMatrixDataset<String, String> outputMatrixSampleSize = new DoubleMatrixDataset<String, String>();

        double[][] metaAnalysisZScores = new double[allCovariates.size()][allEQTL.size()];
        double[][] metaAnalysisZScoresSampleSizes = new double[allCovariates.size()][allEQTL.size()];

        TextFile snpOutputFile = new TextFile(output + "MetaAnalysisSNPSummary.txt", TextFile.W);

        // snp + "\t" + probe + "\t" + dsNames + "\t" + mafStr + "\t" + crStr + "\t" + hwepStr + "\t" + nrCalledStr + "\t" + firstSNPAlleles + "\t" + firstSNPAssessedAllele
        snpOutputFile.writeln("SNP\tProbe\tDatasets\tMAF\tCR\tHWEP\tNrCalled\tAlleles\tAlleleAssessed");

        TextFile bonferroniSig = new TextFile(output + "BonferroniSignificant.txt", TextFile.W);
        // eQTL + "\t" + allCovariateArr[c] + "\t" + ht12v4ToHugo.get(allCovariateArr[c]) + "\t" + metaZ + "\t" + P + "\t" + totalSample

        double bonferroni = 0.05 / (allCovariateArr.length * allEQTLArr.length);
        System.out.println("Bonferroni threshold: " + bonferroni + "\t(" + allCovariateArr.length + " x " + allEQTLArr.length + ")");
        bonferroniSig.writeln("Bonferroni threshold: " + bonferroni + "\t(" + allCovariateArr.length + " x " + allEQTLArr.length + ")");
        bonferroniSig.writeln("eQTL\tCovariate\tHUGO\tMetaZ\tPVal\tSampleSize");

        ArrayList<String> colNames = new ArrayList<String>();
        int nrTests = 0;
        int nrNaN = 0;
        ProgressBar pb = new ProgressBar(allEQTLArr.length, "Running meta");
        boolean[] hasNaN = new boolean[allEQTLArr.length];
        boolean NaNsPresent = false;

        for (int e = 0; e < allEQTLArr.length; e++) {

            String eQTL = allEQTLArr[e];
            String snp = eQTL.split("-")[0];
            String probe = eQTL.split("-")[1];
            Integer snpId = snpToId.get(snp);

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

            snpOutputFile.writeln(snp + "\t" + probe + "\t" + dsNames + "\t" + mafStr + "\t" + crStr + "\t" + hwepStr + "\t" + nrCalledStr + "\t" + firstSNPAlleles + "\t" + firstSNPAssessedAllele);

            colNames.add(eQTL);

            for (int c = 0; c < allCovariateArr.length; c++) {
                double[] zScores = new double[datasets.size()];
                int totalSample = 0;
                int[] sampleSize = new int[datasets.size()];
                int nrDatasetsWithNan = 0;
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
                            nrDatasetsWithNan++;
                            zScores[d] = Double.NaN;
                            sampleSize[d] = 0;
                        }
                    } else {
                        nrDatasetsWithNan++;
                        zScores[d] = Double.NaN;
                        sampleSize[d] = 0;
                    }
                }

                double metaZ = Double.NaN;
                if (forcePresenceOfEQTL) {
                    if (nrDatasetsWithNan == 0) {
                        metaZ = ZScores.getWeightedZ(zScores, sampleSize);
                    }
                } else {
                    metaZ = ZScores.getWeightedZ(zScores, sampleSize);
                }

                nrTests++;

                if (Double.isNaN(metaZ)) {
                    hasNaN[e] = true;
                    NaNsPresent = true;
//                    System.out.println(Strings.concat(zScores, Strings.semicolon) + "\t" + Strings.concat(sampleSize, Strings.semicolon) + "\tresults in NaN for " + allEQTLArr[e] + "\t" + allCovariateArr[c]);
                    for (int d = 0; d < datasets.size(); d++) {
                        Integer snpIdInDataset = snpIndex[d][snpId];
                        Integer eQTLIdInDataset = eQTLIndex[d][e];
                        Integer covariateIdInDataset = covariateIndex[d][c];
//                        System.out.println(d + "\t" + snpIdInDataset + "\t" + eQTLIdInDataset + "\t" + covariateIdInDataset);
                    }
                    nrNaN++;
                } else {
                    metaAnalysisZScores[c][e] = metaZ;
                    metaAnalysisZScoresSampleSizes[c][e] = (double) totalSample;
                    double P = ZScores.zToP(metaZ);
                    if (P < bonferroni) {
                        bonferroniSig.writeln(eQTL + "\t" + allCovariateArr[c] + "\t" + ht12v3ToHugo.get(allCovariateArr[c]) + "\t" + metaZ + "\t" + Strings.concat(zScores, Strings.semicolon) + "\t" + Strings.concat(sampleSize, Strings.semicolon) + "\t" + P + "\t" + totalSample);
                    }
                }
            }
//            pb.iterate();
        }
        pb.close();

        System.out.println(nrNaN + " NaN out of: " + allEQTLArr.length);
        System.out.println(nrTests + " performed / " + (allCovariateArr.length * allEQTLArr.length) + " expected");

        outputMatrix.rawData = metaAnalysisZScores;
        outputMatrix.colObjects = colNames;
        outputMatrix.rowObjects = Arrays.asList(allCovariateArr);
        outputMatrix.recalculateHashMaps();

        outputMatrixSampleSize.rawData = metaAnalysisZScoresSampleSizes;
        outputMatrixSampleSize.colObjects = colNames;
        outputMatrixSampleSize.rowObjects = Arrays.asList(allCovariateArr);
        outputMatrixSampleSize.recalculateHashMaps();

        if (NaNsPresent && !allowNaNs) {
            // remove those stupid NaN values...
            ArrayList<String> newColNames = new ArrayList<String>();
            int nanctr = 0;
            for (int c = 0; c < outputMatrix.nrCols; c++) {
                if (!hasNaN[c]) {
                    nanctr++;
                    newColNames.add(outputMatrix.colObjects.get(c));
                }
            }
            System.out.println(nanctr + " eQTLs without NaN's. Will remove them now.");

            double[][] newRawData = new double[outputMatrix.nrRows][newColNames.size()];
            double[][] newRawDataSampleSize = new double[outputMatrix.nrRows][newColNames.size()];
            for (int r = 0; r < outputMatrix.nrRows; r++) {
                int colCtr = 0;
                for (int c = 0; c < outputMatrix.nrCols; c++) {
                    if (!hasNaN[c]) {
                        newRawData[r][colCtr] = outputMatrix.rawData[r][c];
                        newRawDataSampleSize[r][colCtr] = outputMatrixSampleSize.rawData[r][c];
                        colCtr++;
                    }
                }
            }

            // replace the previous data
            outputMatrix.rawData = newRawData;
            outputMatrix.colObjects = newColNames;
            outputMatrix.recalculateHashMaps();

            outputMatrixSampleSize.rawData = newRawDataSampleSize;
            outputMatrixSampleSize.colObjects = newColNames;
            outputMatrixSampleSize.recalculateHashMaps();

        }

        // flip the Z-scores if main effect Z-score < 0
        Integer mainEffectZScoreRow = outputMatrix.hashRows.get("MainEffectZScore");
        if (mainEffectZScoreRow == null) {
            System.err.println("Apparently there is no main effect zscore. " + mainEffectZScoreRow);
        } else {
            for (int col = 0; col < outputMatrix.nrCols; col++) {
                double maineffectZscore = outputMatrix.rawData[mainEffectZScoreRow][col];
                for (int row = 0; row < outputMatrix.nrRows; row++) {
                    if (!mainEffectZScoreRow.equals(row)) {
                        if (maineffectZscore < 0 && flipeffects) {
                            outputMatrix.rawData[row][col] *= -1;
                        }
                    }
                }
            }
        }
        // check whether the new matrix is NaN free...
        System.out.println("Checking whether NaN's are still present...");
        CheckForNaNs.check(outputMatrix);

        // write proxy interaction Z-score
        Integer id = outputMatrix.hashRows.get("CellTypeInteractionZScore");
        if (id != null) {
            TextFile out = new TextFile(output + "Vector-CellTypeInteractionZScore.txt", TextFile.W);
            out.writeln("eQTL\tCellTypeInteractionZScore");

            for (int col = 0; col < outputMatrix.nrCols; col++) {
                out.writeln(outputMatrix.colObjects.get(col) + "\t" + outputMatrix.rawData[id][col]);
            }
            out.close();
        }

        id = outputMatrix.hashRows.get("MainEffectZScore");
        if (id != null) {
            TextFile out = new TextFile(output + "Vector-MainEffectZScore.txt", TextFile.W);
            out.writeln("eQTL\tMainEffectZScore");
            for (int col = 0; col < outputMatrix.nrCols; col++) {
                out.writeln(outputMatrix.colObjects.get(col) + "\t" + outputMatrix.rawData[id][col]);
            }
            out.close();
        }

        outputMatrix.save(output + "MetaAnalysisZScoreMatrix.txt");
        outputMatrixSampleSize.save(output + "MetaAnalysisZScoreMatrixSampleSize.txt");

        snpOutputFile.close();
        bonferroniSig.close();
    }
}
