/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
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
    public static void main(String[] args) {
        // TODO code application logic here

        try {
            Meta m = new Meta();

            String output = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/MetaAnalysisWithoutEGCUT/";
            Gpio.createDir(output);
            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            String[] datasetFileDirs = new String[]{
//                "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/EGCUT/",
                "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/Groningen/",
                "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-12-12KEQTLs/Rotterdam/"
            };

            String[] datasetNames = new String[]{
//                "EGCUT",
                "BloodHT12", "Rotterdam"
            };

//            boolean[] platformIsHT12v3 = new boolean[]{true, true, false};

            boolean[] platformIsHT12v3 = new boolean[]{
//                true, 
                true, false};


            m.run(output, probetranslationfile, datasetFileDirs, datasetNames, platformIsHT12v3, false);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String output, String probetranslationfile, String[] datasetFileDirs, String[] datasetNames, boolean[] platformIsHT12v3, boolean allowNaNs) throws IOException {

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
        SNP[][] snps = new SNP[datasetFileDirs.length][0]; // JAGGED ARRAY FOR SNP STORAGE
        ArrayList<DoubleMatrixDataset<String, String>> datasets = new ArrayList<DoubleMatrixDataset<String, String>>();
        HashSet<String> allSNPs = new HashSet<String>();
        HashSet<String> allEQTL = new HashSet<String>();
        for (int dsN = 0; dsN < datasetFileDirs.length; dsN++) {
            String dir = datasetFileDirs[dsN];
            // load the dataset....
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(dir + "CellTypeSpecificityMatrix.binary");
            datasets.add(ds);

            List<String> eqtls = ds.colObjects;

            for (int i = 0; i < eqtls.size(); i++) {
                String eQTL = eqtls.get(i);
                if (platformIsHT12v3[dsN]) {
                    allEQTL.add(eQTL);
                } else {
                    String[] eqtlelems = eQTL.split("-");
                    eqtlelems[1] = ht12v4ToHT12v3.get(eqtlelems[1]);
                    allEQTL.add(Strings.concat(eqtlelems, Strings.dash));
                }
            }

            TextFile tfIn = new TextFile(dir + "SNPSummaryStatistics.txt", TextFile.R);
            tfIn.readLine(); // skip the header..
            ArrayList<SNP> allSNPsInDs = new ArrayList<SNP>();
            String[] elems = tfIn.readLineElems(TextFile.tab);
            while (elems != null) {
                // SNP	Chr	ChrPos	Alleles	MinorAllele	MAF	CallRate	HWE	GenotypesCalled
                String snp = elems[0];
                String alleles = elems[3];
                String alleleAssessed = elems[4];
                String maf = elems[5];
                String cr = elems[6];
                String hwe = elems[7];
                String nrCalled = elems[8];

                SNP s = new SNP();
                s.setName(snp);
                byte[] allelesB = new byte[2];
                String[] alleleElems = alleles.split("/");
                allelesB[0] = BaseAnnot.toByte(alleleElems[0]);
                allelesB[1] = BaseAnnot.toByte(alleleElems[1]);
                s.setAlleleCodes(allelesB);
                byte alleleAssessedB = BaseAnnot.toByte(alleleAssessed);
                s.setMinorAllele(alleleAssessedB);

                s.setMAF(Double.parseDouble(maf));
                s.setCR(Double.parseDouble(cr));
                s.setHWEP(Double.parseDouble(hwe));

                s.setNrCalled(Integer.parseInt(nrCalled));
                allSNPs.add(snp);
                allSNPsInDs.add(s);
                elems = tfIn.readLineElems(TextFile.tab);
            }
            tfIn.close();
            snps[dsN] = allSNPsInDs.toArray(new SNP[0]);
        }

        // hash the SNPs for easy lookup..
        Integer[][] snpIndex = new Integer[datasets.size()][allSNPs.size()]; // for each SNP in the meta analysis, this index points to a SNP in a dataset (in the jagged allSNPs array)
        HashMap<String, Integer> snpToId = new HashMap<String, Integer>();
        int ctr = 0;
        for (String snp : allSNPs) {
            snpToId.put(snp, ctr);
            ctr++;
        }
        System.out.println(snpToId.size() + " snps loaded");
        for (int d = 0; d < datasets.size(); d++) {
            SNP[] snpsInDs = snps[d];
            for (int s = 0; s < snpsInDs.length; s++) {
                Integer id = snpToId.get(snps[d][s].getName());
                snpIndex[d][id] = s;
            }
        }

        // hash the eQTLs as well..
        Integer[][] eQTLIndex = new Integer[datasets.size()][allEQTL.size()];
        HashMap<String, Integer> eQTLToId = new HashMap<String, Integer>();
        String[] allEQTLArr = allEQTL.toArray(new String[0]);
        for (int e = 0; e < allEQTLArr.length; e++) {
            String eQTL = allEQTLArr[e];
            eQTLToId.put(eQTL, e);
        }
        for (int d = 0; d < datasets.size(); d++) {
            List<String> eqtls = datasets.get(d).colObjects;
            for (int s = 0; s < eqtls.size(); s++) {

                String eQTL = eqtls.get(s);

                Integer id = null;
                if (platformIsHT12v3[d]) {
                    id = eQTLToId.get(eQTL);
                } else {
                    String[] eqtlelems = eQTL.split("-");
                    eqtlelems[1] = ht12v4ToHT12v3.get(eqtlelems[1]);
                    eQTL = Strings.concat(eqtlelems, Strings.dash);
                    id = eQTLToId.get(eQTL);
                }
                if (id == null) {
                    System.err.println("ERROR: " + id + " is null for " + eQTL + " in dataset " + datasetNames[d]);
                } else {
                    eQTLIndex[d][id] = s;
                }
            }
        }


        /*
         * rowNames.add("CellTypeSNPZScore");
         rowNames.add("CellTypeZScore");
         rowNames.add("CellTypeInteractionZScore");
         rowNames.add("MainEffectZScore");
         */

        // hash the probes / covariates...
        HashSet<String> allCovariates = new HashSet<String>();
        for (int d = 0; d < datasets.size(); d++) {
            DoubleMatrixDataset<String, String> ds = datasets.get(d);
            for (int r = 0; r < ds.nrRows; r++) {

                // each matrix has two 'special' covariates
                if (ds.rowObjects.get(r).equals("CellTypeInteractionZScore") || ds.rowObjects.get(r).equals("MainEffectZScore") || ds.rowObjects.get(r).equals("CellTypeZScore") || ds.rowObjects.get(r).equals("CellTypeSNPZScore")) {
                    allCovariates.add(ds.rowObjects.get(r));
                } else {
                    // parse the 'normal' covariates
                    if (!platformIsHT12v3[d]) {
                        String cov = ht12v4ToHT12v3.get(ds.rowObjects.get(r));
                        if (cov == null) {
                            System.err.println("ERROR: HT12v4 probe not included in translation file? " + cov + "\t" + ds.rowObjects.get(r));
                            System.exit(-1);
                        }
                        if (cov.equals("-")) {
                            System.err.println("ERROR: covariate equals '-' " + ds.rowObjects.get(r));
                        } else {
                            allCovariates.add(cov);
                        }
                    } else {
                        String cov = ds.rowObjects.get(r);
                        if (cov.equals("-")) {
                            System.err.println("ERROR: covariate equals '-' " + ds.rowObjects.get(r));
                        } else {
                            allCovariates.add(cov);
                        }
                    }
                }
            }
        }


        String[] allCovariateArr = allCovariates.toArray(new String[0]);
        Integer[][] covariateIndex = new Integer[datasets.size()][allCovariates.size()];
        HashMap<String, Integer> covariateToId = new HashMap<String, Integer>();
        for (int c = 0; c < allCovariateArr.length; c++) {
            covariateToId.put(allCovariateArr[c], c);
        }
        for (int d = 0; d < datasets.size(); d++) {
            DoubleMatrixDataset<String, String> ds = datasets.get(d);
            for (int r = 0; r < ds.nrRows; r++) {
                String cov = null;

                if (ds.rowObjects.get(r).equals("CellTypeInteractionZScore") || ds.rowObjects.get(r).equals("MainEffectZScore") || ds.rowObjects.get(r).equals("CellTypeZScore") || ds.rowObjects.get(r).equals("CellTypeSNPZScore")) {
                    cov = ds.rowObjects.get(r);
                } else {
                    if (!platformIsHT12v3[d]) {
                        cov = ht12v4ToHT12v3.get(ds.rowObjects.get(r));
                        if (cov == null) {
                            System.out.println("Covatiate not found: " + ds.rowObjects.get(r));
                            System.exit(-1);
                        }
                    } else {
                        cov = ds.rowObjects.get(r);
                    }
                }


                Integer id = covariateToId.get(cov);
                covariateIndex[d][id] = r;
            }
        }




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

                nrTests++;

                double metaZ = ZScores.getWeightedZ(zScores, sampleSize);
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

            // check whether the new matrix is NaN free...
            System.out.println("Checking whether NaN's are still present...");
            CheckForNaNs.check(outputMatrix);

        }

        outputMatrix.save(output + "MetaAnalysisZScoreMatrix.txt");
        outputMatrixSampleSize.save(output + "MetaAnalysisZScoreMatrixSampleSize.txt");

        snpOutputFile.close();
        bonferroniSig.close();
    }
}
