/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.distribution.BinomialDistribution;
import umcg.genetica.containers.Pair;
import umcg.genetica.gwas.Independifier;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class BinomialDiseaseTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String metaAnalysisMatrix = "/Volumes/iSnackHD/AeroFS/2013-11-20-MetaAnalysisForPlots/MetaAnalysisZScoreMatrix.txt"; // /Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/MetaAnalysisZScoreMatrix.binary";
        String metaAnalysisFDRMatrix = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-09-27-FDREstimates/MetaAnalysisZScoreMatrix.binary-FDR.binary";
        String gwasCatalog = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-09-23-gwascatalog-WithAllIBDLoci.txt";
        double gwasSignificanceCutoffPValue = 5E-8;
        int minNrEQTLSNPs = 10;
        int minNrEQTLSNPsCellTypeSpecific = 1;
        double cellTypeSpecificityCutoffZScore = 0;
        String genotypeReferenceDataset = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap2r24-CEU/";
        double ldThreshold = 0.2;
        boolean buildBGDistributionUsingGWASSNPsOnly = false;
        // double fdrThreshold = 1.66E-4; // bonferroni
        double fdrThreshold = 0.05;
//        String selectCovariate = "CellTypeInteractionZScore";

        HashSet<String> selectCovariate = new HashSet<String>();
        selectCovariate.add("CellTypeInteractionZScore");
        BinomialDiseaseTest test = new BinomialDiseaseTest();
        String probeTranslationFile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-18-ProbeAnnotationFile.txt";

        boolean alsoPrintAffectedGenes = true;
        boolean useOnlyTopEffectPerProbe = false;
        // String[] traits = new String[]{"platelet", "corpuscular", "crohn", "neutrophil"};
        String outdir = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/20114-01-15-BinomialDiseaseTestFDR" + fdrThreshold + "/";
        String[] traits = null;
        try {
//        
//            DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaAnalysisMatrix);
//            metaMatrix.save("/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/MetaAnalysisZScoreMatrix.binary");
//            

            Gpio.createDir(outdir);

            String outfile = outdir + "FisherExactOutput-Pruning" + ldThreshold + "-BgBuildWithGWASSNPsOnly" + buildBGDistributionUsingGWASSNPsOnly + ".txt";

            if (selectCovariate != null) {
                outfile = outdir + "FisherExactOutput-" + selectCovariate + "-Pruning" + ldThreshold + "-BgBuildWithGWASSNPsOnly" + buildBGDistributionUsingGWASSNPsOnly + ".txt";
            }

            test.runOriginal(metaAnalysisMatrix, metaAnalysisFDRMatrix, fdrThreshold,
                    gwasCatalog, gwasSignificanceCutoffPValue, minNrEQTLSNPs,
                    minNrEQTLSNPsCellTypeSpecific, cellTypeSpecificityCutoffZScore,
                    genotypeReferenceDataset, ldThreshold, buildBGDistributionUsingGWASSNPsOnly,
                    selectCovariate, probeTranslationFile, outfile, traits, alsoPrintAffectedGenes, useOnlyTopEffectPerProbe);
//            test.runSecondVersion(metaAnalysisMatrix, metaAnalysisFDRMatrix, cellTypeSpecificityCutoffZScore, fdrThreshold, probeTranslationFile, gwasCatalog, gwasSignificanceCutoffPValue, genotypeReferenceDataset, ldThreshold, minNrEQTLSNPs, selectCovariate, outfile, traits);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void runOriginal(String metaAnalysisMatrix, String metaAnalysisFDRMatrix, double fdrThreshold, String gwasCatalog, double gwasSignificanceCutoffPValue, int minNrEQTLSNPs, int minNrEQTLSNPsCellTypeSpecific,
            double cellTypeSpecificityCutoffZScore, String genotypeReferenceDataset, double ldthreshold, boolean buildBGDistributionUsingGWASSNPsOnly, HashSet<String> selectCovariate, String probeTranslation, String outfileName,
            String[] testSpecificTraits, boolean alsoPrintAffectedenes, boolean useOnlyTopEffectPerProbe) throws IOException {

        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> probeToHUGO = pb.getProbeTranslation(probeTranslation, "HT12v3.txt", "Gene");

        GWASCatalog catalog = new GWASCatalog(gwasCatalog, gwasSignificanceCutoffPValue);
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaAnalysisMatrix, selectCovariate);//new DoubleMatrixDataset<String, String>(metaAnalysisMatrix, selectCovariate, null);
        DoubleMatrixDataset<String, String> metaFDRMatrix = new DoubleMatrixDataset<String, String>(metaAnalysisFDRMatrix, selectCovariate);

        // iterate the covariates..
        GWASTrait[] allTraits = catalog.getTraits();
        HashSet<String> allDiseaseSNPs = new HashSet<String>();
        HashMap<String, HashSet<String>> allSNPsPerTrait = new HashMap<String, HashSet<String>>();
        for (GWASTrait trait : allTraits) {
            GWASSNP[] snps = trait.getSNPs(gwasSignificanceCutoffPValue);
            HashSet<String> snpsForTrait = allSNPsPerTrait.get(trait.getName());
            if(snpsForTrait == null){
                snpsForTrait = new HashSet<String>();
            }
            for (GWASSNP snp : snps) {
                String name = snp.getName();
                allDiseaseSNPs.add(name);
                snpsForTrait.add(name);
            }
            allSNPsPerTrait.put(trait.getName(), snpsForTrait);
        }

        HashSet<String> allDiseaseSNPsHavingCisEQTLs = new HashSet<String>();

        for (int eqtl = 0; eqtl < metaMatrix.nrCols; eqtl++) {
            String eqtlName = metaMatrix.colObjects.get(eqtl);
            String snp = eqtlName.split("-")[0];

            if (allDiseaseSNPs.contains(snp)) {
                allDiseaseSNPsHavingCisEQTLs.add(snp);
            }
        }

        /*
         * covariateName + "\t" + probeToHUGO.get(covariateName) + "\t" + eQTLSNPsBeyondThreshold.size() + "\t" + allEQTLSNPs.size() + "\t" + proportionCellTypeSpecificOverAll + "\t" + trait.getName()
         + "\t" + snps.length + "\t" + nrTotalRemainingTraitSNPsAfterPruning + "\t" + nrTraitSNPsThatAreEQTL + "\t" + nrTraitSNPsCellTypeSpecific + "\t" + oneTailedP;
         */
        HashMap<Integer, HashSet<String>> cellTypeSpecificEQTLSNPsPerCovariate = new HashMap<Integer, HashSet<String>>();
        HashMap<Integer, HashSet<String>> cellTypeSpecificEQTLPerCovariate = new HashMap<Integer, HashSet<String>>();
        HashSet<String> allEQTLSNPs = new HashSet<String>();
        for (int covariate = 0; covariate < metaMatrix.nrRows; covariate++) {
            HashSet<String> eQTLSNPsBeyondThreshold = new HashSet<String>();
            HashSet<String> eQTLsBeyondThreshold = new HashSet<String>();

            if (useOnlyTopEffectPerProbe) {
                // find the top effect per probe

                for (int eqtl = 0; eqtl < metaMatrix.nrCols; eqtl++) {
                    double z = metaMatrix.rawData[covariate][eqtl];
                    String eqtlname = metaMatrix.colObjects.get(eqtl);
                    String eqtlsnp = eqtlname.split("-")[0];
                    String eqtlprobe = eqtlname.split("-")[1];
                    double absZ = Math.abs(z);

                }
            }

            for (int eqtl = 0; eqtl < metaMatrix.nrCols; eqtl++) {
                double z = metaMatrix.rawData[covariate][eqtl];
                String eqtlname = metaMatrix.colObjects.get(eqtl);
                String eqtlsnp = eqtlname.split("-")[0];
                String eqtlprobe = eqtlname.split("-")[1];
                if (!buildBGDistributionUsingGWASSNPsOnly || (buildBGDistributionUsingGWASSNPsOnly && allDiseaseSNPsHavingCisEQTLs.contains(eqtlsnp))) {
                    allEQTLSNPs.add(eqtlsnp);
                    double fdr = metaFDRMatrix.rawData[covariate][eqtl];
                    if (fdr < fdrThreshold && z > cellTypeSpecificityCutoffZScore) {
                        eQTLSNPsBeyondThreshold.add(eqtlsnp);
                        eQTLsBeyondThreshold.add(eqtlname);
                    }
                }
            }
            System.out.println(eQTLsBeyondThreshold.size() + " eQTLs beyond threshold for " + metaMatrix.rowObjects.get(covariate));
            cellTypeSpecificEQTLSNPsPerCovariate.put(covariate, eQTLSNPsBeyondThreshold);
            cellTypeSpecificEQTLPerCovariate.put(covariate, eQTLsBeyondThreshold);
        }

        TriTyperGenotypeData genotypeData = new TriTyperGenotypeData(genotypeReferenceDataset);
        SNPLoader loader = genotypeData.createSNPLoader();
        DetermineLD ld = new DetermineLD();
        TextFile outfile = new TextFile(outfileName, TextFile.W);

        String outputheader = "Covariate\tGene\tNrBackGroundCellTypeSpecificeQTLSNPs\tNrEQTLsSNPsTotal\tRatioBackGround\tDisease\t"
                + "NrTraitSNPs\tNrTraitSNPsAfterPrunint\tNrTraitSNPsThatAreEQTL\tNrTraitSNPsThatAreCellTypeSpecificEQTL\tRatio\tOneTailedP\tFisherExactTestP\tFisherExactTestP2";
        if (alsoPrintAffectedenes) {
            outputheader += "\tCellTypeSpecificSNPs\tCellTypeSpecificEQTLs\tCellTypeSpecificEQTLGenes\tPrunedSNPs";
        }
        outfile.writeln(outputheader);

        for (GWASTrait trait : allTraits) {

            boolean testTrait = false;
            if (testSpecificTraits == null) {
                testTrait = true;
            } else {
                String traitName = trait.getName().toLowerCase();
                for (String traitElems : testSpecificTraits) {
                    if (traitName.contains(traitElems)) {
                        testTrait = true;
                    }
                }
            }

            if (testTrait) {

                GWASSNP[] traitSNPs = trait.getSNPs(gwasSignificanceCutoffPValue);
                HashSet<String> uniqueDiseaseSNPs = new HashSet<String>();
                for (GWASSNP snp : traitSNPs) {
                    String snpName = snp.getName();
                    if (allDiseaseSNPsHavingCisEQTLs.contains(snpName)) {
                        uniqueDiseaseSNPs.add(snpName);
                    }
                }

                // prune the SNPs
                String[] uniqueDiseaseSNPArray = uniqueDiseaseSNPs.toArray(new String[0]);
                boolean[] snpExclude = new boolean[uniqueDiseaseSNPArray.length];
                for (int s = 0; s < uniqueDiseaseSNPArray.length; s++) {
                    if (ldthreshold <= 1) {
                        if (!snpExclude[s]) {
                            Integer snpID1 = genotypeData.getSnpToSNPId().get(uniqueDiseaseSNPArray[s]);
                            if (snpID1 == null) {
                                snpExclude[s] = true;
                            } else {
                                SNP snpObj1 = genotypeData.getSNPObject(snpID1);
                                for (int t = s + 1; t < uniqueDiseaseSNPArray.length; t++) {
                                    loader.loadGenotypes(snpObj1);
                                    if (!snpExclude[t]) {
                                        Integer snpID2 = genotypeData.getSnpToSNPId().get(uniqueDiseaseSNPArray[t]);
                                        if (snpID2 == null) {
                                            snpExclude[s] = true;
                                        } else {
                                            SNP snpObj2 = genotypeData.getSNPObject(snpID2);
                                            if (snpObj1.getChr() == snpObj2.getChr() && snpObj1.getChr() > 0) {
                                                loader.loadGenotypes(snpObj2);
                                                double rSquared = ld.getRSquared(snpObj1, snpObj2, genotypeData, ld.RETURN_R_SQUARED, ld.INCLUDE_CASES_AND_CONTROLS, false);
                                                if (rSquared > ldthreshold) {
                                                    snpExclude[t] = true;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                int nrTraitSNPsThatAreEQTL = 0;
                ArrayList<String> snpsAfterPruning = new ArrayList<String>();
                for (int s = 0; s < uniqueDiseaseSNPArray.length; s++) {
                    if (!snpExclude[s]) {
                        String rsName = uniqueDiseaseSNPArray[s];
                        if (allEQTLSNPs.contains(rsName)) {
                            nrTraitSNPsThatAreEQTL++;
                            snpsAfterPruning.add(rsName);
                        }
                    }
                }

                if (nrTraitSNPsThatAreEQTL >= minNrEQTLSNPs) {
                    System.out.println("Testing: " + trait.getName() + " nrEQTL SNPs: " + nrTraitSNPsThatAreEQTL);
                    for (int covariate = 0; covariate < metaMatrix.nrRows; covariate++) {
                        String covariateName = metaMatrix.rowObjects.get(covariate);

                        HashSet<String> eQTLSNPsBeyondThreshold = cellTypeSpecificEQTLSNPsPerCovariate.get(covariate);
                        HashSet<String> eQTLsBeyondThreshold = cellTypeSpecificEQTLPerCovariate.get(covariate);

                        if (useOnlyTopEffectPerProbe) {

                            HashSet<String> probesBeyondThreshold = new HashSet<String>();
                            HashSet<String> allEQTLProbes = new HashSet<String>();
                            for (int eqtl = 0; eqtl < metaMatrix.nrCols; eqtl++) {
                                String eqtlname = metaMatrix.colObjects.get(eqtl);
                                String probe = eqtlname.split("-")[1];
                                if (eQTLsBeyondThreshold.contains(eqtlname)) {
                                    probesBeyondThreshold.add(probe);
                                }
                                allEQTLProbes.add(probe);
                            }

                            double proportionCellTypeSpecificOverAll = (double) probesBeyondThreshold.size() / allEQTLProbes.size(); // this is not pruned.. because..?

                            HashSet<String> allTraitProbesAfterPruningBeyondThreshold = new HashSet<String>();
                            HashSet<String> allTraitProbesAfterPruning = new HashSet<String>();
                            HashSet<String> allTraitProbes = new HashSet<String>();

                            // get pruned snps
                            for (int s = 0; s < uniqueDiseaseSNPArray.length; s++) {
                                String rsName = uniqueDiseaseSNPArray[s];
                                for (int eqtl = 0; eqtl < metaMatrix.nrCols; eqtl++) {
                                    String[] eqtlelems = metaMatrix.colObjects.get(eqtl).split("-");

                                    if (eqtlelems[0].equals(rsName)) {
                                        allTraitProbes.add(eqtlelems[1]);
                                        if (!snpExclude[s]) {
                                            allTraitProbesAfterPruning.add(eqtlelems[1]);
                                            if (probesBeyondThreshold.contains(eqtlelems[1])) {
                                                allTraitProbesAfterPruningBeyondThreshold.add(eqtlelems[1]);
                                            }
                                        }
                                    }
                                }
                            }

                            BinomialDistribution dist = new BinomialDistribution(allTraitProbesAfterPruning.size(), proportionCellTypeSpecificOverAll);
//                    if (nrTraitSNPsThatAreEQTL >= minNrEQTLSNPs) { // && nrTraitSNPsCellTypeSpecific >= minNrEQTLSNPsCellTypeSpecific) {
                            double oneTailedP = 0;

                            for (int a = allTraitProbesAfterPruningBeyondThreshold.size(); a <= allTraitProbesAfterPruning.size(); a++) {
                                oneTailedP += dist.probability(a);
                            }
                            FisherExactTest fet = new FisherExactTest();
                            double fetp = fet.getFisherPValue(
                                    allTraitProbesAfterPruningBeyondThreshold.size(),
                                    allTraitProbesAfterPruning.size() - allTraitProbesAfterPruningBeyondThreshold.size(),
                                    probesBeyondThreshold.size(),
                                    allEQTLProbes.size() - probesBeyondThreshold.size());

//                            double fetp2 = 1d;fet.getFisherPValue(nrTraitSNPsCellTypeSpecific,
//                                    nrTraitSNPsThatAreEQTL - nrTraitSNPsCellTypeSpecific,
//                                    eQTLSNPsBeyondThreshold.size() - nrTraitSNPsCellTypeSpecific,
//                                    allEQTLSNPs.size() - eQTLSNPsBeyondThreshold.size() - nrTraitSNPsCellTypeSpecific - (nrTotalRemainingTraitSNPsAfterPruning - nrTraitSNPsCellTypeSpecific));
                            double ratio = (double) allTraitProbesAfterPruningBeyondThreshold.size() / (allTraitProbesAfterPruning.size() - allTraitProbesAfterPruningBeyondThreshold.size());
                            String outputStr
                                    = covariateName + "\t"
                                    + probeToHUGO.get(covariateName) + "\t"
                                    + probesBeyondThreshold.size() + "\t"
                                    + allEQTLProbes.size() + "\t"
                                    + proportionCellTypeSpecificOverAll + "\t"
                                    + trait.getName() + "\t"
                                    + allTraitProbes.size() + "\t"
                                    + allTraitProbesAfterPruning.size() + "\t"
                                    + allTraitProbesAfterPruningBeyondThreshold.size() + "\t"
                                    + allTraitProbesAfterPruningBeyondThreshold.size() + "\t"
                                    + ratio + "\t" + oneTailedP + "\t" + fetp + "\t" + 1d;
                            outfile.writeln(outputStr);
                            System.out.println(outputStr);
                        } else {
                            double proportionCellTypeSpecificOverAll = (double) eQTLSNPsBeyondThreshold.size() / allEQTLSNPs.size(); // this is not pruned.. because..?
                            int nrTotalRemainingTraitSNPsAfterPruning = 0;

                            int nrTraitSNPsCellTypeSpecific = 0;
                            HashSet<String> traitSNPsWithSignificantInteractionTerms = new HashSet<String>();
                            for (int s = 0; s < uniqueDiseaseSNPArray.length; s++) {
                                if (!snpExclude[s]) {
                                    String rsName = uniqueDiseaseSNPArray[s];
                                    if (eQTLSNPsBeyondThreshold.contains(rsName)) {
                                        nrTraitSNPsCellTypeSpecific++;
                                        traitSNPsWithSignificantInteractionTerms.add(rsName);
                                    }
                                    nrTotalRemainingTraitSNPsAfterPruning++;
                                }
                            }

                            BinomialDistribution dist = new BinomialDistribution(nrTraitSNPsThatAreEQTL, proportionCellTypeSpecificOverAll);
//                    if (nrTraitSNPsThatAreEQTL >= minNrEQTLSNPs) { // && nrTraitSNPsCellTypeSpecific >= minNrEQTLSNPsCellTypeSpecific) {
                            double oneTailedP = 0;

                            for (int a = nrTraitSNPsCellTypeSpecific; a <= nrTraitSNPsThatAreEQTL; a++) {
                                oneTailedP += dist.probability(a);
                            }
                            FisherExactTest fet = new FisherExactTest();
                            double fetp = fet.getFisherPValue(nrTraitSNPsCellTypeSpecific,
                                    nrTraitSNPsThatAreEQTL - nrTraitSNPsCellTypeSpecific,
                                    eQTLSNPsBeyondThreshold.size(),
                                    allEQTLSNPs.size() - eQTLSNPsBeyondThreshold.size());

                            double fetp2 = fet.getFisherPValue(nrTraitSNPsCellTypeSpecific,
                                    nrTraitSNPsThatAreEQTL - nrTraitSNPsCellTypeSpecific,
                                    eQTLSNPsBeyondThreshold.size() - nrTraitSNPsCellTypeSpecific,
                                    allEQTLSNPs.size() - eQTLSNPsBeyondThreshold.size() - nrTraitSNPsCellTypeSpecific - (nrTotalRemainingTraitSNPsAfterPruning - nrTraitSNPsCellTypeSpecific));

                            double ratio = (double) nrTraitSNPsCellTypeSpecific / (nrTraitSNPsThatAreEQTL - nrTraitSNPsCellTypeSpecific);
                            String outputStr = covariateName 
                                    + "\t" + probeToHUGO.get(covariateName) 
                                    + "\t" + eQTLSNPsBeyondThreshold.size() 
                                    + "\t" + allEQTLSNPs.size() 
                                    + "\t" + proportionCellTypeSpecificOverAll 
                                    + "\t" + trait.getName()
                                    + "\t" + allSNPsPerTrait.get(trait.getName()).size() 
                                    + "\t" + nrTotalRemainingTraitSNPsAfterPruning 
                                    + "\t" + nrTraitSNPsThatAreEQTL 
                                    + "\t" + nrTraitSNPsCellTypeSpecific 
                                    + "\t" + ratio 
                                    + "\t" + oneTailedP 
                                    + "\t" + fetp 
                                    + "\t" + fetp2;

                            if (alsoPrintAffectedenes) {

                                outputStr += "\t" + Strings.concat(traitSNPsWithSignificantInteractionTerms.toArray(new String[0]), Strings.comma);
                                HashSet<String> eqtls = new HashSet<String>();
                                HashSet<String> probes = new HashSet<String>();
                                for (String rs : traitSNPsWithSignificantInteractionTerms) {
                                    for (String e : eQTLsBeyondThreshold) {
                                        if (e.contains(rs)) {
                                            eqtls.add(e);
                                            probes.add(e.split("-")[1]);
                                        }
                                    }
                                }
                                outputStr += "\t" + Strings.concat(eqtls.toArray(new String[0]), Strings.comma);

                                outputStr += "\t";
                                for (String probe : probes) {
                                    outputStr += probeToHUGO.get(probe) + ",";
                                }

                                outputStr += "\t" + Strings.concat(snpsAfterPruning.toArray(new String[0]), Strings.comma);

                            }

                            outfile.writeln(outputStr);
                            System.out.println(outputStr);
                        }
//                    }

                    }
                } else {
                    // System.out.println("Not testing trait: " + trait.getName() + " nrEQTL SNPs: " + nrTraitSNPsThatAreEQTL);
                }
            }
        }

        outfile.close();
        loader.close();
    }

    public void runSecondVersion(String metaAnalysisMatrix, String metaAnalysisFDRMatrix, double minimalCellTypeSpecificityZScore, double fdrThreshold, String probeTranslation, String gwasCatalog, double gwasSignificanceCutoffPValue,
            String genotypeReferenceDataset, double ldthreshold, int minimumNumberOfCisEQTLSNPs, String queryCovariate, String outputFileName, String[] testSpecificTraits) throws IOException {
        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> probeToHUGO = pb.getProbeTranslation(probeTranslation, "HT12v3.txt", "Gene");

        System.out.println("Writing to: " + outputFileName);

        GWASCatalog catalog = new GWASCatalog(gwasCatalog, gwasSignificanceCutoffPValue);
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaAnalysisMatrix);
        DoubleMatrixDataset<String, String> metaFDRMatrix = new DoubleMatrixDataset<String, String>(metaAnalysisFDRMatrix);

        TriTyperGenotypeData genotypeData = new TriTyperGenotypeData(genotypeReferenceDataset);
        SNPLoader loader = genotypeData.createSNPLoader();
        Independifier independifier = new Independifier(genotypeData, loader);

        // determine the cis-eQTLs that are cell type specific
        TextFile output = new TextFile(outputFileName, TextFile.W);
        output.writeln("Covariate\tGene\tTrait\t"
                + "\tTraitSpecificBeforePruning\tTraitNonSpecificBeforePruning\tOtherSpecificBeforePruning\tOtherNonSpecificBeforePruning\tFetPUnpruned"
                + "\tTraitSpecificAfterPruning\tTraitNonSpecificAfterPruning\tOtherSpecificAfterPruning\tOtherNonSpecificAfterPruning\tFetPPruned");
        System.out.println("Selecting SNPs with P-val<" + gwasSignificanceCutoffPValue);
        for (int covariate = 0; covariate < metaMatrix.nrRows; covariate++) {

            String covariateName = metaMatrix.rowObjects.get(covariate);
            if (queryCovariate == null || queryCovariate.equals(covariateName)) {
                System.out.println("Testing: " + covariateName);
                HashSet<String> cellTypeSpecificCisEQTLSNPs = new HashSet<String>();
                HashSet<String> allOtherCisEQTLSNPsTMP = new HashSet<String>();
                HashSet<String> allEQTLSNPs = new HashSet<String>();
                for (int eqtl = 0; eqtl < metaMatrix.nrCols; eqtl++) {
                    String snp = metaMatrix.colObjects.get(eqtl).split("-")[0];
                    if (metaMatrix.rawData[covariate][eqtl] > minimalCellTypeSpecificityZScore && metaFDRMatrix.rawData[covariate][eqtl] < fdrThreshold) {
                        cellTypeSpecificCisEQTLSNPs.add(snp);
                    } else {
                        allOtherCisEQTLSNPsTMP.add(snp);
                    }
                }

                HashSet<String> allOtherCisEQTLSNPs = new HashSet<String>();
                for (String snp : allOtherCisEQTLSNPsTMP) {
                    if (!cellTypeSpecificCisEQTLSNPs.contains(snp)) {
                        allOtherCisEQTLSNPs.add(snp);
                    }
                }

                System.out.println("Trait Specific/Aspecific: " + cellTypeSpecificCisEQTLSNPs.size() + "/" + allOtherCisEQTLSNPs.size());

                allEQTLSNPs.addAll(cellTypeSpecificCisEQTLSNPs);
                allEQTLSNPs.addAll(allOtherCisEQTLSNPs);

                // prune both for the background
                String[] remainingSNPsCellTypeSpecific = cellTypeSpecificCisEQTLSNPs.toArray(new String[0]);
                String[] remainingSNPsNonSpecific = allOtherCisEQTLSNPs.toArray(new String[0]);
                int nrBeforePruningSpecific = remainingSNPsCellTypeSpecific.length;
                int nrBeforePruningNonSpecific = remainingSNPsNonSpecific.length;
                if (ldthreshold <= 1) {
                    remainingSNPsCellTypeSpecific = independifier.independify(cellTypeSpecificCisEQTLSNPs.toArray(new String[0]), ldthreshold, 1000000);
                    remainingSNPsNonSpecific = independifier.independify(allOtherCisEQTLSNPs.toArray(new String[0]), ldthreshold, 1000000);
                }
                int nrAfterPruningSpecific = remainingSNPsCellTypeSpecific.length;
                int nrAfterPruningNonSpecific = remainingSNPsNonSpecific.length;

                // now for the traits
                for (GWASTrait t : catalog.getTraits()) {
                    boolean testTrait = false;
                    if (testSpecificTraits == null) {
                        testTrait = true;
                    } else {
                        String traitName = t.getName().toLowerCase();
                        for (String traitElems : testSpecificTraits) {
                            if (traitName.contains(traitElems)) {
                                testTrait = true;
                            }
                        }
                    }

                    if (testTrait) {

                        GWASSNP[] traitSNPs = t.getSNPs(gwasSignificanceCutoffPValue);
                        HashSet<String> allTraitSNPsThatAreEQTL = new HashSet<String>();
                        for (GWASSNP s : traitSNPs) {
                            String snp = s.getName();
                            if (allEQTLSNPs.contains(snp)) {
                                allTraitSNPsThatAreEQTL.add(snp);
                            }
                        }

                        if (allTraitSNPsThatAreEQTL.size() >= minimumNumberOfCisEQTLSNPs) {

                            // split between cell type specific and non specific.
                            HashSet<String> traitCellTypeSpecificCisEQTLSNPs = new HashSet<String>();
                            HashSet<String> traitAllOtherCisEQTLSNPs = new HashSet<String>();
                            for (String snp : allTraitSNPsThatAreEQTL) {
                                if (cellTypeSpecificCisEQTLSNPs.contains(snp)) {
                                    traitCellTypeSpecificCisEQTLSNPs.add(snp);
                                } else {
                                    traitAllOtherCisEQTLSNPs.add(snp);
                                }
                            }

                            String[] remainingTraitSNPsCellTypeSpecific = traitCellTypeSpecificCisEQTLSNPs.toArray(new String[0]);
                            String[] remainingTraitSNPsNonSpecific = traitAllOtherCisEQTLSNPs.toArray(new String[0]);
                            int nrBeforePruningTraitSpecific = remainingTraitSNPsCellTypeSpecific.length;
                            int nrBeforePruningTraitNonSpecific = remainingTraitSNPsNonSpecific.length;
                            if (ldthreshold <= 1) {
                                remainingTraitSNPsCellTypeSpecific = independifier.independify(traitCellTypeSpecificCisEQTLSNPs.toArray(new String[0]), ldthreshold, 1000000);
                                remainingTraitSNPsNonSpecific = independifier.independify(traitAllOtherCisEQTLSNPs.toArray(new String[0]), ldthreshold, 1000000);
                            }
                            int nrAfterPruningTraitSpecific = remainingTraitSNPsCellTypeSpecific.length;
                            int nrAfterPruningTraitNonSpecific = remainingTraitSNPsNonSpecific.length;

                            FisherExactTest fet = new FisherExactTest();
                            double fetp = fet.getFisherPValue(nrBeforePruningTraitSpecific, nrBeforePruningTraitNonSpecific, nrBeforePruningSpecific, nrBeforePruningNonSpecific);
                            double fetpAfterPruning = fet.getFisherPValue(nrAfterPruningTraitSpecific, nrAfterPruningTraitNonSpecific, nrAfterPruningSpecific, nrAfterPruningNonSpecific);

                            output.writeln(covariateName + "\t" + probeToHUGO.get(covariateName) + "\t" + t.getName()
                                    + "\t" + nrBeforePruningTraitSpecific
                                    + "\t" + nrBeforePruningTraitNonSpecific
                                    + "\t" + nrBeforePruningSpecific
                                    + "\t" + nrBeforePruningNonSpecific
                                    + "\t" + fetp
                                    + "\t" + nrAfterPruningTraitSpecific
                                    + "\t" + nrAfterPruningTraitNonSpecific
                                    + "\t" + nrAfterPruningSpecific
                                    + "\t" + nrAfterPruningNonSpecific
                                    + "\t" + fetpAfterPruning);
                            System.out.println(covariate + "\t" + t.getName() + "\t" + fetp + "\t" + fetpAfterPruning);
                        }
                    }
                }
            }

        }

        output.close();
        loader.close();
    }
}
