/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package yanglicausalinferencedataextractor;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class YangLiCausalInferenceDataExtractor {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        System.out.println("");
        System.out.println("This is YangLiCausalInferenceDataExtractor v1.");
        System.out.println("");

        if (args.length < 6) {
            System.out.println("Usage: YangLiCausalInferenceDataExtractor.jar snpcistranspairs.txt trityperdir cisExpressionFile.txt.gz transExpressionFile.txt.gz probetranslation.txt.gz outfile.txt.gz [genotypeToExpressionCoupling.txt]");
        } else {
            try {
                System.out.println("Dear collaborator, how are you today?");
                YangLiCausalInferenceDataExtractor y = new YangLiCausalInferenceDataExtractor();
                if (args.length == 6) {
                    y.run(args[0], args[1], args[2], args[3], args[4], args[5], null);
                } else {
                    y.run(args[0], args[1], args[2], args[3], args[4], args[5], args[6]);
                }
            } catch (IOException e) {
                System.err.println(e.getMessage());
                e.printStackTrace();
            } catch (Exception e) {
                System.err.println(e.getMessage());
                e.printStackTrace();
            }
        }
    }

    public void run(String snpCisTransPairs, String genotypeData, String cisExpressionData, String transExpressionData, String probeTranslation, String outfileName, String genotypeToExpressionCoupling) throws IOException, Exception {

        if (!Gpio.exists(snpCisTransPairs)) {
            throw new IOException("File does not exist: " + snpCisTransPairs);
        }
        if (!Gpio.isDir(genotypeData)) {
            throw new IOException("Location is not a directory: " + genotypeData);
        }
        if (!Gpio.exists(cisExpressionData)) {
            throw new IOException("File does not exist: " + cisExpressionData);
        }
        if (!Gpio.exists(transExpressionData)) {
            throw new IOException("File does not exist: " + transExpressionData);
        }
        if (genotypeToExpressionCoupling == null) {
            System.out.println("Will not use explicit genotype to expression coupling file");
        } else if (!Gpio.exists(genotypeToExpressionCoupling)) {
            throw new IOException("File does not exist: " + genotypeToExpressionCoupling);
        }

        TriTyperGeneticalGenomicsDatasetSettings cisSettings = new TriTyperGeneticalGenomicsDatasetSettings();
        cisSettings.genotypeLocation = genotypeData;
        cisSettings.expressionLocation = cisExpressionData;
        cisSettings.genotypeToExpressionCoupling = genotypeToExpressionCoupling;
        cisSettings.cisAnalysis = true;
        cisSettings.transAnalysis = true;
        cisSettings.name = "CIS";

        TriTyperGeneticalGenomicsDatasetSettings transSettings = new TriTyperGeneticalGenomicsDatasetSettings();
        transSettings.genotypeLocation = genotypeData;
        transSettings.expressionLocation = transExpressionData;
        transSettings.genotypeToExpressionCoupling = genotypeToExpressionCoupling;
        transSettings.cisAnalysis = true;
        transSettings.transAnalysis = true;
        transSettings.name = "TRANS";

        TriTyperGeneticalGenomicsDataset cisData = new TriTyperGeneticalGenomicsDataset(cisSettings);
        TriTyperGeneticalGenomicsDataset transData = new TriTyperGeneticalGenomicsDataset(transSettings);

        if (cisData.getTotalGGSamples() != transData.getTotalGGSamples()) {
            throw new Exception("Error: cis and trans data do not contain the same amount of samples!");
        }

        System.out.println("");
        System.out.println("");
        // maps probes to metaProbes
        HashMap<String, String> probeToMetaMap = null;
//        HashMap<String, String> metaToProbeMap = null;
        TextFile probeMapIn = new TextFile(probeTranslation, TextFile.R);
        probeToMetaMap = (HashMap<String, String>) probeMapIn.readAsHashMap(1, 0);
        probeMapIn.close();
//        probeMapIn.open();
//        metaToProbeMap = (HashMap<String, String>) probeMapIn.readAsHashMap(1, 0);
//        probeMapIn.close();

        // map genotype samples across datasets.
        TriTyperExpressionData cisExpData = cisData.getExpressionData();
        TriTyperExpressionData transExpData = transData.getExpressionData();
        TriTyperGenotypeData cisGenotypeData = cisData.getGenotypeData();

        double[][] rawCisExp = cisExpData.getMatrix();
        double[][] rawTransExp = transExpData.getMatrix();

        String[] cisSamples = cisExpData.getIndividuals();

        int[] cisExpToGenotype = cisData.getExpressionToGenotypeIdArray();
        int[] cisGenotypeToExp = new int[cisGenotypeData.getIndividuals().length];

        HashMap<Integer, Integer> cisGenotypeToExpressionId = cisData.getGenotypeToExpressionIdHash();
        for (int i = 0; i < cisGenotypeToExp.length; i++) {
            if (cisGenotypeToExpressionId.get(i) == null) {
                cisGenotypeToExp[i] = -1;
            } else {
                cisGenotypeToExp[i] = cisGenotypeToExpressionId.get(i);
            }
        }
        int[] cisExpToTransExp = new int[cisExpToGenotype.length];

        boolean error = false;
        int sharedSamples = 0;
        for (int i = 0; i < cisSamples.length; i++) {
            String sample = cisSamples[i];
            Integer transSampleId = transExpData.getIndividualId(sample);
            if (transSampleId == null) {
                System.err.println("ERROR: " + sample + "\tdoes not exist in trans gene expression data");
                error = true;
            } else {
                cisExpToTransExp[i] = transSampleId;
                sharedSamples++;
            }
        }

        if (error) {
            System.err.println("---");
            throw new Exception("Errors encountered during processing.");
        }
        System.out.println("Shared samples: " + sharedSamples);

        System.out.println("Everything looks fine. Will now load: " + snpCisTransPairs);

        HashSet<Triple<String, String, String>> combos = new HashSet<Triple<String, String, String>>();

        TextFile snpIn = new TextFile(snpCisTransPairs, TextFile.R);
        String[] elems = snpIn.readLineElems(TextFile.tab);
        while (elems != null) {

            combos.add(new Triple<String, String, String>(elems[0], elems[1], elems[2]));
            elems = snpIn.readLineElems(TextFile.tab);
        }
        snpIn.close();


        System.out.println("Hold on while I convert your data....");
        // write header
        TextFile tfOut = new TextFile(outfileName, TextFile.W);
        TextFile tfOutFailed = new TextFile(outfileName + "-FailedCombos.txt", TextFile.W);
        String sampleOutput = "";
        for (int i = 0; i < cisData.getTotalGGSamples(); i++) {
            sampleOutput += "\tSample" + i + "\tCisExpSample" + i + "\tTransExpSample" + i;
        }
        tfOut.writeln("SNP\tCisProbe\tTransProbe" + sampleOutput);

        // write combos

        SNPLoader loader = cisGenotypeData.createSNPLoader();
        for (Triple<String, String, String> t : combos) {

            String snpName = t.getLeft();
            String cisProbe = t.getMiddle();
            String transProbe = t.getRight();

            Integer cisProbeId = cisExpData.getProbeToId().get(cisProbe);
            Integer transProbeId = transExpData.getProbeToId().get(transProbe);
            Integer snpId = cisGenotypeData.getSnpToSNPId().get(snpName);

            if (cisProbeId == null || transProbeId == null || snpId == null) {
                System.err.println("ERROR: snp, cis or trans probe not found in expression data: snp:\t" + snpName + " - ID:" + snpId
                        + "\tcis: " + cisProbe + " - ID: " + cisProbeId + "\ttrans: " + transProbe + " - ID: " + transProbeId);
                tfOutFailed.writeln("ERROR: snp, cis or trans probe not found in expression data: snp:\t" + snpName + " - ID:" + snpId
                        + "\tcis: " + cisProbe + " - ID: " + cisProbeId + "\ttrans: " + transProbe + " - ID: " + transProbeId);
            } else {

                SNP snpObj = cisGenotypeData.getSNPObject(snpId);
                loader.loadGenotypes(snpObj);

                if (loader.hasDosageInformation()) {
                    loader.loadDosage(snpObj);
                }

                if (snpObj.getMAF() < 0.05 || snpObj.getHWEP() < 0.0001 || snpObj.getCR() < 0.95) {
                    tfOutFailed.writeln("ERROR: snp failed QC thresholds:\t" + snpName + " - ID:" + snpId
                            + "\tcis: " + cisProbe + " - ID: " + cisProbeId + "\ttrans: " + transProbe + " - ID: " + transProbeId + "\t" + snpObj.getMAF() + "\t" + snpObj.getHWEP() + "\t" + snpObj.getCR());
                } else {
                    if (probeToMetaMap.get(cisProbe) == null) {
                        System.err.println("Error: probe " + cisProbe + " has no META ID?");
                    }
                    if (probeToMetaMap.get(transProbe) == null) {
                        System.err.println("Error: probe " + transProbe + " has no META ID?");
                    }
                    
                    
                    String output = snpName + "\t" + probeToMetaMap.get(cisProbe) + "\t" + probeToMetaMap.get(transProbe);
                    byte[] alleles1 = snpObj.getAllele1();
                    byte[] alleles2 = snpObj.getAllele2();
                    for (int i = 0; i < alleles1.length; i++) {
                        if (cisGenotypeToExp[i] != -1 && cisGenotypeData.getIsIncluded()[i]) {
                            int cisId = cisGenotypeToExp[i];
                            int transId = cisExpToTransExp[cisId];
//                            System.out.println(cisId + "\t" + transId);
                            output += "\t" + BaseAnnot.toString(alleles1[i]) + BaseAnnot.toString(alleles2[i]) + "\t" + rawCisExp[cisProbeId][cisId] + "\t" + rawTransExp[transProbeId][transId];
                        }
                    }
                    tfOut.writeln(output);
                }
                snpObj.clearGenotypes();
            }
        }

        tfOut.close();
        tfOutFailed.close();

        loader.close();

        // check the output
        // 3 cols per sample + 3 cols headers
        System.out.println("Done converting!");
        System.out.println("Now testing output file...");
        int nrColsReq = (cisData.getTotalGGSamples() * 3) + 3;
        TextFile tfIn = new TextFile(outfileName, TextFile.R);
        String header = tfIn.readLine();
        elems = tfIn.readLineElems(TextFile.tab);
        boolean containserrors = false;
        while (elems != null) {
            if (elems.length != nrColsReq) {
                System.err.println("Error: expecting " + nrColsReq + "\tbut found: " + elems.length);
                containserrors = true;
            }
            elems = tfIn.readLineElems(TextFile.tab);
        }
        tfIn.close();

        if (!containserrors) {
            System.out.println("Output file contains no errors! Have a nice day!...\n"
                    + "Now all your base are belong to us\nLock phasers on target!\nLet's get some cake.");
        }

    }
}
