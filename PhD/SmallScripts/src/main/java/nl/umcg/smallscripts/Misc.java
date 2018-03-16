/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Log2Transform;
import umcg.genetica.math.stats.QuantileNormalization;
import umcg.genetica.math.stats.TTest;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class Misc {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here
//            Misc.excludeIndividuals("/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/AllCEUSamplesWithoutFinnish.txt", "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/Merged/PhenotypeInformation.txt");
//            Misc.probeMapForSinaGharib();
//            Misc.SashaSNPData();
//            Misc.compareExpressionLevelsKarin();
//            Misc.determineNumberOfEQTLsWithoutGeneName(
//                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PostQC/eQTLProbesFDR0.05.txt",
//                    //                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz",
//                    "/Volumes/iSnackHD/SkyDrive/MetaAnalysisAnnotationFiles/2012-08-08-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt");

//            Misc.correlationForKarin(
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/Individuals.txt",
//                    "/Volumes/iSnackHD/Data/Projects/KarinFransen/2012-07-19-ProbesForGenesForEQTLs.txt",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/SatVatLiverMuscleBloodHT12/Blood-SAT-VAT-Liver-MuscleHT12CombinedExpressionData/BloodSATVATLiverMuscleHT12ProbesCentered.txt.50PCAsOverSamplesRemoved.txt.TriTyperFormat.txt");
//            Misc.uniquelncount("/Volumes/iSnackHD/SkyDrive/Cis-SNPs/RealData/AllSNPs.txt");
//            Misc.convertToSNPNexusFormat("/Volumes/iSnackHD/SkyDrive/Cis-SNPs/RealData/AllSNPs.txt-WithProxies.txt");
//            Misc.probeOverlap("/Volumes/iSnackHD/Datasets/CisDatasets/sortedForVinod/Vinod-eQTLsFDR0.05.txt", "/Volumes/iSnackHD/Datasets/CisDatasets/Sorted/eQTLs.txt.gz_HumanHT-12_V3_0_R2_11283641_A.txt.txt.gz");
//            Misc.determineASESNPsInEQTLData();
//            Misc.filterEQTLsGivenThreshold("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt", 1.3141052405861965E-4);

            Misc.ensemblify("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt-FilteredForProbeLevelFDR.txt.gz", "/Volumes/iSnackHD/SkyDrive/MetaAnalysisAnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-OnlyEnsemblAnnotation.txt");
        } catch (IOException ex) {
            Logger.getLogger(Misc.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public static void excludeIndividuals(String individualsToInclude, String phenotypeFile) throws IOException {
        TextFile tf = new TextFile(individualsToInclude, TextFile.R);
        HashSet<String> samplesToInclude = new HashSet<String>();
        samplesToInclude.addAll(tf.readAsArrayList());
        tf.close();

        TextFile pheno = new TextFile(phenotypeFile, TextFile.R);
        TextFile phenoOut = new TextFile(phenotypeFile + "-SelectedInds.txt", TextFile.W);
        String[] elems = pheno.readLineElems(TextFile.tab);

        while (elems != null) {
            String sample = elems[0];
            if (samplesToInclude.contains(sample)) {
                elems[2] = "include";
            } else {
                elems[2] = "exclude";
            }

            phenoOut.writeln(Strings.concat(elems, Strings.tab));

            elems = pheno.readLineElems(TextFile.tab);
        }
        phenoOut.close();
        pheno.close();
    }

    public static void probeMapForSinaGharib() throws IOException {
        TextFile in = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt", TextFile.R);
        TextFile pt = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", TextFile.R);
        TextFile out = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt-WithHT12Identifiers.txt", TextFile.W);

        String[] ids = in.readAsArray();

        HashMap<String, String> probeMap = new HashMap<String, String>();

        String[] elems = pt.readLineElems(TextFile.tab);

        elems = pt.readLineElems(TextFile.tab);
        while (elems != null) {

            String probe = elems[0];
            String id = elems[5];

            probeMap.put(probe, id);

            elems = pt.readLineElems(TextFile.tab);
        }

        in.close();
        pt.close();

        for (String probe : ids) {
            out.writeln(probe + "\t" + probeMap.get(probe));
        }

        out.close();


    }

    private static void convertHT12v3ToMetaProbeID() throws IOException {
        TextFile pt = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", TextFile.R);

        String[] header = pt.readLineElems(TextFile.tab);
        String[] elems = pt.readLineElems(TextFile.tab);

        HashMap<String, String> metaProbeMap = new HashMap<String, String>();
        int len = header.length;
        int ctr = 1;
        while (elems != null) {

            if (elems.length > 5) {
                metaProbeMap.put(elems[5], elems[0]);
            }
            elems = pt.readLineElems(TextFile.tab);
            ctr++;
        }

        pt.close();

        TextFile tf = new TextFile("/Volumes/iSnackHD/SkyDrive/PhD/2012-07-30-Vinod_lincRNA_eQTL/2012-08-02-eQTLProbesFDR0.05.txt", TextFile.R);
        String line = tf.readLine();
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];
            String metaProbe = metaProbeMap.get(probe);
            System.out.println(snp + "\t" + metaProbe);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
    }

    public static void probeOverlap(String file1, String file2) throws IOException {
        TextFile tf1 = new TextFile(file1, TextFile.R);
        HashSet<String> probes1 = new HashSet<String>();
        tf1.readLine();
        String[] elems = tf1.readLineElems(TextFile.tab);
        while (elems != null) {
            probes1.add(elems[4]);
            elems = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();



        TextFile tf2 = new TextFile(file2, TextFile.R);
        HashSet<String> probes2 = new HashSet<String>();
        tf2.readLine();
        elems = tf2.readLineElems(TextFile.tab);
        while (elems != null) {
            probes2.add(elems[4]);
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();


        String[] probes1Arr = probes1.toArray(new String[0]);
        String[] probes2Arr = probes2.toArray(new String[0]);


        int ctr1 = 0;
        for (String s : probes1Arr) {
            if (!probes2.contains(s)) {
                System.out.println(s + " unique for ds1");
                ctr1++;
            }
        }

        int ctr2 = 0;
        for (String s : probes2Arr) {
            if (!probes1.contains(s)) {
                System.out.println(s + " unique for ds2");
                ctr2++;
            }
        }

        System.out.println(ctr1);
        System.out.println(ctr2);
    }

    public static void SashaSNPData() throws IOException {
        HashSet<String> genotypedSNPs = new HashSet<String>();
        for (int i = 1; i < 23; i++) {
            TextFile tf = new TextFile("/Volumes/Data/LifeLinesRelease3/LL3Hap2PedAndMap/chr" + i + ".map", TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[1];
                genotypedSNPs.add(snp);
                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();
        }

        TextFile snpmap = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-27-SNPMappings-dbSNP130.txt.gz", TextFile.R);
        String[] snpelems = snpmap.readLineElemsReturnObjects(TextFile.tab);

        HashMap<String, String> snpToChr = new HashMap<String, String>();
        HashMap<String, String> snpToChrPos = new HashMap<String, String>();

        while (snpelems != null) {

            String chr = snpelems[0];
            String chrpos = snpelems[1];
            String snp = snpelems[2];

            snpToChr.put(snp, chr);
            snpToChrPos.put(snp, chrpos);
            snpelems = snpmap.readLineElemsReturnObjects(TextFile.tab);
        }

        snpmap.close();

        TextFile tf2 = new TextFile("/Volumes/iSnackHD/Data/GWAS/LifeLines3/BeagleQualityScores.txt", TextFile.R);
        tf2.readLineElems(TextFile.tab);
        String[] elems = tf2.readLineElems(TextFile.tab);
        TextFile out = new TextFile("/Volumes/iSnackHD/LL3SNPImputationQuality.txt.gz", TextFile.W);
        out.writeln("SNP\tChr\tChrPos\tImputed?\tImputationQuality");
        while (elems != null) {
            if (elems.length > 2) {
                String snp = elems[1];
                boolean snpImputed = true;
                if (genotypedSNPs.contains(snp)) {
                    snpImputed = false;
                }
                out.writeln(snp + "\t" + snpToChr.get(snp) + "\t" + snpToChrPos.get(snp) + "\t" + snpImputed + "\t" + elems[elems.length - 2]);
                System.out.println(snp + "\t" + snpToChr.get(snp) + "\t" + snpToChrPos.get(snp) + "\t" + snpImputed + "\t" + elems[elems.length - 2]);
            }
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();
        out.close();

    }

    public static void compareExpressionLevelsKarin() throws IOException {
//      DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/SatVatLiverMuscleBloodHT12/Blood-SAT-VAT-Liver-MuscleHT12CombinedExpressionData/BloodSATVATLiverMuscleHT12ProbesCentered.txt.50PCAsOverSamplesRemoved.txt.TriTyperFormat.txt.gz");

        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData/ExpressionData.txt");

        QuantileNormalization.quantilenormalize(ds.rawData);
        Log2Transform.log2transform(ds.rawData);
        TextFile probeFile = new TextFile("/Volumes/iSnackHD/Data/Projects/KarinFransen/2012-07-25-ProbesForGenesForEQTLs.txt", TextFile.R);

        int offset = 0;
        ArrayList<String> probes = probeFile.readAsArrayList();
        probeFile.close();

        // load eQTL file
        HashSet<String> probesHavingCisEffect = new HashSet<String>();
        TextFile eqtlfile = new TextFile("/Volumes/iSnackHD/Data/Projects/KarinFransen/eQTLResults/2012-07-25-Blood-CIS-1M-100P/eQTLProbesFDR0.05.txt", TextFile.R);
        String[] elems = eqtlfile.readLineElems(TextFile.tab);
        elems = eqtlfile.readLineElems(TextFile.tab);
        while (elems != null) {
            probesHavingCisEffect.add(elems[4]);
            elems = eqtlfile.readLineElems(TextFile.tab);
        }
        eqtlfile.close();


        HashSet<String> probesWithoutCisEffect = new HashSet<String>();

        for (String p : probes) {
            if (!probesHavingCisEffect.contains(p)) {
                probesWithoutCisEffect.add(p);
            }
        }

        // average expression for the probes showing cis-effect
        double[] expressionForProbesWithCisEffect = new double[ds.colObjects.size() - offset];

        System.out.println("");
        System.out.println("--------------------");
        System.out.println("Probes Having Cis-FX");
        System.out.println("--------------------");
        System.out.println("");
        int nrNan = 0;
        for (String p : probesHavingCisEffect) {
            Integer probeId = ds.hashRows.get(p);


            for (int i = offset; i < ds.colObjects.size(); i++) {
                if (Double.isNaN(ds.rawData[probeId][i])) {
                    nrNan++;
                    System.out.println(probeId + "\t" + i + "\tisNAN");
                }
                expressionForProbesWithCisEffect[i - offset] += ds.rawData[probeId][i];
            }
            double mean = Descriptives.mean(expressionForProbesWithCisEffect);

            System.out.println(p + "\t" + probeId + "\t" + mean);
        }
        System.out.println("Total NaN: " + nrNan);

        for (int i = 0; i < expressionForProbesWithCisEffect.length; i++) {
            expressionForProbesWithCisEffect[i] /= (double) probesHavingCisEffect.size();
        }

        System.out.println("");
        System.out.println("------------------");
        System.out.println("Now testing probes");
        System.out.println("------------------");
        System.out.println("");
        for (String p : probesWithoutCisEffect) {
            Integer probeId = ds.hashRows.get(p);
            if (probeId != null) {
                double[] expressionForProbeWithoutCisEffect = new double[expressionForProbesWithCisEffect.length];
                for (int i = offset; i < ds.rawData[probeId].length; i++) {
                    expressionForProbeWithoutCisEffect[i - offset] = ds.rawData[probeId][i];
                }
                double mean = Descriptives.mean(expressionForProbeWithoutCisEffect);
                double tTestPval = TTest.test(expressionForProbesWithCisEffect, expressionForProbeWithoutCisEffect);
                WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();
                double wilcoxonPval = mwm.returnWilcoxonMannWhitneyPValue(expressionForProbesWithCisEffect, expressionForProbeWithoutCisEffect);
                System.out.println(p + "\t" + probeId + "\t" + mean + "\t" + tTestPval + "\t" + wilcoxonPval);
            }
        }
    }

    public static void determineWhicheQTLsAreMissing(String fehrmannFile, String newMetaFile) throws IOException {

        TextFile tf = new TextFile(fehrmannFile, TextFile.R);

        tf.close();

        TextFile newFile = new TextFile(newMetaFile, TextFile.R);

        newFile.close();
    }

    public static void determineNumberOfEQTLsWithoutGeneName(String eqtlfile, String annotation) throws IOException {
        HashMap<String, String> probeToGene = new HashMap<String, String>();
        TextFile tf = new TextFile(annotation, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        elems = tf.readLineElems(TextFile.tab);
        int nrProbes = 0;
        int nrProbesWithGenes = 0;
        HashSet<String> uniqueGenesInAnnot = new HashSet<String>();
        while (elems != null) {

            if (elems.length >= 5) {
                String probe = elems[0];
                String gene = elems[4];
                if (gene.trim().length() > 0 && !gene.equals("-")) {
                    probeToGene.put(probe, gene);
                    nrProbesWithGenes++;

//                    System.out.println(gene);
                    uniqueGenesInAnnot.add(gene);
                }
            }
            nrProbes++;
            elems = tf.readLineElems(TextFile.tab);
        }

        System.out.println("Annotation file");
        System.out.println("--------------------------");
        System.out.println("Nr Probes: " + nrProbes);
        System.out.println("Nr Probes With Ensembl Genes: " + nrProbesWithGenes);
        System.out.println("Unique Genes in Annotation: " + uniqueGenesInAnnot.size());
        System.out.println("--------------------------");
        System.out.println("");
        tf.close();

        System.out.println("eQTL file: ");
        System.out.println(eqtlfile);
        System.out.println("--------------------------");
        TextFile eFile = new TextFile(eqtlfile, TextFile.R);
        elems = eFile.readLineElems(TextFile.tab); // header
        elems = eFile.readLineElems(TextFile.tab);
        HashSet<String> uniqueGenes = new HashSet<String>();
        int nrWithGenes = 0;
        int nrWithoutGenes = 0;
        int nrEQTLs = 0;
        while (elems != null) {
            String probe = elems[4];
            String gene = probeToGene.get(probe);
            if (gene == null) {
                if (!elems[eQTLTextFile.HUGO].equals("-")) {
                    System.out.println(Strings.concat(elems, Strings.tab));
                }
                nrWithoutGenes++;
//                System.out.println(Strings.concat(elems, Strings.tab));
            } else {
                nrWithGenes++;
                uniqueGenes.add(gene);
            }
            nrEQTLs++;
            elems = eFile.readLineElems(TextFile.tab);
        }
        eFile.close();
        System.out.println("Nr EQTLs: " + nrEQTLs);
        System.out.println("Unique genes: " + uniqueGenes.size());
        System.out.println("eQTLs with gene name: " + nrWithGenes);
        System.out.println("eQTLs without gene name: " + nrWithoutGenes);
    }

    public static void correlationForKarin(String indsInBloodDataset, String probeFile, String infile) throws IOException {

        TextFile tf1 = new TextFile(indsInBloodDataset, TextFile.R);
        HashSet<String> indsHT12 = new HashSet<String>();
        indsHT12.addAll(tf1.readAsArrayList());
        tf1.close();

        System.out.println(indsHT12.size() + " individuals");

        TextFile tf = new TextFile(probeFile, TextFile.R);
        String[] probesToTest = tf.readAsArray();
        tf.close();

        System.out.println(probesToTest.length + " probes");

        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(infile);

        int ctr = 0;
        ArrayList<Integer> indsArl = new ArrayList<Integer>();
        for (int c = 0; c < ds.nrCols; c++) {
            if (indsHT12.contains(ds.colObjects.get(c))) {
                indsArl.add(c);
            }
        }

        double[][] correlationMatrix = new double[probesToTest.length][probesToTest.length];

        Integer[] probeIds = new Integer[probesToTest.length];
        for (int p = 0; p < probesToTest.length; p++) {
            probeIds[p] = ds.hashRows.get(probesToTest[p]);
        }

        System.out.println("Running correlations");
        int[] binsOrig = new int[1000];
        for (int i = 0; i < probesToTest.length; i++) {
            int probeId1 = probeIds[i];
            double[] probe1vals = getVals(ds.rawData[probeId1], indsArl);
            correlationMatrix[i][i] = 1d;
            for (int j = i + 1; j < probesToTest.length; j++) {
                int probeId2 = probeIds[j];
                double[] probe2vals = getVals(ds.rawData[probeId2], indsArl);
                correlationMatrix[i][j] = JSci.maths.ArrayMath.correlation(probe1vals, probe2vals);
                correlationMatrix[j][i] = correlationMatrix[i][j] * correlationMatrix[i][j];
                correlationMatrix[i][j] = correlationMatrix[j][i];
                int binNo = (int) Math.floor((correlationMatrix[i][j] * binsOrig.length));
                binsOrig[binNo]++;
            }
        }


        System.out.println("Done. Saving file");
        DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>();
        ds2.rawData = correlationMatrix;
        ds2.nrCols = probesToTest.length;
        ds2.nrRows = probesToTest.length;

        ds2.colObjects = Arrays.asList(probesToTest);
        ds2.rowObjects = Arrays.asList(probesToTest);
        ds2.save("/Volumes/iSnackHD/Data/Projects/KarinFransen/2012-08-13-CorrelationMatrix.txt");


        HashSet<Integer> probesInRealData = new HashSet<Integer>();
        probesInRealData.addAll(Arrays.asList(probeIds));

        int[] bins = new int[1000];
        for (int perm = 0; perm < 1000; perm++) {

            System.out.println("Perm: " + perm);
            Integer[] randomProbeIds = null;

            HashSet<Integer> selectedIds = new HashSet<Integer>();
            int probesSelected = 0;
            while (probesSelected < probesToTest.length) {
                int randomId = (int) Math.floor(Math.random() * ds.nrRows);
                if (!selectedIds.contains(randomId) && !probesInRealData.contains(randomId)) {
                    selectedIds.add(randomId);
                    probesSelected++;
                }
            }

            randomProbeIds = selectedIds.toArray(new Integer[0]);

            for (int i = 0; i < randomProbeIds.length; i++) {
                int probeId1 = randomProbeIds[i];
                double[] probe1vals = getVals(ds.rawData[probeId1], indsArl);
                correlationMatrix[i][i] = 1d;
                for (int j = i + 1; j < randomProbeIds.length; j++) {
                    int probeId2 = randomProbeIds[j];
                    double[] probe2vals = getVals(ds.rawData[probeId2], indsArl);
                    double corr = JSci.maths.ArrayMath.correlation(probe1vals, probe2vals);
                    int binNo = (int) Math.floor(((corr * corr) * bins.length));
                    bins[binNo]++;
                }
            }

        }

        TextFile randomCorr = new TextFile("/Volumes/iSnackHD/Data/Projects/KarinFransen/2012-08-13-RandomGeneCorrelations.txt", TextFile.W);
        for (int i = 0; i < bins.length; i++) {
            randomCorr.writeln(i + "\t" + ((double) i / bins.length) + "\t" + bins[i] + "\t" + binsOrig[i]);
        }
        randomCorr.close();


    }

    private static double[] getVals(double[] data, ArrayList<Integer> indsArl) {
        double[] out = new double[indsArl.size()];
        for (int i = 0; i < indsArl.size(); i++) {
            out[i] = data[indsArl.get(i)];
        }
        return out;

    }

    private static void uniquelncount(String f) throws IOException {
        HashSet<String> unique = new HashSet<String>();
        TextFile tf = new TextFile(f, TextFile.R);
        unique.addAll(tf.readAsArrayList());
        tf.close();
        System.out.println(unique.size() + " unique elements.");
    }

    private static void convertToSNPNexusFormat(String string) throws IOException {
        HashSet<String> snps = new HashSet<String>();
        TextFile tf = new TextFile(string, TextFile.R);


        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            snps.add(elems[0]);
            snps.add(elems[1]);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tfout = new TextFile(string + "-SNPNexusFormat.txt", TextFile.W);
        for (String snp : snps) {
            tfout.writeln("dbsnp\t" + snp);
        }
        tfout.close();

    }

    private static void determineASESNPsInEQTLData() throws IOException {
        HashSet<String> snpsASE = new HashSet<String>();
        TextFile in = new TextFile("/Volumes/iSnackHD/Data/Projects/MarijkeVanDerSijde/snps.txt", TextFile.R);
        String[] elems = in.readLineElems(TextFile.tab);
        while (elems != null) {
            snpsASE.add(elems[0]);
            elems = in.readLineElems(TextFile.tab);
        }

        in.close();

        System.out.println(snpsASE.size() + " unique SNPs with ASE");

        HashSet<String> snpsNotASE = new HashSet<String>();
        TextFile in3 = new TextFile("/Volumes/iSnackHD/Data/Projects/MarijkeVanDerSijde/ase_results_metaAnalysis_noThreshold.txt", TextFile.R);
        String[] elems3 = in3.readLineElems(TextFile.tab);
        while (elems3 != null) {
            if (!snpsASE.contains(elems3[0])) {
                snpsNotASE.add(elems3[0]);
            }
            elems3 = in3.readLineElems(TextFile.tab);
        }
        in3.close();

        System.out.println(snpsNotASE.size() + " unique SNPs without ASE");

        HashSet<String> snpsASETestedForEQTLs = new HashSet<String>();
        HashSet<String> snpsNotASETestedForEQTLs = new HashSet<String>();
        TextFile in4 = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/QC/TruePositives.txt.gz", TextFile.R);
        System.out.println("Reading: /Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/QC/TruePositives.txt.gz");
        String[] elems4 = in4.readLineElems(TextFile.tab);
        while (elems4 != null) {
            String snp = elems4[1];

            if (snpsASE.contains(snp)) {
                snpsASETestedForEQTLs.add(snp);
            }
            if (snpsNotASE.contains(snp)) {
                snpsNotASETestedForEQTLs.add(snp);
            }
            elems4 = in4.readLineElems(TextFile.tab);
        }

        in4.close();
        System.out.println(snpsASETestedForEQTLs.size() + " unique SNPs with ASE tested for eQTLs");
        System.out.println(snpsNotASETestedForEQTLs.size() + " unique SNPs without ASE  tested for eQTLs");


        TextFile in2 = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PostQC/eQTLsFDR0.05.txt", TextFile.R);
        String[] elems2 = in2.readLineElems(TextFile.tab);
        int nrSNPsWitheQTL = 0;
        HashSet<String> uniqueSNPs = new HashSet<String>();

        int nrSNPsASEWeQTL = 0;
        int nrSNPsNotASEWeQTL = 0;


        while (elems2 != null) {
            String snp = elems2[1];
            if (!uniqueSNPs.contains(snp)) {
                if (snpsASETestedForEQTLs.contains(snp)) {
                    nrSNPsASEWeQTL++;
                }
                if (snpsNotASETestedForEQTLs.contains(snp)) {
                    nrSNPsNotASEWeQTL++;
                }
                uniqueSNPs.add(snp);
            }
            elems2 = in2.readLineElems(TextFile.tab);
        }


        int nrSNPsASEWoQTL = snpsASETestedForEQTLs.size() - nrSNPsASEWeQTL;
        int nrSNPsNotASEW0QTL = snpsNotASETestedForEQTLs.size() - nrSNPsNotASEWeQTL;

        System.out.println("2x2:\teQTL\tnoEQTL\tTotal");
        System.out.println("ASE\t" + nrSNPsASEWeQTL + "\t" + nrSNPsASEWoQTL + "\t" + (nrSNPsASEWeQTL + nrSNPsASEWoQTL));
        System.out.println("NASE\t" + nrSNPsNotASEWeQTL + "\t" + nrSNPsNotASEW0QTL + "\t" + (nrSNPsNotASEWeQTL + nrSNPsNotASEW0QTL));
        System.out.println("Tot:\t" + (nrSNPsASEWeQTL + nrSNPsNotASEWeQTL) + "\t" + (nrSNPsASEWoQTL + nrSNPsNotASEW0QTL) + "\t" + ((nrSNPsASEWeQTL + nrSNPsASEWoQTL) + (nrSNPsNotASEWeQTL + nrSNPsNotASEW0QTL)));

        double x = ChiSquare.getX(nrSNPsASEWeQTL, nrSNPsASEWoQTL, nrSNPsNotASEWeQTL, nrSNPsNotASEW0QTL);
        double p = ChiSquare.getP(1, x);
        System.out.println("X:" + x);
        System.out.println("X:" + p);

        in2.close();
    }

    private static void filterEQTLsGivenThreshold(String string, double d) throws IOException {
        TextFile out = new TextFile(string + "-FilteredForProbeLevelFDR.txt.gz", TextFile.W);
        TextFile in = new TextFile(string, TextFile.R);
        out.writeln(in.readLine());
        String[] elems = in.readLineElems(TextFile.tab);
        while (elems != null) {

            double p = Double.parseDouble(elems[0]);
            if (p <= d) {
                out.writeln(Strings.concat(elems, Strings.tab));
            }


            elems = in.readLineElems(TextFile.tab);
        }

        in.close();





        out.close();
    }

    private static void ensemblify(String volumesiSnackHDeQTLMetaAnalysisMetaAna, String volumesiSnackHDSkyDriveMetaAnalysisAnn) throws IOException {
        HashMap<String, String> probeMap = new HashMap<String, String>();

        TextFile in = new TextFile(volumesiSnackHDSkyDriveMetaAnalysisAnn, TextFile.R);
        String[] elems = in.readLineElems(TextFile.tab);
        elems = in.readLineElems(TextFile.tab);
        while (elems != null) {
            String probe = elems[0];
            if (elems.length > 4) {
                String ensembl = elems[4];
                probeMap.put(probe, ensembl);
            }
            elems = in.readLineElems(TextFile.tab);
        }

        in.close();

        TextFile out = new TextFile(volumesiSnackHDeQTLMetaAnalysisMetaAna + "-Ensembl.txt", TextFile.W);
        TextFile in2 = new TextFile(volumesiSnackHDeQTLMetaAnalysisMetaAna, TextFile.R);
        out.writeln(in2.readLine());
        elems = in2.readLineElems(TextFile.tab);
        while (elems != null) {
            String probe = elems[4];
            String ensembl = probeMap.get(probe);
            if (ensembl != null) {
                elems[eQTLTextFile.HUGO] = ensembl;
            }
            out.writeln(Strings.concat(elems, Strings.tab));
            elems = in2.readLineElems(TextFile.tab);
        }

        in2.close();
        out.close();
    }
}
