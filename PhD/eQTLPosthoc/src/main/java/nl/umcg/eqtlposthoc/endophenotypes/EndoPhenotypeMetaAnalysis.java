/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.endophenotypes;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.TTest;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class EndoPhenotypeMetaAnalysis {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            
            String probeTranslation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            String[] desc = new String[]{
                "DILGOM",
                "EGCUT",
                "Rotterdam Study",
                "SHIP",
                "GroningenHT12v3",
                "InCHIANTI",
                "HVH-HT12v3",
                "HVH-HT12v4"};
//            String[] files = new String[]{
//                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/DILGOM/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt.gz_Correlations.txt",
//                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/EGCUT/QuantLog-Correlations.txt",
//                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Rotterdam/CorrelationOutput-FDR0.5-0PCs-adj/Correlations.txt",
//                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/SHIP/2013-02-01_PhenoGXcorrelation/Correlations.txt",
//                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Groningen/0PCs/Correlations.txt.gz"
//            };

            String[] phenodescriptions = new String[]{
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/DILGOM/EndophenotypeDescriptions.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/EGCUT/EndophenotypeDescriptions.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Rotterdam/EndophenotypeDescriptions.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/SHIP/EndophenotypeDescriptions.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Groningen/EndophenotypeDescriptions.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/InChianti/EndophenotypeDescriptions.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/HVH/EndophenotypeDescriptions.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/HVH/EndophenotypeDescriptions.txt"
            };

            String[] files = new String[]{
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/DILGOM/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt.gz_Correlations.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/EGCUT/QuantLog-Correlations.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Rotterdam/CorrelationOutput-FDR0.5-0PCs-adj/Correlations.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/SHIP/2013-02-01_PhenoGXcorrelation/Correlations.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Groningen/0PCs/Correlations.txt.gz",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/InChianti/Correlations_CovariatesRemoved.txt.gz",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/HVH/v3z/Correlations.txt.gz",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/HVH/v4Fz/Correlations.txt.gz"
            };
//            String output = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Meta/0PCRun2/0PCs";
            String output = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Meta0PCTransEQTL/0PCs";

            EndoPhenotypeMetaAnalysis m = new EndoPhenotypeMetaAnalysis();
            String[] phenotypes1 = m.run(files, output, desc, phenodescriptions, probeTranslation);

            String[] files2 = new String[]{
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/DILGOM/ExpressionData.txt.gz.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt_Correlations.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/EGCUT/40PCs-Correlations.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Rotterdam/CorrelationOutput-FDR0.5-40PCs-adj/Correlations.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/SHIP/2013-02-01_PhenoGXcorrelation-40PCs/Correlations.txt",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Groningen/40PCs/Correlations.txt.gz",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/InChianti/Correlations_GeneticVectorsNotRemoved.txt.gz",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/HVH/v3gvnr/Correlations.txt.gz",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/HVH/v4gvnr/Correlations.txt.gz"
            };
//
            output = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Meta0PCTransEQTL/40PCs";

            String[] phenotypes2 = m.run(files2, output, desc, phenodescriptions, probeTranslation);

            HashSet<String> allphenotypes = new HashSet<String>();
            allphenotypes.addAll(Arrays.asList(phenotypes1));
            allphenotypes.addAll(Arrays.asList(phenotypes2));

            output = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/Meta0PCTransEQTL/";
//            String eQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.5.txt.gz";

            String gwasCatalog = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt";
            String eQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-11-16-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-TRANS-0PCs-0GWASPCs/eQTLsFDR.txt.gz";
            m.combine(output, new String[]{"0PCs", "40PCs"}, allphenotypes.toArray(new String[0]), output + "MetaMerged.txt", eQTLFile, gwasCatalog);


        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    /*
     * Platform: HT12v4
     Phenotype	Probe	PlatformProbe	NumSamples	SpearmanR	PearsonR
     */
    int probeCol = 1;
    int phenoCol = 0;
    int pearsoncol = 5;
    int spearmancol = 4;
    int samplecol = 3;
    private int minNrSamples = 1500;

    public String[] run(String[] files, String output, String[] desc, String[] phenodescriptions, String probeTranslation) throws IOException {

        // reading endophenotype description table

        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> probeToGeneMap = pb.getProbeTranslation(probeTranslation, "Probe", "HUGO");



        // iterate through files, determine available phenotypes
        HashSet<String> uniquePhenotypes = new HashSet<String>();
        HashSet<String> uniqueProbes = new HashSet<String>();
        HashMap<String, Integer> maxSamplesPerPhenotype = new HashMap<String, Integer>();
        for (int q = 0; q < files.length; q++) {
            String f = files[q];
            HashMap<String, String> phenotypeTranslationTable = new HashMap<String, String>();
            TextFile phe = new TextFile(phenodescriptions[q], TextFile.R);
            String[] phelems = phe.readLineElems(TextFile.tab);
            while (phelems != null) {
                String metapheno = phelems[1].toLowerCase().trim().replaceAll(" ", "");
                phenotypeTranslationTable.put(phelems[0].toLowerCase().trim(), metapheno);
                phelems = phe.readLineElems(TextFile.tab);
            }
            phe.close();

            TextFile tf = new TextFile(f, TextFile.R);
            System.out.println("Platform: " + tf.readLine() + "\t" + f);
            String[] elems = tf.readLineElems(TextFile.tab);
            elems = tf.readLineElems(TextFile.tab);

            HashMap<String, Integer> maxSamplesPerPhenotypeForDataset = new HashMap<String, Integer>();
            while (elems != null) {

                String pheno = elems[phenoCol].toLowerCase().trim();
                String metapheno = phenotypeTranslationTable.get(pheno);
                if (metapheno == null) {
                    System.err.println("ERROR PHENO NOT FOUND: " + pheno);
                } else {
                    Integer nrSamples = Integer.parseInt(elems[samplecol]);
                    Integer currentMax = maxSamplesPerPhenotypeForDataset.get(metapheno);
                    if (currentMax == null || nrSamples > currentMax) {
                        maxSamplesPerPhenotypeForDataset.put(metapheno, nrSamples);
                    }
                    uniquePhenotypes.add(metapheno);
                }
                uniqueProbes.add(elems[probeCol]);

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            for (String pheno : uniquePhenotypes) {
                Integer maxSamplesForDS = maxSamplesPerPhenotypeForDataset.get(pheno);
                if (maxSamplesForDS != null) {
                    Integer currentSampleCount = maxSamplesPerPhenotype.get(pheno);
                    if (currentSampleCount == null) {
                        currentSampleCount = maxSamplesForDS;
                    } else {
                        currentSampleCount += maxSamplesForDS;
                    }
                    maxSamplesPerPhenotype.put(pheno, currentSampleCount);
                }
            }
        }

        System.out.println("Detected phenotypes:");
        for (String s : uniquePhenotypes) {
            System.out.println(s);
        }

        System.out.println(".. and " + uniqueProbes.size() + " probes.");

        HashMap<String, Integer> probeToId = new HashMap<String, Integer>();
        HashMap<Integer, String> idToProbe = new HashMap<Integer, String>();
        int pctr = 0;
        for (String p : uniqueProbes) {
            probeToId.put(p, pctr);
            idToProbe.put(pctr, p);
            pctr++;
        }

        System.out.println("");





        Gpio.createDir(Gpio.getParentDir(output));


        for (String queryPheno : uniquePhenotypes) {
//            if (queryPheno.equals("age") || queryPheno.equals("sex")) {
            TextFile outfile = new TextFile(output + "-" + queryPheno + ".txt", TextFile.W);
            System.out.println("DatasetOrder: " + Strings.concat(desc, Strings.comma));
            outfile.writeln("DatasetOrder: " + Strings.concat(desc, Strings.comma));
            outfile.writeln("Phenotype\tProbe\tGene\tPSpearman\tPPearson\tMetaZSpearman\tMetaZPearson\tNrTotalSamples\tDatasetCorrSpearman\tDatasetCorrPearson\tDatasetZScoresSpearman\tDatasetZScoresPearson\tSampleSizesPerDataset");
            System.out.println(queryPheno);

            int maxNrSamples = maxSamplesPerPhenotype.get(queryPheno);
            Correlation.correlationToZScore(maxNrSamples);
            int[][] samples = new int[uniqueProbes.size()][files.length];
            double[][] correlationSpearman = new double[uniqueProbes.size()][files.length];
            double[][] correlationPearson = new double[uniqueProbes.size()][files.length];

            for (int ds = 0; ds < files.length; ds++) {
                HashMap<String, String> phenotypeTranslationTable = new HashMap<String, String>();
                TextFile phe = new TextFile(phenodescriptions[ds], TextFile.R);
                String[] phelems = phe.readLineElems(TextFile.tab);
                while (phelems != null) {
                    String metapheno = phelems[1].toLowerCase().trim().replaceAll(" ", "");
                    phenotypeTranslationTable.put(phelems[0].toLowerCase().trim(), metapheno);
                    phelems = phe.readLineElems(TextFile.tab);
                }
                phe.close();

                String filename = files[ds];
                TextFile tf = new TextFile(filename, TextFile.R);
                tf.readLine();
                String[] elems = tf.readLineElems(TextFile.tab);
                elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {

                    String pheno = elems[phenoCol].toLowerCase().trim();
                    String metapheno = phenotypeTranslationTable.get(pheno);
                    if (metapheno != null && metapheno.equals(queryPheno)) {
                        String probe = elems[probeCol];
                        Integer probeId = probeToId.get(probe);
                        Integer nrSamples = Integer.parseInt(elems[samplecol]);
                        Double rho = Double.parseDouble(elems[spearmancol]);
                        Double r = Double.parseDouble(elems[pearsoncol]);
                        correlationPearson[probeId][ds] = r;
                        correlationSpearman[probeId][ds] = rho;
                        samples[probeId][ds] = nrSamples;
                    }


                    elems = tf.readLineElems(TextFile.tab);
                }

                tf.close();
            }

            // meta-analysis magic.
            double[][] dsZScoresPearson = new double[uniqueProbes.size()][files.length];
            double[][] dsZScoresSpearman = new double[uniqueProbes.size()][files.length];
            for (int ds = 0; ds < files.length; ds++) {
                for (int probe = 0; probe < dsZScoresPearson.length; probe++) {
                    int nrSamples = samples[probe][ds];

                    if (nrSamples == 0) {
                        dsZScoresPearson[probe][ds] = Double.NaN;
                        dsZScoresSpearman[probe][ds] = Double.NaN;
                    } else {
                        dsZScoresPearson[probe][ds] = Correlation.convertCorrelationToZScore(nrSamples, correlationPearson[probe][ds]);
                        dsZScoresSpearman[probe][ds] = Correlation.convertCorrelationToZScore(nrSamples, correlationSpearman[probe][ds]);
                    }

                }
            }

            for (int probe = 0; probe < dsZScoresPearson.length; probe++) {
                int nrDs = 0;
                int nrSamples = 0;
                for (int q = 0; q < files.length; q++) {
                    if (samples[probe][q] > 0) {
                        nrDs++;
                        nrSamples += samples[probe][q];
                    }
                }



                double[] zP = new double[nrDs];
                double[] zS = new double[nrDs];
                int[] nrS = new int[nrDs];
                int ctr = 0;
                for (int q = 0; q < files.length; q++) {
                    if (samples[probe][q] > 0) {

                        zP[ctr] = dsZScoresPearson[probe][q];
                        zS[ctr] = dsZScoresSpearman[probe][q];
                        nrS[ctr] = samples[probe][q];
                        ctr++;
                    }
                }

                double wzPearson = ZScores.getWeightedZ(zP, nrS);
                double wzSpearman = ZScores.getWeightedZ(zS, nrS);
                double pPearson = ZScores.zToP(wzPearson);
                double pSpearman = ZScores.zToP(wzSpearman);

                String dsCorrSpearmanStr = Strings.concat(correlationSpearman[probe], Strings.comma);
                String dsCorrPearsonStr = Strings.concat(correlationPearson[probe], Strings.comma);
                String dsZScoresSpearmanStr = Strings.concat(dsZScoresSpearman[probe], Strings.comma);
                String dsZScoresPearsonStr = Strings.concat(dsZScoresPearson[probe], Strings.comma);
                String dsSamplesStr = Strings.concat(samples[probe], Strings.comma);

                String probeName = idToProbe.get(probe);
                outfile.writeln(queryPheno + "\t" + probeName + "\t" + probeToGeneMap.get(probeName) + "\t" + pSpearman + "\t" + pPearson + "\t" + wzSpearman + "\t" + wzPearson + "\t" + nrSamples + "\t" + dsCorrSpearmanStr + "\t" + dsCorrPearsonStr + "\t" + dsZScoresSpearmanStr + "\t" + dsZScoresPearsonStr + "\t" + dsSamplesStr);
            }
            outfile.close();
//            }

        }


        return uniquePhenotypes.toArray(new String[0]);
    }

    // parse files with header:
    // Phenotype	Probe	Gene	PSpearman	PPearson	MetaZSpearman	MetaZPearson	NrTotalSamples	DatasetCorrSpearman	DatasetCorrPearson	DatasetZScoresSpearman	DatasetZScoresPearson	SampleSizesPerDataset
    public void combine(String dir, String[] prefixes, String[] phenotypes, String out, String eQTLFile, String gwasCatalog) throws IOException {

        GWASCatalog gw = new GWASCatalog();
        gw.read(gwasCatalog);

        HashMap<String, HashSet<String>> probeToSNP = new HashMap<String, HashSet<String>>();
        HashMap<String, Double> bestFDR = new HashMap<String, Double>();
        HashMap<String, Double> eqtlzscores = new HashMap<String, Double>();
        HashSet<String> snps = new HashSet<String>();

        TextFile t = new TextFile(eQTLFile, TextFile.R);
        String[] header = t.readLineElems(TextFile.tab);
        int fdrcol = -1;
        for (int q = 0; q < header.length; q++) {
            if (header[q].contains("FDR")) {
                fdrcol = q;
            }
        }

        System.out.println("FDR COLUMN: " + fdrcol);
        String[] telems = t.readLineElems(TextFile.tab);
        while (telems != null) {
            String snp = telems[1];
            snps.add(snp);
            String probe = telems[4];

            Double fdr = 0d;
            try {
                Double.parseDouble(telems[fdrcol]);
            } catch (NumberFormatException e) {
            }

            bestFDR.put(snp + "-" + probe, fdr);


            Double metaZ = Double.parseDouble(telems[eQTLTextFile.METAZ]);
            eqtlzscores.put(snp + "-" + probe, metaZ);

            HashSet<String> psnps = probeToSNP.get(probe);
            if (psnps == null) {
                psnps = new HashSet<String>();
            }
            psnps.add(snp);
            probeToSNP.put(probe, psnps);

            telems = t.readLineElems(TextFile.tab);
        }
        t.close();

        HashMap<String, String> probeToGene = new HashMap<String, String>();
        HashSet<String> uniqueProbes = new HashSet<String>();

        for (String prefix : prefixes) {
            for (String phenotype : phenotypes) {
                String file = dir + prefix + "-" + phenotype + ".txt";
                TextFile tf = new TextFile(file, TextFile.R);
                tf.readLine();// read pre-header
                tf.readLine();// read header
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {
                    uniqueProbes.add(elems[1]);
                    probeToGene.put(elems[1], elems[2]);
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
            }
        }

        System.out.println("Detected: " + uniqueProbes.size() + " probes.");

        double[][][] pvalues = new double[prefixes.length][phenotypes.length][uniqueProbes.size()];
        double[][][] zscores = new double[prefixes.length][phenotypes.length][uniqueProbes.size()];

        int[][][] samples = new int[prefixes.length][phenotypes.length][uniqueProbes.size()];

        HashMap<String, Integer> probeToIndex = new HashMap<String, Integer>();
        int pbctr = 0;
        for (String probe : uniqueProbes) {
            probeToIndex.put(probe, pbctr);
            pbctr++;
        }


        for (int p = 0; p < prefixes.length; p++) {
            String prefix = prefixes[p];
            for (int q = 0; q < phenotypes.length; q++) {
                String phenotype = phenotypes[q];

                String file = dir + prefix + "-" + phenotype + ".txt";
                TextFile tf = new TextFile(file, TextFile.R);
                tf.readLine();// read pre-header
                tf.readLine();// read header
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {

                    double pvalue = Double.parseDouble(elems[3]);
                    double zscore = Double.parseDouble(elems[5]);
                    Integer probeId = probeToIndex.get(elems[1]);
                    int sampleSize = Integer.parseInt(elems[7]);
                    pvalues[p][q][probeId] = pvalue;
                    zscores[p][q][probeId] = zscore;
                    samples[p][q][probeId] = sampleSize;
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();

            }
        }

        TextFile outFile = new TextFile(out, TextFile.W);
        int[] avgsamples = new int[phenotypes.length];
        String outheader = "Probe\tGene\tSNP(s)\teQTLZScores\tFDR\tTraits";
        for (int q = 0; q < phenotypes.length; q++) {
            for (int p = 0; p < prefixes.length; p++) {
                int nrSamplesTotal = 0;
                int nrProbesWData = 0;
                for (int probe = 0; probe < samples[p][q].length; probe++) {
                    int nr = samples[p][q][probe];
                    if (nr > 0) {
                        nrSamplesTotal += nr;
                        nrProbesWData++;
                    }
                }
                int avSamples = nrSamplesTotal / nrProbesWData;
                avgsamples[q] = avSamples;
                if (avSamples > minNrSamples) {
                    String prefix = prefixes[p];
                    String phenotype = phenotypes[q];
                    outheader += "\tP-" + prefix + "-" + phenotype;
                    outheader += "\tZ-" + prefix + "-" + phenotype;
                    outheader += "\tN-" + prefix + "-" + phenotype;
                }
            }
        }
        outFile.writeln(outheader);



        for (String probe : uniqueProbes) {
            Integer probeId = probeToIndex.get(probe);
            HashSet<String> psnps = probeToSNP.get(probe);
            if (psnps != null) {
                String[] snparr = psnps.toArray(new String[0]);
                String snpStr = Strings.concat(snparr, Strings.semicolon);
                HashSet<GWASTrait> traits = new HashSet<GWASTrait>();
                String traitStr = "";
                String eqtlStr = "";
                String fdrStr = "";
                for (String snp : snparr) {
                    GWASSNP snpObj = gw.getSnpToObj().get(snp);
                    GWASTrait[] traitsForSNP = snpObj.getAssociatedTraits().toArray(new GWASTrait[0]);
                    traitStr += Strings.concat(traitsForSNP, Strings.comma) + ";";

                    Double eqtl = eqtlzscores.get(snp + "-" + probe);
                    Double fdr = bestFDR.get(snp + "-" + probe);
                    eqtlStr += eqtl + ";";
                    fdrStr += fdr + ";";
                }

                String output = probe + "\t" + probeToGene.get(probe) + "\t" + snpStr + "\t" + eqtlStr + "\t" + fdrStr + "\t" + traitStr;
                for (int q = 0; q < phenotypes.length; q++) {
                    int avSamples = avgsamples[q];
                    if (avSamples > minNrSamples) {
                        for (int p = 0; p < prefixes.length; p++) {
                            output += "\t" + pvalues[p][q][probeId];
                            output += "\t" + zscores[p][q][probeId];
                            output += "\t" + samples[p][q][probeId];
                        }
                    }
                }
                outFile.writeln(output);
            }

        }
        outFile.close();

        System.out.println("Done combining.. Now running eQTL correlations");

        // correlate eQTL ZScores with endophenotype zscores
        outFile = new TextFile(out + "CorrelationsWithEQTLZScoresAllVsAll.txt", TextFile.W);
        TextFile toutFile = new TextFile(out + "TTestZScoresWilCoxon.txt", TextFile.W);


        outheader = "SNP\tNrProbes\tTraits";
        String toutheader = "SNP\tNrProbes\tTraits\tNrPos\tNrNeg";
        // write header

        for (int q = 0; q < phenotypes.length; q++) {
            for (int p = 0; p < prefixes.length; p++) {
                int avSamples = avgsamples[q];
                if (avSamples > minNrSamples) {
                    String prefix = prefixes[p];
                    String phenotype = phenotypes[q];
                    outheader += "\tR-" + prefix + "-" + phenotype;
                    toutheader += "\tPval-" + prefix + "-" + phenotype;
                }
            }

        }
        outFile.writeln(outheader);
        toutFile.writeln(toutheader);


        // snps vs probes
        for (String snp : snps) {
            ArrayList<String> probes = new ArrayList<String>();
            ArrayList<String> posprobes = new ArrayList<String>();
            ArrayList<String> negprobes = new ArrayList<String>();
            for (String probe : uniqueProbes) {
                Double ez = eqtlzscores.get(snp + "-" + probe);
                if (ez != null) {
                    if (ez >= 0) {
                        posprobes.add(probe);
                    } else {
                        negprobes.add(probe);
                    }
                    probes.add(probe);
                }
            }
            if (probes.size() > 5) {

                // correlate
                GWASSNP snpObj = gw.getSnpToObj().get(snp);
                GWASTrait[] traitsForSNP = snpObj.getAssociatedTraits().toArray(new GWASTrait[0]);
                String traitStr = Strings.concat(traitsForSNP, Strings.comma);
                String outln = snp + "\t" + probes.size() + "\t" + traitStr;
                for (int q = 0; q < phenotypes.length; q++) {

                    int avSamples = avgsamples[q];
                    if (avSamples > minNrSamples) {
                        for (int p = 0; p < prefixes.length; p++) {
                            double[] snpeqtlzscores = new double[probes.size()];
                            double[] snpendozscores = new double[probes.size()];
                            for (int pid = 0; pid < probes.size(); pid++) {
                                int probeId = probeToIndex.get(probes.get(pid));
                                snpeqtlzscores[pid] = eqtlzscores.get(snp + "-" + probes.get(pid));
                                snpendozscores[pid] = zscores[p][q][probeId];
                            }
                            double r = JSci.maths.ArrayMath.correlation(snpeqtlzscores, snpendozscores);
                            outln += "\t" + r;
                        }
                    }
                }
                outFile.writeln(outln);

                // ttest
                if (posprobes.size() > 5 && negprobes.size() > 5) {
                    outln = snp + "\t" + probes.size() + "\t" + traitStr + "\t" + posprobes.size() + "\t" + negprobes.size();
                    for (int q = 0; q < phenotypes.length; q++) {
                        int avSamples = avgsamples[q];
                        if (avSamples > minNrSamples) {

                            for (int p = 0; p < prefixes.length; p++) {

                                double[] snpendozscorespos = new double[probes.size()];
                                double[] snpendozscoresneg = new double[probes.size()];
                                for (int pid = 0; pid < posprobes.size(); pid++) {
                                    int probeId = probeToIndex.get(posprobes.get(pid));
                                    snpendozscorespos[pid] = zscores[p][q][probeId];
                                }
                                for (int pid = 0; pid < negprobes.size(); pid++) {
                                    int probeId = probeToIndex.get(negprobes.get(pid));
                                    snpendozscoresneg[pid] = zscores[p][q][probeId];
                                }

//                                double pval = TTest.test(snpendozscorespos, snpendozscoresneg);

                                WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();
                                double pval = mwm.returnWilcoxonMannWhitneyPValue(snpendozscorespos, snpendozscoresneg);
                                outln += "\t" + pval;
                            }

                        }
                    }
                    toutFile.writeln(outln);
                }
            }
        }

        toutFile.close();
        outFile.close();




    }
}
