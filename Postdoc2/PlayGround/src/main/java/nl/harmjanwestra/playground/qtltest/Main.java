package nl.harmjanwestra.playground.qtltest;

import nl.harmjanwestra.playground.legacy.vcf.VCFTabix;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.HWE;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;
import umcg.genetica.util.RankArray;
import umcg.genetica.util.RunTimer;
import umontreal.iro.lecuyer.probdist.BetaDist;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class Main {

    HashMap<String, String> RNAToDNA = new HashMap<>();
    HashMap<String, String> DNAtoRNA = new HashMap<>();
    HashMap<String, String> RNAToDataset = new HashMap<>();
    HashMap<Object, Object> RNAsPerDataset;

    double mafthreshold = 0.01;
    double callratethreshold = 0.95;
    double hwepthreshold = 0.0001;

    public static void main(String[] args) {

        Main main = new Main();
        if (args.length < 7) {
            System.out.println("Usage:");
            System.out.println("vcffile chromosome linkfile genelimitfile geneexpressionfile geneannotationfile outfile");
        } else {
            String vcfile = args[0];
            int chromosome = Integer.parseInt(args[1]);
            String linkfile = args[2];
            String genelimitfile = args[3];
            String geneexpressionfile = args[4];
            String geneannotationfile = args[5];
            String outfile = args[6];
            try {
                main.validate(vcfile, chromosome, linkfile, genelimitfile, geneexpressionfile, geneannotationfile, outfile);

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void fastqtlclone(String vcfFile,
                             int chromosome,
                             String linkfile,
                             String geneLimitFile,
                             String geneExpressionDataFile,
                             String geneAnnotationFile,
                             String outfile) throws IOException {


        // load genotyped samples
        ArrayList<String> genotypeSamples = getGenotypeSamples(vcfFile);

        HashMap<String, Integer> genotypeSampleHash = Util.hash(genotypeSamples);

        // determine genotype samples with links to an RNA-seq sample
        // load sample links
        loadSampleLinks(linkfile, genotypeSampleHash.keySet());

        // determine RNA samples to load
        HashSet<String> RNASamplesMatchedToDNA = new HashSet<String>();
        for (int i = 0; i < genotypeSamples.size(); i++) {
            String RNASample = DNAtoRNA.get(genotypeSamples.get(i));
            if (RNASample != null) {
                RNASamplesMatchedToDNA.add(RNASample);
            }
        }

        // load set of genes to limit to
        Set<String> geneLimitSet = null;
        if (geneLimitFile != null) {
            System.out.println("Limiting to genes found in : " + geneLimitFile);
            TextFile tf2 = new TextFile(geneLimitFile, TextFile.R);
            geneLimitSet = tf2.readAsSet(0, TextFile.tab);
            tf2.close();
            System.out.println("File contains " + geneLimitSet.size() + " genes.");
        }

        // load gene annotation
        GeneAnnotation geneAnnotation = new GeneAnnotation(geneAnnotationFile, chromosome, geneLimitSet);

        // load expression data
        GeneExpressionData expressionData = new GeneExpressionData(geneExpressionDataFile, geneAnnotation.getAllGenes(), RNASamplesMatchedToDNA);

        // determine which genotype samples to load using the available RNA samples
        HashSet<String> availableExpressionSamples = new HashSet<>();
        for (String s : expressionData.samples) {
            availableExpressionSamples.add(s);
        }

        // make a list of columns to load from the VCF file
        boolean[] genotypeSamplesToInclude = new boolean[genotypeSamples.size()];
        ArrayList<String> loadedGenotypeSamples = new ArrayList<>();
        for (int i = 0; i < genotypeSamples.size(); i++) {
            String DNA = genotypeSamples.get(i);
            String RNA = DNAtoRNA.get(DNA);
            if (availableExpressionSamples.contains(RNA)) {
                genotypeSamplesToInclude[i] = true;
                loadedGenotypeSamples.add(DNA);
            }
        }

        // replace genotype sample list, and rehash
        genotypeSamples = loadedGenotypeSamples;
        genotypeSampleHash = Util.hash(genotypeSamples);
        System.out.println(genotypeSamples.size() + " DNA samples match loaded RNA samples");

        // divide samples up into datasets
        Dataset[] datasets = createDatasets(genotypeSampleHash, expressionData.sampleMap);

        Dataset combinedDataset = createCombinedDataset(genotypeSampleHash, expressionData.sampleMap);
        System.out.println("Combined dataset has: " + combinedDataset.expressionIds.length + " samples.");

        // initialize Z-score lookup table
        Correlation.correlationToZScore(combinedDataset.expressionIds.length);

        // iterate genes
        Chromosome chromosomeObj = Chromosome.parseChr("" + chromosome);
        System.out.println("Processing: " + vcfFile);
        VCFTabix tabix = new VCFTabix(vcfFile);

        TextFile output = new TextFile(outfile + "-TopFxPerGene.txt", TextFile.W);
        TextFile outputAll = new TextFile(outfile + "-AllAssociations.txt.gz", TextFile.W);
        TextFile outputsnplog = new TextFile(outfile + "SNPLog.txt.gz", TextFile.W);
        output.writeln("Gene\tSNP\tAlleles\tEffectallele\tMetaN\tMetaZ\tMetaP\tBDistAlpha\tBDistBeta\tAdjP");
        outputAll.writeln("Gene\tSNP\tAlleles\tEffectallele\tMetaN\tMetaZ\tMetaP");
        String snplogheader = "SNP\tDatasetsPassQC\tNTotal\tNJoint1\tMAFJoint1\tCallRateJoint1\tHWEPJoint1\tNJoint1\tMAFJoint1\tCallRateJoint1\tHWEPJoint1";
        HashSet<String> seenSNP = new HashSet<>();
        for (int d = 0; d < datasets.length; d++) {
            String name = datasets[d].name;
            snplogheader += "\t" + name + "-MAF";
            snplogheader += "\t" + name + "-CallRate";
            snplogheader += "\t" + name + "-HWEP";
        }
        outputsnplog.writeln(snplogheader);

        long[] seed = new long[1000];
        Random rand = new Random(123456789);
        for (int i = 0; i < seed.length; i++) {
            seed[i] = rand.nextLong();
        }

        // iterate genes
        RunTimer timer = new RunTimer();
        timer.start();
        RankArray ranker = new RankArray();
        for (int g = 0; g < expressionData.genes.length; g++) {
            String gene = expressionData.genes[g];

            Integer geneAnnotationId = geneAnnotation.getGeneId(gene);
            if (geneAnnotationId != null) {
                double[] expData = expressionData.data[g];

                // rank once, then center+scale
                double[][] expDataPerDatasetRanked = new double[datasets.length][];
                IntStream.range(0, datasets.length).parallel().forEach(d -> {
                    Dataset thisDataset = datasets[d];
                    double[] datasetExpressionData = thisDataset.select(expData, thisDataset.expressionIds);
                    expDataPerDatasetRanked[d] = ranker.rank(datasetExpressionData, true);
                });

                // get variants, 1mb up and downstream
                int pos = geneAnnotation.getPos(geneAnnotationId);
                int start = pos - 1000001;
                if (start < 0) {
                    start = 0;
                }
                int stop = pos + 1000001;

                Feature region = new Feature(chromosomeObj, start, stop);
                Iterator<VCFVariant> snpIterator = tabix.getVariants(region, genotypeSamplesToInclude);
                int vctr = 0;
                int tctr = 0;
                boolean print = false;

                double[] lowestPermP = new double[1000];
                Arrays.fill(lowestPermP, 1);

                String topsnp = null;
                double topsnpp = 1;
                double topsnpz = 0;
                double topsnpn = 0;
                String topSNPAlleles = null;
                String topSNPEffectAllele = null;


                while (snpIterator.hasNext()) {
                    VCFVariant variant = snpIterator.next();

                    if (variant != null) {
                        String snpid = variant.getId();
//                        if (snpid.equals("22:23987520:rs140245:A_G")) {
//                            System.out.println();
//                            System.out.println("-------");
//                            print = true;
//                        } else {
//                            print = false;
//                        }
                        // collect data per dataset
                        double[][] datasetsGenotypeData = new double[datasets.length][];
                        double[][] datasetsGenotypeDosageData = new double[datasets.length][];
                        QCObj[] datasetsQCObjs = new QCObj[datasets.length];
                        AtomicInteger datasetsPassQC = new AtomicInteger();

                        // System.out.println(variant.getId() + "\t" + gene + "\t" + Strings.concat(variant.getAlleles(), Strings.backwardslash));
                        IntStream.range(0, datasets.length).parallel().forEach(d -> {
                            Dataset thisDataset = datasets[d];

                            double[] datasetGenotypeData = thisDataset.select(variant.getGenotypesAsByteVector(), thisDataset.genotypeIds);

                            // do some QC checks? test MAF, Callrate, HWE-P
                            QCObj qcobj = checkVariant(datasetGenotypeData);
                            datasetsQCObjs[d] = qcobj;
                            datasetsGenotypeData[d] = datasetGenotypeData;
//                            double[] datasetDosageValues = thisDataset.select(variant.getGenotypeDosage(), thisDataset.genotypeIds);
                            double[] datasetDosageValues = thisDataset.select(variant.getDosage(), thisDataset.genotypeIds);
                            datasetsGenotypeDosageData[d] = datasetDosageValues;

//                            System.out.println(variant.getGenotypesAsByteVector().length);
                            if (print) {
                                System.out.println(thisDataset.name + "\t" + datasetGenotypeData.length + "\t" + qcobj.toString());
                            }
//                            for (int q = 0; q < datasetGenotypeData.length; q++) {
//                                System.out.println(datasetGenotypeData[q]);
//                            }
//                            System.exit(-1);


                            if (qcobj.passqc) {
                                datasetsPassQC.getAndIncrement();
                            }
                        });


                        if (datasetsPassQC.get() > 1) {
                            int totalN = 0;
                            // perform different types of analysis:

                            ////////////////////////////////
                            // 1. EMP style meta-analysis //
                            ////////////////////////////////
                            double[] zScores = new double[datasets.length];
                            double[][] zScoresPerm = new double[1000][datasets.length];
                            int[] datasetSampleSizes = new int[datasets.length];

                            int dctr = 0;
                            for (int d = 0; d < datasets.length; d++) {
                                if (datasetsQCObjs[d].passqc) {
                                    // prune the dataset; remove missing values

                                    Pair<double[], double[]> prunedDatasetData = pruneMissing(datasetsGenotypeData[d], datasetsGenotypeDosageData[d], expDataPerDatasetRanked[d]);
                                    if (prunedDatasetData.getRight().length == 0 || prunedDatasetData.getLeft().length == 0) {
                                        System.out.println(datasets[d].name + " has 0 values, but " + datasets[d].expressionIds.length + " samples");
                                        for (int v = 0; v < datasetsGenotypeData[d].length; v++) {
                                            System.out.println(datasetsGenotypeData[d][v] + "\t" + expDataPerDatasetRanked[d][v]);
                                        }
                                        System.out.println(datasetsQCObjs[d].toString());
                                        System.exit(-1);
                                    }
                                    double[] expDataCenterScale = centerScale(prunedDatasetData.getRight());
                                    double[] gtDataCenterScale = centerScale(prunedDatasetData.getLeft());


                                    // perform correlation
                                    double r = Correlation.correlate(expDataCenterScale, gtDataCenterScale);
                                    // convert correlation to z-score
                                    int n = gtDataCenterScale.length;
                                    datasetSampleSizes[d] = n;
                                    totalN += n;
                                    double z = Correlation.convertCorrelationToZScore(n, r);
                                    zScores[d] = z;

                                    int finalD = d;
                                    IntStream.range(0, 1000).parallel().forEach(p -> {
                                        double[] exp = new double[expDataCenterScale.length];
                                        System.arraycopy(exp, 0, expDataCenterScale, 0, expDataCenterScale.length);
                                        shuffleArray(exp, seed[p]);
                                        double rp = Correlation.correlate(exp, gtDataCenterScale);
                                        double zp = Correlation.convertCorrelationToZScore(n, rp);
                                        zScoresPerm[p][finalD] = zp;
                                    });
//                                    System.out.println(datasets[d].name + "\t" + n + "\t" + r + "\t" + z + "\t" + datasetsQCObjs[d]);
                                    dctr++;
                                } else {
                                    zScores[d] = Double.NaN;

                                    int finalD1 = d;
                                    IntStream.range(0, 1000).parallel().forEach(p -> {
                                        zScoresPerm[p][finalD1] = Double.NaN;
                                    });
//                                    System.out.println(datasets[d].name + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + datasetsQCObjs[d]);
                                }
                                if (print) {
                                    System.out.println(datasets[d].name + "\t" + zScores[d]);
                                }
                            }

                            // 1.2 meta-analyze z-scores
                            double metaZ = ZScores.getWeightedZ(zScores, datasetSampleSizes);
                            double metaP = ZScores.zToP(metaZ);
                            String[] alleles = variant.getAlleles();
                            String alleleStr = alleles[0] + "/" + alleles[1];
                            if (metaP < topsnpp) {
                                topsnpp = metaP;
                                topsnpz = metaZ;
                                topsnpn = totalN;
                                topsnp = snpid;
                                topSNPAlleles = alleleStr;
                                topSNPEffectAllele = alleles[1];
                            }

                            // meta-analyze permuted p-values
                            IntStream.range(0, 1000).parallel().forEach(p -> {
                                double metaZP = ZScores.getWeightedZ(zScoresPerm[p], datasetSampleSizes);
                                double metaPP = ZScores.zToP(metaZP);
                                if (metaPP < lowestPermP[p]) {
                                    lowestPermP[p] = metaPP;
                                }
                            });


                            outputAll.writeln(gene + "\t" + snpid + "\t" + alleleStr + "\t" + alleles[1] + "\t" + totalN + "\t" + metaZ + "\t" + metaP);

                        }


                    }
                    vctr++;
                    if (vctr % 1000 == 0) {
                        System.out.print("Gene:" + g + "/" + expressionData.genes.length + "; " + vctr + " variants loaded, " + tctr + " variants tested. T: " + timer.getTimeDesc() + "\r");
                    }
                }

                BetaDist bdist = BetaDist.getInstanceFromMLE(lowestPermP, lowestPermP.length);
                double pperm = bdist.cdf(topsnpp);
                double alpha = bdist.getAlpha();
                double beta = bdist.getBeta();

                // output somewhere
                // // "Gene\tSNP\tAlleles\tEffectallele\tMetaN\tMetaZ\tMetaP\tBDistAlpha\tBDistBeta\tAdjP"
                output.writeln(gene + "\t" + topsnp + "\t" + topSNPAlleles + "\t" + topSNPEffectAllele + "\t" + topsnpn + "\t" + topsnpz + "\t" + topsnpp + "\t" + alpha + "\t" + beta + "\t" + pperm);

                System.out.print("Gene:" + g + "/" + expressionData.genes.length + "; " + vctr + " variants loaded, " + tctr + " variants tested. T: " + timer.getTimeDesc() + "\n");
            }
        }
        outputsnplog.close();
        output.close();
        outputAll.close();
        tabix.close();

    }

    private void shuffleArray(double[] ar, long seed) {
        Random rnd = new Random(seed);
        for (int i = ar.length - 1; i >= 0; i--) {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            double a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }

    public void validate(String vcfFile,
                         int chromosome,
                         String linkfile,
                         String geneLimitFile,
                         String geneExpressionDataFile,
                         String geneAnnotationFile,
                         String outfile) throws IOException {


        // load genotyped samples
        ArrayList<String> genotypeSamples = getGenotypeSamples(vcfFile);

        HashMap<String, Integer> genotypeSampleHash = Util.hash(genotypeSamples);

        // determine genotype samples with links to an RNA-seq sample
        // load sample links
        loadSampleLinks(linkfile, genotypeSampleHash.keySet());

        // determine RNA samples to load
        HashSet<String> RNASamplesMatchedToDNA = new HashSet<String>();
        for (int i = 0; i < genotypeSamples.size(); i++) {
            String RNASample = DNAtoRNA.get(genotypeSamples.get(i));
            if (RNASample != null) {
                RNASamplesMatchedToDNA.add(RNASample);
            }
        }

        // load set of genes to limit to
        Set<String> geneLimitSet = null;
        if (geneLimitFile != null) {
            System.out.println("Limiting to genes found in : " + geneLimitFile);
            TextFile tf2 = new TextFile(geneLimitFile, TextFile.R);
            geneLimitSet = tf2.readAsSet(0, TextFile.tab);
            tf2.close();
            System.out.println("File contains " + geneLimitSet.size() + " genes.");
        }

        // load gene annotation
        GeneAnnotation geneAnnotation = new GeneAnnotation(geneAnnotationFile, chromosome, geneLimitSet);

        // load expression data
        GeneExpressionData expressionData = new GeneExpressionData(geneExpressionDataFile, geneAnnotation.getAllGenes(), RNASamplesMatchedToDNA);

        // determine which genotype samples to load using the available RNA samples
        HashSet<String> availableExpressionSamples = new HashSet<>();
        for (String s : expressionData.samples) {
            availableExpressionSamples.add(s);
        }

        // make a list of columns to load from the VCF file
        boolean[] genotypeSamplesToInclude = new boolean[genotypeSamples.size()];
        ArrayList<String> loadedGenotypeSamples = new ArrayList<>();
        for (int i = 0; i < genotypeSamples.size(); i++) {
            String DNA = genotypeSamples.get(i);
            String RNA = DNAtoRNA.get(DNA);
            if (availableExpressionSamples.contains(RNA)) {
                genotypeSamplesToInclude[i] = true;
                loadedGenotypeSamples.add(DNA);
            }
        }

        // replace genotype sample list, and rehash
        genotypeSamples = loadedGenotypeSamples;
        genotypeSampleHash = Util.hash(genotypeSamples);
        System.out.println(genotypeSamples.size() + " DNA samples match loaded RNA samples");

        // divide samples up into datasets
        Dataset[] datasets = createDatasets(genotypeSampleHash, expressionData.sampleMap);

        Dataset combinedDataset = createCombinedDataset(genotypeSampleHash, expressionData.sampleMap);
        System.out.println("Combined dataset has: " + combinedDataset.expressionIds.length + " samples.");

        // initialize Z-score lookup table
        Correlation.correlationToZScore(combinedDataset.expressionIds.length);

        // iterate genes
        Chromosome chromosomeObj = Chromosome.parseChr("" + chromosome);
        System.out.println("Processing: " + vcfFile);
        VCFTabix tabix = new VCFTabix(vcfFile);

        TextFile output = new TextFile(outfile, TextFile.W);
        TextFile outputsnplog = new TextFile(outfile + "SNPLog.txt.gz", TextFile.W);
        output.writeln("Gene\tSNP\tAlleles\tEffectallele\t" +
                "MetaN\tMetaZ\tMetaP\t" +
                "MetaZPreRanked\tMetaPMPreRanked\t" +
                "JointN1\tJointZ1\tJointP1\t" +
                "JointN2\tJointZ2\tJointP2");

        String snplogheader = "SNP\tDatasetsPassQC\tNTotal\tNJoint1\tMAFJoint1\tCallRateJoint1\tHWEPJoint1\tNJoint1\tMAFJoint1\tCallRateJoint1\tHWEPJoint1";
        HashSet<String> seenSNP = new HashSet<>();
        for (int d = 0; d < datasets.length; d++) {
            String name = datasets[d].name;
            snplogheader += "\t" + name + "-MAF";
            snplogheader += "\t" + name + "-CallRate";
            snplogheader += "\t" + name + "-HWEP";
        }
        outputsnplog.writeln(snplogheader);

        // iterate genes
        RunTimer timer = new RunTimer();
        timer.start();
        RankArray ranker = new RankArray();
        for (int g = 0; g < expressionData.genes.length; g++) {
            String gene = expressionData.genes[g];

            Integer geneAnnotationId = geneAnnotation.getGeneId(gene);
            if (geneAnnotationId != null) {
                double[] expData = expressionData.data[g];

                // rank once, then center+scale
                double[][] expDataPerDatasetRanked = new double[datasets.length][];
                double[][] expDataPerDatasetRankedCenterScaled = new double[datasets.length][];

                for (int d = 0; d < datasets.length; d++) {
                    Dataset thisDataset = datasets[d];
                    double[] datasetExpressionData = thisDataset.select(expData, thisDataset.expressionIds);
                    expDataPerDatasetRanked[d] = ranker.rank(datasetExpressionData, true);
                    expDataPerDatasetRankedCenterScaled[d] = centerScale(expDataPerDatasetRanked[d]);
                }

                // get variants, 1mb up and downstream
                int pos = geneAnnotation.getPos(geneAnnotationId);
                int start = pos - 1000001;
                if (start < 0) {
                    start = 0;
                }
                int stop = pos + 1000001;
                Feature region = new Feature(chromosomeObj, start, stop);
                Iterator<VCFVariant> snpIterator = tabix.getVariants(region, genotypeSamplesToInclude);
                int vctr = 0;
                int tctr = 0;
                boolean print = false;
                while (snpIterator.hasNext()) {
                    VCFVariant variant = snpIterator.next();

                    if (variant != null) {
                        String snpid = variant.getId();
                        if (snpid.equals("22:18908190:rs424512:T_C")) {
                            System.out.println();
                            System.out.println("-------");
                            print = true;
                        } else {
                            print = false;
                        }
                        // collect data per dataset
                        double[][] datasetsGenotypeData = new double[datasets.length][];
                        double[][] datasetsGenotypeDosageData = new double[datasets.length][];
                        QCObj[] datasetsQCObjs = new QCObj[datasets.length];
                        int datasetsPassQC = 0;

//                        ArrayList<Double> genotypeDataForJointAnalysis = new ArrayList<>();
//                        ArrayList<Double> genotypeDosageDataForJointAnalysis = new ArrayList<>();
//                        ArrayList<Double> expressionDataForJointAnalysis = new ArrayList<>();
//
//                        ArrayList<Double> genotypeDataForJointAnalysisDatasetsPassingQC = new ArrayList<>();
//                        ArrayList<Double> genotypeDosageDataForJointAnalysisDatasetsPassingQC = new ArrayList<>();
//                        ArrayList<Double> expressionDataForJointAnalysisDatasetsPassingQC = new ArrayList<>();

                        // System.out.println(variant.getId() + "\t" + gene + "\t" + Strings.concat(variant.getAlleles(), Strings.backwardslash));
                        for (int d = 0; d < datasets.length; d++) {
                            Dataset thisDataset = datasets[d];

                            double[] datasetGenotypeData = thisDataset.select(variant.getGenotypesAsByteVector(), thisDataset.genotypeIds);

                            // do some QC checks? test MAF, Callrate, HWE-P
                            QCObj qcobj = checkVariant(datasetGenotypeData);
                            datasetsQCObjs[d] = qcobj;
                            datasetsGenotypeData[d] = datasetGenotypeData;
//                            double[] datasetDosageValues = thisDataset.select(variant.getGenotypeDosage(), thisDataset.genotypeIds);
                            double[] datasetDosageValues = thisDataset.select(variant.getDosage(), thisDataset.genotypeIds);
                            datasetsGenotypeDosageData[d] = datasetDosageValues;

//                            System.out.println(variant.getGenotypesAsByteVector().length);
                            if (print) {
                                System.out.println(thisDataset.name + "\t" + datasetGenotypeData.length + "\t" + qcobj.toString());
                            }
//                            for (int q = 0; q < datasetGenotypeData.length; q++) {
//                                System.out.println(datasetGenotypeData[q]);
//                            }
//                            System.exit(-1);

                            // System.out.println(datasets[d].name + "\t" + qcobj.toString());
                            for (int v = 0; v < datasetGenotypeData.length; v++) {
//                                genotypeDataForJointAnalysis.add(datasetGenotypeData[v]);
//                                genotypeDosageDataForJointAnalysis.add(datasetDosageValues[v]);
//                                expressionDataForJointAnalysis.add(expDataPerDatasetRankedCenterScaled[d][v]);
                                if (qcobj.passqc) {
//                                    genotypeDataForJointAnalysisDatasetsPassingQC.add(datasetGenotypeData[v]);
//                                    genotypeDosageDataForJointAnalysisDatasetsPassingQC.add(datasetDosageValues[v]);
//                                    expressionDataForJointAnalysisDatasetsPassingQC.add(expDataPerDatasetRankedCenterScaled[d][v]);
                                }
                            }
                            if (qcobj.passqc) {
                                datasetsPassQC++;
                            }
                        }

                        if (datasetsPassQC > 1) {
                            int totalN = 0;
                            // perform different types of analysis:

                            ////////////////////////////////
                            // 1. EMP style meta-analysis //
                            ////////////////////////////////
                            double[] zScores = new double[datasets.length];
                            int[] datasetSampleSizes = new int[datasets.length];

                            int dctr = 0;
                            for (int d = 0; d < datasets.length; d++) {
                                if (datasetsQCObjs[d].passqc) {
                                    // prune the dataset; remove missing values

                                    Pair<double[], double[]> prunedDatasetData = pruneMissing(datasetsGenotypeData[d], datasetsGenotypeDosageData[d], expDataPerDatasetRanked[d]);
                                    if (prunedDatasetData.getRight().length == 0 || prunedDatasetData.getLeft().length == 0) {
                                        System.out.println(datasets[d].name + " has 0 values, but " + datasets[d].expressionIds.length + " samples");
                                        for (int v = 0; v < datasetsGenotypeData[d].length; v++) {
                                            System.out.println(datasetsGenotypeData[d][v] + "\t" + expDataPerDatasetRanked[d][v]);
                                        }
                                        System.out.println(datasetsQCObjs[d].toString());
                                        System.exit(-1);
                                    }
                                    double[] expDataCenterScale = centerScale(prunedDatasetData.getRight());
                                    double[] gtDataCenterScale = centerScale(prunedDatasetData.getLeft());


                                    // perform correlation
                                    double r = Correlation.correlate(expDataCenterScale, gtDataCenterScale);
                                    // convert correlation to z-score
                                    int n = gtDataCenterScale.length;
                                    datasetSampleSizes[d] = n;
                                    totalN += n;
                                    double z = Correlation.convertCorrelationToZScore(n, r);
                                    zScores[d] = z;
//                                    System.out.println(datasets[d].name + "\t" + n + "\t" + r + "\t" + z + "\t" + datasetsQCObjs[d]);
                                    dctr++;
                                } else {
                                    zScores[d] = Double.NaN;
//                                    System.out.println(datasets[d].name + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\t" + datasetsQCObjs[d]);
                                }
                                if (print) {
                                    System.out.println(datasets[d].name + "\t" + zScores[d]);
                                }
                            }


                            // 1.2 meta-analyze z-scores
                            double metaZ = ZScores.getWeightedZ(zScores, datasetSampleSizes);
                            double metaP = ZScores.zToP(metaZ);
                            if (print) {
                                System.out.println(metaZ + "\t" + metaP);
                                System.out.println("-------");
//                                System.exit(0);
                            }
//                            System.out.println(metaZ);
//                            System.out.println(metaP);
//                            System.exit(0);

                            String[] alleles = variant.getAlleles();
                            String alleleStr = alleles[0] + "/" + alleles[1];
                            tctr++;

                            ////////////////////////////////////////////////////////////////
                            // 2. joint analysis including all datasets, regardless of QC //
                            ////////////////////////////////////////////////////////////////

                            // prune the data
//                            Pair<double[], double[]> prunedData = pruneMissing(Primitives.toPrimitiveArr(genotypeDataForJointAnalysis),
//                                    Primitives.toPrimitiveArr(genotypeDosageDataForJointAnalysis),
//                                    Primitives.toPrimitiveArr(expressionDataForJointAnalysis));
//                            // perform correlation
//                            QCObj combinedQCObj1 = checkVariant(Primitives.toPrimitiveArr(genotypeDataForJointAnalysis));
                            double zJoint1 = Double.NaN;
                            double pJoint1 = Double.NaN;
                            int nJoint1 = 0;
//                            //if (combinedQCObj1.passqc) {
//                            double[] genotypeDataCenterScale = centerScale(prunedData.getLeft());
//                            double rJoint1 = Correlation.correlate(prunedData.getRight(), genotypeDataCenterScale);
//                            nJoint1 = genotypeDataCenterScale.length;
//                            zJoint1 = Correlation.convertCorrelationToZScore(nJoint1, rJoint1);
//                            pJoint1 = ZScores.zToP(zJoint1);
                            //}

                            //////////////////////////////////////////////////////////////
                            // 3. joint analysis including only datasets passing SNP QC //
                            //////////////////////////////////////////////////////////////

                            // prune the data
//                            Pair<double[], double[]> prunedData2 = pruneMissing(Primitives.toPrimitiveArr(genotypeDataForJointAnalysisDatasetsPassingQC),
//                                    Primitives.toPrimitiveArr(genotypeDosageDataForJointAnalysisDatasetsPassingQC),
//                                    Primitives.toPrimitiveArr(expressionDataForJointAnalysisDatasetsPassingQC));
//                            // perform correlation
//                            QCObj combinedQCObj2 = checkVariant(Primitives.toPrimitiveArr(genotypeDataForJointAnalysisDatasetsPassingQC));
                            double zJoint2 = Double.NaN;
                            double pJoint2 = Double.NaN;
                            int nJoint2 = 0;
//
//                            genotypeDataCenterScale = centerScale(prunedData2.getLeft());
//                            double rJoint2 = Correlation.correlate(prunedData2.getRight(), genotypeDataCenterScale);
//                            nJoint2 = genotypeDataCenterScale.length;
//                            zJoint2 = Correlation.convertCorrelationToZScore(nJoint2, rJoint2);
//                            pJoint2 = ZScores.zToP(zJoint2);


                            // "Gene\tSNP\tAlleles\tEffectallele\tMetaN\tMetaZ\tMetaP\tRankedOverAllSamplesZ\tRankedOverAllSamplesP\tRankedPerDatasetZ\tRankedPerDatasetP"
                            output.writeln(gene + "\t" + snpid + "\t" + alleleStr + "\t" + alleles[1] + "\t" +
                                    totalN + "\t" + metaZ + "\t" + metaP + "\t" +
                                    nJoint1 + "\t" + zJoint1 + "\t" + pJoint1 + "\t" +
                                    nJoint2 + "\t" + zJoint2 + "\t" + pJoint2);

                            StringBuilder snplogstr = null;
                            if (!seenSNP.contains(snpid)) {
                                // SNP	DatasetsPassQC	NTotal	NJoint MAFJoint	CallRateJoint	HWEPJoint
//                                snplogstr = new StringBuilder().append(snpid)
//                                        .append("\t").append(datasetsPassQC)
//                                        .append("\t").append(totalN)
//                                        .append("\t").append(nJoint1)
//                                        .append("\t").append(combinedQCObj1.maf)
//                                        .append("\t").append(combinedQCObj1.cr)
//                                        .append("\t").append(combinedQCObj1.hwep)
//                                        .append("\t").append(nJoint2)
//                                        .append("\t").append(combinedQCObj2.maf)
//                                        .append("\t").append(combinedQCObj2.cr)
//                                        .append("\t").append(combinedQCObj2.hwep);
                                snplogstr = new StringBuilder().append(snpid)
                                        .append("\t").append(datasetsPassQC)
                                        .append("\t").append(totalN)
                                        .append("\t").append(nJoint1)
                                        .append("\t").append(0)
                                        .append("\t").append(0)
                                        .append("\t").append(1)
                                        .append("\t").append(nJoint2)
                                        .append("\t").append(0)
                                        .append("\t").append(0)
                                        .append("\t").append(1);
                                for (int d = 0; d < datasetsQCObjs.length; d++) {
                                    snplogstr.append("\t").append(datasetsQCObjs[d].maf)
                                            .append("\t").append(datasetsQCObjs[d].cr)
                                            .append("\t").append(datasetsQCObjs[d].hwep);
                                }
                                outputsnplog.writeln(snplogstr.toString());
                                seenSNP.add(snpid);
                            }
                        }


                    }
                    vctr++;
                    if (vctr % 1000 == 0) {
                        System.out.print("Gene:" + g + "/" + expressionData.genes.length + "; " + vctr + " variants loaded, " + tctr + " variants tested. T: " + timer.getTimeDesc() + "\r");
                    }
                }
                System.out.print("Gene:" + g + "/" + expressionData.genes.length + "; " + vctr + " variants loaded, " + tctr + " variants tested. T: " + timer.getTimeDesc() + "\n");
            }
        }
        outputsnplog.close();
        output.close();
        tabix.close();

    }

    private double[] centerScale(double[] data) {
        double mean = JSci.maths.ArrayMath.mean(data);
        double sd = JSci.maths.ArrayMath.standardDeviation(data);
        double[] output = new double[data.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = (data[i] - mean) / sd;
        }
        return output;
    }


    private Pair<double[], double[]> pruneMissing(double[] datasetGenotypeData, double[] datasetGenotypeDosages, double[] datasetExpressionData) {
        int nrmissing = 0;
        for (int i = 0; i < datasetGenotypeData.length; i++) {
            if (datasetGenotypeData[i] == -1) {
                nrmissing++;
            }
        }
        double[] dosages = new double[datasetGenotypeData.length - nrmissing];
        double[] expression = new double[datasetGenotypeData.length - nrmissing];
        int ctr = 0;
        for (int i = 0; i < datasetGenotypeData.length; i++) {
            if (datasetGenotypeData[i] != -1) {
                expression[ctr] = datasetExpressionData[i];
                dosages[ctr] = datasetGenotypeDosages[i];
//                dosages[ctr] = datasetGenotypeData[i];
                ctr++;
            }
        }
        return new Pair<>(dosages, expression);
    }

    public class QCObj {
        double maf;
        double cr;
        double hwep;
        boolean passqc;

        @Override
        public String toString() {
            return "QCObj{" +
                    "maf=" + maf +
                    ", cr=" + cr +
                    ", hwep=" + hwep +
                    ", passqc=" + passqc +
                    '}';
        }
    }

    private QCObj checkVariant(double[] gt) {
        int obsAA = 0;
        int obsAB = 0;
        int obsBB = 0;

        double freqA = 0;
        double called = 0;
        for (int i = 0; i < gt.length; i++) {
            double gti = gt[i];
            if (gti != -1) {
                called++;
                if (gt[i] == 0) {
                    freqA += 2;
                    obsAA++;
                } else if (gt[i] == 1) {
                    freqA += 1;
                    obsAB++;
                } else {
                    obsBB++;
                }
            }
        }
        double maf = 0;

        double hwep = HWE.calculateExactHWEPValue(obsAB, obsAA, obsBB);

        if (called == 0) {
            maf = 0;
        } else {
            maf = freqA / (called * 2);
            called /= gt.length;
            if (maf > 0.5) {
                maf = 1 - maf;
            }
        }

        QCObj obj = new QCObj();
        obj.maf = maf;
        obj.cr = called;
        obj.hwep = hwep;
        obj.passqc = ((obsAA > 0 && obsAB > 0) || (obsAB > 0 && obsBB > 0) || (obsAA > 0 && obsBB > 0))
                && (maf >= mafthreshold && called >= callratethreshold && hwep >= hwepthreshold);
        return obj;
    }

    private Dataset[] createDatasets(HashMap<String, Integer> DNASampleMap, HashMap<String, Integer> RNASampleMap) {
        System.out.println("Creating datasets.");
        HashMap<String, Dataset> datasetMap = new HashMap<>();
        ArrayList<Dataset> datasets = new ArrayList<>();
        for (String RNA : RNASampleMap.keySet()) {
            String DNA = RNAToDNA.get(RNA);

            if (DNA != null) {
                String datasetName = RNAToDataset.get(RNA);
                Dataset dataset = datasetMap.get(datasetName);
                if (dataset == null) {
                    dataset = new Dataset();
                    dataset.name = datasetName;
                    datasets.add(dataset);
                    datasetMap.put(datasetName, dataset);
                }
                Integer RNAId = RNASampleMap.get(RNA);
                Integer DNAId = DNASampleMap.get(DNA);
                if (RNAId != null && DNAId != null) {
                    dataset.append(DNAId, RNAId);
                }
            }
        }
        Dataset[] output = datasets.toArray(new Dataset[0]);
        Arrays.sort(output);
        System.out.println(datasets.size() + " datasets defined.");
        for (Dataset d : output) {
            System.out.println(d.name + "\t" + d.RNAIds.size() + " samples.");
            d.toArr();
        }
        return output;
    }

    private Dataset createCombinedDataset(HashMap<String, Integer> DNASampleMap, HashMap<String, Integer> RNASampleMap) {
        Dataset dataset = new Dataset();
        dataset.name = "Combined";
        for (String RNA : RNASampleMap.keySet()) {
            String DNA = RNAToDNA.get(RNA);

            if (DNA != null) {
                Integer RNAId = RNASampleMap.get(RNA);
                Integer DNAId = DNASampleMap.get(DNA);
                if (RNAId != null && DNAId != null) {
                    dataset.append(DNAId, RNAId);
                }
            }
        }
        dataset.toArr();
        return dataset;
    }

    private void loadSampleLinks(String linkfile, Set<String> includedDNAs) throws IOException {
        System.out.println("Reading: " + linkfile);
        TextFile tf = new TextFile(linkfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        RNAsPerDataset = new HashMap<>();
        while (elems != null) {
            String gt = elems[0];
            if (includedDNAs.contains(gt)) {
                String rna = elems[1];
                String dataset = elems[2];
                RNAToDNA.put(rna, gt);
                DNAtoRNA.put(gt, rna);
                RNAToDataset.put(rna, dataset);
            }
            elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();
        System.out.println(RNAToDataset.size() + " links between RNA and datasets.");
        System.out.println(RNAToDNA.size() + " links between RNA and DNA");
        System.out.println(DNAtoRNA.size() + " links between DNA and RNA");
    }

    private ArrayList<String> getGenotypeSamples(String vcfFile) throws IOException {
        TextFile tf = new TextFile(vcfFile, TextFile.R);
        String ln = tf.readLine();
        ArrayList<String> samples = new ArrayList<>();
        while (ln != null) {
            if (ln.startsWith("#")) {
                if (ln.startsWith("#CHROM")) {
                    String[] elems = ln.split("\t");
                    for (int i = 9; i < elems.length; i++) {
                        samples.add(elems[i]);
                    }
                }
            } else {
                break;
            }
            ln = tf.readLine();
        }
        tf.close();
        System.out.println(samples.size() + " samples in " + vcfFile);
        return samples;
    }

}
