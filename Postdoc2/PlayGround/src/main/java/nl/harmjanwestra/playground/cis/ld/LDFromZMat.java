package nl.harmjanwestra.playground.cis.ld;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix.SymmetricFloatDistanceMatrix;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class LDFromZMat {


    public static void main(String[] args) {
        LDFromZMat z = new LDFromZMat();
        try {
//            z.createCombos("D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\sortedGeneSNPCombos.txt.gz",
//                    "D:\\snpgenecombos\\");
//
            if (args.length < 3) {

                z.test(args[0], args[1]);

            } else if (args.length >= 7) {
                String genecovariancematrixloc = args[0];
                String cohortselectionfile = args[1];
                String zmatloc = args[2];
                String genedefloc = args[3];
                String genelistfile = args[4];
                String referenceallelefile = args[5];
                String outloc = args[6];
                Integer minnrvals = null;
                boolean providelist = false;
                boolean binaryoutput = true;

                if (args.length >= 8) {
                    minnrvals = Integer.parseInt(args[7]);
                }

                if (args.length > 9) {
                    providelist = Boolean.parseBoolean(args[8]);
                }
                if (args.length > 10) {
                    binaryoutput = Boolean.parseBoolean(args[9]);
                }
                z.calculateLD(genecovariancematrixloc, cohortselectionfile, zmatloc, genedefloc, genelistfile, referenceallelefile, outloc, minnrvals, providelist, binaryoutput);
            } else {
                System.out.println("Arguments:\ngenecovariancematrix cohortselectionfile zmatloc genedefloc genelistfile referenceallelefile outploc [minnrvals] [providelist] [binaryoutput]");

            }


        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.exit(0);
    }

    private void test(String input, String output) throws IOException {

//        outputmatrix.save(new File(outloc + querygene + "-ldmatrix.dat.gz"));
//        TextFile tfsnp = new TextFile(outloc + querygene + "-ldmatrix.snp.gz", TextFile.W);

        TextFile tf = new TextFile(input + "-ldmatrix.snp.gz", TextFile.R);
        ArrayList<String> snps = tf.readAsArrayList();
        tf.close();

        SymmetricFloatDistanceMatrix m = new SymmetricFloatDistanceMatrix(snps.size());

        m.load(new File(input + "-ldmatrix.dat.gz"));

        TextFile tfo = new TextFile(output, TextFile.W);
        String header = "-\t" + Strings.concat(snps, Strings.tab);
        tfo.writeln(header);
        for (int i = 0; i < snps.size(); i++) {
            String ln = snps.get(i);
            for (int j = 0; j < snps.size(); j++) {
                ln += "\t" + m.get(i, j);
            }
            tfo.writeln(ln);
        }
        tfo.close();

    }

    public void createCombos(String genesnpcombos, String snpgeneout) throws IOException {


        TextFile tf = new TextFile(genesnpcombos, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);


        HashMap<String, Integer> geneMap = new HashMap<String, Integer>();
        ArrayList<String> geneList = new ArrayList<>();

        HashMap<String, ArrayList<Integer>> snpgenecombos = new HashMap<>();


        int gctr = 0;
        int lctr = 0;
        while (elems != null) {

            String gene = elems[0];
            String snp = elems[1];

            if (!geneMap.containsKey(gene)) {
                geneMap.put(gene, gctr);
                geneList.add(gene);
                gctr++;
                if (gctr < 0) {
                    System.exit(-1);
                }
            }

            Integer geneId = geneMap.get(gene);


            ArrayList<Integer> genes = snpgenecombos.get(snp);
            if (genes == null) {
                genes = new ArrayList<>();
            }
            genes.add(geneId);
            snpgenecombos.put(snp, genes);

            lctr++;
            if (lctr % 100000 == 0) {
                System.out.print("\r" + lctr + " lines parsed ");
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        System.out.println(snpgenecombos.size() + " snps loaded. ");
        tf.open();

        elems = tf.readLineElems(TextFile.tab);
        String prevgene = null;
        TextFile out = null;
        lctr = 0;
        while (elems != null) {
            String gene = elems[0];
            String snp = elems[1];

            if (prevgene == null || !gene.equals(prevgene)) {
                if (out != null) {
                    out.close();
                }
                out = new TextFile(snpgeneout + gene + ".txt.gz", TextFile.W);
            }

            prevgene = gene;

            ArrayList<Integer> combos = snpgenecombos.get(snp);
            String[] genes = new String[combos.size()];
            for (int i = 0; i < genes.length; i++) {
                genes[i] = geneList.get(combos.get(i));
            }
            out.writeln(snp + "\t" + Strings.concat(genes, Strings.semicolon));

            lctr++;
            if (lctr % 100000 == 0) {
                System.out.println("Lines parsed: " + lctr + " current gene: " + prevgene);
            }
            elems = tf.readLineElems(TextFile.tab);
        }

        if (out != null) {
            out.close();
        }

    }


    public void calculateLD(String genecovariancematrixloc,
                            String cohortSelectionFile,
                            String zmatloc,
                            String geneDefLoc,
                            String geneListFile,
                            String referenceAlleleFile,
                            String outloc,
                            Integer minNrVals,
                            boolean writelist,
                            boolean writebinary) throws Exception {

        boolean run = true;

        if (!Gpio.exists(cohortSelectionFile)) {
            System.out.println("Could not find cohort selection file: " + cohortSelectionFile);
            run = false;
        }
        if (!Gpio.exists(zmatloc)) {
            System.out.println("Could not find z-score matrix location: " + zmatloc);
            run = false;
        }
        if (!Gpio.exists(geneDefLoc)) {
            System.out.println("Could not find gene definition location: " + geneDefLoc);
            run = false;
        }

        if (!Gpio.exists(referenceAlleleFile)) {
            System.out.println("Could not find reference allele file: " + referenceAlleleFile);
            run = false;
        }

        if (!Gpio.exists(outloc)) {
            Gpio.createDir(outloc);
        }

        //		String[] cohortsEuropean = {"inCHIANTI", "HVH_HT12v3", "EGCUT_RNAseq", "LIFE_Adult_plus", "EGCUT_HT12v4", "Fehrmann_HT12v3", "NTR_NESDA", "LIFE_Heart_plus", "BSGS", "CODAM", "LLS_660Q", "DILGOM", "Fehrmann_H8v2", "GoNL_WGS", "NTR_GoNL", "SHIP_TREND", "CARTaGENE_freeze2", "Rotterdam_HT12v4", "KORA_F4", "Cardiology", "NTR_AFFY", "GTEx", "EGCUT_HT12v3", "DGN", "CHDWB", "CARTaGENE_freeze1", "PAN", "LLS_OmniExpr", "Sorbs", "FHS", "HVH_HT12v4", "Rotterdam_RNASeq", "LL"};
//		String[] genes = {"ENSG00000002726", "ENSG00000002933", "ENSG00000055118", "ENSG00000106565", "ENSG00000164867", "ENSG00000177590", "ENSG00000241134"};
//		String genecovariancematrixloc = "/Users/lude/Documents/Genetica/eQTLGen/ZScoreMarices-20180125//GeneCorrelation/ProbeCovariance-Perms1To10.binary";
//
//		String cohortSelectionFile = null;
//		String zmatloc = null;
//		String geneDefLoc = null;
//		String geneListFile = null;
//		String outloc;

        if (minNrVals == null) {
            minNrVals = 0;
        }

        ArrayList<String> genesToRun = new ArrayList<>();

        if (run && geneListFile.startsWith("ENSG")) {
            genesToRun.add(geneListFile);
        } else {
            if (!Gpio.exists(geneListFile)) {
                System.out.println("Could not find gene list file: " + geneListFile);
                run = false;
            }
            if (!run) {
                System.out.println("Some gtfToProbeAnnotationFile conditions not met.");
                System.exit(-1);
            }
            System.out.println("Loading genes from: " + geneListFile);
            TextFile tfc = new TextFile(geneListFile, TextFile.R);
            genesToRun = tfc.readAsArrayList();
            tfc.close();
        }


        System.out.println("Selected " + genesToRun.size() + " genes to calculate LD for.");
        if (genesToRun.isEmpty()) {
            System.out.println("No genes found to gtfToProbeAnnotationFile");
            System.exit(-1);
        }

        HashSet<String> hashCohortsEuropean = new HashSet<String>();
        {
            System.out.println("Loading cohorts from file: " + cohortSelectionFile);
            TextFile tfcs = new TextFile(cohortSelectionFile, TextFile.R);
            ArrayList<String> cohortsToInclude = tfcs.readAsArrayList();
            tfcs.close();
            System.out.println("Cohorts to include: " + cohortsToInclude.size());

            for (int c = 0; c < cohortsToInclude.size(); c++) {
                for (int perm = 1; perm <= 10; perm++) {
                    hashCohortsEuropean.add(cohortsToInclude.get(c) + "-perm-" + perm);
                }
            }
            System.out.println(hashCohortsEuropean.size() + " cohorts/permutation names selected");
        }


        if (!referenceAlleleFile.endsWith(".bin") && !Gpio.exists(referenceAlleleFile + ".bin")) {
            System.out.println("Converting reference allele map: " + referenceAlleleFile);
            {
                TextFile tf = new TextFile(referenceAlleleFile, TextFile.R);
                BinaryFile bf = new BinaryFile(referenceAlleleFile + ".bin", BinaryFile.W);
                String[] elems = tf.readLineElems(TextFile.tab);
                int ctr = 0;
                while (elems != null) {
                    String snp = elems[0];
                    String alleles = elems[1] + ";" + elems[2];
//                    referenceAlleleMap.put(snp, Strings.cache(alleles));
                    bf.writeString(snp);
                    bf.writeString(alleles);

                    ctr++;
                    if (ctr % 10000 == 0) {
                        System.out.print(ctr + " lines parsed sofar\r");
                    }
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
                bf.close();
                System.out.println();
            }
        }

        System.out.println("");


        if (hashCohortsEuropean.isEmpty()) {
            System.out.println("No cohorts specified");
            System.exit(-1);
        }

        // load reference alleles
        for (String querygene : genesToRun) {
            // determine which snps and genes to read..

            String geneDef = geneDefLoc + querygene + ".txt.gz";
            if (!Gpio.exists(geneDef)) {
                System.out.println(querygene + "\tCould not find gene definition: " + geneDef);
                break;
            } else {

                HashSet<String> snpSet = new HashSet<String>();
                ArrayList<String> snpList = new ArrayList<>();
                LinkedHashSet<String> geneSet = new LinkedHashSet<String>();


                System.out.println(querygene + "\t" + "Loading gene definition: " + geneDef);

                TextFile tf = new TextFile(geneDef, TextFile.R);
                String[] elems = tf.readLineElems(TextFile.tab);
                int sctr = 0;
                while (elems != null) {
                    String snp = elems[0];
                    if (sctr < 10) {
                        System.out.println(snp);
                    }
                    sctr++;

                    String[] genes = elems[1].split(";");
                    snpSet.add(snp);

                    for (String gene : genes) {
                        geneSet.add(new String(gene.getBytes("UTF-8")));
                    }
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();

                System.out.println(querygene + "\t" + geneSet.size() + " genes loaded for gene");
                System.out.println(querygene + "\t" + snpSet.size() + " snps for gene");

                snpList.addAll(snpSet);

                double[] weights = null;
                String[] genesWithWeights = null;
                {
                    DoubleMatrixDataset<String, String> datasetCorr = DoubleMatrixDataset.loadSubsetOfBinaryDoubleData(
                            genecovariancematrixloc, geneSet, geneSet);

                    System.out.println(querygene + "\t" + datasetCorr.rows() + " genes present in gene covariance matrix");

                    //Calculate the factorloadings:
                    Jama.EigenvalueDecomposition eig = eigenValueDecomposition(datasetCorr.getMatrix().toArray());
                    double[] eigenValues = eig.getRealEigenvalues();

                    DoubleMatrix2D factorLoadings = new DenseDoubleMatrix2D(datasetCorr.rows(), datasetCorr.rows());
                    for (int comp = 0; comp < datasetCorr.rows(); comp++) {
                        double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
                        if (eigenvalue < 0) {
                            eigenvalue = 0;
                        }
                        double sqrtEigenvalue = Math.sqrt(eigenvalue);
                        double[] eigenvector = getEigenVector(eig, comp);

                        for (int a = 0; a < datasetCorr.rows(); a++) {
                            double v = sqrtEigenvalue * eigenvector[a];
                            factorLoadings.setQuick(comp, a, v);
                        }
                    }

                    //Calculate the weights of the individual genes, to be used for the weighed correlation, to account for co-expression between genes:
                    weights = new double[datasetCorr.rows()];
                    genesWithWeights = new String[datasetCorr.rows()];
                    for (int p = 0; p < datasetCorr.rows(); p++) {
                        double weight = 0;
                        genesWithWeights[p] = datasetCorr.getRowObjects().get(p);
                        for (int comp = 0; comp < datasetCorr.rows(); comp++) {
                            double eigenvalue = eigenValues[eigenValues.length - 1 - comp];
                            if (eigenvalue < 1) {
                                eigenvalue = 1;
                            }
                            double factorLoadingCompP = factorLoadings.getQuick(comp, p);
                            weight += factorLoadingCompP * factorLoadingCompP / eigenvalue;
                        }

                        weights[p] = weight;
                        System.out.println(querygene + "\t" + p + "\t" + datasetCorr.getRowObjects().get(p) + "\t" + weights[p]);
                    }
                }

                HashMap<String, String> referenceAlleleMap = new HashMap<String, String>();

                System.out.println("Reading alleles from: " + referenceAlleleFile + ".bin");
                BinaryFile bf = new BinaryFile(referenceAlleleFile + ".bin", BinaryFile.R);

//                        String snp = bf.readString();
//                        String alleles = bf.readString();

                int actr = 0;
                while (bf.available() > 0) {
                    String nextSNP = bf.readString();
                    String alleles = bf.readString();
                    if (snpSet.contains(nextSNP)) {
                        referenceAlleleMap.put(nextSNP, alleles);
                    }

                    actr++;
                    if (actr % 100000 == 0) {
                        System.out.print("\r" + actr + " snps loaded.");
                    }

                }
                bf.close();
                System.out.println();


                ArrayList<DoubleMatrixDataset<String, String>> zmats = new ArrayList<>();
                // load relevant z-score matrices..
                for (int g = 0; g < genesWithWeights.length; g++) {
                    String geneToLoad = zmatloc + genesWithWeights[g] + "-zmat.txt.gz";
                    String allelefile = zmatloc + genesWithWeights[g] + "-referenceAlleles.txt.gz";

                    if (Gpio.exists(geneToLoad)) {
                        System.out.println(querygene + "\tLoading gene z-score matrix: " + g + "/" + genesWithWeights.length + "\t" + geneToLoad);
                        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(geneToLoad, '\t', hashCohortsEuropean, snpSet);
                        zmats.add(ds);

                        System.out.println("Reference alleles loaded for: " + referenceAlleleMap.size() + " snps.");

                        // flip alleles
                        TextFile tfa = new TextFile(allelefile, TextFile.R);
                        elems = tfa.readLineElems(TextFile.tab);
                        while (elems != null) {

                            String snp = elems[0];
                            Integer snpid = ds.getHashCols().get(snp);
                            if (snpid != null) {
                                String alleles = elems[1];
                                String assessed = elems[2];
                                String referenceAlleles = referenceAlleleMap.get(snp);
                                if (referenceAlleles != null) {
                                    String[] refAlleleElems = referenceAlleles.split(";");
                                    Boolean flip = BaseAnnot.flipalleles(refAlleleElems[0], refAlleleElems[1], alleles, assessed);
                                    if (flip == null) {
                                        // kill snp
                                        for (int r = 0; r < ds.rows(); r++) {
                                            ds.getMatrix().setQuick(r, snpid, Double.NaN);
                                        }
                                    } else if (flip) {
                                        // kill snp
                                        for (int r = 0; r < ds.rows(); r++) {
                                            ds.getMatrix().setQuick(r, snpid, ds.getMatrix().getQuick(r, snpid) * -1);
                                        }
                                    }
                                } else {
                                    // kill snp
                                    for (int r = 0; r < ds.rows(); r++) {
                                        ds.getMatrix().setQuick(r, snpid, Double.NaN);
                                    }
                                }
                            }

                            elems = tfa.readLineElems(TextFile.tab);
                        }


                        tfa.close();
                    } else {
                        zmats.add(null);
                    }
                }


                int nrSNPs = snpSet.size();


//                DoubleMatrixDataset<String, String> ldmatrix = new DoubleMatrixDataset<String, String>(nrSNPs, nrSNPs);
//                ldmatrix.getMatrix().assign(Double.NaN);
//                ldmatrix.setRowObjects(snpList);
//                ldmatrix.setColObjects(snpList);
                SymmetricFloatDistanceMatrix outputmatrix = new SymmetricFloatDistanceMatrix(nrSNPs);
//                SymmetricShortDistanceMatrix outputnmatrix = new SymmetricShortDistanceMatrix(nrSNPs);

                System.out.println(querygene + "\tSaving list file here: " + outloc + querygene + "-list.txt.gz");

                TextFile tfout = null;
                if (writelist) {
                    tfout = new TextFile(outloc + querygene + "-list.txt.gz", TextFile.W);
                    String header = "SNP1\tSNP2\tNrValues\tNrSamples\tR\tRSq";
                    tfout.writeln(header);
                }
//SNP1    SNP2    NrValues        NrSamples       R(PEARSON)      Rsq


//                DoubleMatrix2D matrixObj = ldmatrix.getMatrix();
                ProgressBar pb = new ProgressBar(nrSNPs, querygene + " - Calculating correlations:");

                double[] finalWeights = weights;
                Integer finalMinNrVals = minNrVals;
                AtomicInteger ctr = new AtomicInteger();
                TextFile finalTfout = tfout;
                IntStream.range(0, nrSNPs).parallel().forEach(v -> {
                    try {
                        run(snpList, v, nrSNPs, zmats, finalWeights, finalMinNrVals, outputmatrix, ctr, pb, finalTfout);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                });


                pb.close();
                if (writelist) {
                    tfout.close();
                }


                System.out.println(querygene + "\tSaving LD matrix here: " + outloc + querygene + "-ldmatrix.dat.gz");

                // save upper triangle in a binary format
                outputmatrix.save(new File(outloc + querygene + "-ldmatrix.dat.gz"));
                TextFile tfsnp = new TextFile(outloc + querygene + "-ldmatrix.snp.gz", TextFile.W);
                tfsnp.writeList(snpList);
                tfsnp.close();



//				if (writebinary) {
//                BinaryFile bf = new BinaryFile(outloc + querygene + "-ldmatrix.bin.gz", BinaryFile.W);
//                bf.writeInt(snpList.size());
//                for (String s : snpList) {
//                    bf.writeString(s);
//                }
//
//                for (int i = 0; i < snpList.size(); i++) {
//                    for (int j = i + 1; j < snpList.size(); j++) {
//                        float f = (float) ldmatrix.getMatrix().getQuick(i, j);
//                        bf.writeFloat(f);
//                    }
//                }
//                bf.close();
//                ldmatrix.save(outloc + querygene + "-ldmatrix.txt.gz");
//				} else {
//
//				}
            }
        }


    }


    public void run(ArrayList<String> snpList, int snp1, int nrSNPs,
                    ArrayList<DoubleMatrixDataset<String, String>> zmats, double[] weights,
                    int minNrVals, SymmetricFloatDistanceMatrix matrixObj, AtomicInteger ctr, ProgressBar pb, TextFile tfout) throws IOException {
        String snp1Name = snpList.get(snp1);

        ArrayList<double[]> snp1dd = new ArrayList<>();
        for (int g = 0; g < zmats.size(); g++) {
            DoubleMatrixDataset<String, String> geneZmat = zmats.get(g);
            Integer snp1id = geneZmat.getHashCols().get(snp1Name);
            if (snp1id == null) {
                snp1dd.add(null);
            } else {
                double[] col1 = geneZmat.getMatrix().viewColumn(snp1id).toArray();
                snp1dd.add(col1);
            }
        }

        for (int snp2 = snp1 + 1; snp2 < nrSNPs; snp2++) {
            ArrayList<Double> snp1d = new ArrayList<>();
            ArrayList<Double> snp2d = new ArrayList<>();
            ArrayList<Double> weightsd = new ArrayList<>();
            String snp2Name = snpList.get(snp2);

            // accumulate data over genes
            for (int g = 0; g < zmats.size(); g++) {
                DoubleMatrixDataset<String, String> geneZmat = zmats.get(g);
                if (geneZmat != null) {
                    double[] col1 = snp1dd.get(g);
                    if (col1 != null) {

                        LinkedHashMap<String, Integer> zmatcolhash = geneZmat.getHashCols();
                        Integer snp2id = zmatcolhash.get(snp2Name);
                        if (snp2id != null) {
                            double[] col2 = geneZmat.getMatrix().viewColumn(snp2id).toArray();

                            for (int row = 0; row < geneZmat.rows(); row++) {
                                double v1 = col1[row]; //col1.get(row); //geneZmat.getElementQuick(row, snp1id);
                                double v2 = col2[row]; //col2.get(row); // geneZmat.getElementQuick(row, snp2id);

                                if (!Double.isNaN(v1) && !Double.isNaN(v2)) {
                                    snp1d.add(v1);
                                    snp2d.add(v2);
                                    weightsd.add(weights[g]);
                                }
                            }
                        }
                    }
                }
            }

            // calculate weighted correlation
            if (snp1d.size() >= minNrVals) {
                double[] vals1 = Primitives.toPrimitiveArr(snp1d);
                double[] vals2 = Primitives.toPrimitiveArr(snp2d);
                double[] valWeights = Primitives.toPrimitiveArr(weightsd);

//				double corr = JSci.maths.ArrayMath.correlation(vals1, vals2);
                double weightedCorr = weightedCorrelation(vals1, vals2, valWeights);
//                            System.out.println(snp1 + "\t" + snp2 + "\t" + snp1Name + "\t" + snp2Name + "\t" + snp1d.size() + "\t" + corr + "\t" + weightedCorr);
                if (tfout != null) {
                    tfout.writelnsynced(snp1Name + "\t" + snp2Name + "\t" + snp1d.size() + "\t" + snp1d.size() + "\t" + weightedCorr + "\t" + (weightedCorr * weightedCorr));
                }
                matrixObj.set(snp1, snp2, (float) weightedCorr);

//				matrixObj.setQuick(snp2, snp1, weightedCorr);
            }
        }

        matrixObj.set(snp1, snp1, 1f);

        ctr.getAndIncrement();
        pb.set(ctr.get());
    }

    public double weightedCorrelation(double[] x, double[] y, double[] weights) {
        double wmX = weightedMean(x, weights);
        double wmY = weightedMean(y, weights);
        double sumWeights = JSci.maths.ArrayMath.mass(weights);
        double covXX = 0;
        double covXY = 0;
        double covYY = 0;
        for (int s = 0; s < x.length; s++) {
            covXX += weights[s] * (x[s] - wmX) * (x[s] - wmX);
            covXY += weights[s] * (x[s] - wmX) * (y[s] - wmY);
            covYY += weights[s] * (y[s] - wmY) * (y[s] - wmY);
        }
        covXX /= sumWeights;
        covXY /= sumWeights;
        covYY /= sumWeights;
        double corr = covXY / (Math.sqrt(covXX * covYY));
        return corr;
    }

    public double weightedMean(double[] x, double[] weights) {
        double m = 0;
        double sumWeights = 0;
        for (int s = 0; s < x.length; s++) {
            m += x[s] * weights[s];
            sumWeights += weights[s];
        }
        return m / sumWeights;
    }


    private Jama.EigenvalueDecomposition eigenValueDecomposition(double[][] data) {
        Jama.Matrix m = new Jama.Matrix(data);
        Jama.EigenvalueDecomposition eig = m.eig();
        return eig;
    }

    private double[] getEigenVector(Jama.EigenvalueDecomposition eig, double[] eigenValues, int pca) {
        Jama.Matrix eigenValueMatrix = eig.getV();
        double[][] eigenValueMat = eigenValueMatrix.getArray();
        double[] eigenVector = new double[eigenValueMat.length];
        for (int i = 0; i < eigenValueMat.length; i++) {
            eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
        }
        return eigenVector;
    }

    private double[] getEigenVector(Jama.EigenvalueDecomposition eig, int pca) {
        Jama.Matrix eigenValueMatrix = eig.getV();
        double[][] eigenValueMat = eigenValueMatrix.getArray();
        double[] eigenVector = new double[eigenValueMat.length];
        for (int i = 0; i < eigenValueMat.length; i++) {
            eigenVector[i] = eigenValueMat[i][eigenValueMat.length - 1 - pca]; // * Math.sqrt(eigenValues[eigenValues.length - 1 - pca]);
        }
        return eigenVector;
    }


    private double getEigenValueVar(double[] eigenValues, int pca) {
        double sumEigenvalues = 0.0;
        for (Double d : eigenValues) {
            sumEigenvalues += Math.abs(d);
        }
        double result = eigenValues[eigenValues.length - 1 - pca] / sumEigenvalues;
        return result;
    }

    private double[] getEigenVectorSVD(Jama.SingularValueDecomposition svd, double[] singularValues, int pca) {
        Jama.Matrix eigenValueMatrix = svd.getV();
        double[][] eigenValueMat = eigenValueMatrix.getArray();
        double[] eigenVector = new double[eigenValueMat.length];
        for (int i = 0; i < eigenValueMat.length; i++) {
            eigenVector[i] = eigenValueMat[i][pca] * Math.sqrt(singularValues[pca]);
        }
        return eigenVector;
    }
}
