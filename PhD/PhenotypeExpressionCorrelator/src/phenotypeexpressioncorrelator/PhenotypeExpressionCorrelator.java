/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypeexpressioncorrelator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.RankDoubleArray;

/**
 *
 * @author harmjan
 */
public class PhenotypeExpressionCorrelator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            PhenotypeExpressionCorrelator p = new PhenotypeExpressionCorrelator(args);
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("");
        System.exit(0);
    }
    private String fte;
    private String inexp;
    private String inexpplatform;
    private String phenotypes;
    private String transprobes;
    private String out;
    private String annotation;

    public PhenotypeExpressionCorrelator(String[] args) throws IOException {
        // TODO code application logic here
        System.out.println("");
        System.out.println("------------------------------------");
        System.out.println("Phenotype - Expression correlator v1");
        System.out.println("------------------------------------");
        fte = null;
        inexp = null;
        inexpplatform = null;
        phenotypes = null;
        transprobes = null;
        out = null;

        if (args.length < 2) {
            printUsage();
            System.out.println("FATAL ERROR: some required settings were not specified in your command line");
        } else {
            for (int i = 0; i < args.length; i++) {
                String val = null;
                if (i + 1 < args.length) {
                    val = args[i + 1];
                }
                if (args[i].equals("--inexp")) {
                    inexp = val;
                } else if (args[i].equals("--inexpplatform")) {
                    inexpplatform = val;
                } else if (args[i].equals("--transprobes")) {
                    transprobes = val;
                } else if (args[i].equals("--phenotypes")) {
                    phenotypes = val;
                } else if (args[i].equals("--fte")) {
                    fte = val;
                } else if (args[i].equals("--out")) {
                    out = val;
                } else if (args[i].equals("--annotation")) {
                    annotation = val;
                }
            }
            if (inexp == null || inexpplatform == null || phenotypes == null || transprobes == null || annotation == null) {
                printUsage();
                printSettings();
                System.out.println("FATAL ERROR: some required settings were not specified in your command line");
            } else {
                printSettings();
                if (checkFiles()) {
                    System.out.println("All files are present.");
                    run();
                } else {
                    System.out.println("");
                    System.out.println("FATAL ERROR: There were some problems with the files specified. Please review the messages above");
                }
            }
        }
    }

    private void printUsage() {
        System.out.println("");
        System.out.println("Usage: ");
        System.out.println("------");
        System.out.println("--inexp\t\t\tFull path to expression file");
        System.out.println("--inexpplatform\t\tExpression platform (see file specified at --annotation for your platform identifier)");
        System.out.println("--transprobes\t\tList of trans probes to test");
        System.out.println("--phenotypes\t\tFull path to phenotype file");
        System.out.println("--fte\t\t\tOPTIONAL: path to file that links phenotypes to expression samples");
        System.out.println("--out\t\t\tDirectory to write results to");
        System.out.println("--annotation\t\tPath to probe annotation file");
        System.out.println("");

    }

    private void printSettings() {
        System.out.println("");
        System.out.println("Currently used settings: ");
        System.out.println("-------------------------");
        System.out.println("--inexp:\t\t" + inexp);
        System.out.println("--inexpplatform:\t" + inexpplatform);
        System.out.println("--transprobes:\t\t" + transprobes);
        System.out.println("--phenotypes:\t\t" + phenotypes);
        System.out.println("--fte:\t\t\t" + fte);
        System.out.println("--out:\t\t\t" + out);
        System.out.println("--annotation:\t\t" + annotation);
        System.out.println("");
    }

    private boolean checkFiles() throws IOException {
        System.out.println("");
        System.out.println("Checking whether all required files are there...");
        System.out.println("------------------------------------------------");
        boolean filesOk = true;

        if (fte != null && !Gpio.exists(fte)) {
            System.err.println("fte file specified, but not found: " + fte);
            filesOk = false;
        }

        if (!Gpio.exists(inexp)) {
            System.err.println("inexp file specified, but not found: " + inexp);
            filesOk = false;
        }



        if (!Gpio.exists(transprobes)) {
            System.err.println("transprobes file specified, but not found: " + transprobes);
            filesOk = false;
        }

        if (!Gpio.exists(phenotypes)) {
            System.err.println("phenotypes file specified, but not found: " + phenotypes);
            filesOk = false;
        }

        if (!Gpio.exists(annotation)) {
            System.err.println("annotation file specified, but not found: " + annotation);
            filesOk = false;


        } else {
            // check whether inexpplatform is available in the annottion file.
            System.out.println("Checking whether the " + inexpplatform + " platform is present in the probe annotation file.");
            ProbeTranslation pb = new ProbeTranslation();
            if (pb.getProbeTranslation(annotation, "Probe", inexpplatform) == null) {
                filesOk = false;
            } else {
                System.out.println("Your platform seems to be present..");
            }
        }


        out = Gpio.formatAsDirectory(out);
        if (!Gpio.exists(out)) {
            Gpio.createDir(out);
        }

        if (!Gpio.exists(out)) {
            System.err.println("Tried to create output directory, but somehow failed! " + out);
            filesOk = false;
        }

        return filesOk;
    }

    private void run() throws IOException {

        System.out.println("");
        System.out.println("Now running main program...");
        System.out.println("---------------------------");

        System.out.println("Loading your expression data:");
        DoubleMatrixDataset<String, String> exp = new DoubleMatrixDataset<String, String>(inexp);
        System.out.println("Loading your phenotype data:");
        DoubleMatrixDataset<String, String> phe = new DoubleMatrixDataset<String, String>(phenotypes);

        HashMap<String, String> fteHash = null;
        HashMap<Integer, Integer> fteIndexHash = new HashMap<Integer, Integer>();
        boolean createFTEHash = true;
        if (fte != null) {
            System.out.println("Reading phenotype to expression couplings: " + fte);
            TextFile tf = new TextFile(fte, TextFile.R);
            fteHash = (HashMap<String, String>) tf.readAsHashMap(0, 1);
            tf.close();
            System.out.println(fteHash.size() + " couplings loaded from file. We will check whether these couplings are present in the data at a later stage.");
            createFTEHash = false;
        } else {
            System.out.println("Determining phenotype to expression couplings from loaded matrices..");
        }

        if (createFTEHash) {
            fteHash = new HashMap<String, String>();
        }
        List<String> pheSamples = phe.colObjects;
        for (int i = 0; i < pheSamples.size(); i++) {
            String pheSample = pheSamples.get(i);
            String expSample = pheSample;
            if (!createFTEHash) {
                expSample = fteHash.get(pheSample);
            }
            if (pheSample != null) {
                Integer expSampleIndex = exp.hashCols.get(expSample);
                if (expSampleIndex != null) {
                    if (createFTEHash) {
                        fteHash.put(pheSample, expSample);
                    }
                    fteIndexHash.put(i, expSampleIndex);
                }

            }
        }

        System.out.println("Final number of coupled samples: " + fteIndexHash.size());

        // check whether there are any probes to test
        ProbeTranslation pb = new ProbeTranslation();


        System.out.println("Loading query probes: " + transprobes);
        TextFile tf = new TextFile(transprobes, TextFile.R);
        HashSet<String> queryProbes = new HashSet<String>();
        queryProbes.addAll(tf.readAsArrayList());
        tf.close();
        System.out.println("Maximum number of probes to test: " + queryProbes.size());



        HashMap<String, String> metaToPlatform = pb.getProbeTranslation(annotation, "Probe", inexpplatform);
        HashMap<String, String> platformToMeta = pb.getProbeTranslation(annotation, inexpplatform, "Probe");

        HashSet<Integer> finalProbeList = new HashSet<Integer>();
        for (String probe : queryProbes) {
            String platformProbe = metaToPlatform.get(probe);
            Integer platformId = exp.hashRows.get(platformProbe);
            if (platformId != null) {
                finalProbeList.add(platformId);

            }
        }

        System.out.println("Final number of coupled samples is: " + fteIndexHash.size());
        if (fteIndexHash.isEmpty()) {
            System.err.println("Error: no coupled samples found. Is there something wrong with your sample name? Possibly you should provide --fte ?");
        } else {
            // we can start correlating.

            System.out.println("Everything loaded correctly. Will now start correlating...");
            TextFile outFile = new TextFile(out + "Correlations.txt.gz", TextFile.W);

            outFile.writeln("Platform: " + inexpplatform);
            outFile.writeln("Phenotype\tProbe\tPlatformProbe\tNumSamples\tSpearmanR\tPearsonR");

            // check which samples should be included
            boolean[] pheSampleIsToBeIncluded = new boolean[phe.nrCols];
            for (int i = 0; i < phe.nrCols; i++) {
                if (fteIndexHash.containsKey(i)) {
                    pheSampleIsToBeIncluded[i] = true;
                }
            }

            for (int phenoId = 0; phenoId < phe.nrRows; phenoId++) {
                ArrayList<Double> pheValues = new ArrayList<Double>();

                ArrayList<Integer> queryPheIndexes = new ArrayList<Integer>();
                for (int j = 0; j < phe.nrCols; j++) {
                    if (pheSampleIsToBeIncluded[j] && !Double.isNaN(phe.rawData[phenoId][j])) {
                        pheValues.add(phe.rawData[phenoId][j]);
                        queryPheIndexes.add(j);
                    }
                }

                Double[] x = pheValues.toArray(new Double[0]);
                double[] xpri = new double[x.length];
                for (int q = 0; q < xpri.length; q++) {
                    xpri[q] = x[q];
                }

                RankDoubleArray rda = new RankDoubleArray();
                double[] xpriranked = rda.rank(xpri);

                for (Integer platformProbeId : finalProbeList) {
                    ArrayList<Double> expValues = new ArrayList<Double>();
                    for (int q = 0; q < queryPheIndexes.size(); q++) {
                        Integer expIdIndex = fteIndexHash.get(queryPheIndexes.get(q));
                        expValues.add(exp.rawData[platformProbeId][expIdIndex]);
                    }
                    Double[] y = expValues.toArray(new Double[0]);

                    double[] ypri = new double[y.length];
                    for (int q = 0; q < ypri.length; q++) {
                        ypri[q] = y[q];
                    }


                    double[] ypriranked = rda.rank(ypri);

                    double rpearson = JSci.maths.ArrayMath.correlation(xpri, ypri);
                    double rspearman = JSci.maths.ArrayMath.correlation(xpriranked, ypriranked);

                    // pheno\tprobe\tnrSamples\tr
                    String platformProbe = platformToMeta.get(exp.rowObjects.get(platformProbeId));
                    outFile.writeln(phe.rowObjects.get(phenoId) + "\t" + platformProbe + "\t" + exp.rowObjects.get(platformProbeId) + "\t" + xpri.length + "\t" + rspearman + "\t" + rpearson);


                }
            }

            outFile.close();
            System.out.println("Done! Have a nice day!");
        }
    }
}
