/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package csvtott;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CSVToTT {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
//        if (args.length == 1 && args[0].equals("1240")) {
//            try {
//                CSVToTT t = new CSVToTT();
//                for (int chr = 1; chr < 23; chr++) {
//                    String infile = "/Volumes/iSnackHD/ImpData/set_1240_chr" + chr + "_imputed.csv";
//                    String outdir = "/Volumes/iSnackSSD/Data/ImpData/TriTyper/Set1240-Chr" + chr + "/";
//                    String annot = "/Volumes/iSnackHD/ImpData/set_1240_chr" + chr + "_imputed.support.txt";
//                    t.run(infile, outdir, annot);
//                }
//
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//
//        } else {
//
//            try {
//                CSVToTT t = new CSVToTT();
//                for (int chr = 1; chr < 23; chr++) {
//                    String infile = "/Volumes/iSnackHD/ImpData/set_229_chr" + chr + "_imputed.csv";
//                    String outdir = "/Volumes/iSnackSSD/Data/ImpData/TriTyper/Set229-Chr" + chr + "/";
//                    String annot = "/Volumes/iSnackHD/ImpData/set_229_chr" + chr + "_imputed.support.txt";
//                    t.run(infile, outdir, annot);
//                }
//
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }

//        try {
//            if (true) {
//                CSVToTT t = new CSVToTT();
//                int chr = 1;
//                t.verify("/Volumes/iSnackHD/ImpData/set_1240_chr" + chr + "_imputed.support.txt", "/Volumes/iSnackSSD/Data/ImpData/TriTyper/Set1240-Chr" + chr + "/stats.txt");
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
        try {
            if (true) {
                CSVToTT t = new CSVToTT();
                int chr = 1;
                String infile = "/Volumes/iSnackHD/ImpData/set_229_chr" + chr + "_imputed.csv";
                String ttdir = "/Volumes/iSnackSSD/Data/ImpData/TriTyper/Set229-Chr" + chr + "/";
                t.verifyTTFile(infile, ttdir);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void check(String inDir, String pattern) throws IOException {
        String[] lsof = Gpio.getListOfFiles(inDir);
        System.out.println(lsof.length + " files match pattern.");
        HashSet<String> indsInFirstFile = null;
        HashSet<String> snps = new HashSet<String>();
        boolean firstfile = true;
        Arrays.sort(lsof);
        for (String file : lsof) {

            if (file.startsWith(pattern) && file.endsWith("csv")) {
                System.out.println(file);
                file = inDir + file;
                if (indsInFirstFile == null) {
                    indsInFirstFile = new HashSet<String>();
                } else {
                    firstfile = false;
                }
                TextFile tf = new TextFile(file, TextFile.R);
                String[] headerLineElems = tf.readLineElems(Strings.comma);
                for (int i = 1; i < headerLineElems.length; i++) {
                    String snp = headerLineElems[i];
                    if (snps.contains(snp)) {
                        System.out.println("ERROR: snp " + snp + " is present multiple times");
                    } else {
                        snps.add(snp);
                    }
                }


                String[] lineElems = tf.readLineElems(Strings.comma);
                int nrcols = 0;
                double max = Double.MIN_VALUE;
                double min = Double.MAX_VALUE;
                while (lineElems != null) {
                    String ind = lineElems[0];

                    if (nrcols == 0) {
                        nrcols = lineElems.length;
                    } else {
                        if (nrcols != lineElems.length) {
                            System.out.println("Different column length. Was: " + nrcols + ". Now is: " + lineElems.length);
                        }
                    }
                    if (firstfile) {
                        if (indsInFirstFile.contains(ind)) {
                            System.out.println("Duplicate ind in first file: " + ind);
                        } else {
                            indsInFirstFile.add(ind);
                        }
                    } else {
                        if (indsInFirstFile.contains(ind)) {
                        } else {
                            System.out.println("Could not find ind in first file: " + ind);
                        }
                    }
                    for (int i = 1; i < lineElems.length; i++) {
                        double d = Double.parseDouble(lineElems[i]);
                        if (d > max) {
                            max = d;
                        }
                        if (d < min) {
                            min = d;
                        }
                    }

                    lineElems = tf.readLineElems(Strings.comma);
                }
                tf.close();
                System.out.println(file + "\t" + snps.size() + "\tsnps\t x " + indsInFirstFile.size() + " inds\t nrcols:\t" + nrcols + " max: " + max + " min " + min);
            }
        }

        System.out.println("Final size: " + snps.size() + " x " + indsInFirstFile.size());
    }

    public void run(String file, String outdir, String annot) throws IOException {


        System.out.println("Running: " + file);

        Gpio.createDir(outdir);


        ArrayList<String> inds = new ArrayList<String>();

        ArrayList<String> snps = new ArrayList<String>();

        TextFile tf = new TextFile(file, TextFile.R);
        String[] headerLineElems = tf.readLineElems(Strings.comma);
        TextFile snpFile = new TextFile(outdir + "SNPs.txt", TextFile.W);
        TextFile snpMapFile = new TextFile(outdir + "SNPMappings.txt.gz", TextFile.W);
        for (int i = 1; i < headerLineElems.length; i++) {
            String snp = headerLineElems[i];
            snp = snp.replaceAll("\"", "");
            snps.add(snp);
            snpMapFile.writeln("1\t1\t" + snp);
            snpFile.writeln(snp);
        }
        snpMapFile.close();
        snpFile.close();
        System.out.println("Written info for: " + snps.size() + " SNPs");

        String[] lineElems = tf.readLineElems(Strings.comma);
        int nrcols = 0;
        double max = Double.MIN_VALUE;
        double min = Double.MAX_VALUE;

        TextFile indFile = new TextFile(outdir + "Individuals.txt", TextFile.W);
        TextFile phenoFile = new TextFile(outdir + "PhenotypeInformation.txt", TextFile.W);
        while (lineElems != null) {
            String ind = lineElems[0];
            ind = ind.replaceAll("\"", "");
            indFile.writeln(ind);
            phenoFile.writeln(ind + "\tcontrol\tinclude\tmale");
            inds.add(ind);
            lineElems = tf.readLineElems(Strings.comma);
        }
        phenoFile.close();
        indFile.close();
        tf.close();
        System.out.println("Written info for: " + inds.size() + " Individuals");

        // now parse the actual file
        tf.open();
        tf.readLine(); // header
        int indId = 0;



        // load allele annotation for snps
        // temporary lookup table:
        HashSet<String> snpsInFile = new HashSet<String>();
        snpsInFile.addAll(snps);

        System.out.println("Loading: " + annot);
        HashMap<String, byte[]> snpAlleles = new HashMap<String, byte[]>();
        System.out.println("SNPs in file: " + snpsInFile.size());
        TextFile alleleAnnotFile = new TextFile(annot, TextFile.R);
        String[] elems = alleleAnnotFile.readLineElems(TextFile.tab);
        while (elems != null) {

            if (snpsInFile.contains(elems[1])) {
                byte[] alleles = new byte[2];
                alleles[0] = BaseAnnot.toByte(elems[2]);
                alleles[1] = BaseAnnot.toByte(elems[3]);
                snpAlleles.put(elems[1], alleles);
            }
            elems = alleleAnnotFile.readLineElems(TextFile.tab);
        }
        alleleAnnotFile.close();

        System.out.println("Loaded annotation for: " + snpAlleles.size() + " SNPs.");
        if (snpAlleles.size() == snps.size()) {
            System.out.println("This is in accordance with the SNPs in the genotype data");
        } else {
            System.out.println("WARNING: " + (snpAlleles.size() - snps.size()) + " difference in nr of SNPs with annotation!!!!@#!@#!");
        }

        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(snps.size(), inds.size(), new File(outdir + "GenotypeMatrix.dat"), false);
        WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(snps.size(), inds.size(), new File(outdir + "ImputedDosageMatrix.dat"), false);

        ProgressBar pb = new ProgressBar(inds.size(), "Converting genotypes....");
        int lnctr = 0;

        int nrprocs = Runtime.getRuntime().availableProcessors();
        System.out.println("Using " + nrprocs + " threads for parsing the file");
        ExecutorService threadPool = Executors.newFixedThreadPool(nrprocs);
        CompletionService<Triple<Integer, String, byte[]>> pool = new ExecutorCompletionService<Triple<Integer, String, byte[]>>(threadPool);



        int bufferSize = 50;

        int nrSubmitted = 0;
        int bufferNr = 0;
        boolean finished = false;
        while (!finished) {

            String line = tf.readLine();
            if (line == null) {
                // write buffer (if any)
                if (nrSubmitted > 0) {
                    processResults(pool, fileMatrixGenotype, matrixImputedDosage, nrSubmitted, snpAlleles, snps, bufferNr, bufferSize, inds.size());
                }
                finished = true;
            } else {
                ByteParseTask task = new ByteParseTask(line, 1, indId, null, null, Strings.comma);
                pool.submit(task);
                nrSubmitted++;
            }

            if (nrSubmitted == bufferSize) {
                processResults(pool, fileMatrixGenotype, matrixImputedDosage, nrSubmitted, snpAlleles, snps, bufferNr, bufferSize, inds.size());
                bufferNr++;
                nrSubmitted = 0;
            }

            indId++;
            //lineElems = tf.readLineElems(Strings.comma);

            pb.set(indId);
        }
        pb.close();

        fileMatrixGenotype.close();
        matrixImputedDosage.close();

        tf.close();

        System.out.println("------------");
        System.out.println("");
        threadPool.shutdown();
    }

    private void processResults(CompletionService<Triple<Integer, String, byte[]>> pool,
            WGAFileMatrixGenotype fileMatrixGenotype, WGAFileMatrixImputedDosage matrixImputedDosage,
            int nrResultsToProcess, HashMap<String, byte[]> snpAlleles, ArrayList<String> snps, int bufferNr, int bufferSize, int nrIndsTotal) throws IOException {
        int returnedResults = 0;

        byte[][] dosageBuffer = new byte[nrResultsToProcess][0];

        while (returnedResults < nrResultsToProcess) {
            try {
                Triple<Integer, String, byte[]> result = pool.take().get();
                if (result != null) {

                    int actualIndId = result.getLeft();
                    int bufferIndId = actualIndId - (bufferSize * bufferNr);
                    dosageBuffer[bufferIndId] = result.getRight();
                    returnedResults++;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }




        // transpose
        int bufferStart = (bufferSize * bufferNr);
//        System.out.println(bufferSize+"\t"+bufferNr+"\t"+bufferStart);
        for (int snp = 0; snp < snps.size(); snp++) {
            byte[] alleles1 = new byte[nrResultsToProcess];
            byte[] alleles2 = new byte[nrResultsToProcess];
            byte[] dosages = new byte[nrResultsToProcess];
            byte[] alleles = snpAlleles.get(snps.get(snp));
            for (int ind = 0; ind < nrResultsToProcess; ind++) {

//            int actualInd = ind + (bufferSize * bufferNr);
                byte d = dosageBuffer[ind][snp];
                double dosageValue = ((double) (-Byte.MIN_VALUE + d)) / 100;
                byte allele1 = 0;
                byte allele2 = 0;
                if (dosageValue < 0.5) {
                    allele1 = alleles[0];
                    allele2 = alleles[0];
                } else {
                    if (dosageValue > 1.5) {
                        allele1 = alleles[1];
                        allele2 = alleles[1];
                    } else {
                        allele1 = alleles[0];
                        allele2 = alleles[1];
                    }
                }
                alleles1[ind] = allele1;
                alleles2[ind] = allele2;
                dosages[ind] = d;
            }

            // write the results
            fileMatrixGenotype.setAllele1(snp, bufferStart, alleles1);
            fileMatrixGenotype.setAllele1(snp, bufferStart + nrIndsTotal, alleles2);
            matrixImputedDosage.setDosage(snp, bufferStart, dosages);
        }



        /*
         * for (int snpId = 0; snpId < doubles.length; snpId++) {

         byte[] alleles = snpAlleles.get(snps.get(snpId));

         double d = doubles[snpId];


         
         int dosageInt = (int) Math.round(dosageValue * 100d);
         byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
         byte allele1 = 0;
         byte allele2 = 0;
         if (dosageValue < 0.5) {
         allele1 = alleles[0];
         allele2 = alleles[0];
         } else {
         if (dosageValue > 1.5) {
         allele1 = alleles[1];
         allele2 = alleles[1];
         } else {
         allele1 = alleles[0];
         allele2 = alleles[1];
         }
         }

         //                        System.out.println(indId + "\t" + snpId + "t" + snps.get(snpId) + "\t" + d + "\t" + dosageInt + "\t" + BaseAnnot.toString(alleles[0]) + "/" + BaseAnnot.toString(alleles[1]));
         //
         //                        if (snpId > 3) {
         //                            System.exit(0);
         //                        }

         
         }
         */
    }

    private void verify(String supportfile, String statfile) throws IOException {
        HashMap<String, String> minorAllele = new HashMap<String, String>();
        HashMap<String, String> otherAllele = new HashMap<String, String>();
        HashMap<String, Double> minorAlleleFreq = new HashMap<String, Double>();

        TextFile tf1 = new TextFile(supportfile, TextFile.R);
        tf1.readLine();
        String[] elems = tf1.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String allele1 = elems[2];
            String allele2 = elems[3];
            double freq1 = Double.parseDouble(elems[4]);
            double maf = Double.parseDouble(elems[5]);
            minorAlleleFreq.put(snp, maf);
            if (freq1 == maf) {
                minorAllele.put(snp, allele1);
                otherAllele.put(snp, allele2);
            } else {
                minorAllele.put(snp, allele2);
                otherAllele.put(snp, allele1);
            }
            elems = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();

        TextFile tf2 = new TextFile(statfile, TextFile.R);
        tf2.readLine();
        elems = tf2.readLineElems(TextFile.tab);
        int concordant = 0;
        int discordant = 0;
        double maxDiff = Double.MIN_VALUE;
        double minDiff = Double.MAX_VALUE;
        while (elems != null) {
            if (elems.length > 5) {
                String snp = elems[0];
                String[] alleles = elems[1].split("/");
                String allele1 = alleles[0];
                String allele2 = alleles[1];
                String minorallele = elems[2];
                String otherallele = "";
                if (allele1.equals(minorallele)) {
                    otherallele = allele2;
                } else {
                    otherallele = allele1;
                }

                Double freq = Double.parseDouble(elems[4]);

                String origMinorAllele = minorAllele.get(snp);
                String origOtherAllele = otherAllele.get(snp);
                Double origAlleleFreq = minorAlleleFreq.get(snp);
                if (otherallele.equals(origOtherAllele) && minorallele.equals(origMinorAllele)) {
                    // nothing aan het handje
                    double f = Math.abs(origAlleleFreq - freq);
                    if (f > maxDiff) {
                        maxDiff = f;
                    }
                    if (f < minDiff) {
                        minDiff = f;
                    }

                    if (f > 0.005) {
                        System.out.println("Discrepancy with converted data: " + snp + "\t" + minorallele + "/" + otherallele + "\t" + origMinorAllele + "/" + origOtherAllele + "\t" + freq + "\t" + origAlleleFreq + "\t" + f);
                    }

                    concordant++;
                } else {
//                    System.out.println("Discrepancy with converted data: "+snp+"\t"+minorallele+"/"+otherallele+"\t"+origMinorAllele+"/"+origOtherAllele+"\t"+freq+"\t"+origAlleleFreq);
                    discordant++;
                }

            }
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();

        System.out.println(concordant + "\t" + discordant);
        System.out.println(maxDiff + "\t" + minDiff);

    }

    private void verifyTTFile(String infile, String ttdir) throws IOException {
        System.out.println("Comparing: " + infile + "\t" + ttdir);
        TriTyperGenotypeData ds = new TriTyperGenotypeData();
        ds.load(ttdir);
        SNPLoader loader = ds.createSNPLoader();

        TextFile cs = new TextFile(infile, TextFile.R);
        String ln = cs.readLine();
        String[] snp = Strings.comma.split(ln);
        cs.close();

        int nrTested = 0;
        int nrFullyConcordant = 0;
        for (int i = 1; i < snp.length; i++) {


//            snp[i] = snp[i].replaceAll("\"", "");
            Integer snpId = ds.getSnpToSNPId().get(snp[i].replaceAll("\"", ""));
            if (snpId == null) {
                System.err.println("ERROR: " + snp[i] + " not converted????? ?");
            } else {
                SNP snpObj = ds.getSNPObject(snpId);

                loader.loadGenotypes(snpObj);
                loader.loadDosage(snpObj);

                byte[] alleles = snpObj.getAlleles();
                if (alleles[0] == 0 || alleles[1] == 0) {
                    // don't compare
                } else {
                    cs.open();
                    cs.readLine();
                    String[] elems = cs.readLineElems(Strings.comma);
                    double[] dosages = new double[ds.getIndividuals().length];
                    while (elems != null) {
                        String value = elems[i].replaceAll("\"", "");
                        String sample = elems[0].replaceAll("\"", "");
                        Integer sampleId = ds.getIndividualId(sample);
                        if (sampleId == null) {
                            System.err.println("ERROR: " + sample + " not converted???? ");
                        }
                        dosages[sampleId] = Double.parseDouble(value);
                        elems = cs.readLineElems(Strings.comma);
                    }
                    cs.close();

                    double[] dosagesInTT = snpObj.getDosageValues();
                    int nrEqual = 0;
                    for (int d = 0; d < dosagesInTT.length; d++) {
                        if (dosagesInTT[d] == dosages[d]) {
                            nrEqual++;
                        } else {
                            System.out.println(ds.getIndividuals()[d] + "\t" + dosagesInTT[d] + "\t" + dosages[d] + "\t" + Math.abs(dosagesInTT[d] - dosages[d]));
                        }
                    }
                    if (nrEqual == dosages.length) {
                        nrFullyConcordant++;
                    }
                    nrTested++;
                }
                snpObj.clearGenotypes();
            }

            System.out.println(i + "\t" + nrTested + "\t" + nrFullyConcordant);
        }
        loader.close();
    }
}
