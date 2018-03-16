/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package regulomedb;

import JSci.maths.statistics.ChiSqrDistribution;
import cern.jet.stat.Probability;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import jsc.distributions.ChiSquared;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author harmjan
 */
public class RegulomeDB {

    /**
     * @param args the command line arguments //
     */
//    public static void main(String[] args) {
//        // TODO code application logic here
//        // parse args
//
//        String snps = null;
//        String permsnps = null;
//        String regdb = null;
//        String snpmaffile = null;
//        int nrRandomperms = 10;
//
//        System.out.println("RegulomeDB Parser v1");
//        System.out.println("Written by Harm-Jan Westra");
//        System.out.println("--------------------------");
//
//
//        if (args.length < 2) {
//            System.out.println("Usage: RegulomeDB.jar --snps snplist --regdb regulomedb [--perm permsnplist] [--permrandom snpmaffile [nrperms]]");
//            System.out.println("--perm and --permrandom are optional");
//            System.out.println("--permrandom takes preference above --perm");
//            System.out.println("Nrperms is optional for --permrandom. Default permutations is 10");
//        } else {
//            for (int i = 0; i < args.length - 1; i++) {
//                String arg = args[i];
//                String val = args[i + 1];
//
//                if (arg.equals("--snps")) {
//                    snps = val;
//                } else if (arg.equals("--regdb")) {
//                    regdb = val;
//                } else if (arg.equals("--perm")) {
//                    permsnps = val;
//                } else if (arg.equals("--permrandom")) {
//                    snpmaffile = val;
//                    nrRandomperms = Integer.parseInt(args[i+2]);
//                    
//                }
//            }
//
//            if (permsnps != null && snpmaffile != null) {
//                System.out.println("Warning: --permrandom overrides --perm permsnplist");
//            }
//
//            try {
//                RegulomeDB rdb = new RegulomeDB();
//                rdb.run(snps, regdb, permsnps, snpmaffile, nrRandomperms);
//            } catch (IOException e) {
//                System.out.println("ERROR: There was an error in your files.");
//                e.printStackTrace();
//            }
//        }
//
//
//
//    }
    public static void main(String[] args) {



        RegulomeDB r = new RegulomeDB();

        try {

            r.run("/Volumes/iSnackHD/SkyDrive/SNPFunctionalAnnotation/Trans-SNPs/Real/AllSNPs.txt-WithProxies.txt",
                    "/Volumes/iSnackHD/RegulomeDB/",
                    "/Volumes/iSnackHD/SkyDrive/SNPFunctionalAnnotation/Cis-SNPs/RealData/AllSNPs.txt-WithProxies.txt",
                    "", 1,
                    "/Volumes/iSnackHD/SkyDrive/SNPFunctionalAnnotation/Trans-SNPs/RegulomeDB.txt");
//            r.run("/Volumes/iSnackHD/SkyDrive/Cis-SNPs/RealData/AllSNPs.txt-WithProxies.txt", 
//                    "/Volumes/iSnackHD/RegulomeDB/", 
//                    "/Volumes/iSnackHD/SkyDrive/Cis-SNPs/AllSNPsPerm.txt-WithProxies.txt", 
//                    "", 1, 
//                    "/Volumes/iSnackHD/SkyDrive/Cis-SNPs/RegulomeDB.txt");
        } catch (IOException ex) {
            Logger.getLogger(RegulomeDB.class.getName()).log(Level.SEVERE, null, ex);
        }


    }

    private RegulomeDB() {
    }

    private void run(String snpfilename, String regdbdir, String permsnpfilename, String takerandomsnps, int nrperm, String out) throws IOException {

        HashSet<String> querySNPs = new HashSet<String>();
        TextFile snpfile = new TextFile(snpfilename, TextFile.R);
        String[] snpelems = snpfile.readLineElems(TextFile.tab);
        while (snpelems != null) {
            querySNPs.add(snpelems[0]);
            querySNPs.add(snpelems[1]);
            snpelems = snpfile.readLineElems(TextFile.tab);
        }
        snpfile.close();

        System.out.println("Real SNPs read: " + querySNPs.size());

        HashSet<String> permSNPs = new HashSet<String>();
        snpfile = new TextFile(permsnpfilename, TextFile.R);
        snpelems = snpfile.readLineElems(TextFile.tab);
        while (snpelems != null) {
            permSNPs.add(snpelems[0]);
            permSNPs.add(snpelems[1]);
            snpelems = snpfile.readLineElems(TextFile.tab);
        }
        snpfile.close();

        System.out.println("Perm SNPs read: " + permSNPs.size());

        HashSet<String> uniqueClasses = new HashSet<String>();
        String[] regdbfiles = Gpio.getListOfFiles(regdbdir, "gz");
        HashSet<Pair<String, String>> eqtls = new HashSet<Pair<String, String>>();

        HashSet<String> uniqueSNPs = new HashSet<String>();
        for (String f : regdbfiles) {
            String filename = f;
            System.out.println("Parsing file: " + filename);
            TextFile tf = new TextFile(filename, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);

            while (elems != null) {

                String chr = elems[0];
                String snp = elems[2];
                String pos = elems[1];
                String annot = elems[3];

                if (querySNPs.contains(snp)) {
                    
                    String[] annotelems = annot.split(", ");
                    for (String s : annotelems) {
//                    System.out.println(s);
                        if (s.contains("eQTL")) {
                            uniqueClasses.add("eQTL");
                            // add eQTL
//                        System.out.println(s);
                            String gene = s.split("|")[2];
                            Pair<String, String> eqtl = new Pair<String, String>(snp, gene);
                            eqtls.add(eqtl);

                        } else {
                            uniqueClasses.add(s);
                        }

                    }

                    if (uniqueSNPs.contains(snp)) {
                        System.out.println("Double SNP: " + snp);
                    }
                    uniqueSNPs.add(snp);
                }

                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();
        }

        int ctr = 0;
        HashMap<String, Integer> classIndex = new HashMap<String, Integer>();
        for (String func : uniqueClasses) {
            classIndex.put(func, ctr);
            ctr++;
        }

        int[] realDataCtr = new int[uniqueClasses.size()];
        int[] permDataCtr = new int[uniqueClasses.size()];

        int ctdRealSNPs = 0;
        int ctdPermSNPs = 0;

        System.out.println("Now counting data..");
        for (int perm = 0; perm < 2; perm++) {
            System.out.println("Perm: " + perm);
            HashMap<String, byte[]> snpToAnnot = new HashMap<String, byte[]>();
            HashSet<String> q = querySNPs;
            if (perm == 1) {
                q = permSNPs;
            }

            for (String f : regdbfiles) {
                String filename = f;
                System.out.println("Parsing file: " + filename);
                TextFile tf = new TextFile(filename, TextFile.R);
                String[] elems = tf.readLineElems(TextFile.tab);

                while (elems != null) {

                    String chr = elems[0];
                    String snp = elems[2];
                    String pos = elems[1];
                    String annot = elems[3];

                    if (q.contains(snp)) {
                        byte[] byteannot = snpToAnnot.get(snp);

                        if (byteannot == null) {
                            byteannot = new byte[uniqueClasses.size()];
                        }

                        String[] annotelems = annot.split(", ");
                        for (String s : annotelems) {
                            Integer index = classIndex.get(s);
                            if (index != null) {
                                byteannot[index] = 1;
                            }
                        }
                        snpToAnnot.put(snp, byteannot);
                    }
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
            }



            int[] countarr = realDataCtr;

            if (perm == 1) {
                countarr = permDataCtr;
                ctdPermSNPs = snpToAnnot.size();
            } else {
                ctdRealSNPs = snpToAnnot.size();
            }

            for (String s : q) {
                byte[] annotation = snpToAnnot.get(s);
                if (annotation != null) {
                    for (int i = 0; i < annotation.length; i++) {
                        if (annotation[i] > 0) {
                            countarr[i]++;
                        }
                    }
                }
            }





        }

        System.out.println("Real SNPs: " + ctdRealSNPs);
        System.out.println("Perm SNPs: " + ctdPermSNPs);
        TextFile outfile = new TextFile(out, TextFile.W);
        ChiSqrDistribution d = new ChiSqrDistribution(1);
        for (String s : uniqueClasses) {
            Integer index = classIndex.get(s);

            int real = realDataCtr[index];
            int perm = permDataCtr[index];
            int remaining = ctdRealSNPs - real;
            int remainingperm = ctdPermSNPs - perm;

            double xsq = 0;
            double pval = 1;
            double fet = 1;
            double pval2 = 1;
            if (real > 1 && remaining > 1 && perm > 1 && remainingperm > 1) {
                xsq = ChiSquare.getX(real, remaining, perm, remainingperm);
                pval2 = Probability.chiSquareComplemented(1, xsq);
                fet = new FisherExactTest().getFisherPValue(real, remaining, perm, remainingperm);
            }
            outfile.writeln(s + "\t" + real + "\t" + remaining + "\t" + perm + "\t" + remainingperm + "\t" + xsq + "\t" + pval2 + "\t" + fet);

        }


        outfile.close();
        jsc.distributions.ChiSquared c = new jsc.distributions.ChiSquared(1);


        System.out.println("Chi-test: JSCi: " + d.probability(3.180818876));
        System.out.println("Chi-test: Colt: " + (1 - Probability.chiSquare(1, 3.180818876)));
        System.out.println("Chi-test: jsc: " + (1 - c.cdf(3.180818876)));

    }
}
