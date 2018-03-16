/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.cis.figure1;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author harmjan
 */
public class SNPNexusParser {

    /**
     * @param args the command line arguments
     *
     * perm + "\t" + (iter+1000) + "\t" + cnp + "\t" + cosmic + "\t" + cpg +
     * "\t" + indel + "\t" + inversion + "\t" + mirna + "\t" + tfbs + "\t" +
     * vista + "\t" + a.numSNPs + "\t" + a.coding + "\t" + a.syn + "\t" +
     * a.nonsyn + "\t" + a.frameshift + "\t" + a.noncoding + "\t" +
     * a.noncodingintronic + "\t" + a.downstream3 + "\t" + a.upstream5 + "\t" +
     * a.intronic
     *
     *
     * + "\t" + a.utr3 + "\t" + a.utr5
     */
    public static void main(String[] args) {
        // TODO code application logic here
        SNPNexusParser p = new SNPNexusParser();
//        p.parseCis();
        p.parseTrans();
    }
    HashSet<String> countedSNPs = null;

    public void parseTrans() {
        try {
//            String dir = "";
////            String outfile = "d:\\Skydrive\\Cis-SNPs\\SNPNexusOutput.txt";
//            String outfile = "D:\\SkyDrive\\latesteQTLs\\transAnnot\\SNPSelectionsForFuncAnalysis\\SNPNexusOutput.txt";
//            TextFile out = new TextFile(outfile, TextFile.W);
//            out.writeln("perm\titer\tSnpsInBin\tcnp\tcosmic\tcpg\tindel\tinversion\tmirna\ttfbs\tvista\tensemblnrsnps\tcoding\tsyn\tnonsyn\tframeshift\tnoncoding\tdownstream\tupstream\tintronic\tutr3\tutr5");
//            for (int perm = 0; perm < 11; perm++) {
//                String proxyfile = "";
//                String binfile = "";
//
////                String indir = "d:\\Skydrive\\Cis-SNPs\\";
//                String indir = "D:\\SkyDrive\\latesteQTLs\\transAnnot\\SNPSelectionsForFuncAnalysis\\";
//                String annotDir = "D:\\SkyDrive\\latesteQTLs\\transAnnot\\SNPSelectionsForFuncAnalysis\\";
//                String snpdir =  indir;  //"/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/SNPSelectionsForFuncAnalysis/";
//                if (perm == 0) {
//                    snpdir += "RealData/";
//                    indir += "RealData/";
//                    proxyfile = indir + "AllSNPs.txt-WithProxies.txt";
//                    
//                    annotDir += "RealData/";
//                } else {
//                    
//
//                    snpdir += "PermutationRound" + perm + "/";
//                    indir += "PermutationRound" + perm + "/";
//                    proxyfile = indir + "AllSNPs.txt-WithProxies.txt";
//                    annotDir += "PermutationRound" + perm + "/";
//                }
//
//                if (perm == 0 || perm == 1) {
//                    countedSNPs = new HashSet<String>();
//                }
//                for (int iter = 0; iter < 1000; iter += 1000) {
//                    binfile = snpdir + "Bin" + iter + "-" + (iter + 1000) + ".txt";
//                    String output = run(perm, iter, binfile, proxyfile, annotDir);
//                    out.writeln(output);
//                }
//            }
//            out.close();
//
//            chisquare(outfile);


            String dir = "";
//            String outfile = "d:\\Skydrive\\Cis-SNPs\\SNPNexusOutput.txt";
            String outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/RandomlySNPSelectionsForFuncAnalysis/SNPNexusOutput.txt";
            TextFile out = new TextFile(outfile, TextFile.W);
            out.writeln("perm\titer\tSnpsInBin\tcnp\tcosmic\tcpg\tindel\tinversion\tmirna\ttfbs\tvista\tensemblnrsnps\tcoding\tsyn\tnonsyn\tframeshift\tnoncoding\tdownstream\tupstream\tintronic\tutr3\tutr5");
            for (int perm = 0; perm < 21; perm++) {
                String proxyfile = "";
                String binfile = "";

//                String indir = "d:\\Skydrive\\Cis-SNPs\\";
                String indir = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/RandomlySNPSelectionsForFuncAnalysis/";
                String annotDir = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/RandomlySNPSelectionsForFuncAnalysis/";
                String snpdir = indir;  //"/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/SNPSelectionsForFuncAnalysis/";
                if (perm == 0) {
                    snpdir += "RealData/";
                    indir += "RealData/";
                    proxyfile = indir + "AllSNPs.txt-WithProxies.txt";
                    annotDir += "RealData/";
                } else {
                    snpdir += "Set" + (perm - 1) + "/";
                    indir += "Set" + (perm - 1) + "/";
                    proxyfile = indir + "Set-" + (perm - 1) + ".txt-WithProxies.txt";
                    annotDir += "Set" + (perm - 1) + "/";
                }

                if (perm == 0 || perm == 1) {
                    countedSNPs = new HashSet<String>();
                }
                for (int iter = 0; iter < 1000; iter += 1000) {
                    if (perm == 0) {
                        binfile = snpdir + "Bin" + iter + "-" + (iter + 1000) + ".txt";
                    } else {
                        binfile = snpdir + "Set-" + (perm - 1) + ".txt";
                    }

                    String output = run(perm, iter, binfile, proxyfile, annotDir);
                    out.writeln(output);
                }
            }
            out.close();

            chisquare(outfile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void parseCis() {
        try {
            String dir = "";
//            String outfile = "d:\\Skydrive\\Cis-SNPs\\SNPNexusOutput.txt";
            String outfile = "D:\\SkyDrive\\latesteQTLs\\cisAnnot\\2012-10-23-PostQC-SNPsForFuncAnnot\\SNPNexusOutput.txt";
            TextFile out = new TextFile(outfile, TextFile.W);
            out.writeln("perm\titer\tSnpsInBin\tcnp\tcosmic\tcpg\tindel\tinversion\tmirna\ttfbs\tvista\tensemblnrsnps\tcoding\tsyn\tnonsyn\tframeshift\tnoncoding\tdownstream\tupstream\tintronic\tutr3\tutr5");
            for (int perm = 0; perm < 11; perm++) {
                String proxyfile = "";
                String binfile = "";

//                String indir = "d:\\Skydrive\\Cis-SNPs\\";
                String indir = "D:\\SkyDrive\\latesteQTLs\\cisAnnot\\2012-10-23-PostQC-SNPsForFuncAnnot\\";
                String annotDir = "D:\\SkyDrive\\latesteQTLs\\cisAnnot\\OldAnnot\\Cis-SNPs\\";
                String snpdir = indir;  //"/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/SNPSelectionsForFuncAnalysis/";
                if (perm == 0) {
                    snpdir += "RealData/";
                    indir += "RealData/";
                    proxyfile = indir + "AllSNPs.txt-WithProxies.txt";

                    annotDir += "RealData/";
                } else {
                    proxyfile = indir + "AllSNPsPerm.txt-WithProxies.txt";

                    snpdir += "PermutationRound" + perm + "/";
                    indir += "PermutationRound" + perm + "/";
                    annotDir += "PermutationRound" + perm + "/";
                }

                if (perm == 0 || perm == 1) {
                    countedSNPs = new HashSet<String>();
                }
                for (int iter = 0; iter < 8000; iter += 1000) {
                    binfile = snpdir + "Bin" + iter + "-" + (iter + 1000) + ".txt";
                    String output = run(perm, iter, binfile, proxyfile, annotDir);
                    out.writeln(output);
                }
            }
            out.close();

            chisquare(outfile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    private HashSet<String> snpsInBin = null;

    public String run(int perm, int iter, String binfile, String proxyfile, String inDir) throws IOException {

        // read snps in bin
        TextFile binSNPs = new TextFile(binfile, TextFile.R);
        HashSet<String> querySNPs = new HashSet<String>();

        querySNPs.addAll(binSNPs.readAsArrayList());
        binSNPs.close();

        // now read proxies
        TextFile pf = new TextFile(proxyfile, TextFile.R);
        String[] elems = pf.readLineElems(TextFile.tab);
        ArrayList<String> proxies = new ArrayList<String>();
        while (elems != null) {
            String query = elems[0];
            if (querySNPs.contains(query)) {
                proxies.add(elems[1]);
            }
            elems = pf.readLineElems(TextFile.tab);
        }
        pf.close();

        querySNPs.addAll(proxies);

        snpsInBin = new HashSet<String>();
        for (String snp : querySNPs) {
            if (!countedSNPs.contains(snp)) {
                snpsInBin.add(snp);
                countedSNPs.add(snp);
            }
        }

        String[] allTextFilesInDir = Gpio.getListOfFiles(inDir);
        int cnp = 0;
        int cosmic = 0;
        int cpg = 0;
        int indel = 0;
        int inversion = 0;
        int mirna = 0;
        int tfbs = 0;
        int vista = 0;
        EnsemblAnnotation a = null;

        for (String annotationFile : allTextFilesInDir) {

            if (annotationFile.startsWith("cnp")) {
                cnp = parseFile(inDir + annotationFile);
            }
            if (annotationFile.startsWith("consv")) {
            }
            if (annotationFile.startsWith("cosmic")) {
                cosmic = parseFile(inDir + annotationFile);
            }
            if (annotationFile.startsWith("cpg")) {
                cpg = parseFile(inDir + annotationFile);
            }
            if (annotationFile.startsWith("firstef")) {
            }
            if (annotationFile.startsWith("gad")) {
            }
            if (annotationFile.startsWith("gwas")) {
            }
            if (annotationFile.startsWith("indel")) {
                indel = parseFile(inDir + annotationFile);
            }
            if (annotationFile.startsWith("invb")) {
            }
            if (annotationFile.startsWith("inversion")) {
                inversion = parseFile(inDir + annotationFile);
            }
            if (annotationFile.startsWith("mirbase")) {
            }
            if (annotationFile.startsWith("mirna1")) {
                mirna = parseFile(inDir + annotationFile);
            } else if (annotationFile.startsWith("mirna1")) {
            }
            if (annotationFile.startsWith("tfbs")) {
                tfbs = parseTFBSFile(inDir + annotationFile);
            }
            if (annotationFile.startsWith("vista")) {
                vista = parseFile(inDir + annotationFile);
            }
            if (annotationFile.startsWith("ensembl")) {

                a = parseEnsemblFile(inDir + annotationFile);

            }
        }

        if (a == null) {
            System.out.println("NO ENSEMBL FOR " + perm + "\t" + iter);
        }

        return (perm + "\t" + (iter + 1000)
                + "\t" + snpsInBin.size()
                + "\t" + cnp
                + "\t" + cosmic
                + "\t" + cpg
                + "\t" + indel
                + "\t" + inversion
                + "\t" + mirna
                + "\t" + tfbs
                + "\t" + vista
                + "\t" + a.numSNPs
                + "\t" + a.coding
                + "\t" + a.syn
                + "\t" + a.nonsyn
                + "\t" + a.frameshift
                + "\t" + a.noncoding
                + "\t" + a.downstream3
                + "\t" + a.upstream5
                + "\t" + a.intronic
                + "\t" + a.utr3
                + "\t" + a.utr5);

    }

    /*
     * public int intronic = 0;
     public int utr3 = 0;
     public int utr5 = 0;
     public int coding = 0;
     public int downstream3 = 0;
     public int upstream5 = 0;
     public int noncoding = 0;
     private int noncodingintronic;
     private int frameshift;
     private int numSNPs;
     */
    public int parseTFBSFile(String filename) throws IOException {

        int ctr = 0;
        HashSet<String> visitedSNPs = new HashSet<String>();

        TextFile tf = new TextFile(filename, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab); // skip header
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length > 6) {
                String snp = elems[0];
                String type = elems[6];
                if (type.contains("human") && snpsInBin.contains(snp) && !visitedSNPs.contains(snp)) {
                    ctr++;
                    visitedSNPs.add(snp);
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return ctr;
    }

    public int parseFile(String filename) throws IOException {

        int ctr = 0;
        HashSet<String> visitedSNPs = new HashSet<String>();

        TextFile tf = new TextFile(filename, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab); // skip header
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[0];
            if (snpsInBin.contains(snp) && !visitedSNPs.contains(snp)) {
                ctr++;
                visitedSNPs.add(snp);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return ctr;
    }

    public Triple<Integer, Integer, Integer> parseFirstEfFile(String filename) throws IOException {

        int exonctr = 0;
        int promoterctr = 0;
        int cpgctr = 0;
        HashSet<String> visitedSNPs = new HashSet<String>();

        TextFile tf = new TextFile(filename, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab); // skip header
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[0];
            if (snpsInBin.contains(snp) && !visitedSNPs.contains(snp)) {

                visitedSNPs.add(snp);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return new Triple<Integer, Integer, Integer>(exonctr, promoterctr, cpgctr);
    }

    private EnsemblAnnotation parseEnsemblFile(String annotationFile) throws IOException {
        TextFile in = new TextFile(annotationFile, TextFile.R);
        String[] elems = in.readLineElems(TextFile.tab);
        elems = in.readLineElems(TextFile.tab);

        int intronic = 0;
        int utr3 = 0;
        int utr5 = 0;
        int coding = 0;
        int downstream3 = 0;
        int upstream5 = 0;
        int noncoding = 0;

        int frameshift = 0;
        int noncodingintronic = 0;

        int nonsyn = 0;
        int syn = 0;

        HashSet<String> visitedSNPs = new HashSet<String>();
        HashMap<String, HashSet<String>> functionsPerSNP = new HashMap<String, HashSet<String>>();


        while (elems != null) {
            if (elems.length > 5) {
                String snp = elems[0];

                if (snpsInBin.contains(snp)) {
                    HashSet<String> functionsOfSNP = functionsPerSNP.get(snp);
                    if (functionsOfSNP == null) {
                        functionsOfSNP = new HashSet<String>();
                    }

                    visitedSNPs.add(snp);
                    String func = elems[5];
                    if (func.equals("coding")) {
                        if (!functionsOfSNP.contains("coding")) {
                            coding++;
                            if (elems.length > 10) {
                                String addFunc = elems[10];
                                if (addFunc.equals("nonsyn") || addFunc.equals("*nonsyn") || addFunc.equals("nonsyn:stop-gain") || addFunc.equals("*nonsyn:stop-gain") || addFunc.equals("nonsyn:stop-loss")) {
                                    if (!functionsOfSNP.contains("nonsyn")) {
                                        nonsyn++;
                                        functionsOfSNP.add("nonsyn");
                                    }
                                } else if (addFunc.equals("syn") || addFunc.equals("*syn")) {
                                    if (!functionsOfSNP.contains("syn")) {
                                        syn++;
                                        functionsOfSNP.add("syn");
                                    }
                                } else if (addFunc.equals("frameshift:stop-gain")) {
                                    if (!functionsOfSNP.contains("frameshift:stop-gain")) {
                                        frameshift++;
                                        functionsOfSNP.add("frameshift:stop-gain");
                                    }
                                } else {
                                    System.out.println("AddFunc: " + addFunc);
                                }
                            }
                            functionsOfSNP.add("coding");
                        }

                    } else if (func.equals("3utr")) {
                        if (!functionsOfSNP.contains("utr3")) {
                            utr3++;
                            functionsOfSNP.add("utr3");
                        }

                    } else if (func.equals("5utr")) {
                        if (!functionsOfSNP.contains("utr5")) {
                            utr5++;
                            functionsOfSNP.add("utr5");
                        }

                    } else if (func.equals("intronic")) {
                        if (!functionsOfSNP.contains("intronic")) {
                            intronic++;
                            functionsOfSNP.add("intronic");
                        }

                    } else if (func.equals("non-coding")) {
                        if (!functionsOfSNP.contains("non-coding")) {
                            noncoding++;
                            functionsOfSNP.add("non-coding");
                        }

                    } else if (func.equals("5upstream")) {
                        if (!functionsOfSNP.contains("5upstream")) {
                            upstream5++;
                            functionsOfSNP.add("5upstream");
                        }

                    } else if (func.equals("3downstream")) {
                        if (!functionsOfSNP.contains("3downstream")) {
                            downstream3++;
                            functionsOfSNP.add("3downstream");
                        }

                    } else if (func.equals("non-coding intronic")) {
                        if (!functionsOfSNP.contains("intronic")) {
                            intronic++;
                            functionsOfSNP.add("intronic");
                        }

                    } else {
                        System.out.println(func);
                    }

                    functionsPerSNP.put(snp, functionsOfSNP);
                }
            }
            elems = in.readLineElems(TextFile.tab);
        }


        in.close();


        EnsemblAnnotation a = new EnsemblAnnotation();
        a.numSNPs = functionsPerSNP.size();
        a.coding = coding;
        a.downstream3 = downstream3;
        a.intronic = intronic;
        a.noncoding = noncoding;
        a.upstream5 = upstream5;
        a.utr3 = utr3;
        a.utr5 = utr5;

        a.frameshift = frameshift;

        a.syn = syn;
        a.nonsyn = nonsyn;
        return a;

    }

    private void chisquare(String outfile) throws IOException {
        // read the collected data
        System.out.println("Performing chi-squared");
        TextFile in = new TextFile(outfile, TextFile.R);
        String[] header = in.readLineElems(TextFile.tab);
        int numElements = header.length;

        String[] elems = in.readLineElems(TextFile.tab);
        int[][] elemsPerPerm = new int[2][22];
        while (elems != null) {

            int perm = Integer.parseInt(elems[0]);
            if (perm > 0) {
                perm = 1;
            }
            for (int e = 1; e < elems.length; e++) {
                elemsPerPerm[perm][e] += Integer.parseInt(elems[e]);
            }
            elems = in.readLineElems(TextFile.tab);
        }

        // done loading data, perform chi-squared
        in.close();
        TextFile out = new TextFile(outfile + "-ChiSquared.txt", TextFile.W);
        // 
        out.writeln("Class\tChi-Squared\tP-val\tFisherExact\tRealSNPsInClass\tRemainingRealSNPs\tPermSNPsInClass\tPermSNPsRemaining");
        System.out.println(elemsPerPerm[0].length);
        for (int e = 1; e < elemsPerPerm[0].length; e++) {
//            System.out.println(header[e]);
            int snpsRealTotal = 0;
            int snpsRealInFunc = 0;
            int snpsPermTotal = 0;
            int snpsPermInFunc = 0;
            if (e == 1 || e == 2 || e == 11) {
                // do nothing
            } else {
                if (e >= 11) {
                    // ensembl annotation data should be treated differently
                    snpsRealTotal = elemsPerPerm[0][11];
                    snpsRealInFunc = elemsPerPerm[0][e];
                    snpsPermTotal = elemsPerPerm[1][11];
                    snpsPermInFunc = elemsPerPerm[1][e];
                } else {
                    snpsRealTotal = elemsPerPerm[0][2];
                    snpsRealInFunc = elemsPerPerm[0][e];
                    snpsPermTotal = elemsPerPerm[1][2];
                    snpsPermInFunc = elemsPerPerm[1][e];
                }

                int snpsRealRemaining = snpsRealTotal - snpsRealInFunc;
                int snpsPermRemaining = snpsPermTotal - snpsPermInFunc;

//                if (snpsRealRemaining > 0 && snpsPermRemaining > 0 && snpsRealInFunc > 0 && snpsPermInFunc > 0) {
//                    double xsq = ChiSquare.getX(snpsRealInFunc, snpsRealRemaining, snpsPermInFunc, snpsPermRemaining);
//
//                    double p = 1;
//                    try {
//                        p = ChiSquare.getP(1, xsq);
//                    } catch (JSci.maths.statistics.OutOfRangeException ex) {
//                        System.out.println("Error with xsq: " + xsq);
//                        System.out.println("Problem with datasource: " + header[e] + "\t" + snpsRealInFunc + "\t" + snpsRealRemaining + "\t" + snpsPermInFunc + "\t" + snpsPermRemaining);
////                            ex.printStackTrace();
//
//                    }
                FisherExactTest fet = new FisherExactTest();
                double pfisher = fet.getFisherPValue(snpsRealInFunc, snpsRealRemaining, snpsPermInFunc, snpsPermRemaining);
                System.out.println(e + "\t" + header[e] + "\t" + 0 + "\t" + 1 + "\t" + pfisher + "\t" + snpsRealInFunc + "\t" + snpsRealRemaining + "\t" + snpsPermInFunc + "\t" + snpsPermRemaining);
                out.writeln(header[e] + "\t" + 0 + "\t" + 1 + "\t" + pfisher + "\t" + snpsRealInFunc + "\t" + snpsRealRemaining + "\t" + snpsPermInFunc + "\t" + snpsPermRemaining);
//                } else {
//                    System.out.println("Problem with datasource: " + header[e] + "\t" + snpsRealInFunc + "\t" + snpsRealRemaining + "\t" + snpsPermInFunc + "\t" + snpsPermRemaining);
//                }

            }
        }
        out.close();
    }

    public class EnsemblAnnotation {

        public int intronic = 0;
        public int utr3 = 0;
        public int utr5 = 0;
        public int coding = 0;
        public int downstream3 = 0;
        public int upstream5 = 0;
        public int noncoding = 0;
        private int frameshift;
        private int numSNPs;
        private int syn;
        private int nonsyn;
    }
}
