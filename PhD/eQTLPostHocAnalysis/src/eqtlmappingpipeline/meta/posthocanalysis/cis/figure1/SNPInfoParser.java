/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.cis.figure1;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class SNPInfoParser {

    HashSet<String> snpsAlreadyCounted = new HashSet<String>();

    public static void main(String[] args) {
        SNPInfoParser p = new SNPInfoParser();
        p.run();
    }

    public void run() {
        try {

//            TextFile out = new TextFile("/Volumes/iSnackHD/Skydrive/Cis-SNPs/FuncAnnotFractions.txt", TextFile.W);
//            out.writeln("Perm\tBin\tsnps\ttfbs\tsplice1\tsplice2\tsplice3\tmirna1\tmirna2\tnssnp\tstop\treg\tcons");
//
//
//            SNPInfoParser p = new SNPInfoParser();
//            int iterCtr = 0;
//            for (int iter = 0; iter < 15000; iter += 1000) {
//                String iterStr = "Bin" + iter + "-" + (iter + 1000) + ".txt";
//                String result = p.collect(0, (iter + 1000),
//                        "/Volumes/iSnackHD/Skydrive/Cis-SNPs/RealData/AllSNPs.txt-WithProxies.txt",
//                        "/Volumes/iSnackHD/Skydrive/Cis-SNPs/RealData/" + iterStr,
//                        "/Volumes/iSnackHD/Skydrive/Cis-SNPs/RealData/snpfunc.csv", true);
//                iterCtr++;
//                out.writeln(result);
//            }
//
//            for (int perm = 1; perm < 11; perm++) {
//                iterCtr = 0;
//                for (int iter = 0; iter < 15000; iter += 1000) {
//                    String iterStr = "Bin" + iter + "-" + (iter + 1000) + ".txt";
//                    String result = p.collect(perm, (iter + 1000),
//                            "/Volumes/iSnackHD/Skydrive/Cis-SNPs/AllSNPsPerm.txt-WithProxies.txt",
//                            "/Volumes/iSnackHD/Skydrive/Cis-SNPs/PermutationRound" + perm + "/" + iterStr,
//                            "/Volumes/iSnackHD/Skydrive/Cis-SNPs/PermutationRound" + perm + "/snpfunc.csv", true);
//                    iterCtr++;
//                    out.writeln(result);
//                }
//            }
//
//            out.close();

//            String dir = "d:\\SkyDrive\\Cis-SNPs\\";

//            String dir = "D:\\SkyDrive\\latesteQTLs\\transAnnot\\SNPSelectionsForFuncAnalysis\\";
//            String snpdir = dir;
//
//            TextFile out = new TextFile(snpdir + "/FuncAnnot-SNPInfo.txt", TextFile.W);
//            out.writeln("Perm\tBin\tsnps\ttfbs\tsplice1\tsplice2\tsplice3\tmirna1\tmirna2\tnssnp\tstop\treg\tcons");
//
//
//            int max = 1000;
//
//            int iterCtr = 0;
//            snpsAlreadyCounted = new HashSet<String>();
//
//            for (int iter = 0; iter < max; iter += 1000) {
//                String iterStr = "Bin" + iter + "-" + (iter + 1000) + ".txt";
//                String result = collect(0, (iter + 1000),
//                        dir + "/RealData/AllSNPs.txt-WithProxies.txt",
//                        snpdir + "/RealData/" + iterStr,
//                        dir + "/RealData/snpfunc.csv", false);
//                iterCtr++;
//
//                out.writeln(result);
//            }
//
//            System.out.println("Actually counted: " + ctr);
//            System.out.println("In HashMap: " + snpsAlreadyCounted.size());
//
//            snpsAlreadyCounted = new HashSet<String>();
//            for (int perm = 1; perm < 11; perm++) {
//                iterCtr = 0;
//                for (int iter = 0; iter < max; iter += 1000) {
//                    String iterStr = "Bin" + iter + "-" + (iter + 1000) + ".txt";
//                    String result = collect(perm, (iter + 1000),
//                            snpdir + "/PermutationRound" + perm + "/AllSNPs.txt-WithProxies.txt",
//                            snpdir + "/PermutationRound" + perm + "/" + iterStr,
//                            dir + "/PermutationRound" + perm + "/snpfunc.csv", false);
//                    iterCtr++;
//                    out.writeln(result);
//                }
//            }
//
//            out.close();


//            String dir = "D:\\SkyDrive\\latesteQTLs\\cisAnnot\\2012-10-23-PostQC-SNPsForFuncAnnot\\";
//            String annotDir = "D:\\SkyDrive\\latesteQTLs\\cisAnnot\\OldAnnot\\Cis-SNPs\\";
//            
//            String snpdir = dir;
//
//            TextFile out = new TextFile(snpdir + "/FuncAnnot-SNPInfo.txt", TextFile.W);
//            out.writeln("Perm\tBin\tsnps\ttfbs\tsplice1\tsplice2\tsplice3\tmirna1\tmirna2\tnssnp\tstop\treg\tcons");
//
//
//            int max = 1000;
//
//            int iterCtr = 0;
//            snpsAlreadyCounted = new HashSet<String>();
//
//            for (int iter = 0; iter < max; iter += 1000) {
//                String iterStr = "Bin" + iter + "-" + (iter + 1000) + ".txt";
//                String proxyStr = "AllSNPs.txt-WithProxies.txt";
//                String result = collect(0, (iter + 1000),
//                        annotDir + "/RealData/"+proxyStr,
//                        snpdir + "/RealData/" + iterStr,
//                        annotDir + "/RealData/snpfunc.csv", false);
//                iterCtr++;
//
//                out.writeln(result);
//            }
//
//            System.out.println("Actually counted: " + ctr);
//            System.out.println("In HashMap: " + snpsAlreadyCounted.size());
//
//            snpsAlreadyCounted = new HashSet<String>();
//            for (int perm = 1; perm < 11; perm++) {
//                iterCtr = 0;
//                for (int iter = 0; iter < max; iter += 1000) {
//                    String iterStr = "Bin" + iter + "-" + (iter + 1000) + ".txt";
//                    String proxyStr = "AllSNPsPerm.txt-WithProxies.txt";
//                    String result = collect(perm, (iter + 1000),
//                            annotDir + "/"+proxyStr,
//                            snpdir + "/PermutationRound" + perm + "/" + iterStr,
//                            annotDir + "/PermutationRound" + perm + "/snpfunc.csv", false);
//                    iterCtr++;
//                    out.writeln(result);
//                }
//            }
//
//            out.close();
//            // average accross bins
//            average(snpdir + "/FuncAnnot-SNPInfo.txt", snpdir + "/FuncAnnot-SNPInfo-Summary.txt");



            String dir = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/RandomlySNPSelectionsForFuncAnalysis/";
            String annotDir = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/RandomlySNPSelectionsForFuncAnalysis/";

            String snpdir = dir;

            TextFile out = new TextFile(snpdir + "/FuncAnnot-SNPInfo.txt", TextFile.W);
            out.writeln("Perm\tBin\tsnps\ttfbs\tsplice1\tsplice2\tsplice3\tmirna1\tmirna2\tnssnp\tstop\treg\tcons");


            int max = 1000;
//
//            int iterCtr = 0;
            snpsAlreadyCounted = new HashSet<String>();

            String proxyfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/eQTLSNPFunctionAnnotation/transAnnot/SNPSelectionsForFuncAnalysis/RealData/AllSNPs.txt-WithProxies.txt";
            String queryfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/eQTLSNPFunctionAnnotation/transAnnot/SNPSelectionsForFuncAnalysis/RealData/Bin0-1000.txt";
            String annotationfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/eQTLSNPFunctionAnnotation/transAnnot/SNPSelectionsForFuncAnalysis/RealData/snpfunc.csv";
            String result = collect(0, 1000, proxyfile, queryfile, annotationfile, false);

            out.writeln(result);

            System.out.println("Actually counted: " + ctr);
            System.out.println("In HashMap: " + snpsAlreadyCounted.size());

            snpsAlreadyCounted = new HashSet<String>();
            for (int perm = 1; perm < 21; perm++) {
                proxyfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/RandomlySNPSelectionsForFuncAnalysis/Set"+(perm-1)+"/Set-"+(perm-1)+".txt-WithProxies.txt";
                queryfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/RandomlySNPSelectionsForFuncAnalysis/Set"+(perm-1)+"/Set-"+(perm-1)+".txt";
                annotationfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/RandomlySNPSelectionsForFuncAnalysis/Set"+(perm-1)+"/snpfunc.csv";
                result = collect(perm, 1000, proxyfile, queryfile, annotationfile, false);
                out.writeln(result);
            }

            out.close();
            // average accross bins
            average(snpdir + "/FuncAnnot-SNPInfo.txt", snpdir + "/FuncAnnot-SNPInfo-Summary.txt");

            // perform Chi-squared

//            chisquared("/Volumes/iSnackHD/Skydrive/Cis-SNPs/FuncAnnot.txt", "/Volumes/iSnackHD/Skydrive/Cis-SNPs/FuncAnnot-Summary.txt");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    int ctr = 0;

    public String collect(int perm, int iter, String proxyFile, String queryFile, String annotationFile, boolean returnFractions) throws IOException {
        TextFile queryf = new TextFile(queryFile, TextFile.R);
        HashSet<String> initialquerySNPs = new HashSet<String>();
        initialquerySNPs.addAll(queryf.readAsArrayList());
        queryf.close();

        // now read in all proxies, add them to querySNPs..
        // format query\tproxy
        TextFile proxyfile = new TextFile(proxyFile, TextFile.R);
        String[] proxyElems = proxyfile.readLineElems(TextFile.tab);
        HashSet<String> allUniqueSNPs = new HashSet<String>();
        while (proxyElems != null) {
            String q = proxyElems[0];
            String p = proxyElems[1];
            if (initialquerySNPs.contains(q)) {
                initialquerySNPs.add(p);
            }
            allUniqueSNPs.add(p);
            allUniqueSNPs.add(q);
            proxyElems = proxyfile.readLineElems(TextFile.tab);
        }
        proxyfile.close();

        System.out.println("Total Unique SNPs: " + allUniqueSNPs.size() + "\tsnps");

//        System.out.println("Annotating: "+querySNPs.size()+" snps");
        TextFile tf = new TextFile(annotationFile, TextFile.R);
        String[] header = tf.readLineElems(TextFile.comma);

//        ctr += querySNPs.size();

        HashSet<String> querySNPs = new HashSet<String>();
        for (String s : initialquerySNPs) {
            if (!snpsAlreadyCounted.contains(s)) {
                querySNPs.add(s);
            }
        }

        int snpscol = -1;
        int tfbscol = -1;
        int splicecol = -1;
        int splicesscol = -1;
        int spliceabolishdomain = -1;
        int mirnanda = -1;
        int mirnasanger = -1;
        int nssnpcol = -1;
        int stopcodon = -1;
        int regpotential = -1;
        int conservation = -1;
        /*
         * 
         * No.,
         * rs,
         * Chromosome,
         * Position,
         * Allele,
         * LDsnp,
         * Pop/LD,
         * TFBS,
         * Splicing(site),
         * Splicing(ESE or ESS),
         * Splicing(abolish domain),
         * miRNA(miRanda),
         * miRNA(Sanger),
         * nsSNP,
         * Stop Codon,
         * Polyphen,
         * SNPs3D(svm profile),
         * SNPs3D(svm structure),
         * RegPotential,
         * Conservation,
         * Nearby Gene,
         * Distance (bp)
         */
        for (int h = 0; h < header.length; h++) {
            String colname = header[h];
            if (colname.equals("rs")) {
                snpscol = h;
            } else if (colname.equals("TFBS")) {
                tfbscol = h;
            } else if (colname.equals("Splicing(site)")) {
                splicecol = h;
            } else if (colname.equals("Splicing(ESE or ESS)")) {
                splicesscol = h;
            } else if (colname.equals("Splicing(abolish domain)")) {
                spliceabolishdomain = h;
            } else if (colname.equals("miRNA(miRanda)")) {
                mirnanda = h;
            } else if (colname.equals("miRNA(Sanger)")) {
                mirnasanger = h;
            } else if (colname.equals("nsSNP")) {
                nssnpcol = h;
            } else if (colname.equals("Stop Codon")) {
                stopcodon = h;
            } else if (colname.equals("RegPotential")) {
                regpotential = h;
            } else if (colname.equals("Conservation")) {
                conservation = h;
            }
        }

        String[] elems = tf.readLineElems(TextFile.comma);
        HashSet<String> snpsVisited = new HashSet<String>();
        int snpsCounter = 0;
        int tfbsctr = 0;
        int splice1ctr = 0;
        int splice2ctr = 0;
        int splice3ctr = 0;
        int mirna1ctr = 0;
        int mirna2ctr = 0;
        int nssnpctr = 0;
        int stopctr = 0;
        double regsum = 0;
        double conssum = 0;
        while (elems != null) {

            String snp = elems[snpscol];

            if (querySNPs.contains(snp) && !snpsVisited.contains(snp)) {
//                if (!snpsAlreadyCounted.contains(snp)) {
                String tfbs = elems[tfbscol];
                String splice1 = elems[splicecol];
                String splice2 = elems[splicesscol];
                String splice3 = elems[spliceabolishdomain];
                String mirna1 = elems[mirnanda];
                String mirna2 = elems[mirnasanger];
                String nssnp = elems[nssnpcol];
                String stop = elems[stopcodon];
                String reg = elems[regpotential];
                String cons = elems[conservation];

//                System.out.println(snp + "\t" + reg + "\t" + cons);

                if (tfbs.equals("Y")) {
                    tfbsctr++;
                }
                if (splice1.equals("Y")) {
                    splice1ctr++;
                }
                if (splice2.equals("Y")) {
                    splice2ctr++;
                }
                if (splice3.equals("Y")) {
                    splice3ctr++;
                }

                if (mirna1.equals("Y")) {
                    mirna1ctr++;
                }
                if (mirna2.equals("Y")) {
                    mirna2ctr++;
                }

                if (nssnp.equals("Y")) {
                    nssnpctr++;
                }

                if (stop.equals("Y")) {
                    stopctr++;
                }

                if (!reg.equals("NA")) {
                    double dreg = Double.parseDouble(reg);
                    regsum += dreg;
                }


                if (!cons.equals("NA")) {
                    double dcons = Double.parseDouble(cons);
                    conssum += dcons;
                }

                snpsCounter++;
                snpsAlreadyCounted.add(snp);
//                }
                snpsVisited.add(snp);
            }




            elems = tf.readLineElems(TextFile.comma);
        }

        tf.close();
//        snpsCounter = querySNPs.size();
        String output = "";
        if (returnFractions) {
            output = (perm
                    + "\t" + iter
                    + "\t" + snpsCounter
                    + "\t" + (querySNPs.size() - snpsCounter)
                    + "\t" + ((double) tfbsctr / snpsCounter)
                    + "\t" + ((double) splice1ctr / snpsCounter)
                    + "\t" + ((double) splice2ctr / snpsCounter)
                    + "\t" + ((double) splice3ctr / snpsCounter)
                    + "\t" + ((double) mirna1ctr / snpsCounter)
                    + "\t" + ((double) mirna2ctr / snpsCounter)
                    + "\t" + ((double) nssnpctr / snpsCounter)
                    + "\t" + ((double) stopctr / snpsCounter)
                    + "\t" + ((double) regsum / snpsCounter)
                    + "\t" + ((double) conssum / snpsCounter));
        } else {
            output = (perm
                    + "\t" + iter
                    + "\t" + snpsCounter
                    + "\t" + (querySNPs.size() - snpsCounter)
                    + "\t" + tfbsctr
                    + "\t" + splice1ctr
                    + "\t" + splice2ctr
                    + "\t" + splice3ctr
                    + "\t" + mirna1ctr
                    + "\t" + mirna2ctr
                    + "\t" + nssnpctr
                    + "\t" + stopctr
                    + "\t" + regsum
                    + "\t" + conssum);
        }
//        System.out.println(output);
        return output;
    }

    private void average(String annotation, String outfilename) throws IOException {
        TextFile tf = new TextFile(annotation, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab); // header

        String[] elems = tf.readLineElems(TextFile.tab);

        // [bin][col]
        int[][] realData = new int[16][10];
        int[][] permutedData = new int[16][10];
        while (elems != null) {

            int perm = Integer.parseInt(elems[0]);

            int bin = (Integer.parseInt(elems[1]) / 1000) - 1;
            int nrSNPs = Integer.parseInt(elems[2]);
            int nrSNPsWoData = Integer.parseInt(elems[3]);
            int tfbs = Integer.parseInt(elems[4]);
            int splice1 = Integer.parseInt(elems[5]);
            int splice2 = Integer.parseInt(elems[6]);
            int splice3 = Integer.parseInt(elems[7]);
            int mirna1 = Integer.parseInt(elems[8]);
            int mirna2 = Integer.parseInt(elems[9]);
            int nssnp = Integer.parseInt(elems[10]);
            int stopctr = Integer.parseInt(elems[11]);
            double regsum = Double.parseDouble(elems[12]);
            double cons = Double.parseDouble(elems[13]);

            if (bin < 15) {
                if (perm > 0) {
//                perm -= 1;
                    permutedData[bin][0] += nrSNPs;
                    permutedData[bin][1] += nrSNPsWoData;
                    permutedData[bin][2] += tfbs;
                    permutedData[bin][3] += splice1;
                    permutedData[bin][4] += splice2;
                    permutedData[bin][5] += splice3;
                    permutedData[bin][6] += mirna1;
                    permutedData[bin][7] += mirna2;
                    permutedData[bin][8] += nssnp;
                    permutedData[bin][9] += stopctr;
                } else {
                    realData[bin][0] = nrSNPs;
                    realData[bin][1] = nrSNPsWoData;
                    realData[bin][2] = tfbs;
                    realData[bin][3] = splice1;
                    realData[bin][4] = splice2;
                    realData[bin][5] = splice3;
                    realData[bin][6] = mirna1;
                    realData[bin][7] = mirna2;
                    realData[bin][8] = nssnp;
                    realData[bin][9] = stopctr;
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile out = new TextFile(outfilename + "-Real.txt", TextFile.W);
        out.writeln(Strings.concat(header, Strings.tab, 0, header.length - 2));

        for (int bin = 0; bin < realData.length; bin++) {
            String real = Strings.concat(realData[bin], Strings.tab);
            out.writeln("RealData\t" + bin + "\t" + real);
        }


        out.writeln();

        for (int bin = 0; bin < realData.length; bin++) {
            String perm = Strings.concat(permutedData[bin], Strings.tab);
            out.writeln("PermData\t" + bin + "\t" + perm);
        }

        out.close();


        TextFile outchi = new TextFile(outfilename + "-ChiSquared.txt", TextFile.W);

        outchi.writeln("Class\tChi-Squared\tP-val\tFisherExact\tRealSNPsInClass\tRemainingRealSNPs\tPermSNPsInClass\tPermSNPsRemaining\tRealSNPsNotAnnotated\tPermSNPsNotAnnotated\tChi-Sauared(+unannotatedSNPs)\tPval(+unannotatedSNPs)");
        String[] funcAnnot = new String[]{
            "NrSNPs",
            "NrSNPsWoData",
            "TFBS",
            "Splice1",
            "Splice2",
            "Splice3",
            "MiRNA1",
            "MiRNA2",
            "nsSNP",
            "stop"
        };

        for (int func = 2; func < realData[0].length; func++) {
            double realSNPsInFunc = 0;
            double realTotalSNPs = 0;
            double permSNPsInFunc = 0;
            double permTotalSNPs = 0;

            int snpsRealNotAnnotated = 0;
            int snpsPermNotAnnotated = 0;

            for (int bin = 0; bin < realData.length; bin++) {
                snpsRealNotAnnotated += realData[bin][1];
                snpsPermNotAnnotated += permutedData[bin][1];
                realSNPsInFunc += realData[bin][func];
                realTotalSNPs += realData[bin][0];
                permSNPsInFunc += permutedData[bin][func];
                permTotalSNPs += permutedData[bin][0];

//                int realNrSNPs = ;
//                int permNrSNPs = ;
//                int realSNPsInFunc = ;
//                int permSNPsInFunc = ;
//
//                int remainingReal = realNrSNPs - realSNPsInFunc;
//                int remainingPerm = permNrSNPs - permSNPsInFunc;
//
//                if (realSNPsInFunc > 0 && permSNPsInFunc > 0 && remainingReal > 0 && remainingPerm > 0) {
//                    
//                    

//                    if (p < 0.001) {
//                        System.out.println(bin + "\t" + funcAnnot[func] + "\t" + p + "\t" + realSNPsInFunc + "\t" + remainingReal + "\t" + permSNPsInFunc + "\t" + remainingPerm);
//                    }
//                } else {
////                    System.out.println(bin + "\t" + funcAnnot[func] + "\tNA\t" + realSNPsInFunc + "\t" + remainingReal + "\t" + permSNPsInFunc + "\t" + remainingPerm);
//                }
            }

            double remainingReal = (realTotalSNPs - realSNPsInFunc);
            double remainingPerm = (permTotalSNPs - permSNPsInFunc);
//            if (realSNPsInFunc > 0 && permSNPsInFunc > 0 && remainingPerm > 0 && remainingReal > 0) {

            double xsq = ChiSquare.getX(realSNPsInFunc, remainingReal, permSNPsInFunc, remainingPerm);
            double xsq2 = ChiSquare.getX(realSNPsInFunc, remainingReal + snpsRealNotAnnotated, permSNPsInFunc, remainingPerm + snpsPermNotAnnotated);
            double p = ChiSquare.getP(1, xsq);
            double p2 = ChiSquare.getP(1, xsq2);
            FisherExactTest fet = new FisherExactTest();
            double pfisher = fet.getFisherPValue((int) realSNPsInFunc, (int) remainingReal, (int) permSNPsInFunc, (int) remainingPerm);
            outchi.writeln(funcAnnot[func] + "\t" + xsq + "\t" + p + "\t" + pfisher + "\t" + realSNPsInFunc + "\t" + remainingReal + "\t" + permSNPsInFunc + "\t" + remainingPerm + "\t" + snpsRealNotAnnotated + "\t" + snpsPermNotAnnotated + "\t" + xsq2 + "\t" + p2);
            System.out.println(func + "\t" + funcAnnot[func] + "\t" + xsq + "\t" + p + "\t" + pfisher);
        }
        outchi.close();
    }
//}
}

/*
 * (perm + "\t" + iter
 + "\t" + snpsCounter
 + "\t" + ((double)tfbsctr/snpsCounter)
 + "\t" + ((double)splice1ctr/snpsCounter)
 + "\t" + ((double)splice2ctr/snpsCounter)
 + "\t" + ((double)splice3ctr/snpsCounter)
 + "\t" + ((double)mirna1ctr/snpsCounter)
 + "\t" + ((double)mirna2ctr/snpsCounter)
 + "\t" + ((double)nssnpctr/snpsCounter)
 + "\t" + ((double)stopctr/snpsCounter)
 + "\t" + ((double)regsum/snpsCounter)
 + "\t" + ((double)conssum/snpsCounter))
 * 
 */
