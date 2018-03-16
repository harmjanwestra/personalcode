/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.cis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

/**
 *
 * @author harmjan
 */
public class GetPerfectProxiesForSNPs {

    private static TriTyperGenotypeData d;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String reference = "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/Merged/";
        try {
//            String indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC-SNPsForFuncAnnot/RealData/";
//            for (int i = 0; i < 9000; i += 1000) {
//                String fileIn = indir + "Bin" + i + "-" + (i + 1000) + ".txt";
//                GetPerfectProxiesForSNPs.collect(fileIn,
//                        reference,
//                        fileIn + "-WithProxies.txt");
//            }
//
//            for (int perm = 1; perm < 11; perm++) {
//                indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC-SNPsForFuncAnnot/PermutationRound" + perm + "/";
//                for (int i = 0; i < 9000; i += 1000) {
//                    String fileIn = indir + "Bin" + i + "-" + (i + 1000) + ".txt";
//                    GetPerfectProxiesForSNPs.collect(fileIn,
//                            reference,
//                            fileIn + "-WithProxies.txt");
//                }
//
//            }

//            for (int perm = 0; perm < 11; perm++) {
//
//                String indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/SNPSelectionsForFuncAnalysis/";
//                String fileIn = indir + "PermutationRound" + perm + "/PermutationRound" + perm + ".txt";
//                if (perm == 0) {
//                    fileIn = indir + "RealData/RealData.txt";
//                }
////                GetPerfectProxiesForSNPs.collect(fileIn,
////                        reference,
////                        fileIn + "-WithProxies.txt");
//
//                GetPerfectProxiesForSNPs.convertToList(fileIn + "-WithProxies.txt", fileIn + "-WithProxies-List.txt");
//
//            }

            String fileIn = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/SNPsForFuncAnalysisGWASSignificantUnpruned/querySNPs.txt";
            GetPerfectProxiesForSNPs.collect(fileIn,
                        reference,
                        fileIn + "-WithProxies.txt");

                GetPerfectProxiesForSNPs.convertToList(fileIn + "-WithProxies.txt", fileIn + "-WithProxies-List.txt");
            
//             for (int perm = 0; perm < 10; perm++) {
//
//                String indir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/RandomlySNPSelectionsForFuncAnalysisSet2";
//                String fileIn = indir + "/Set-" + perm + ".txt";
////                if (perm == 0) {
////                    fileIn = indir + "RealData/RealData.txt";
////                }
//                GetPerfectProxiesForSNPs.collect(fileIn,
//                        reference,
//                        fileIn + "-WithProxies.txt");
//
//                GetPerfectProxiesForSNPs.convertToList(fileIn + "-WithProxies.txt", fileIn + "-WithProxies-List.txt");
//
//            }
//            
//            
//            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/SNPSelectionsForFuncAnalysis/AllSNPsPerm.txt", TextFile.W);
//            HashSet<String> allSNPs = new HashSet<String>();
//            for (int perm = 1; perm < 11; perm++) {
//                indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/SNPSelectionsForFuncAnalysis/PermutationRound" + perm + "/";
//                fileIn = indir + "AllSNPs.txt";
//                TextFile snpin = new TextFile(fileIn, TextFile.R);
//                allSNPs.addAll(snpin.readAsArrayList());
//                snpin.close();
//            }
//            ArrayList<String> list = new ArrayList<String>();
//            list.addAll(allSNPs);
//            tf.writeList(list);
//            tf.close();

//            String fileIn = "/Volumes/iSnackHD/SkyDrive/SNPFunctionalAnnotation/Trans-SNPs/Real/AllSNPs.txt";
//            GetPerfectProxiesForSNPs.collect(fileIn,
//                        reference,
//                        fileIn + "-WithProxies.txt");

//            for (int perm = 5; perm < 11; perm++) {
//                String indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/SNPSelectionsForFuncAnalysis/PermutationRound" + perm + "/";
//                String fileIn = indir + "AllSNPs.txt";
//                
//
//            }

        } catch (IOException es) {
            es.printStackTrace();
        }


    }

    public static void collect(String snplist, String reference, String out) throws IOException {
        System.out.println("Collecting proxies for: " + snplist);
        System.out.println("Writing to: " + out);
        System.out.println("Reference: " + reference);
        // read SNPs from file
        TextFile tf = new TextFile(snplist, TextFile.R);
        HashSet<String> inputSNPs = new HashSet<String>();
        inputSNPs.addAll(tf.readAsArrayList());
        tf.close();

        System.out.println(inputSNPs.size() + " unique SNPs to query");
        TextFile tfout = new TextFile(out, TextFile.W);

        // load the genotype data
        if (d == null) {
            d = new TriTyperGenotypeData();
            d.load(reference);
        }
        String[] snps = d.getSNPs();

        SortableSNP[] sortedSNPs = new SortableSNP[snps.length];

        for (int s = 0; s < snps.length; s++) {
            sortedSNPs[s] = new SortableSNP(snps[s], s, d.getChr(s), d.getChrPos(s), SortableSNP.SORTBY.CHRPOS);
        }

        // sort the snps
        System.out.println("Sorting SNPs");
        java.util.Arrays.sort(sortedSNPs);
        String[] inputSNPArr = inputSNPs.toArray(new String[0]);



        int snpsWithProxies = 0;
        int snpsNotInRef = 0;
        int snpsWithoutProxies = 0;

        ExecutorService executor = Executors.newFixedThreadPool(16);
        CompletionService<Pair<String, ArrayList<String>>> compService = new ExecutorCompletionService<Pair<String, ArrayList<String>>>(executor);

        // now find the proxies for the inputSNPs
        System.out.println("Submitting tasks");
        for (int insnp = 0; insnp < inputSNPArr.length; insnp++) {
            String query = inputSNPArr[insnp];
            CalculateLDTask t = new CalculateLDTask(d, query, sortedSNPs);
            compService.submit(t);
        }
        System.out.println("Done submitting tasks");

        System.out.println("Retrieving results");
        for (int insnp = 0; insnp < inputSNPArr.length; insnp++) {
            try {
                Pair<String, ArrayList<String>> future = compService.take().get();
                String query = future.getLeft();
                ArrayList<String> result = future.getRight();
                tfout.writeln(query + "\t" + query);
                if (result == null) {
                    System.out.println(insnp + "/" + inputSNPArr.length + "\t" + query + "\thas\t" + 0 + "\tproxies");
                } else {
                    System.out.println(insnp + "/" + inputSNPArr.length + "\t" + query + "\thas\t" + result.size() + "\tproxies");
                    for (String s : result) {
                        tfout.writeln(query + "\t" + s);
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        executor.shutdown();
        tfout.close();

        System.out.println("SNPs with proxies: " + snpsWithProxies);
        System.out.println("SNPs without proxies: " + snpsWithoutProxies);
        System.out.println("SNPs not in ref: " + snpsNotInRef);
        System.out.println("Done.");
        System.out.println("");
    }

    private static void convertToList(String fileIn, String fileOut) throws IOException {
        TextFile in = new TextFile(fileIn, TextFile.R);
        String[] elems = in.readLineElems(TextFile.tab);
        HashSet<String> snps = new HashSet<String>();
        while (elems != null) {
            snps.add(elems[0]);
            snps.add(elems[1]);
            elems = in.readLineElems(TextFile.tab);
        }
        in.close();

        TextFile out = new TextFile(fileOut, TextFile.W);
        out.writeList(Arrays.asList(snps.toArray(new String[0])));
        out.close();

        TextFile out2 = new TextFile(fileOut + "-SNPNexus.txt", TextFile.W);
        for (String snp : snps) {
            out2.writeln("dbsnp\t" + snp);
        }
        out2.close();
    }
}
