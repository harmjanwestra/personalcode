/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class ProxyFinderForASE {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        String[] referenceDS = new String[]{"/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap2r24-CEU/", "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/Merged/"};
        String[] referenceDSNames = new String[]{"HapMap2CEU", "100GenomesCEU"};
        String vectorFile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/Vector-CellTypeInteractionZScore.txt";
        String output = "/Volumes/iSnackHD/AeroFS/2013-12-09-ProxiesForASE/";

        ProxyFinderForASE f = new ProxyFinderForASE(referenceDS, referenceDSNames, vectorFile, output);
        try {
            f.run();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
    private final String[] referenceDS;
    private final String vectorFile;
    private final String output;
    private final String[] referenceDSNames;

    private ProxyFinderForASE(String[] referenceDS, String[] referenceDSNames, String vectorFile, String output) {
        this.referenceDS = referenceDS;
        this.referenceDSNames = referenceDSNames;
        this.vectorFile = vectorFile;
        this.output = output;
    }

    private void run() throws IOException {

        HashSet<String> querySNPs = new HashSet<String>();
        TextFile tf = new TextFile(vectorFile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String eqtl = elems[0];
            String snp = eqtl.split("-")[0];
            querySNPs.add(snp);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        System.out.println(querySNPs.size() + " unique SNPs in " + vectorFile);

        DetermineLD ldcalc = new DetermineLD();

        for (int d = 0; d < referenceDS.length; d++) {
            TriTyperGenotypeData ds = new TriTyperGenotypeData(referenceDS[d]);
            SNPLoader loader = ds.createSNPLoader();

            double threshold = 0.99;
            int distance = 10000000;

            TextFile outfile = new TextFile(output + referenceDSNames[d] + "-Proxies-" + threshold + ".txt", TextFile.W);
            System.out.println("Writing output to: " + output + referenceDSNames[d] + "-Proxies-" + threshold + ".txt");
            outfile.writeln("QuerySNP\tChr\tChrPos\tProxies\tR-squared\tProxyPos");

            HashSet<String> doNotTestTheseSNPsAgain = new HashSet<String>();
            for (int chr = 1; chr < 26; chr++) {

                // put all possible SNPs in a hashMap.
                HashSet<Integer> possibleProxies = new HashSet<Integer>();
                for (int s = 0; s < ds.getSNPs().length; s++) {
                    if (ds.getChr(s) == chr) {
                        possibleProxies.add(s);
                    }
                }
                System.out.println(possibleProxies.size() + "\tSNPs on chr:\t" + chr + "\tin reference:\t" + referenceDSNames[d]);

                // get the chromosomes of the query snps..
                for (String snp : querySNPs) {
                    if (!doNotTestTheseSNPsAgain.contains(snp)) {
                        Integer SNPId = ds.getSnpToSNPId().get(snp);
                        if (SNPId != null) {
                            if (ds.getChr(SNPId) == chr) {
                                SNP snpObj1 = ds.getSNPObject(SNPId);
                                loader.loadGenotypes(snpObj1);
                                ArrayList<Integer> proxies = new ArrayList<Integer>();
                                ArrayList<Double> proxiesLD = new ArrayList<Double>();
                                for (Integer s : possibleProxies) {
                                    SNP snpObj2 = ds.getSNPObject(s);
                                    if (Math.abs(snpObj2.getChrPos() - snpObj1.getChrPos()) < distance) {
                                        loader.loadGenotypes(snpObj2);
                                        double r2 = ldcalc.getRSquared(snpObj1, snpObj2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                        if (r2 >= threshold) {
                                            proxies.add(s);
                                            proxiesLD.add(r2);
                                        }
                                        snpObj2.clearGenotypes();
                                    }
                                }
                                snpObj1.clearGenotypes();

                                String outputline = snp + "\t"
                                        + ChrAnnotation.parseByte(snpObj1.getChr()) + "\t"
                                        + +snpObj1.getChrPos();

                                if (proxies.isEmpty()) {
                                    outputline += "\tNo proxies found\t-\t-";
                                } else {
                                    String output1 = "";
                                    String output2 = "";
                                    String output3 = "";
                                    for (int result = 0; result < proxies.size(); result++) {
                                        if (result == 0) {
                                            output1 += ds.getSNPs()[proxies.get(result)];
                                            output2 += ds.getChrPos(proxies.get(result));
                                            output3 += proxiesLD.get(result);
                                        } else {
                                            output1 += ";" + ds.getSNPs()[proxies.get(result)];
                                            output2 += ";" + ds.getChrPos(proxies.get(result));
                                            output3 += ";" + proxiesLD.get(result);
                                        }
                                    }
                                    outputline += "\t" + output1 + "\t" + output3 + "\t" + output2;
                                }
                                doNotTestTheseSNPsAgain.add(snp);
                                outfile.writeln(outputline);
                            }

//                            else {
//                                System.out.println(chr + "\t" + ds.getChr(SNPId) + "\t" + (ds.getChr(SNPId) == chr));
//                            }
                        } else {
                            doNotTestTheseSNPsAgain.add(snp);
                            outfile.writeln(snp + "\tSNP not in reference dataset\t-\t-");
                        }
                    }
                }
            }
            outfile.close();
        }
    }

}
