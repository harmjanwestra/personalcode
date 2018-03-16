/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class DetermineNumberOfCisFXForSignificantSNPs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String cisFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt-FilteredForProbeLevelFDR.txt.gz";
            String dir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPTransEQTLEnrichment/";
            String reference = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
            String dontcountthesesnpsinrealdata = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPTransEQTLEnrichment/NotInReference.txt";
            TriTyperGenotypeData ds = new TriTyperGenotypeData(reference);
            SNPLoader loader = ds.createSNPLoader();
            DetermineLD ldcalc = new DetermineLD();

            HashSet<String> cisSNPs = new HashSet<String>();
            HashSet<String> topCisSNPs = new HashSet<String>();
            HashSet<String> visitedProbes = new HashSet<String>();
            HashMap<String, String> topSNPPerProbe = new HashMap<String, String>();

            HashMap<String, HashSet<String>> snpsPerProbe = new HashMap<String, HashSet<String>>();
            HashMap<String, HashSet<String>> probesPerSNP = new HashMap<String, HashSet<String>>();

            TextFile tf1 = new TextFile(cisFile, TextFile.R);
            tf1.readLine();
            String[] ciselems = tf1.readLineElems(TextFile.tab);

            // determine top effect per probe
            while (ciselems != null) {
                String probe = ciselems[4];
                HashSet<String> snpsForProbe = snpsPerProbe.get(probe);
                if (snpsForProbe == null) {
                    snpsForProbe = new HashSet<String>();
                }
                snpsForProbe.add(ciselems[1]);
                snpsPerProbe.put(probe, snpsForProbe);
                if (!visitedProbes.contains(probe)) {
                    topCisSNPs.add(ciselems[1]);
                    topSNPPerProbe.put(probe, ciselems[1]);
                    visitedProbes.add(probe);
                }

                HashSet<String> probesSNP = probesPerSNP.get(ciselems[1]);
                if (probesSNP == null) {
                    probesSNP = new HashSet<String>();
                }
                probesSNP.add(probe);
                probesPerSNP.put(ciselems[1], probesSNP);

                cisSNPs.add(ciselems[1]);
                ciselems = tf1.readLineElems(TextFile.tab);
            }
            tf1.close();

            System.out.println(cisSNPs.size());
            System.out.println("");

            // determine LD between cis-eQTL SNPs.
            HashMap<String, HashSet<String>> independentCisFx = new HashMap<String, HashSet<String>>();
            for (String probe : visitedProbes) {
                HashSet<String> cisEQTLSNPs = snpsPerProbe.get(probe);
                HashSet<String> independentFx = new HashSet<String>();
                String[] snparr = cisEQTLSNPs.toArray(new String[0]);

                String topSNP = topSNPPerProbe.get(probe);
                Integer snp1Id = ds.getSnpToSNPId().get(topSNP);
                if (snp1Id != null) {
                    SNP snp1Obj = ds.getSNPObject(snp1Id);
                    loader.loadGenotypes(snp1Obj);
                    for (int j = 0; j < snparr.length; j++) {
                        Integer snp2Id = ds.getSnpToSNPId().get(snparr[j]);
                        if (snp2Id != null) {
                            SNP snp2Obj = ds.getSNPObject(snp2Id);
                            loader.loadGenotypes(snp2Obj);

                            double r2 = ldcalc.getRSquared(snp1Obj, snp2Obj, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                            if (!Double.isNaN(r2)) {
                                if (r2 >= 0.95) {
                                    independentFx.add(snparr[j]);
                                }
                            }

                            snp2Obj.clearGenotypes();
                        }
                    }
                    snp1Obj.clearGenotypes();
                }
                independentCisFx.put(probe, independentFx);
            }


            System.out.println(cisSNPs.size() + " cis eQTLs");
            System.out.println(topCisSNPs.size() + " cis eQTL top SNPs for " + visitedProbes.size() + " probes.");

            for (int iter = 0; iter < 101; iter++) {
                String transfile = "";
                HashSet<String> dontcounttheseSNPs = null;
                if (iter == 0) {
                    dontcounttheseSNPs = new HashSet<String>();
                    TextFile tfq = new TextFile(dontcountthesesnpsinrealdata, TextFile.R);
                    dontcounttheseSNPs.addAll(tfq.readAsArrayList());
                    tfq.close();
                    transfile = dir + "querySNPs.txt";
                } else {
                    transfile = dir + "Set-" + (iter - 1) + ".txt";
                }

                HashSet<String> snps = new HashSet<String>();
                TextFile tf = new TextFile(transfile, TextFile.R);
                tf.readLine();
                String[] elems = tf.readLineElems(TextFile.tab);
                int nrTop = 0;
                int nrCis = 0;
                int nrTot = 0;
                HashSet<String> uniqueCisSNPs = new HashSet<String>();
                HashSet<String> uniqueTopSNPs = new HashSet<String>();
                while (elems != null) {
                    String snp = elems[0];
                    if (dontcounttheseSNPs == null || !dontcounttheseSNPs.contains(snp)) {
                        boolean isTopSNP = false;
                        if (cisSNPs.contains(snp)) {
                            nrCis++;
                            uniqueCisSNPs.add(snp);
                            // determine probes to which this snp is associated
                            HashSet<String> probes = probesPerSNP.get(snp);
                            for (String probe : probes) {
                                // get independent effects per probe
                                HashSet<String> fx = independentCisFx.get(probe);
                                if (fx.contains(snp)) {
                                    isTopSNP = true;
                                    uniqueTopSNPs.add(snp);
                                }
                            }
                        }

                        if (isTopSNP) {
                            nrTop++;
                        }
                        nrTot++;
                    }

                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
                System.out.println(iter + "\t" + uniqueCisSNPs.size() + "\t" + uniqueTopSNPs.size());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
