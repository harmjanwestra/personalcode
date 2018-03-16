/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.replication;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ReplicationCheck {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            HashSet<String> snps = new HashSet<String>();
            HashSet<String> probes = new HashSet<String>();
            HashSet<Pair<String, String>> pairs = new HashSet<Pair<String, String>>();


            TextFile tf1 = new TextFile("/Volumes/iSnackHD/SNPProbe-HT12v4.txt", TextFile.R); // snp/probe file
            String[] elems = tf1.readLineElems(TextFile.tab); // header
            elems = tf1.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[0];
                String probe = elems[1];
                snps.add(snp);
                probes.add(probe);
                Pair<String, String> p = new Pair<String, String>(snp, probe);
                elems = tf1.readLineElems(TextFile.tab);
            }
            tf1.close();

            TextFile tf2 = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/MonocyteReplication/SubsetOfProbesUsedDuringNormalization/eQTLs.txt", TextFile.R);
            String[] elemstf2 = tf2.readLineElems(TextFile.tab); // header
            elemstf2 = tf2.readLineElems(TextFile.tab);
            HashSet<String> snpsIneQTLFile = new HashSet<String>();
            HashSet<String> probesIneQTLFile = new HashSet<String>();
            HashSet<Pair<String, String>> pairsIneQTLFile = new HashSet<Pair<String, String>>();
            while (elemstf2 != null) {

                String snp = elemstf2[1];
                String probe = elemstf2[4];
                Pair<String, String> p = new Pair<String, String>(snp, probe);
                snpsIneQTLFile.add(snp);
                probesIneQTLFile.add(probe);
                pairsIneQTLFile.add(p);
                elemstf2 = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            int snpsshared = 0;
            int probesshared = 0;
            int pairsshared = 0;
            HashSet<String> snpsThatAreNotSharedBetweenBothSets = new HashSet<String>();
            for (String snp : snpsIneQTLFile) {
                if (snps.contains(snp)) {
                    snpsshared++;
                } else {
                    snpsThatAreNotSharedBetweenBothSets.add(snp);
                }
            }

            for (String probe : probesIneQTLFile) {
                if (probes.contains(probe)) {
                    probesshared++;
                } else {
                    System.out.println(probe);
                }
            }

            for (Pair<String, String> p : pairsIneQTLFile) {
                if (pairs.contains(p)) {
                    pairsshared++;
                }
            }

            System.out.println("Shared probes: " + probesshared + " out of: " + probes.size() + " - " + probesIneQTLFile.size());
            System.out.println("Shared snps: " + snpsshared + " out of: " + snps.size() + " - " + snpsIneQTLFile.size());
            System.out.println("Shared pairs: " + pairsshared + " out of: " + pairs.size() + " - " + pairsIneQTLFile.size());

            
            // determine what happened to the snp
            TextFile tfsnp = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/MonocyteReplication/SubsetOfProbesUsedDuringNormalization/excludedSNPsBySNPProbeCombinationFilter.txt", TextFile.R);
            elems = tfsnp.readLineElems(TextFile.tab);
            int excludedByFilter = 0;
            while(elems!=null){
                String snp = elems[0];
                if(snpsThatAreNotSharedBetweenBothSets.contains(snp)){
                    excludedByFilter++;
                }
                elems = tfsnp.readLineElems(TextFile.tab);
            }
            tfsnp.close();
            
            System.out.println("Excluded by snp-probe combination filter: "+excludedByFilter);
            
        } catch (IOException e) {
            e.printStackTrace();

        }
    }
}
