/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harmjan
 */
public class eQTLFileDetermineNumberOfIndependentEffects {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            eQTLFileDetermineNumberOfIndependentEffects.determineNumberOfIndependentEffects(
                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/trans/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.5.txt",
                    "");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void determineNumberOfIndependentEffects(String eQTLFile, String snpReference) throws IOException {
        HashMap<String, ArrayList<String>> eQTLSNPsPerProbe = new HashMap<String, ArrayList<String>>();
        HashMap<String, ArrayList<String>> eQTLSNPsPerGene = new HashMap<String, ArrayList<String>>();
        HashMap<String, ArrayList<String>> eQTLProbesPerSNP = new HashMap<String, ArrayList<String>>();


        TextFile tf = new TextFile(eQTLFile, TextFile.R);
        tf.readLineElems(TextFile.tab);
        String[] elems = tf.readLineElemsReturnObjects(TextFile.tab);

        ArrayList<String> genes = new ArrayList<String>();
        ArrayList<String> snps = new ArrayList<String>();
        ArrayList<String> probes = new ArrayList<String>();

        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];
            String gene = elems[eQTLTextFile.HUGO];

            if (!gene.equals("-")) {
                ArrayList<String> eqtlSNPsForGene = eQTLSNPsPerGene.get(gene);
                if (eqtlSNPsForGene == null) {
                    eqtlSNPsForGene = new ArrayList<String>();
                    genes.add(gene);
                }
                eqtlSNPsForGene.add(snp);
                eQTLSNPsPerGene.put(gene, eqtlSNPsForGene);
            }

            ArrayList<String> eqtlSNPsForProbe = eQTLSNPsPerProbe.get(probe);
            if (eqtlSNPsForProbe == null) {
                eqtlSNPsForProbe = new ArrayList<String>();
                probes.add(probe);
            }
            eqtlSNPsForProbe.add(snp);
            eQTLSNPsPerProbe.put(probe, eqtlSNPsForProbe);

            ArrayList<String> probesForSNP = eQTLProbesPerSNP.get(snp);
            if (probesForSNP == null) {
                probesForSNP = new ArrayList<String>();
                snps.add(snp);
            }
            probesForSNP.add(probe);
            eQTLProbesPerSNP.put(snp, probesForSNP);

            elems = tf.readLineElemsReturnObjects(TextFile.tab);
        }
        tf.close();

//        TriTyperGenotypeData ds = new TriTyperGenotypeData();
//        ds.load(snpReference);

        int snpsWithMultipleEffects = 0;
        int probesWithMultipleEffects = 0;
        int genesWithMultipleEffects = 0;

        for (String s : snps) {
            ArrayList<String> al = eQTLProbesPerSNP.get(s);
            if (al != null && al.size() > 1) {
                snpsWithMultipleEffects++;
            }
        }

        for (String g : genes) {
            ArrayList<String> al = eQTLSNPsPerGene.get(g);
            if (al != null && al.size() > 1) {
                genesWithMultipleEffects++;
            }
        }

        for (String p : probes) {
            ArrayList<String> al = eQTLSNPsPerProbe.get(p);
            if (al != null && al.size() > 1) {
                probesWithMultipleEffects++;
            }
        }



        System.out.println("SNPs: " + snpsWithMultipleEffects);

        System.out.println("Probes: " + probesWithMultipleEffects);
        System.out.println("Genes: " + genesWithMultipleEffects);
    }
}
