/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.qc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harm-jan
 */
public class TransEQTLFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {
            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/trans/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/QC/alignmentsummary.txt", TextFile.R);

            String[] elems = tf.readLineElems(TextFile.tab);

            HashSet<Pair<String, String>> falsepositives = new HashSet<Pair<String, String>>();

            while (elems != null) {
                String probe = elems[0].replace(">", "");
                String snp = elems[1];
                String id = elems[11];
                id = id.replace("Id: ", "");
                Integer idint = Integer.parseInt(id);
                if (idint >= 15) {
                    falsepositives.add(new Pair<String, String>(snp, probe));
                }
                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();

            TextFile eQTLfile = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/trans/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/eQTLsFDR0.5.txt", TextFile.R);

            TextFile eQTLfileOut = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/trans/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/QC/eQTLsFDR0.5.txt", TextFile.W);
            elems = eQTLfile.readLineElems(TextFile.tab);
            eQTLfileOut.writeln(Strings.concat(elems, Strings.tab));
            elems = eQTLfile.readLineElems(TextFile.tab);
            HashSet<String> snps = new HashSet<String>();
            HashSet<String> probes = new HashSet<String>();
            HashSet<String> genes = new HashSet<String>();

            HashSet<String> fpprobes = new HashSet<String>();
            HashSet<String> fpsnps = new HashSet<String>();
            while (elems != null) {

                String probe = elems[4];
                String snp = elems[1];
                String gene = elems[16];

                if (falsepositives.contains(new Pair<String, String>(snp, probe))) {
                    System.out.println("FALSE POSITIVE: " + Strings.concat(elems, Strings.tab));
                    fpprobes.add(probe);
                    fpsnps.add(snp);
                } else {
                    snps.add(snp);
                    probes.add(probe);
                    genes.add(gene);
                    eQTLfileOut.writeln(Strings.concat(elems, Strings.tab));
                }

                elems = eQTLfile.readLineElems(TextFile.tab);
            }

            eQTLfile.close();
            eQTLfileOut.close();
            System.out.println(probes.size());
            System.out.println(snps.size());
            System.out.println(genes.size());

            System.out.println("");
            System.out.println(fpprobes.size());
            System.out.println(fpsnps.size());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
