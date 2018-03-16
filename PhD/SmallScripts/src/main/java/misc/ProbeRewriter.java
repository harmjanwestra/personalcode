/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ProbeRewriter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {


            HashMap<String, String> probeToEns = new HashMap<String, String>();
            TextFile ptbt = new TextFile("/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/All.0.96%Identity-Merged-PerProbe-UniqueMappings-Ensembl70.txt", TextFile.R);
            String[] elems = ptbt.readLineElems(TextFile.tab);

            while (elems != null) {
                String probe = elems[0];
                String ens = elems[elems.length - 2];
                probeToEns.put(probe, ens);
                elems = ptbt.readLineElems(TextFile.tab);
            }
            ptbt.close();


            TextFile tfOut = new TextFile("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05-EnsB70.txt", TextFile.W);
            TextFile tfOut1lpg = new TextFile("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05-EnsB70-OneLinePerGene.txt", TextFile.W);
            TextFile tf = new TextFile("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt", TextFile.R);
            tfOut.writeln(tf.readLine());
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String probe = elems[4];
                String ens = probeToEns.get(probe);
                if (ens == null) {
                    ens = elems[4];
                }
                String[] enselems = ens.split(";");

                for (int i = 0; i < enselems.length; i++) {
                    elems[4] = enselems[i];
                    String output = Strings.concat(elems, Strings.tab);
                    tfOut.writeln(output);
                }

                elems[4] = ens;
                String output = Strings.concat(elems, Strings.tab);
                tfOut1lpg.writeln(output);


                elems = tf.readLineElems(TextFile.tab);
            }
            tf.readLine();
            tf.close();
            tfOut1lpg.close();
            tfOut.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
