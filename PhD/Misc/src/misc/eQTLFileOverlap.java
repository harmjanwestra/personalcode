/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class eQTLFileOverlap {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String f1 = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05-HT12v3.txt.gz";
            String f2 = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGSKoraMeta/MetaOut/eQTLs.txt";

            TextFile tf = new TextFile(f1, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            HashSet<Pair<String, String>> e1 = new HashSet<Pair<String, String>>();
            HashSet<Pair<String, String>> e2 = new HashSet<Pair<String, String>>();
            while (elems != null) {
                e1.add(new Pair<String, String>(elems[1], elems[4]));
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();


            TextFile tf2 = new TextFile(f2, TextFile.R);
            tf2.readLine();
            elems = tf2.readLineElems(TextFile.tab);
            while (elems != null) {
                e2.add(new Pair<String, String>(elems[1], elems[4]));
                elems = tf2.readLineElems(TextFile.tab);
            }
            tf.close();

            for(Pair<String, String> p: e1){
                if(!e2.contains(p)){
                    System.out.println(p.toString()+" not found in f2");
                }
            }
            
        } catch (IOException e) {
            e.printStackTrace();

        }
    }
}
