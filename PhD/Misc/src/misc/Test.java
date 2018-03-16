/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class Test {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            HashSet<String> allowedSNPs = new HashSet<String>();
            TextFile tf3 = new TextFile("/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLude.txt", TextFile.R);
            allowedSNPs.addAll(tf3.readAsArrayList());
            tf3.close();

//            TextFile tf1 = new TextFile("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGSFixedKoraMeta/Compare/SNPsForPruning.txt", TextFile.W);
//            tf1.writeln("P\tSNP");

            for (int i = 0; i < 100; i++) {
                TextFile tf = new TextFile("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPCisEQTLEnrichment/Set-"+i+".txt", TextFile.R);
                String ln = tf.readLine();
                int ctr = 0;
                while (ln != null) {
                    if (allowedSNPs.contains(ln)) {
//                    tf1.writeln("1E-9\t"+ln);
                        ctr++;
                    }
                    ln = tf.readLine();
                }
                tf.close();
//            tf1.close();

                System.out.println(ctr);
            }
        } catch (IOException e) {
        }
    }
}
