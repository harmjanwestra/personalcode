/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class DetermineSNPsHavingLymphocyteAndNeutrophilSpecificEffects {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String file = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-09-27-FDREstimates/CellTypeInteractionZScore-null.txt";
        try {
            TextFile tf = new TextFile(file, TextFile.R);
            String header = tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            HashSet<String> signNeg = new HashSet<String>();
            HashSet<String> signPos = new HashSet<String>();
            while (elems != null) {
                String snp = elems[0].split("-")[0];
                Double Z = Double.parseDouble(elems[2]);
                Double f = Double.parseDouble(elems[4]);
                if (f < 0.05) {
                    if (Z > 0) {
                        signPos.add(snp);
                    } else {
                        signNeg.add(snp);
                    }
                }
                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();

            System.out.println("SNPs in both:");
            int ctr = 0;
            for (String snp : signNeg) {
                if (signPos.contains(snp)) {
                    System.out.println(snp);
                    ctr++;
                }
            }
            System.out.println(ctr + " in total");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
