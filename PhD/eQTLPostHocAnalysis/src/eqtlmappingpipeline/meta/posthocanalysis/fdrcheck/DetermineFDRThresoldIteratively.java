/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.fdrcheck;

import eqtlmappingpipeline.metaqtl3.FDR;
import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harmjan
 */
public class DetermineFDRThresoldIteratively {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
    }

    public void determineFDR(int nrMaxPermutations, String indir, int nrMaxEQTLs, double threshold) throws IOException {
        for (int perm = 1; perm < nrMaxPermutations; perm++) {
            FDR f = new FDR();
            FDR.calculateFDR(indir, perm, nrMaxEQTLs, 0.05, false, null, null);
            // determine lowest P-value from eQTLsFDR.txt.gz
            TextFile tf = new TextFile(indir + "eQTLsFDR.txt.gz", TextFile.R);
            String header = tf.readLine();

            String[] elems = tf.readLineElems(TextFile.tab);
            double pvalue = 1;
            double z = 0;
            while (elems != null) {
                String fdrstr = elems[elems.length - 1];
                double fdrdbl = Double.parseDouble(fdrstr);
                if (fdrdbl > threshold) {
                    break;
                } else {
                    pvalue = Double.parseDouble(elems[0]);
                    z = Double.parseDouble(elems[eQTLTextFile.METAZ]);
                }

                elems = tf.readLineElems(TextFile.tab);
            }
            System.out.println("Perm:\t" + perm + "\t" + pvalue + "\t" + z);
            tf.close();

        }
    }
}
