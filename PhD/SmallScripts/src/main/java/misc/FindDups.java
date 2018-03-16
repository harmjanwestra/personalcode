/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class FindDups {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            TextFile tf = new TextFile("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/EGCUT1000GImputed/SNPs.txt", TextFile.R);
            HashMap<String, Integer> snptoPos = new HashMap<String, Integer>();
            String ln = tf.readLine();
            int ctr = 0;
            while (ln != null) {
                String snp = ln;
                Integer pos = snptoPos.get(snp);
                if (pos == null) {
                    snptoPos.put(snp, ctr);
                } else {
                    System.out.println("Found this before: " + snp + "\t" + pos);
                }

                ctr++;
                ln = tf.readLine();
            }
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
