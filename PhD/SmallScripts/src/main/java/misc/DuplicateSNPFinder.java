/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
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
public class DuplicateSNPFinder {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            TextFile tf = new TextFile("/Volumes/iSnackHD/Data/Projects/PathwayQTL/2013-11-20-PathwayQTLInitialTestRun-Reactome+-GWASSNPs-EGCUTMeta/-MetaAnalysis-RowNames.txt.gz", TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            HashSet<String> uniqueSNPs = new HashSet<String>();
            int nrlines = 0;
            while (elems != null) {
                String snp = elems[0];
                if (uniqueSNPs.contains(snp)) {
                    System.err.println("Found snp twice: " + snp);
                }
                uniqueSNPs.add(snp);
                nrlines++;
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            System.out.println(nrlines + "\t" + uniqueSNPs.size());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
