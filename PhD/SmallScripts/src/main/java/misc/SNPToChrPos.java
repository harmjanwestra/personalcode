/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
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
public class SNPToChrPos {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String efile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-11-HT12v3SNPProbeCombos.txt";
        String efileOut = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-11-HT12v3SNPProbeCombos-SNPPositionsB36.txt";
        String snpmapfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-27-SNPMappings-dbSNP130.txt.gz";
        try {
            HashMap<String, String> snpMap = new HashMap<String, String>();
            TextFile tf = new TextFile(efile, TextFile.R);
            String[] elems1 = tf.readLineElems(TextFile.tab);
            while (elems1 != null) {
                String snp = elems1[0];
                snpMap.put(snp, null);
                elems1 = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile tf2 = new TextFile(snpmapfile, TextFile.R);
            String[] elems2 = tf2.readLineElems(TextFile.tab);
            while (elems2 != null) {
                if (snpMap.containsKey(elems2[2])) {
                    snpMap.put(elems2[2], elems2[0] + "\t" + elems2[1]);
                }
                elems2 = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            tf = new TextFile(efile, TextFile.R);
            TextFile tfOut = new TextFile(efileOut, TextFile.W);
            elems1 = tf.readLineElems(TextFile.tab);
            while (elems1 != null) {
                String map = snpMap.get(elems1[0]);
                if (map == null) {
                    tfOut.writeln(elems1[0] + "\t" + null + "\t" + null);
                } else {
                    tfOut.writeln(elems1[0] + "\t" + map);
                }

                elems1 = tf.readLineElems(TextFile.tab);
            }
            tfOut.close();
            tf.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
