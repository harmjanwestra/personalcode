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
public class eQTLFileProbeCounter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        // 

        try {
            //
            /*
             * 19708
             4533
             */
            TextFile tf = new TextFile("/Volumes/iSnackHD/Datasets/TestRum3/eQTLs.txt", TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            HashSet<String> probes = new HashSet<String>();
            HashSet<String> snps = new HashSet<String>();
            while (elems != null) {
                probes.add(elems[4]);
                snps.add(elems[1]);
                elems = tf.readLineElems(TextFile.tab);
            }

            System.out.println(probes.size());
            System.out.println(snps.size());
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
