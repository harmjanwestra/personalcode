/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ProbeMapCtr {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            TextFile tf = new TextFile("/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-HT12v4.txt", TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            int ctr = 0;
            int ctr2 = 0;
            while (elems != null) {
                String chr = elems[3];
                if (chr.equals("-")) {
                    ctr++;
                    System.out.println(Strings.concat(elems, Strings.tab));
                } else {
                    ctr2++;
                }
                elems = tf.readLineElems(TextFile.tab);
            }

            System.out.println(ctr + "\t" + ctr2);
            tf.close();
        } catch (IOException e) {
        }
    }
}
