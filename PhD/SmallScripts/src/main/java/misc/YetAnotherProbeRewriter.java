/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;

/**
 *
 * @author harmjan
 */
public class YetAnotherProbeRewriter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String pbAnnot = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";
        String inAnnot = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-11-HT12v3SNPProbeCombos.txt";
        String outAnnot = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-11-HT12v4SNPProbeCombos.txt";

        String inplatform = "HT12v3.txt";
        String outplatform = "HT12v4.txt";

        try {
            ProbeTranslation pbt = new ProbeTranslation();

            HashMap<String, String> map = pbt.getProbeTranslation(pbAnnot, inplatform, outplatform);

            TextFile tf = new TextFile(inAnnot, TextFile.R);
            TextFile tfout = new TextFile(outAnnot, TextFile.W);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                String probe = elems[1];
                String other = map.get(probe);
                if (other != null) {
                    tfout.writeln(elems[0] + "\t" + other);
                }

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            tfout.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
