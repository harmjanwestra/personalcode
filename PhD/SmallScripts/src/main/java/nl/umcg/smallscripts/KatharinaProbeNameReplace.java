/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class KatharinaProbeNameReplace {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            TextFile tf = new TextFile("/Volumes/iSnackHD/SkyDrive/MetaAnalysisAnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", TextFile.R);
            HashMap<String, String> probeSeqToHT12v3 = new HashMap<String, String>();
            HashMap<String, String> probeSeqToHT12v4 = new HashMap<String, String>();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems.length >= 6) {
                    String seq = elems[1];
                    String ht12v3 = elems[5];
                    String ht12v4 = elems[7];


                    probeSeqToHT12v3.put(seq, ht12v3);
                    probeSeqToHT12v4.put(seq, ht12v4);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile tf2 = new TextFile("/Users/harmjan/Downloads/eQTLsFDR0.5ForReplication.txt", TextFile.R);
            TextFile tf3 = new TextFile("/Users/harmjan/Downloads/eQTLsFDR0.5ForReplication-HT12v3Identifiers.txt", TextFile.W);

            elems = tf2.readLineElems(TextFile.tab);
            tf3.writeln(Strings.concat(elems, Strings.tab));

            elems = tf2.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems.length >= 8) {
                    String seq = elems[6];
                    String ht12v3 = probeSeqToHT12v3.get(seq);
                    elems[7] = ht12v3;
                    String ht12v4 = probeSeqToHT12v4.get(seq);
                    if (!ht12v4.equals(ht12v3)) {
                        System.out.println(ht12v3 + "\t" + ht12v4);
                    }
                    tf3.writeln(Strings.concat(elems, Strings.tab));
                }
                elems = tf2.readLineElems(TextFile.tab);
            }

            tf3.close();
            tf2.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
