/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class CombineSeq {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String out = "/Volumes/iSnackHD/Data/ProbeAnnotation/IlluminaProbeAnnotation/AllProbesFasta/AllIlluminaProbes.fa";

            TextFile tf = new TextFile("/Volumes/iSnackHD/SkyDrive2/SkyDrive/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt", TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            HashMap<String, String> seqToNewId = new HashMap<String, String>();
            TextFile fastaOut = new TextFile(out, TextFile.W);
            while (elems != null) {
                if (elems.length > 1) {
                    seqToNewId.put(elems[1].toUpperCase(), elems[0]);
                    fastaOut.writeln(">" + elems[0] + "\n" + elems[1].toUpperCase());
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            String dir = "/Volumes/iSnackHD/Data/ProbeAnnotation/IlluminaProbeAnnotation/HJ/";
            String[] filesInDir = Gpio.getListOfFiles(dir);
            HashSet<String> sequences = new HashSet<String>();
            for (String s : filesInDir) {
                TextFile tf2 = new TextFile(dir + s, TextFile.R);
                String[] elems2 = tf2.readLineElems(TextFile.tab);
                while (elems2 != null) {
                    if (elems2.length > 1) {
                        if (!seqToNewId.containsKey(elems2[1].toUpperCase())) {
                            sequences.add(elems2[1].toUpperCase());
                        }
                    }
                    elems2 = tf2.readLineElems(TextFile.tab);
                }
                tf2.close();
            }
            System.out.println(sequences.size() + "new sequences");
            int ctr = 0;
            for (String s : sequences) {
                fastaOut.writeln(">" + (seqToNewId.size() + ctr) + "\n" + s);
                ctr++;
            }

            fastaOut.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
