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
public class RenumberProbeSNPAnnotation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String in = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-SNPsPerProbe.txt";
        String out = "/Volumes/iSnackHD/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-SNPsPerProbe.txt";
        String annotFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";

        try {
            HashMap<String, String> probeToHT12 = new HashMap<String, String>();
            TextFile annot = new TextFile(annotFile, TextFile.R);

            String[] elems = annot.readLineElems(TextFile.tab);
            while (elems != null) {
                String metaId = elems[0];
                String ht12Id = elems[5];
                probeToHT12.put(metaId, ht12Id);
                elems = annot.readLineElems(TextFile.tab);
            }

            annot.close();

            TextFile tf = new TextFile(in, TextFile.R);
            TextFile tfOut = new TextFile(out, TextFile.W);

            String[] elems2 = tf.readLineElems(TextFile.tab); // ik heb geen inspiratie vandaag
            String output = Strings.concat(elems2, Strings.tab);
            tfOut.writeln(output);
            elems2 = tf.readLineElems(TextFile.tab);
            while (elems2 != null) {
                String metaProbe = elems2[0];
                String ht12Probe = probeToHT12.get(metaProbe);

                
                elems2[0] = ht12Probe;

                output = Strings.concat(elems2, Strings.tab);
                tfOut.writeln(output);
                elems2 = tf.readLineElems(TextFile.tab);
            }

            tf.close();
            tfOut.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
