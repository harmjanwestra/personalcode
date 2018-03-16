/*
 * To change this template, choose Tools | Templates
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
public class SNPProbeCombinationTranslation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String probetranslationfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt";
            String file2 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix-Meta.txt";
            String file3 = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix-Meta-HUGO.txt";

            String src = "HumanHT-12_V4_0_R1_15002873_B.txt";
            String dst = "HUGO";

            ProbeTranslation pb = new ProbeTranslation();
            HashMap<String, String> probeToProbe = pb.getProbeTranslation(probetranslationfile, src, dst);

            TextFile tf = new TextFile(file2, TextFile.R);
            TextFile tfout = new TextFile(file3, TextFile.W);
            tfout.writeln(tf.readLine());
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String probe = elems[0];
                String otherprobe = probeToProbe.get(probe);
                if (otherprobe != null) {
                    tfout.writeln(otherprobe + "\t" + elems[1]);
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
