/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;

/**
 *
 * @author harmjan
 */
public class ProbeListCreator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String probetranslationfile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";
        String fileIn = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Unsorted/eQTLs.txt.gz";
        String fileOut = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Unsorted/eQTLs-AllTestedProbes.txt";


        try {
            String ht12v3String = "HT12v3.txt";

            //       PROBE LOOKUP TABLES..
            ProbeTranslation pbt = new ProbeTranslation();
            HashMap<String, String> probeToHT12v3 = new HashMap<String, String>();
            probeToHT12v3 = pbt.getProbeTranslation(probetranslationfile, "Probe", ht12v3String);

            TextFile tf = new TextFile(fileIn, TextFile.R);
            TextFile tfOut = new TextFile(fileOut, TextFile.W);

            String[] elems = tf.readLineElems(TextFile.tab);
            tf.readLine();
            HashSet<String> visitedProbe = new HashSet<String>();
            while (elems != null) {
                String probe = elems[4];
                if (!visitedProbe.contains(probe)) {
                    tfOut.writeln(probeToHT12v3.get(probe));
                    visitedProbe.add(probe);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            tfOut.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
