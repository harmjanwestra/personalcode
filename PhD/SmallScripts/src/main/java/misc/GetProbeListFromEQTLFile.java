/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;

/**
 *
 * @author harmjan
 */
public class GetProbeListFromEQTLFile {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String eQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Unsorted/eQTLs.txt.gz";
        String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
        GetProbeListFromEQTLFile o = new GetProbeListFromEQTLFile();
        try {
            o.run(eQTLFile, probetranslationfile);
        } catch (IOException ex) {
            Logger.getLogger(GetProbeListFromEQTLFile.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void run(String eQTLFile, String probetranslationfile) throws IOException {
        TextFile tf = new TextFile(eQTLFile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        HashSet<String> probes = new HashSet<String>();
        while (elems != null) {
            String probe = elems[4];
            probes.add(probe);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        System.out.println(probes.size());

        ProbeTranslation pbt = new ProbeTranslation();

        String probeStr = "Probe";
        String ht12v3String = "HumanHT-12_V3_0_R2_11283641_A.txt";
        String ht12v4String = "HumanHT-12_V4_0_R1_15002873_B.txt";
        HashMap<String, String> probeToHT12v3 = pbt.getProbeTranslation(probetranslationfile, probeStr, ht12v3String);
        HashMap<String, String> probeToHT12v4 = pbt.getProbeTranslation(probetranslationfile, probeStr, ht12v4String);
        
        for (String probe : probes) {
            System.out.println(probe + "\t" + probeToHT12v3.get(probe) + "\t" + probeToHT12v4.get(probe));
        }
    }
}
