/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;

/**
 *
 * @author harmjan
 */
public class ConvertListOfArrayIdsToHUGO {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here

            String ht12v3String = "HumanHT-12_V3_0_R2_11283641_A.txt";
            String ht12v4String = "HumanHT-12_V4_0_R1_15002873_B.txt";

            String hugoStr = "HUGO";
            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            //       PROBE LOOKUP TABLES..
            ProbeTranslation pbt = new ProbeTranslation();
            HashMap<String, String> ht12v4ToHT12v3 = new HashMap<String, String>();
            HashMap<String, String> ht12v3ToHugo = new HashMap<String, String>();
            ht12v4ToHT12v3 = pbt.getProbeTranslation(probetranslationfile, ht12v4String, ht12v3String);
            ht12v3ToHugo = pbt.getProbeTranslation(probetranslationfile, ht12v3String, hugoStr);
            
            
            TextFile tf = new TextFile("/Volumes/iSnackHD/tmp/probelist.txt", TextFile.R);
            String[] array = tf.readAsArray();
            tf.close();

            for(String probe: array){
                System.out.println(probe+"\t"+ht12v3ToHugo.get(probe));
            }
            
        } catch (IOException ex) {
            Logger.getLogger(ConvertListOfArrayIdsToHUGO.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
