/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.iv;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class IVProbeReplace {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            ProbeTranslation pb = new ProbeTranslation();
            HashMap<String, String> probeMap = pb.getProbeTranslation("/Volumes/iSnackHD/SkyDrive/MetaAnalysisAnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", 
                    "HumanHT-12_V3_0_R2_11283641_A.txt", 
                    "HumanHT-12_V4_0_R1_15002873_B.txt");
            
            TextFile tfIn = new TextFile("/Volumes/iSnackHD/SkyDrive/latesteQTLs/trans_ma_results_31072012-SNPCisTransPairs.txt", TextFile.R);
            TextFile tfOut = new TextFile("/Volumes/iSnackHD/SkyDrive/latesteQTLs/trans_ma_results_31072012-SNPCisTransPairs-HT12v4.txt", TextFile.W);
            String[] elems = tfIn.readLineElems(TextFile.tab);
            while(elems!=null){
                tfOut.writeln(elems[0]+"\t"+probeMap.get(elems[1])+"\t"+probeMap.get(elems[2]));
                elems = tfIn.readLineElems(TextFile.tab);
            }
            tfIn.close();
            tfOut.close();
            
        } catch (IOException e){
            
        }
    }
}
