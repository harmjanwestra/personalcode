/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;

/**
 *
 * @author harmjan
 */
public class CreateEnsemblAnnotationFileForLude {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            HashMap<String, String> ht12ToProbeId = new HashMap<String, String>();

            ProbeTranslation pbt = new ProbeTranslation();
            String ht12v3String = "HumanHT-12_V3_0_R2_11283641_A.txt";
            ht12ToProbeId = pbt.getProbeTranslation(probetranslationfile, ht12v3String, "Probe");

            TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-OnlyEnsemblAnnotation.txt", TextFile.R);
            HashMap<String, String> probeToEnsembl = (HashMap<String, String>) tf.readAsHashMap(0, 4);
            
            TextFile out = new TextFile("/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/2013-06-26-MetaAnalysis/HT12ToProbeToEnsemblId.txt", TextFile.W);
            out.writeln("HT12v3Id\tMetaAnalysisProbeId\tEnsemblId");
            for(String key: ht12ToProbeId.keySet()){
                out.writeln(key+"\t"+ht12ToProbeId.get(key)+"\t"+probeToEnsembl.get(ht12ToProbeId.get(key)));
            }
            
            out.close();
            
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
