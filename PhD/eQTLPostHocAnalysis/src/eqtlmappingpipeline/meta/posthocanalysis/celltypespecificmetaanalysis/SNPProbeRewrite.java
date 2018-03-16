/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class SNPProbeRewrite {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            String ht12v3String = "HumanHT-12_V3_0_R2_11283641_A.txt";
            String ht12v4String = "HumanHT-12_V4_0_R1_15002873_B.txt";
            ProbeTranslation pbt = new ProbeTranslation();
            HashMap<String, String> ht12v3ToHT12v4 = pbt.getProbeTranslation(probetranslationfile, ht12v3String, ht12v4String);

            
            String fileIn = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/2013-07-11-HT12v3SNPProbeCombosFiltered.txt";
            String fileOut = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/2013-07-11-HT12v4SNPProbeCombosFiltered.txt";

            TextFile tfIn = new TextFile(fileIn, TextFile.R);
            TextFile tfOut = new TextFile(fileOut, TextFile.W);
            String[] elems = tfIn.readLineElems(TextFile.tab);
            while (elems != null) {
                String ht12v4 = ht12v3ToHT12v4.get(elems[1]);
                if (!ht12v4.equals("-")) {
//                    if (elems[1].equals("6060674")) {
                        System.out.println(ht12v4);

                        elems[1] = ht12v4;
                        String output = Strings.concat(elems, Strings.tab);
                        System.out.println(output);
                        tfOut.writeln(output);
//                    }


                    
                    
                }
                elems = tfIn.readLineElems(TextFile.tab);
            }
            tfOut.close();
            tfIn.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
