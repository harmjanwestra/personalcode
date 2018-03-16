/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ProbeRewrite {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {

//            ProbeTranslation pb = new ProbeTranslation();
//            HashMap<String, String> ht12ToMetaId = pb.getProbeTranslation("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", "HumanHT-12_V4_0_R1_15002873_B.txt", "Probe");
//
//
//            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/Replication/BSGS/ConcordantTransEQTLsAutstralianDataSameZScoreDistributionAsDiscordantTransEQTLs.txt", TextFile.R);
//            TextFile tfOut = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/Replication/BSGS/ConcordantTransEQTLsAutstralianDataSameZScoreDistributionAsDiscordantTransEQTLs-MetaID.txt", TextFile.W);
//            tf.readLine();
//            String line = tf.readLine();
//            while (line != null) {
//                String[] elems = line.split("-");
//                String metaId = ht12ToMetaId.get(elems[1]);
//                tfOut.writeln(elems[0] + "\t" + metaId);
//                line = tf.readLine();
//            }
//
//            tfOut.close();
//
//            tf.close();
            
//            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/Replication/BSGS/ConcordantTransEQTLsAutstralianDataSameZScoreDistributionAsDiscordantTransEQTL-MetaID.txt", TextFile.R);
            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/Replication/BSGS/DiscordantTransEQTLsAustralianData-MetaID.txt", TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            HashMap hash = new HashMap();
            while(elems!=null){
                hash.put(elems[0] + "\t" + elems[1], null);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            
            //tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/Replication/KatharinaHeim/Kora2FDR0.5/Comp-OppositeEQTLs.txt", TextFile.R);
            tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/Replication/KatharinaHeim/Kora2/eQTLsFDR-MetaIDs.txt", TextFile.R);
            elems = tf.readLineElems(TextFile.tab);
            while(elems!=null){
            
                if (hash.containsKey(elems[1] + "\t" + elems[4])) {
                    System.out.println("Discordant" + "\t" + elems[0] + "\t" + elems[1]);
                    
                } else {
                    //System.out.println("Concordant" + "\t" + elems[0] + "\t" + elems[1]);
                    
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
