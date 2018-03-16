/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.replication;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ReplicationIDConvert {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            HashMap<String, String> ht12v4ids = new HashMap<String, String>();
            HashMap<String, String> hugo = new HashMap<String, String>();
            HashMap<String, String> ensembl = new HashMap<String, String>();
            HashMap<String, String> sequence = new HashMap<String, String>();
            TextFile ensemblfile = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", TextFile.R);

            HashMap<String, String> pos = new HashMap<String, String>();
            HashMap<String, String> chr = new HashMap<String, String>();
            String[] elems = ensemblfile.readLineElems(TextFile.tab);
            while (elems != null) {

                if (elems.length > 7) {
                    String probe = elems[0];
                    String seq = elems[1];
                    String ht12v4 = elems[7];
                    String hugostr = elems[4];


                    ht12v4ids.put(probe, ht12v4);
                    sequence.put(probe, seq);
                    hugo.put(probe, hugostr);

                    chr.put(probe, elems[2]);
                    pos.put(probe, elems[3]);
                }
                elems = ensemblfile.readLineElems(TextFile.tab);
            }



            ensemblfile.close();


// read ensembl
            TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-OnlyEnsemblAnnotation.txt", TextFile.R);
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                if (elems.length > 4) {
                    String probe = elems[0];
                    String ensemblstr = elems[4];

                    ensembl.put(probe, ensemblstr);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile eQTLFileIn = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/trans/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.5.txt", TextFile.R);
            TextFile eQTLFileOut = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/trans/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.5ForReplication.txt", TextFile.W);
            elems = eQTLFileIn.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[1];
                String snpchr = elems[2];
                String snpchrpos = elems[3];
                String probe = elems[4];
                String probechr = chr.get(probe);
                String probechrpos = pos.get(probe);
                String ensemblid = ensembl.get(probe);
                String hugoid = hugo.get(probe);
                String ht12v4 = ht12v4ids.get(probe);
                String seString = sequence.get(probe);

                
                String out = snp + "\t" + snpchr + "\t" + snpchrpos + "\t" + probe + "\t" + probechr + "\t" + probechrpos + "\t" + seString + "\t" + ht12v4 + "\t" + ensemblid + "\t" + hugoid;

                out = out.replaceAll("null", "-");
                
                eQTLFileOut.writeln(out);

                elems = eQTLFileIn.readLineElems(TextFile.tab);
            }
            eQTLFileOut.close();
            eQTLFileIn.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
