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
public class ProbeRewriteHT12v3AndH8v2 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            int h8v2col = elems.length - 1;
            int ht12v3col = 5;


            HashMap<String, String> probeToMetaMap = new HashMap<String, String>();

            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems[h8v2col].contains("Human_RefSeq")) {
                    probeToMetaMap.put(elems[h8v2col], elems[0]);
                    System.out.println(elems[h8v2col] + "\t" + elems[0]);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            tf.open();
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                if (!elems[ht12v3col].equals("-")) {
                    probeToMetaMap.put(elems[ht12v3col], elems[0]);
                    System.out.println(elems[ht12v3col] + "\t" + elems[0]);
                } else {
                    System.out.println(elems[ht12v3col]);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();


            TextFile efilein = new TextFile("/Volumes/iSnackHD/Data/Projects/Giant/eqtlresults/2012-08-15-eQTLs-Cis-conditional/SNP-Initial/eQTLsFDR0.05.txt", TextFile.R);
            TextFile efileout = new TextFile("/Volumes/iSnackHD/Data/Projects/Giant/eqtlresults/2012-08-15-eQTLs-Cis-conditional/SNP-Initial/QC/eQTLsFDR0.05.txt", TextFile.W);
            efileout.writeln(efilein.readLine());
            String[] efileelems = efilein.readLineElems(TextFile.tab);
            while(efileelems!=null){
                efileelems[4] = probeToMetaMap.get(efileelems[4]);
                efileout.writeln(Strings.concat(efileelems, Strings.tab));
                efileelems = efilein.readLineElems(TextFile.tab);
            }
            efilein.close();
            efileout.close();


        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
