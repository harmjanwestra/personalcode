/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class JustALittleCheck {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic her
        try{
            String f = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-GWASSNPTransEQTLEnrichment/eQTLs/QuerySet/eQTLSNPsFDR0.05.txt";
            String f2 = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-GWASSNPTransEQTLEnrichment-Run2/NotInReference.txt";
            TextFile tf = new TextFile(f, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            elems = tf.readLineElems(TextFile.tab);
            HashSet<String> snps = new HashSet<String>();
            while(elems!=null){
                snps.add(elems[1]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            
            TextFile tf2 = new TextFile(f2, TextFile.R);
            HashSet<String> snpsNotInRef = new HashSet<String>();
            snpsNotInRef.addAll(tf2.readAsArrayList());
            tf2.close();
            
            int ctr = 0;
            for(String s: snpsNotInRef){
                if(snps.contains(s)){
                    ctr++;
                }
            }
            System.out.println(ctr);
        } catch (IOException e){
            
        }
    }
}
