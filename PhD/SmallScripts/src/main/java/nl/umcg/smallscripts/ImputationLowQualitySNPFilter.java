/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ImputationLowQualitySNPFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            
            String SNPList = "/Volumes/iSnackHD/tmp/HT12v3/snpsR2gt03.txt";
            String mappingfile = "/Volumes/iSnackHD/tmp/HT12v3/SNPMappings.txt";
            String mappingout = "/Volumes/iSnackHD/tmp/HT12v3/SNPMappingsOut.txt";
            
            
            HashSet<String> allowedSNPs = new HashSet<String>();
            TextFile tf =new TextFile(SNPList, TextFile.R);
            allowedSNPs.addAll(tf.readAsArrayList());
            tf.close();

            TextFile tfout = new TextFile(mappingout, TextFile.W);
            tf =new TextFile(mappingfile, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while(elems!=null){
                String snp = elems[2];
                if(allowedSNPs.contains(snp)){
                    tfout.writeln(Strings.concat(elems, Strings.tab));
                }
                
                elems = tf.readLineElems(TextFile.tab);
            }
            tfout.close();
            tf.close();
        } catch (IOException e){
            
            e.printStackTrace();
        }
    } 
}
