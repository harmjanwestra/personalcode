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
public class SampleRewrite {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            
            String gtefile = "";
            String fileIn = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12Combined-FehrmannHT12v3-PedMap/covariates.txt";
            String fileOut = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12Combined-FehrmannHT12v3-PedMap/covariates-expsampleIds.txt";
            
            HashMap<String, String> gte = new HashMap<String, String>();
            TextFile tf = new TextFile(gtefile, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while(elems!=null){
                gte.put(elems[0], elems[1]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            
            TextFile tfout = new TextFile(fileOut, TextFile.W);
            TextFile tfIn = new TextFile(fileIn, TextFile.R);
            elems = tfIn.readLineElems(TextFile.tab);
            while(elems!=null){
                String sample = elems[0];
                String newSample = gte.get(sample);
                elems[0] = newSample;
                tfout.writeln(Strings.concat(elems, Strings.tab));
                elems = tfIn.readLineElems(TextFile.tab);
            }
            tfIn.close();
            tfout.close();
        } catch (IOException e){
            
        }
    }
}
