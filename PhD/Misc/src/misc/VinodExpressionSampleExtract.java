/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class VinodExpressionSampleExtract {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try{
        VinodExpressionSampleExtract.filter(
                "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodH8v2/BloodH8v2OriginalExpressionData/ExpressionData.txt", 
                "celiac", 
                "/Volumes/iSnackHD/tmp/ExpressionDataForVinod.txt.gz");
        } catch (IOException e){
            e.printStackTrace();
        }
    }
    
    public static void filter(String expressionFile, String sampleStr, String out) throws IOException {
        System.out.println("Reading from: "+expressionFile);
        System.out.println("Writing to: "+out);
        System.out.println("Filtering for: "+sampleStr);
        TextFile tf = new TextFile(expressionFile, TextFile.R);
        TextFile tfout = new TextFile(out, TextFile.W);
        String[] header = tf.readLineElems(TextFile.tab);
        
        boolean[] includeCol = new boolean[header.length];
        includeCol[0] = true;
        ArrayList<String> samplesToInclude = new ArrayList<String>();
        for(int i=1; i<header.length; i++){
            if(header[i].toLowerCase().contains(sampleStr.toLowerCase())){
                includeCol[i] = true;
                samplesToInclude.add(header[i]);
            } else {
                includeCol[i] = false;
            }
        }
        
        System.out.println("Detected samples: "+samplesToInclude.size());
        tfout.writeln("-\t"+Strings.concat(samplesToInclude.toArray(new String[0]), Strings.tab));
        
        String[] elems = tf.readLineElems(TextFile.tab);
        
        while(elems!=null){
            StringBuilder output = new StringBuilder();
            output.append(elems[0]);
            for(int i=1; i<elems.length; i++){
                if(includeCol[i]){
                    output.append("\t").append(elems[i]);
                }
            }
            tfout.writeln(output.toString());
            elems = tf.readLineElems(TextFile.tab);
        }
        tfout.close();
        tf.close();
        
    }
}
