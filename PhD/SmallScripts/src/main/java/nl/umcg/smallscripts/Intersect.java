/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class Intersect {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here
            String f1 = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-03-Trans40PCsCisFXNotRegressedOut/ComparisonTo40PCCisFxRegressedOut/OldVsNewFDR0.05-eQTLComparisonLog.txt";
            String f2 = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-03-Trans40PCsCisFXNotRegressedOut/ComparisonTo40PCCisFxRegressedOut/CisFXForNonOverlappingProbesRun2.txt";

            TextFile tf = new TextFile(f1, TextFile.R);
            TextFile tf2 = new TextFile(f2, TextFile.R);

            String[] elems = tf.readLineElems(TextFile.tab);
            ArrayList<String> queryprobes = new ArrayList<String>();
            while(elems!=null){
                queryprobes.add(elems[2]);
                elems = tf.readLineElems(TextFile.tab);
            }

            elems = tf2.readLineElems(TextFile.tab);
            HashSet<String> cisprobes = new HashSet<String>();
            while(elems!=null){
                cisprobes.add(elems[4]);
                elems = tf2.readLineElems(TextFile.tab);
            }

            
            for(String s: queryprobes){
                if(!cisprobes.contains(s)){
                    System.out.println(s);
                }
            }
            
            
            
            tf.close();
            tf2.close();
        } catch (IOException e) {
        }
    }
}
