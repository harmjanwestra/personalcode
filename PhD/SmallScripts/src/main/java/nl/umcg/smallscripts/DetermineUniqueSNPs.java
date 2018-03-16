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
public class DetermineUniqueSNPs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
//            String f = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-GWAS-SNPsForPruning.txtFilteredForTestedSNPs.txt";
//            TextFile tf = new TextFile(f, TextFile.R);
//            HashSet<String> uniqueValues = new HashSet<String>();
//            uniqueValues.addAll(tf.readAsArrayList());
//            tf.close();
//            
//            System.out.println(uniqueValues.size() + " unique values");
            
            String f = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz";
            TextFile tf = new TextFile(f, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            elems = tf.readLineElems(TextFile.tab);
            HashSet<String> snps = new HashSet<String>();
            while(elems!=null){
                snps.add(elems[1]);
                elems = tf.readLineElems(TextFile.tab);
            }
            
            tf.close();
            
            System.out.println(snps.size() + " unique values");
        } catch (IOException e){
            
        }
    }
}
