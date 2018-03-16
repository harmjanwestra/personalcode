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
public class TransEQTLOverlapWithNonSignificantGWASSNPs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            HashSet<String> insignificantSNPs = new HashSet<String>();
            TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-GWASSNPsWithNonGenomeWideSignificance.txt", TextFile.R);
            insignificantSNPs.addAll(tf.readAsArrayList());
            tf.close();
            
            TextFile tf2 = new TextFile("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLSNPsFDR0.05.txt", TextFile.R);
            
            TextFile tf3 = new TextFile("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-03-R3Q1AdditionalAnalysis/GenomeWideNonSignificantTransEQTLSNPsFDR0.05.txt", TextFile.W);
            String[] elems = tf2.readLineElems(TextFile.tab);
            tf3.writeln("P\tSNP");
            elems = tf2.readLineElems(TextFile.tab);
            int ctr =0;
            while(elems!=null){
                if(insignificantSNPs.contains(elems[1])){
                    ctr++;
                    System.out.println(elems[1]);
                    tf3.writeln(elems[0]+"\t"+elems[1]);
                }
                elems = tf2.readLineElems(TextFile.tab);
            }
            tf3.close();
            System.out.println(ctr);
            tf2.close();
        } catch (IOException e){
            
        }
    }
    
    
}
