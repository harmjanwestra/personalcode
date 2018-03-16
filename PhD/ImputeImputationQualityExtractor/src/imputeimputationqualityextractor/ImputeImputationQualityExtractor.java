/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package imputeimputationqualityextractor;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ImputeImputationQualityExtractor {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {
            String path = "/Volumes/Data/GeneticalGenomicsDataImputedJingyuan/ImputationBloodH8v2Impute/";
            TextFile out = new TextFile("/Data/GeneticalGenomicsDatasets/BloodH8v2ImputeTriTyper/ImputationQualityScores.txt", TextFile.W);
            for (int chr = 1; chr < 23; chr++) {

                System.out.println("Parsing Chr"+chr);
                
                String filename = path+"info_chr"+chr+"_imputed.txt.gz";
                TextFile tf = new TextFile(filename, TextFile.R);
                
                String[] elems = tf.readLineElems(Strings.whitespace);
                while(elems!=null){
                    // format: 0chr 1snp 2pos 3freq 4qual 5est
                    String snp = elems[1];
                    String info = elems[4];
//                    System.out.println(snp+"\t"+info);
                    out.writeln(snp+"\t"+info);
                    elems = tf.readLineElems(Strings.whitespace);
                }
                
                tf.close();
            }
            
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
