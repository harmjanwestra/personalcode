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
public class FilterFile {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
//            String f1 = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/TestedSNPs.txt";
//            String f2 = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLude.txt";
//            String fout = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLudeAndTestedForTrans.txt";
//            
//            HashSet<String> snps1 = new HashSet<String>();
//            TextFile tf1 = new TextFile(f1, TextFile.R);
//            snps1.addAll(tf1.readAsArrayList());
//            tf1.close();
//            
//            HashSet<String> snps2 = new HashSet<String>();
//            TextFile tf2 = new TextFile(f2, TextFile.R);
//            snps2.addAll(tf2.readAsArrayList());
//            tf2.close();
//            
//            TextFile tfout = new TextFile(fout, TextFile.W);
//            int ctr = 0;
//            for(String s: snps1){
//                if(snps2.contains(s)){
//                    System.out.println(s);
//                    tfout.writeln(s);
//                    ctr++;
//                }
//            }
//            tfout.close();
//            System.out.println(ctr);


            String f1 = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/Repruning/GwasSignificant.txt";
            String f2 = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLudeAndTestedForTrans.txt";
//            String fout = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLudeAndTestedForTrans.txt";

            HashSet<String> snps1 = new HashSet<String>();

            TextFile tf1 = new TextFile(f1, TextFile.R);
            String[] elems = tf1.readLineElems(Strings.whitespace);
            while (elems != null) {
                if (elems.length > 1) {
                    snps1.add(elems[1]);
                }
                elems = tf1.readLineElems(Strings.whitespace);
            }
            tf1.close();

            HashSet<String> snps2 = new HashSet<String>();
            TextFile tf2 = new TextFile(f2, TextFile.R);
            snps2.addAll(tf2.readAsArrayList());
            tf2.close();

//            TextFile tfout = new TextFile(fout, TextFile.W);
            int ctr = 0;
            for (String s : snps2) {
                if (!snps1.contains(s)) {
                    System.out.println(s);
//                    tfout.writeln(s);
                    ctr++;
                }
            }
//            tfout.close();
            System.out.println(ctr);

        } catch (IOException e) {
        }
    }
}
