/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class TransSNPIsolate {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            HashSet<String> uniqueSNPs = new HashSet<String>();
            for (int i = 0; i < 11; i++) {
                String file = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/eQTLs.txt";
                if (i > 0) {
                    file = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/PermutedEQTLsPermutationRound"+i+".txt.gz";
                }
                System.out.println(file);
                TextFile in = new TextFile(file, TextFile.R);
                String[] elems = in.readLineElems(TextFile.tab);
                elems = in.readLineElems(TextFile.tab);
                while (elems != null) {
                    uniqueSNPs.add(elems[1]);
                    elems = in.readLineElems(TextFile.tab);
                }
                in.close();
            }
            System.out.println(uniqueSNPs.size());
            TextFile out = new TextFile("/Volumes/iSnackHD/SkyDrive/AllTransSNPs.txt", TextFile.W);
            for(String s: uniqueSNPs){
                out.writeln(s);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
