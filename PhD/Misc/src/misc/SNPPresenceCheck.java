/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class SNPPresenceCheck {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            // check whether there are new SNPs after filtering...
            for (int perm = 0; perm < 11; perm++) {
                String indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/SNPSelectionsForFuncAnalysis/";
                String prevdir = "/Volumes/iSnackHD/SkyDrive/SNPFunctionalAnnotation/Cis-SNPs/";
                if (perm == 0) {
                    indir = indir + "RealData/";
                    prevdir = prevdir + "RealData/";
                } else {
                    indir = indir + "PermutationRound" + perm + "/";
                    prevdir = prevdir + "PermutationRound" + perm + "/";
                }
                HashSet<String> snpsInNewAnalysis = new HashSet<String>();
                TextFile tf = new TextFile(indir + "AllSNPs.txt", TextFile.R);
                snpsInNewAnalysis.addAll(tf.readAsArrayList());
                tf.close();

                HashSet<String> snpsInOldAnalysis = new HashSet<String>();
                TextFile tf2 = new TextFile(prevdir + "AllSNPs.txt", TextFile.R);
                snpsInOldAnalysis.addAll(tf2.readAsArrayList());
                tf2.close();

                int ctr = 0;
                for (String s : snpsInNewAnalysis) {
                    if (snpsInOldAnalysis.contains(s)) {
                        ctr++;
                    }
                }
                System.out.println(perm + "\t" + snpsInNewAnalysis.size() + "\t" + snpsInOldAnalysis.size() + "\t" + ctr);
            }


        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
