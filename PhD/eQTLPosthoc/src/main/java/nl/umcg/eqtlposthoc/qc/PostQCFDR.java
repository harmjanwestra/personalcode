/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.qc;

import eqtlmappingpipeline.metaqtl3.FDR;
import java.io.IOException;
import umcg.genetica.io.Gpio;

/**
 *
 * @author harmjan
 */
public class PostQCFDR {

    public static void main(String[] args) {
//        try {
//            FDR.permutationDir = "/Volumes/iSnackHD/Data/Projects/Isis/eQTLResults/2012-10-30-eQTLs-Cis/SNP-Initial/";
//
//            int nrPerm = 100;
//            int nrSignificant = 5756;
//            double fdrcutoff = 0.05;
//            String baseDir = "/Volumes/iSnackHD/Data/Projects/Isis/eQTLResults/2012-10-30-eQTLs-Cis/SNP-Initial/QC/";
//            
//            FDR.calculateFDR(baseDir, nrPerm, nrSignificant, fdrcutoff, true);
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//        

        for (int i = 1; i < 21; i++) {
            try {
                String outputDir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2013-03-03-IterativeFDR/Iteration"+i+"/";
                String permutationDir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/";

                Gpio.createDir(outputDir);
                int nrPerm = i;
                int nrSignificant = 2000000;
                double fdrcutoff = 0.05;
                String baseDir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/QC/";

                FDR.calculateFDR(baseDir, nrPerm, nrSignificant, fdrcutoff, false, outputDir, permutationDir);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

    }
}
