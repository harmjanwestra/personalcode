/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class eQTLFileProbeNameReplace {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
//            String ptf = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
//
//            //String dir = "/Volumes/iSnackHD/Data/Projects/Vinod/2012-08-01-EGCUT-ReplicationForLincRNAStudy/2012-08-02-ProbesWithSNPsFilteredOut/";
//            String dir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/2012-08-08-GroningenOnly/Sorted/";
//            // /Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/trans/2012-08-08-GroningenOnly/eQTLsFDR0.05.txt
//            String requiredAnnotation = "HumanHT-12_V3_0_R2_11283641_A";
//
//            eQTLFileProbeNameReplace.run(ptf, dir + "eQTLs.txt.gz", requiredAnnotation);
//            eQTLFileProbeNameReplace.run(ptf, dir + "eQTLsFDR0.05.txt", requiredAnnotation);


//            for (int i = 0; i < 11; i++) {
//                if (i == 0) {
//                    eQTLFileProbeNameReplace.rewriteFehrmannResultsToMetaAnalysisIds(
//                            "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz",
//                            "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt",
//                            "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt",
//                            "/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/Monocytes/Trans-40PCAsRemoved-GeneticVectorsNotRemoved/eQTLs.txt",
//                            "/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/Monocytes/Trans-40PCAsRemoved-GeneticVectorsNotRemoved/filtered/eQTLs.txt",
//                            "HumanHT-12_V4_0_R1_15002873_B.txt");
//
//                } else {
//                    eQTLFileProbeNameReplace.rewriteFehrmannResultsToMetaAnalysisIds(
//                            "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz",
//                            "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt",
//                            "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt",
//                            "/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/Monocytes/Trans-40PCAsRemoved-GeneticVectorsNotRemoved/PermutedEQTLsPermutationRound" + i + ".txt.gz",
//                            "/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/Monocytes/Trans-40PCAsRemoved-GeneticVectorsNotRemoved/filtered/PermutedEQTLsPermutationRound" + i + ".txt.gz",
//                            "HumanHT-12_V4_0_R1_15002873_B.txt");
//
//                }
//            }
//
//            eQTLFileProbeNameReplace.rewriteFehrmannResultsToMetaAnalysisIds(
//                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz",
//                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt",
//                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt",
//                    "/Volumes/iSnackHD/Dropbox/eQTLMeta/Fehrmann/journal.pgen.1002197-trans.s013.txt",
//                    "/Volumes/iSnackHD/Dropbox/eQTLMeta/Fehrmann/journal.pgen.1002197-trans.s013.txt-MetaID.txt",
//                    "HumanHT-12_V3_0_R2_11283641_A");

//            eQTLFileProbeNameReplace.rewriteProbeIds(
//                    "/Volumes/iSnackHD/Data/Projects/Vinod/eQTLResults/2012-10-16-LincRNARerun/QC/TruePositives.txt",
//                    "/Volumes/iSnackHD/Data/Projects/Vinod/eQTLResults/2012-10-16-LincRNARerun/QC/eQTLsFDR0.05.txt",
//                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt",
//                    "Probe",
//                    "HumanHT-12_V3_0_R2_11283641_A.txt");

//            for (int i = 1; i < 11; i++) {
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/iSnackHD/SkyDrive/latesteQTLs/Replication/Kora2/PermutedEQTLsPermutationRound" + i + ".txt.gz",
//                        "/Volumes/iSnackHD/SkyDrive/latesteQTLs/Replication/Kora2MetaIDs/PermutedEQTLsPermutationRound" + i + ".txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "HumanHT-12_V3_0_R2_11283641_A.txt",
//                        "Probe");
//            }

//            eQTLFileProbeNameReplace.rewriteProbeIds(
//                    "/Volumes/iSnackHD/tmp/BinaryTest/MetaQTLOutput/eQTLs.txt.gz",
//                    "/Volumes/iSnackHD/tmp/BinaryTest/MetaQTLOutput/eQTLs-MetaID.txt.gz",
//                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                    "H8v2ConvToHT12",
//                    "Probe");
//            

//            int[] arrr = new int[]{20,30,35,40,50};
//            for(int q: arrr){
//            Gpio.createDir("d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\monometa\\"+q+"PC\\");
//            for (int i = 1; i < 11; i++) {
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\mono\\"+q+"PC\\PermutedEQTLsPermutationRound" + i + ".txt.gz",
//                        "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\monometa\\"+q+"PC\\PermutedEQTLsPermutationRound" + i + ".txt.gz",
//                        "d:\\SkyDrive\\MetaAnalysisAnnotationFiles\\2012-04-23-Annot.txt",
//                        "HumanHT-12_V4_0_R1_15002873_B.txt",
//                        "Probe");
//            }
//            
//            eQTLFileProbeNameReplace.rewriteProbeIds(
//                    "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\mono\\"+q+"PC\\eQTLs.txt",
//                    "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\monometa\\"+q+"PC\\eQTLs.txt.gz",
//                    "d:\\SkyDrive\\MetaAnalysisAnnotationFiles\\2012-04-23-Annot.txt",
//                    "HumanHT-12_V4_0_R1_15002873_B.txt",
//                    "Probe");
//            
//            }

            /*
             * /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMapping/2013-02-18-DILGOM-40PCsRemoved-CisEffectsRegressedOut-GeneticVectorsNotRemoved
             * /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMapping/2013-02-19_GroningenH8v2_40PCsRemoved-4GWASPCs-CisEffectsRegressedOut-GeneticVectorsNotRemoved
             * /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMapping/2013-02-19_GroningenHT12v3_40PCsRemoved-4GWASPCs-CisEffectsRegressedOut-GeneticVectorsNotRemoved
             * /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMapping/2013-02-19_SHIP-TREND_40PCsRemoved-4GWASPCs-CisEffectsRegressedOut-GeneticVectorsNotRemoved
             * /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMapping/2013-02-19-EGCUT-40PCsRemoved-CisEffectsRegressedOut-GeneticVectorsNotRemoved
             * /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMapping/2013-02-19-INCHIANTI-40PCsRemoved-CisEffectsRegressedOut-GeneticVectorsNotRemoved
             * /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMapping/2013-02-19-RotterdamStudy-40PCsRemoved-CisEffectsRegressedOut-GeneticVectorsNotRemoved
             * /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMapping/2013-02-20-HVH-HT12v3-40PCsRemoved-CisEffectsRegressedOut-GeneticVectorsNotRemoved
             * /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMapping/2013-02-20-HVH-HT12v4-40PCsRemoved-CisEffectsRegressedOut-GeneticVectorsNotRemoved
             * 
             */
//            
//            for (int i = 0; i < 11; i++) {
//
//                String outdir = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/MetaFolder/PermutationRound" + i + "/";
//                Gpio.createDir(outdir);
//                String infile = "eQTLs.txt.gz";
//                if (i > 0) {
//                    infile = "PermutedEQTLsPermutationRound" + i + ".txt.gz";
//                }
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/2013_02_26_HVH_v3_step6out_4gwaspc_5pc_GVNR_cis_trans_residual_R3Q5/" + infile,
//                        outdir + "HVH-v3.txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "HumanHT-12_V3_0_R2_11283641_A.txt",
//                        "Probe");
//
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/2013_02_26_HVH_v4_step6out_4gwaspc_10pc_GVNR_cis_trans_residual_R3Q5/" + infile,
//                        outdir + "HVH-v4.txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "HumanHT-12_V4_0_R1_15002873_B.txt",
//                        "Probe");
//
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/2013-02-26-EGCUT-KLF14TransEQTLs/" + infile,
//                        outdir + "EGCUT.txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "HumanHT-12_V3_0_R2_11283641_A.txt",
//                        "Probe");
//
//
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/2013-02-27_SHIP-TREND_KLF14SNP_40PCsRemoved-4GWASPCs-CisEffectsRegressedOut-GeneticVectorsNotRemoved/" + infile,
//                        outdir + "SHIP.txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "HumanHT-12_V3_0_R2_11283641_A.txt",
//                        "Probe");
//
//
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/2013-02-27-DILGOM-KLF14TransEQTLs/" + infile,
//                        outdir + "DILGOM.txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "HumanHT-12_V3_0_R2_11283641_A.txt",
//                        "Probe");
//
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/2013-02-27-INCHIANTI-STUDY-KLF14-SNP/" + infile,
//                        outdir + "INCHIANTI.txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "HumanHT-12_V3_0_R2_11283641_A.txt",
//                        "Probe");
//
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/2013-02-27-ROTTERDAM-STUDY-KLF14-SNP/" + infile,
//                        outdir + "RS.txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "HumanHT-12_V4_0_R1_15002873_B.txt",
//                        "Probe");
//
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/2013-02-KLF14TransEQTLs-GRNGH8v2/" + infile,
//                        outdir + "GRNGH8v2.txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "H8v2ConvToHT12",
//                        "Probe");
//
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/2013-02-KLF14TransEQTLs-GRNGHT12v3/" + infile,
//                        outdir + "GRNGHT12v3.txt.gz",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "HumanHT-12_V3_0_R2_11283641_A.txt",
//                        "Probe");
//            }


            // /Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt.gz
//            eQTLFileProbeNameReplace.rewriteProbeIds(
//                    "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt.gz",
//                    "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05-HT12v3.txt.gz",
//                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                    "Probe",
//                    "HumanHT-12_V3_0_R2_11283641_A.txt");


//            eQTLFileProbeNameReplace.rewriteProbeIds(
//                    "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSortedFilteredForFDR0.05InMeta/eQTLsFDR.txt.gz",
//                    "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSortedFilteredForFDR0.05InMeta/eQTLsFDR.txt.gz-HT12v3.txt.gz",
//                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                    "Probe",
//                    "HumanHT-12_V3_0_R2_11283641_A.txt");



//            for (int i = 1; i < 11; i++) {
//                eQTLFileProbeNameReplace.rewriteProbeIds(
//                        "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSortedFilteredForFDR0.05InMetaMinorAlleleInconsistenciesCorrected/PermutedEQTLsPermutationRound"+i+".txt.gz",
//                        "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGSFixedKoraMeta/PermutationRounds/Round"+i+"/BSGS.txt",
//                        "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                        "Probe",
//                        "HumanHT-12_V3_0_R2_11283641_A.txt");
//            }
//            // 

            // 
//            eQTLFileProbeNameReplace.rewriteProbeIds(
//                    "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSorted/eQTLs.txt",
//                    "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSorted/eQTLs-HT12v3.txt",
//                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
//                    "Probe",
//                    "HumanHT-12_V3_0_R2_11283641_A.txt", false);
            ///Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGS/ParsedSortedFilteredForFDR0.05InMeta
            eQTLFileProbeNameReplace.rewriteProbeIds(
                    "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPCisEQTLEnrichment/eQTLs/QuerySet/eQTLsFDR0.05.txt",
                    "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPCisEQTLEnrichment/eQTLs/QuerySet/eQTLsFDR0.05.txt-Meta.txt",
                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
                    "HumanHT-12_V3_0_R2_11283641_A.txt", "Probe", false);

            // 
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void run(String probetranslation, String eqtlin, String selectedAnnotation) throws IOException {
        TextFile tf2 = new TextFile(probetranslation, TextFile.R);

        String[] header = tf2.readLineElems(TextFile.tab);


        ArrayList<String> availableAnnotations = new ArrayList<String>();
        for (int i = 5; i < header.length; i++) {
            availableAnnotations.add(header[i]);
        }
        tf2.close();


        int col = 5;
        for (String annotation : availableAnnotations) {

            if (selectedAnnotation == null || annotation.contains(selectedAnnotation)) {
                tf2.open();
                tf2.readLine(); // skip header
                HashMap<String, String> newProbeToOldProbe = new HashMap<String, String>();
                String[] probeElems = tf2.readLineElems(TextFile.tab);
                while (probeElems != null) {

                    if (probeElems.length > col) {
                        String newProbe = probeElems[0];
                        String oldProbe = probeElems[col];
                        newProbeToOldProbe.put(newProbe, oldProbe);
                    }
                    probeElems = tf2.readLineElems(TextFile.tab);
                }

                tf2.close();

                TextFile tfout = new TextFile(eqtlin + "_" + annotation + ".txt.gz", TextFile.W);
                TextFile tf = new TextFile(eqtlin, TextFile.R);
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {
                    String probe = elems[4];
                    String newprobe = newProbeToOldProbe.get(probe);
                    StringBuilder output = new StringBuilder();
                    output.append(elems[0]);
                    for (int i = 1; i < elems.length; i++) {
                        output.append("\t");
                        if (i == 4) {
                            output.append(newprobe);
                        } else {
                            output.append(elems[i]);
                        }
                    }
                    tfout.writeln(output.toString());
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();

                tfout.close();
            }
            col++;
        }
    }

    public static void rewriteFehrmannResultsToMetaAnalysisIds(String cisEQTLfile, String probeFile, String probeAnnotation, String fileToProcess, String fileOut, String annotationToConvert) throws IOException {  // and filter out unwanted SNPsssssss.
        // read the cis-eQTL file of the meta-analysis (full), as it contains all SNPs passing QC
        HashSet<String> snpsPasssingQC = new HashSet<String>();
        HashSet<String> probesAllowed = new HashSet<String>();
        TextFile efile = new TextFile(cisEQTLfile, TextFile.R);
        String[] efileElems = efile.readLineElems(TextFile.tab);
        efileElems = efile.readLineElems(TextFile.tab);
        while (efileElems != null) {
            String snp = efileElems[1];
            String probe = efileElems[4];
            probesAllowed.add(probe);
            snpsPasssingQC.add(snp);
            efileElems = efile.readLineElems(TextFile.tab);
        }
        efile.close();

        // read the list of allowed probes


        // build a reverse probe lookup table
        HashMap<String, String> probesToConvertToMetaID = new HashMap<String, String>();
        TextFile pb = new TextFile(probeAnnotation, TextFile.R);
        String[] pbElems = pb.readLineElems(TextFile.tab);
        int selectedCol = -1;
        for (int col = 0; col < pbElems.length; col++) {
            System.out.print(pbElems[col]);
            if (pbElems[col].contains(annotationToConvert)) {
                selectedCol = col;
                System.out.print("\tFound it :D");
            }
            System.out.print("\n");
        }
        if (selectedCol == -1) {
            throw new IllegalArgumentException("Selected annotation: " + annotationToConvert + " not found in probe annotation file.");
        }
        pbElems = pb.readLineElems(TextFile.tab);
        while (pbElems != null) {
            String metaprobe = pbElems[0];
            String probeNameToConvert = pbElems[selectedCol];
            if (probesAllowed.contains(metaprobe)) {
                probesToConvertToMetaID.put(probeNameToConvert, metaprobe);
            }
            pbElems = pb.readLineElems(TextFile.tab);
        }
        pb.close();

        // now iterate through the file to process...
        TextFile tf = new TextFile(fileToProcess, TextFile.R);
        TextFile tfout = new TextFile(fileOut, TextFile.W);
        String[] elems = tf.readLineElems(TextFile.tab);
        tfout.writeln(Strings.concat(elems, Strings.tab));
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];
            String metaprobe = probesToConvertToMetaID.get(probe);
            if (metaprobe != null && snpsPasssingQC.contains(snp)) {
                elems[4] = metaprobe;
                tfout.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        tfout.close();

    }

    public static void rewriteProbeIds(String in, String out, String probeTranslation, String source, String destination, boolean keepnulltranslations) throws IOException {
        HashMap<String, String> sourceToDestination = new HashMap<String, String>();
        TextFile tf = new TextFile(probeTranslation, TextFile.R);

        String[] header = tf.readLineElems(TextFile.tab);
        int probesourcecolumn = -1;
        int probedestcolumn = -1;
        for (int i = 0; i < header.length; i++) {
            if (header[i].equals(source)) {
                probesourcecolumn = i;
            }
            if (header[i].equals(destination)) {
                probedestcolumn = i;
            }
        }
        if (probedestcolumn < 0 || probesourcecolumn < 0) {
            System.out.println("Error: could not find destination annotation: " + destination + " - " + probedestcolumn + " or source: " + source + " - " + probesourcecolumn);
            System.exit(0);
        } else {
            System.out.println("Rewriting: " + in);
            System.out.println("Output: " + out);
            System.out.println("Source is in column " + probesourcecolumn);
            System.out.println("Destination is in column " + probedestcolumn);
        }



        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String src = elems[probesourcecolumn];
            String dst = elems[probedestcolumn];
            sourceToDestination.put(src, dst);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile efilein = new TextFile(in, TextFile.R);
        TextFile efileout = new TextFile(out, TextFile.W);
        String efilheader = efilein.readLine();
        efileout.writeln(efilheader);
        String[] efilelems = efilein.readLineElems(TextFile.tab);
        while (efilelems != null) {
            String src = efilelems[eQTLTextFile.PROBE];
            String dest = sourceToDestination.get(src);
            if (dest == null) {
                System.out.println("WARNING: no translation for probe: " + src);
                if (keepnulltranslations) {
                    efilelems[eQTLTextFile.PROBE] = dest;
                    String outstr = Strings.concat(efilelems, Strings.tab);
                    efileout.writeln(outstr);
                }
            } else {
                efilelems[eQTLTextFile.PROBE] = dest;
                String outstr = Strings.concat(efilelems, Strings.tab);
                efileout.writeln(outstr);
            }

            efilelems = efilein.readLineElems(TextFile.tab);
        }
        efileout.close();
        efilein.close();

    }
}
