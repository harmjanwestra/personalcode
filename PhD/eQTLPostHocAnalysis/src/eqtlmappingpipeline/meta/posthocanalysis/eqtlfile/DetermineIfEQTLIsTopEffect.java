/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class DetermineIfEQTLIsTopEffect {

    private static HashSet<String> snpsAllowed;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {

            snpsAllowed = new HashSet<String>();

            TextFile tfgwas = new TextFile("/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLude.txt", TextFile.R);
            snpsAllowed.addAll(tfgwas.readAsArrayList());
            tfgwas.close();

            DetermineIfEQTLIsTopEffect t = new DetermineIfEQTLIsTopEffect();
//            String ref = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
//            String meta = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR-AllProbeLevelQTLs-OnlyGWASSNPs.txt.gz";
//            String finemap = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR-AllProbeLevelQTLs.txt.gz";
//            String outfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-03-CisEQTLFineMapping/ComparisonToFullMetaFDR0.05.txt";
//            String snpmapfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/SNPMappings.txt";
//            String gwasSNPs = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-GWAS-SNPs.txt";
//            String probetranslation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
//            t.run(ref, meta, finemap, outfile, snpmapfile, gwasSNPs, probetranslation);

//            TRANS META
//            String ref = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
//            String meta = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05-WSampleSize.txt";
//            String finemap = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-TransEQTLFineMapping/MetaOutput/eQTLs.txt_sorted.txt.gz";
//            String outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Comparisons/FineMappingResults/ComparisonToFullMetaFDR0.05.txt";
//            String snpmapfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/SNPMappings.txt";
//            String gwasSNPs = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-GWAS-SNPs.txt";
//            String probetranslation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
//            t.run(ref, meta, finemap, outfile, snpmapfile, gwasSNPs, probetranslation);


            // TRANS META -- REBUTTAL 2
//            String ref = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
//            String meta = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/eQTLsFDR0.05-WSampleSize-GWAS.txt";
//            String finemap = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/DataForRebuttal/2013-02-TransEQTLFineMapping/MetaOutput/eQTLs.txt_sorted.txt.gz";
//            String outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/FineMapping/ComparisonToFullMetaFDR0.05v4.txt";
//            String snpmapfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/SNPMappings.txt";
//            String gwasSNPs = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-GWAS-SNPs.txt";
//            String probetranslation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
//            t.run(ref, meta, finemap, outfile, snpmapfile, gwasSNPs, probetranslation);
            
            String ref = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
            String meta = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPCisEQTLEnrichment/eQTLs/QuerySet/eQTLsFDR0.05.txt-Meta.txt";
            String finemap = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR-AllProbeLevelQTLs.txt.gz";
            String outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/FineMapping/CisComparisonToFullMetaFDR0.05v4.txt";
            String snpmapfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/SNPMappings.txt";
            String gwasSNPs = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-GWAS-SNPs.txt";
            String probetranslation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            t.run(ref, meta, finemap, outfile, snpmapfile, gwasSNPs, probetranslation);

//
//            for (int i = 0; i < 100; i++) {
//                meta = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-GWASSNPTransEQTLEnrichment/eQTLs/Set-" + i + "/eQTLSNPsFDR0.05.txt";
//                finemap = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMappingRandomSets/Set" + i + "/eQTLs.txt.gz";
//                outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/FineMapping/Set" + i + "-ComparisonOutput.txt";
//                t.run(ref, meta, finemap, outfile, snpmapfile, gwasSNPs, probetranslation);
//            }
//
//            for (int i = 0; i < 1; i++) {
//                meta = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-GWASSNPTransEQTLEnrichment/eQTLs/QuerySet/eQTLSNPsFDR0.05.txt";
//                finemap = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-TransEQTLFineMappingRandomSets/QuerySet/eQTLs.txt.gz";
//                outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/FineMapping/QuerySet-ComparisonOutput.txt";
//                t.run(ref, meta, finemap, outfile, snpmapfile, gwasSNPs, probetranslation);
//            }

//            for (int i = 0; i < 100; i++) {
//                meta = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPTransEQTLEnrichment/Set-" + i + "/eQTLSNPsFDR0.05.txt";
//                finemap = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-Rebuttal/2013-02-TransEQTLFineMappingRandomSets/Set" + i + "/eQTLs.txt.gz";
//                outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/FineMapping/Set" + i + "-ComparisonOutput.txt";
//                t.run(ref, meta, finemap, outfile, snpmapfile, gwasSNPs, probetranslation);
//            }
//
//            for (int i = 0; i < 1; i++) {
//                meta = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPTransEQTLEnrichment/QuerySet/eQTLSNPsFDR0.05.txt";
//                finemap = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-Rebuttal/2013-02-TransEQTLFineMappingRandomSets/QuerySet/eQTLs.txt.gz";
//                outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/FineMapping/QuerySet-ComparisonOutput.txt";
//                t.run(ref, meta, finemap, outfile, snpmapfile, gwasSNPs, probetranslation);
//            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    private TriTyperGenotypeData ds;

    public void run(String reference, String meta, String finemap, String outfileloc, String snpmappingfile, String gwasSNPsFile, String probetranslation) throws IOException {
        HashSet<String> gwasSNPs = new HashSet<String>();
        TextFile gwasF = new TextFile(gwasSNPsFile, TextFile.R);
        gwasSNPs.addAll(gwasF.readAsArrayList());
        gwasF.close();

        EQTL[] metaqtls = readEQTLFile(meta);
        EQTL[] fineqtls = readEQTLFile(finemap);

        HashMap<String, Byte> snpToChr = new HashMap<String, Byte>();
        HashMap<String, Integer> snpToChrPos = new HashMap<String, Integer>();
        TextFile tfsnp = new TextFile(snpmappingfile, TextFile.R);
        String[] selems = tfsnp.readLineElems(TextFile.tab);
        while (selems != null) {
            String snp = selems[2];
            byte chr = ChrAnnotation.parseChr(selems[0]);
            Integer pos = Integer.parseInt(selems[1]);
            snpToChr.put(snp, chr);
            snpToChrPos.put(snp, pos);
            selems = tfsnp.readLineElems(TextFile.tab);
        }
        tfsnp.close();

        if (ds == null) {
            ds = new TriTyperGenotypeData(reference);
        }

        HashMap<String, String> probeToHUGO = new HashMap<String, String>();
        TextFile tftmp = new TextFile(probetranslation, TextFile.R);
        String[] tfelems = tftmp.readLineElems(TextFile.tab);
        while (tfelems != null) {
            probeToHUGO.put(tfelems[0], tfelems[4]);
            tfelems = tftmp.readLineElems(TextFile.tab);
        }

        tftmp.close();


        SNPLoader loader = ds.createSNPLoader();
        DetermineLD ldcalc = new DetermineLD();

//        TextFile outfile = new TextFile(outfileloc, TextFile.W);
//        TextFile outfileWeaker = new TextFile(outfileloc + "-WeakerEffects.txt", TextFile.W);
//        TextFile outfileIdent = new TextFile(outfileloc + "-IdenticalToTransMeta.txt", TextFile.W);
        TextFile outfileStronger = new TextFile(outfileloc + "-StrongerThanTransMeta.txt", TextFile.W);
//        TextFile nrAround = new TextFile(outfileloc + "-PairsAroundTransSNP.txt", TextFile.W);
        int nrTopSNPs = 0;
        int nrTestedSNPs = 0;
        for (int e = 0; e < metaqtls.length; e++) {
            EQTL metaAnalysisEQTL = metaqtls[e];
            boolean traitAssociatedSNPisTopSNP = true;
            Integer snpId1 = ds.getSnpToSNPId().get(metaAnalysisEQTL.getRsName());
            SNP snpObj1 = null;
            if (snpId1 != null) {
                snpObj1 = ds.getSNPObject(snpId1);
                loader.loadGenotypes(snpObj1);
            }

            byte chr1 = snpToChr.get(metaAnalysisEQTL.getRsName());
            int chr1pos = snpToChrPos.get(metaAnalysisEQTL.getRsName());

            Integer[] sampleCount = metaAnalysisEQTL.getDatasetsSamples();
            int nrSamples = 0;
            for (int i = 0; i < sampleCount.length; i++) {
                if (sampleCount[i] != null) {
                    nrSamples += sampleCount[i];
                }

            }

            if (nrSamples == 0) {
                System.err.println("Nr samples == 0");
                for (int i = 0; i < sampleCount.length; i++) {
                    System.out.println(i + "\t" + sampleCount[i]);
                }

                System.exit(0);
            }


            double explainedVarianceMetaAnalysis = ZScores.zScoreToCorrelation(metaAnalysisEQTL.getZscore(), nrSamples);
            explainedVarianceMetaAnalysis *= explainedVarianceMetaAnalysis;
//            String out = me.getPvalue() + "\t" + explainedVariance1 + "\t" + nrSamples + "\t" + me.getRsName() + "\t" + chr1 + "\t" + chr1pos + "\t" + me.getProbe() + "\t" + probeToHUGO.get(me.getProbe());
//            String outweakers = me.getPvalue() + "\t" + explainedVariance1 + "\t" + nrSamples + "\t" + me.getRsName() + "\t" + chr1 + "\t" + chr1pos + "\t" + me.getProbe() + "\t" + probeToHUGO.get(me.getProbe());
//            String outident = me.getPvalue() + "\t" + explainedVariance1 + "\t" + nrSamples + "\t" + me.getRsName() + "\t" + chr1 + "\t" + chr1pos + "\t" + me.getProbe() + "\t" + probeToHUGO.get(me.getProbe());

            String strongeroutput = metaAnalysisEQTL.getPvalue() + "\t" + explainedVarianceMetaAnalysis + "\t" + nrSamples + "\t" + metaAnalysisEQTL.getRsName() + "\t" + chr1 + "\t" + chr1pos + "\t" + metaAnalysisEQTL.getProbe() + "\t" + probeToHUGO.get(metaAnalysisEQTL.getProbe());
            int nrTested = 0;
            String topLinkedEffect = "";
            double topLinkedExplainedVariance = 0;

            boolean otherEffectHasBeenFound = false;
            boolean topEffectIsIndependent = false;
            for (int f = 0; f < fineqtls.length; f++) {
                EQTL finemapEQTL = fineqtls[f];
                if (finemapEQTL.getProbe().equals(metaAnalysisEQTL.getProbe()) && !finemapEQTL.getRsName().equals(metaAnalysisEQTL.getRsName()) && !gwasSNPs.contains(finemapEQTL.getRsName())) {

                    Integer snpId2 = ds.getSnpToSNPId().get(finemapEQTL.getRsName());
                    SNP snpObj2 = null;
                    if (snpId2 != null) {
                        snpObj2 = ds.getSNPObject(snpId2);
                    }

                    if (snpToChr.get(finemapEQTL.getRsName()) == null) {
                        System.err.println("ERROR: " + finemapEQTL.getRsName() + " not in annotation..?");
                    }

                    byte chr2 = snpToChr.get(finemapEQTL.getRsName());
                    int chr2pos = snpToChrPos.get(finemapEQTL.getRsName());

                    if (chr1 == chr2) {
                        int distance = Math.abs(chr1pos - chr2pos);
                        if (distance < 250000) {
                            nrTested++;
                            Integer[] sampleCount2 = finemapEQTL.getDatasetsSamples();
                            int nrSamples2 = 0;
                            for (int i = 0; i < sampleCount2.length; i++) {
                                if (sampleCount2[i] != null) {
                                    nrSamples2 += sampleCount2[i];
                                }
                            }

                            double explainedVarianceFineMap = ZScores.zScoreToCorrelation(finemapEQTL.getZscore(), nrSamples2);
                            explainedVarianceFineMap *= explainedVarianceFineMap;

                            if (explainedVarianceFineMap > explainedVarianceMetaAnalysis && explainedVarianceFineMap > topLinkedExplainedVariance) {
                                // check whether SNP is linked to GWAS SNP
                                if (snpObj1 != null && snpObj2 != null) {
                                    loader.loadGenotypes(snpObj2);
                                    double ld = ldcalc.getRSquared(snpObj1, snpObj2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);

//                                    if (ld <= 0.2 || ld >= 0.8) {
//                                        // gwas snp is independent
////                                        topLinkedExplainedVariance = explainedVarianceMetaAnalysis;
//                                        topEffectIsIndependent = true;
                                    
//                                        topLinkedEffect = fe.getPvalue() + "\t" + explainedVarianceFineMap + "\t" + nrSamples2 + "\t" + fe.getRsName() + "\t" + chr2 + "\t" + chr2pos + "\t" + distance + "\t" + ld;
//                                    } else {
//                                        topEffectIsIndependent = false;
                                    otherEffectHasBeenFound = true;
                                    topLinkedExplainedVariance = explainedVarianceFineMap;
                                    topLinkedEffect = finemapEQTL.getPvalue() + "\t" + explainedVarianceFineMap + "\t" + nrSamples2 + "\t" + finemapEQTL.getRsName() + "\t" + chr2 + "\t" + chr2pos + "\t" + distance + "\t" + ld;
//                                    }

                                    snpObj2.clearGenotypes();
                                } else {
//                                    topEffectIsIndependent = false;
                                    otherEffectHasBeenFound = true;
                                    topLinkedExplainedVariance = explainedVarianceFineMap;
                                    topLinkedEffect = finemapEQTL.getPvalue() + "\t" + explainedVarianceFineMap + "\t" + nrSamples2 + "\t" + finemapEQTL.getRsName() + "\t" + chr2 + "\t" + chr2pos + "\t" + distance + "\tnull";
                                }

                            }


//                            // bigger or equal effect
//                            if (snpObj2 != null) {
//                                // calculate the LD
//                                loader.loadGenotypes(snpObj2);
//                                double ld = ldcalc.getRSquared(snpObj1, snpObj2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
//                                if (explainedVariance2 > explainedVariance1) {
//                                    if (ld > 0.2 && ld < 0.8) {
//                                        traitAssociatedSNPisTopSNP = false;
//                                    }
//                                    out += "\t" + fe.getPvalue() + "\t" + explainedVariance2 + "\t" + nrSamples2 + "\t" + fe.getRsName() + "\t" + chr2 + "\t" + chr2pos + "\t" + distance + "\t" + ld;
//                                } else {
//                                    out += "\t" + fe.getPvalue() + "\t" + explainedVariance2 + "\t" + nrSamples2 + "\t" + fe.getRsName() + "\t" + chr2 + "\t" + chr2pos + "\t" + distance + "\t" + ld;
//                                }
//                                if (explainedVariance2 > topLinkedExplainedVariance) {
//                                    boolean independent = true;
//                                    if (ld > 0.2 && ld < 0.8) {
//                                        independent = false;
//                                        topLinkedEffect = fe.getPvalue() + "\t" + explainedVariance2 + "\t" + nrSamples2 + "\t" + fe.getRsName() + "\t" + chr2 + "\t" + chr2pos + "\t" + distance + "\t" + ld + "\t" + independent;
//                                        topLinkedExplainedVariance = explainedVariance2;
//                                    }
//                                }
//                                snpObj2.clearGenotypes();
//
//                            } else {
//                                if (explainedVariance2 > topLinkedExplainedVariance) {
//                                    topLinkedEffect = fe.getPvalue() + "\t" + explainedVariance2 + "\t" + nrSamples2 + "\t" + fe.getRsName() + "\t" + chr2 + "\t" + chr2pos + "\t" + distance + "\tnull\tnull";
//                                    topLinkedExplainedVariance = explainedVariance2;
//                                }
//                                if (explainedVariance2 > explainedVariance1) {
//                                    out += "\t" + fe.getPvalue() + "\t" + explainedVariance2 + "\t" + nrSamples2 + "\t" + fe.getRsName() + "\t" + chr2 + "\t" + chr2pos + "\t" + distance + "\tnull";
//                                } else {
//                                    out += "\t" + fe.getPvalue() + "\t" + explainedVariance2 + "\t" + nrSamples2 + "\t" + fe.getRsName() + "\t" + chr2 + "\t" + chr2pos + "\t" + distance + "\tnull";
//                                }
//                            }



                        }
                    }
                }
//                else if (fe.getProbe().equals(me.getProbe()) && fe.getRsName().equals(me.getRsName())) {
//
//                    Integer snpId2 = ds.getSnpToSNPId().get(fe.getRsName());
//                    SNP snpObj2 = null;
//                    if (snpId2 != null) {
//                        snpObj2 = ds.getSNPObject(snpId2);
//                    }
//
//                    if (snpToChr.get(fe.getRsName()) == null) {
//                        System.err.println("ERROR: " + fe.getRsName() + " not in annotation..?");
//                    }
//
//                    int nrSamples2 = 0;
//                    Integer[] sampleCount2 = fe.getDatasetsSamples();
//                    for (int i = 0; i < sampleCount2.length; i++) {
//                        if (sampleCount2[i] != null) {
//                            nrSamples2 += sampleCount2[i];
//                        }
//
//                    }
//                    double r2 = ZScores.zScoreToCorrelation(fe.getZscore(), nrSamples2);
//                    out += "\t" + fe.getPvalue() + "\t" + (r2 * r2) + "\t" + nrSamples2 + "\t" + fe.getRsName() + "\t" + chr1 + "\t" + chr1pos + "\t0\t1";
//                }




            }

//            nrAround.writeln(me.getRsName() + "\t" + me.getProbe() + "\t" + nrTested);
//
//            nrTestedSNPs++;
//            if (traitAssociatedSNPisTopSNP) {
//                nrTopSNPs++;
//            }
//            out = traitAssociatedSNPisTopSNP + "\t" + out;
//            outident = traitAssociatedSNPisTopSNP + "\t" + outident;
//            outweakers = traitAssociatedSNPisTopSNP + "\t" + outweakers;
//            outfile.writeln(out);
//            outfileIdent.writeln(outident);
//            outfileWeaker.writeln(outweakers);


            if (!otherEffectHasBeenFound) {
                strongeroutput += "\t" + metaAnalysisEQTL.getPvalue() + "\t" + explainedVarianceMetaAnalysis + "\t" + nrSamples + "\t" + metaAnalysisEQTL.getRsName() + "\t" + chr1 + "\t" + chr1pos + "\t0\t1";
                outfileStronger.writeln("TRUE\t" + strongeroutput);
            } else {
                strongeroutput += "\t" + topLinkedEffect;
                outfileStronger.writeln("?\t" + strongeroutput);
            }


//            if (!otherEffectHasBeenFound) {
//                strongeroutput += "\t" + me.getPvalue() + "\t" + explainedVarianceMetaAnalysis + "\t" + nrSamples + "\t" + me.getRsName() + "\t" + chr1 + "\t" + chr1pos + "\t0\t1";
//                outfileStronger.writeln("TRUE\t" + strongeroutput);
//            } else {
//                if (topEffectIsIndependent) {
//                    strongeroutput += "\t" + topLinkedEffect;
//                    outfileStronger.writeln("TRUE\t" + strongeroutput);
//                } else {
//                    strongeroutput += "\t" + topLinkedEffect;
//                    outfileStronger.writeln("FALSE\t" + strongeroutput);
//                }
////                if (explainedVarianceMetaAnalysis > topLinkedExplainedVariance) {
////                strongeroutput += "\t" + me.getPvalue() + "\t" + explainedVarianceMetaAnalysis + "\t" + nrSamples + "\t" + me.getRsName() + "\t" + chr1 + "\t" + chr1pos + "\t0\t1";
//
////                }
//            }

//            if (explainedVarianceMetaAnalysis > topLinkedExplainedVariance) {
////                strongeroutput += "\t" + me.getPvalue() + "\t" + explainedVarianceMetaAnalysis + "\t" + nrSamples + "\t" + me.getRsName() + "\t" + chr1 + "\t" + chr1pos + "\t0\t1";
//                strongeroutput += "\t" + topLinkedEffect;
//                outfileStronger.writeln("TRUE\t" + strongeroutput);
//            } else {
//                strongeroutput += "\t" + topLinkedEffect;
//                outfileStronger.writeln("FALSE\t" + strongeroutput);
//            }


            if (snpObj1 != null) {
                snpObj1.clearGenotypes();
            }
        }

//        nrAround.close();
        loader.close();
        System.out.println(nrTopSNPs + "\t" + nrTestedSNPs + "\t" + ((double) nrTopSNPs / nrTestedSNPs));
//        outfileIdent.close();
//        outfileWeaker.close();
//        outfile.close();
        outfileStronger.close();
    }

    private EQTL[] readEQTLFile(String file) throws IOException {
        TextFile tf = new TextFile(file, TextFile.R);
        ArrayList<EQTL> eqtls = new ArrayList<EQTL>();
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length > 13) {
                double p = Double.parseDouble(elems[0]);
                String snp = elems[1];
                String probe = elems[4];
                Double metaz = Double.parseDouble(elems[eQTLTextFile.METAZ]);
                String samples = elems[eQTLTextFile.DATASETSIZE];
                String[] sampleArr = samples.split(";");
                if (sampleArr.length == 1) {
                    sampleArr = samples.split(",");
                }
                Integer[] sizes = new Integer[sampleArr.length];
                for (int i = 0; i < sampleArr.length; i++) {
//                    System.out.println(sampleArr[i]);
                    try {
                        sizes[i] = Integer.parseInt(sampleArr[i]);
                    } catch (NumberFormatException e) {
                        sizes[i] = null;
                    }
                }
                EQTL e = new EQTL();

                e.setRsName(snp);
//                byte chr = ChrAnnotation.parseChr(elems[eQTLTextFile.SNPCHR]);
//                Integer chrpos = Integer.parseInt(elems[eQTLTextFile.SNPLOC]);
//                e.setRsChr(chr);
//                e.setRsChrPos(chrpos);
                e.setProbe(probe);
                e.setPvalue(p);
                e.setDatasetsSamples(sizes);
                e.setZscore(metaz);
                eqtls.add(e);
            } else {
                System.out.println("Less than 10: " + Strings.concat(elems, Strings.tab));
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return eqtls.toArray(new EQTL[0]);
    }
}
