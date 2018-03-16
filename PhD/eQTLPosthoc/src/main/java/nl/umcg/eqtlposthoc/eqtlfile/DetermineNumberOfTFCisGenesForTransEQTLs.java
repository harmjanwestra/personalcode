/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author harmjan
 */
public class DetermineNumberOfTFCisGenesForTransEQTLs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {



//            if (1==1) {
//            FisherExactTest fet = new FisherExactTest();
////            System.out.println(fet.getFisherPValue(6418 - 293, 293, 68 - 8, 8));
//            System.out.println(fet.getFisherPValue(68 - 8, 8, 428 - 24, 24));
//            if (1==1) System.exit(0);
//            }



            DetermineNumberOfTFCisGenesForTransEQTLs df = new DetermineNumberOfTFCisGenesForTransEQTLs();
            String probeannot = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt";
            String finalprobeannot = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-ProbesThatAreTranscriptionFactorGenes.txt";
            String enslist = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2013-02-26-EnsemblGenesEns70-GO0003700.txt";
            df.rewriteEnsTFList(probeannot, enslist, finalprobeannot);

            String allCisEQTLs = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PreQC/eQTLs.txt.gz";
            String cisEQTLs = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt-UpdatedEnsemblAnnotation.txt.gz";
            String transeqtls = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLSNPsFDR0.05.txt";
            String alltranseqtls = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt.gz";

            String snpsThatAreGenomeWideSignificant = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLude.txt";

            String gwasnps = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLudeAndTestedForTrans.txt";
//            String gwasnps = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/TestedSNPs.txt";
            double threshold = 1.3141052405861965E-4; // FDR 0.05
//             threshold = 0.00352058461224396; // FDR 0.5

            String tfprobeannotOut = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-TranscriptionFactorGenes.txt";
            probeannot = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-ProbeToEnsToHugo.txt";
            df.rewriteTFList(enslist, probeannot, tfprobeannotOut);
            df.determineTFGenes(cisEQTLs, tfprobeannotOut, transeqtls, snpsThatAreGenomeWideSignificant, threshold);

//            df.run(finalprobeannot, cisEQTLs, transeqtls, gwasnps, threshold, alltranseqtls, false, true, true, snpsThatAreGenomeWideSignificant, enslist);

            String funcpredSNPList = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-03-GWASSNPAnnotation/SNPWithFuncPred.txt";
//            df.run2(funcpredSNPList, transeqtls, gwasnps);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    private HashMap<String, String> probeToENS;

    public void rewriteTFList(String tfFile, String ensAnnotation, String outfile) throws IOException {
        HashMap<String, String> ensemblToGene = new HashMap<String, String>();
        TextFile tf = new TextFile(ensAnnotation, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String ensembl = elems[2];

            if (!ensembl.equals("-")) {
                String[] ensemblElems = ensembl.split(",");
                for (String ens : ensemblElems) {
                    ensemblToGene.put(ens, elems[1]);
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile out = new TextFile(outfile, TextFile.W);
        TextFile tf2 = new TextFile(tfFile, TextFile.R);
        String line = tf2.readLine();
        HashSet<String> uniqueGenes = new HashSet<String>();
        while (line != null) {
            String gene = ensemblToGene.get(line);
            if (gene != null && !uniqueGenes.contains(gene)) {
                out.writeln(gene);
                uniqueGenes.add(gene);
            }
            line = tf2.readLine();
        }

        tf2.close();
        out.close();
    }

    public void rewriteEnsTFList(String probeAnnot, String enslist, String outfile) throws IOException {
        TextFile tf = new TextFile(enslist, TextFile.R);
        HashSet<String> ensTFs = new HashSet<String>();
        ensTFs.addAll(tf.readAsArrayList());
        tf.close();



        TextFile pb = new TextFile(probeAnnot, TextFile.R);
        String[] elems = pb.readLineElems(TextFile.tab);
        HashSet<String> selectedProbes = new HashSet<String>();


        probeToENS = new HashMap<String, String>();
        while (elems != null) {
            String probe = elems[0];
            String ens = elems[5];
            String[] ensgenes = ens.split(",");
            for (String gene : ensgenes) {
                if (ensTFs.contains(gene)) {
                    selectedProbes.add(probe);
                    String tfs = probeToENS.get(probe);
                    if (tfs == null) {
                        probeToENS.put(probe, gene);
                    } else {
                        tfs += "," + gene;
                        probeToENS.put(probe, tfs);
                    }
                } else {
                    String tfs = probeToENS.get(probe);
                    if (tfs == null) {
                        probeToENS.put(probe, gene);
                    } else {
                        tfs += "," + gene;
                        probeToENS.put(probe, tfs);
                    }
                }
            }
            elems = pb.readLineElems(TextFile.tab);
        }
        pb.close();

        TextFile out = new TextFile(outfile, TextFile.W);
        String[] probes = selectedProbes.toArray(new String[0]);
        for (String probe : probes) {
            out.writeln(probe);
        }
        out.close();
    }

    public void determineTFGenes(String cisEQTL, String tfAnnotation, String transeQTLs, String genomewidesignificanttranssnps, double threshold) throws IOException {


        HashSet<String> transsnpsThatAreGenomeWideSignifcant = new HashSet<String>();

        TextFile gwasSNPs = new TextFile(genomewidesignificanttranssnps, TextFile.R);
        transsnpsThatAreGenomeWideSignifcant.addAll(gwasSNPs.readAsArrayList());
        gwasSNPs.close();

        HashSet<String> transsnps = new HashSet<String>();

        TextFile transeSNPFile = new TextFile(transeQTLs, TextFile.R);
        transeSNPFile.readLine();
        String[] transelems = transeSNPFile.readLineElems(TextFile.tab);
        while (transelems != null) {
            if (transsnpsThatAreGenomeWideSignifcant.contains(transelems[1])) {
                transsnps.add(transelems[1]);
            }
            transelems = transeSNPFile.readLineElems(TextFile.tab);
        }
        transeSNPFile.close();

        HashSet<String> tfGenes = new HashSet<String>();
        TextFile tfFile = new TextFile(tfAnnotation, TextFile.R);
        tfGenes.addAll(tfFile.readAsArrayList());
        tfFile.close();

        TextFile tf = new TextFile(cisEQTL, TextFile.R);
        tf.readLine();

        String[] elems = tf.readLineElems(TextFile.tab);
        HashSet<String> uniqueGenes = new HashSet<String>();
        HashSet<String> uniqueTfGenes = new HashSet<String>();
        
        HashMap hashSNPsAlreadyAssessed = new HashMap();
        while (elems != null) {
            double p = Double.parseDouble(elems[0]);
//            if (p < threshold && transsnpsThatAreGenomeWideSignifcant.contains(elems[1]) && !transsnps.contains(elems[1])) {
            if (p < threshold) {
                
                if (!hashSNPsAlreadyAssessed.containsKey(elems[1])){
                    hashSNPsAlreadyAssessed.put(elems[1], null);
                    String geneStr = elems[eQTLTextFile.HUGO];
                    if (!geneStr.equals("-")) {
                        String[] geneElems = geneStr.split(",");
                        for (String gene : geneElems) {
                            uniqueGenes.add(gene);

                            if (tfGenes.contains(gene)) {
                                
                                if(!uniqueTfGenes.contains(gene)){
                                    System.out.println(gene);
                                    uniqueTfGenes.add(gene);
                                }
                                
                            }
//                        System.out.println(gene);
                        }
                    }
                }
            }

            elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();

        System.out.println(uniqueGenes.size());
        System.out.println(uniqueTfGenes.size());

    }

    public void run(String transcriptionFactorProbes, String cisEQTLs, String transEQTLs, String gwasSNPFile,
            double threshold, String alltranseqtls, boolean assessSNPsWithMultipleFX, boolean useAllCisEQTLSNPsAsBackground, boolean printAllTFs, String limitSNPsToList, String ensemblAnnotation) throws IOException {


//        HashSet<String> limitSNPs = null;
//        if (limitSNPsToList != null) {
//            limitSNPs = new HashSet<String>();
//            TextFile lTf = new TextFile(limitSNPsToList, TextFile.R);
//            limitSNPs.addAll(lTf.readAsArrayList());
//            lTf.close();
//            System.out.println("Limiting to : " + limitSNPs.size());
//        }


        HashSet<String> tfsEnsemblIds = new HashSet<String>();
        TextFile tfFileEns = new TextFile(ensemblAnnotation, TextFile.R);
        tfsEnsemblIds.addAll(tfFileEns.readAsArrayList());
        tfFileEns.close();

        HashSet<String> tfs = new HashSet<String>();
        TextFile tfFile = new TextFile(transcriptionFactorProbes, TextFile.R);
        tfs.addAll(tfFile.readAsArrayList());
        tfFile.close();

        HashSet<String> gwasSNPs = new HashSet<String>();
        HashSet<String> transSNPs = new HashSet<String>();

        TextFile tf = new TextFile(gwasSNPFile, TextFile.R);
        gwasSNPs.addAll(tf.readAsArrayList());
        tf.close();

        HashMap<String, HashSet<String>> transEQTLGeneCounter = new HashMap<String, HashSet<String>>();
        TextFile tfall = new TextFile(alltranseqtls, TextFile.R);
        tfall.readLine();
        String[] allelems = tfall.readLineElems(TextFile.tab);
        HashSet<String> allSNPs = new HashSet<String>();
        while (allelems != null) {
            String snp = allelems[1];
            String probe = allelems[4];
            String ens = probeToENS.get(probe);
            if (ens != null && ens.length() > 1) {
                HashSet<String> assocEnsGenes = transEQTLGeneCounter.get(snp);
                if (assocEnsGenes == null) {
                    assocEnsGenes = new HashSet<String>();
                }
                assocEnsGenes.add(ens);
                allSNPs.add(snp);
                transEQTLGeneCounter.put(snp, assocEnsGenes);
            }
            allelems = tfall.readLineElems(TextFile.tab);
        }
        tfall.close();

        HashSet<String> snpsWithMultipleEffects = null;
        if (assessSNPsWithMultipleFX) {
            snpsWithMultipleEffects = new HashSet<String>();
            for (String snp : allSNPs) {
                HashSet<String> pbs = transEQTLGeneCounter.get(snp);
                if (pbs != null && pbs.size() > 1) {
                    snpsWithMultipleEffects.add(snp);
                }
            }
            System.out.println("SNPs with multiple effect: " + snpsWithMultipleEffects.size());
        }




        tf = new TextFile(transEQTLs, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            transSNPs.add(snp);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        int nrTransEQTLWTFCis = 0;
        int nrTransEQTLWCis = 0;
        int nrGWASWTFCis = 0;
        int nrGWASWCis = 0;

        // iterate cis-eQTL file
        tf = new TextFile(cisEQTLs, TextFile.R);
        tf.readLine();
        elems = tf.readLineElems(TextFile.tab);
        HashSet<String> tfsThatCauseTransFX = new HashSet<String>();
        HashSet<String> tfsThatDoNotCauseTransFX = new HashSet<String>();
        HashSet<String> nonTFsThatCauseTransFX = new HashSet<String>();
        HashSet<String> nonTFsThatDoNotCauseTransFX = new HashSet<String>();


        HashSet<String> uniqueSNPsWTFCis = new HashSet<String>();
        HashSet<String> uniqueSNPsWoFCis = new HashSet<String>();
        HashSet<String> uniqueGWASSNPsWFCis = new HashSet<String>();
        HashSet<String> uniqueGWASSNPsWoFCis = new HashSet<String>();
        HashSet<String> visitedSNPs = new HashSet<String>();


        HashSet<String> uniqueEnsemblGenesThatareTF = new HashSet<String>();
        HashSet<String> uniqueEnsemblGenesThatareNotTF = new HashSet<String>();

        while (elems != null) {

            double p = Double.parseDouble(elems[0]);
            String snp = elems[1];

//            if (!visitedSNPs.contains(snp)) {
            String probe = elems[4];
            String ensembl = probeToENS.get(probe);
            if (p < threshold && ensembl != null && ensembl.length() > 1) {
                boolean isTF = false;

                String[] ensemblelems = ensembl.split(",");
                for (String s : ensemblelems) {
                    if (tfsEnsemblIds.contains(s)) {
                        uniqueEnsemblGenesThatareTF.add(s);
                    } else {
                        uniqueEnsemblGenesThatareNotTF.add(s);
                    }
                }


//                    if (tfs.contains(probe)) {
//                        isTF = true;
//                        uniqueEnsemblGenesThatareTF.add(ensembl);
//                    } else {
//                        uniqueEnsemblGenesThatareNotTF.add(ensembl);
//                    }
                if (transSNPs.contains(snp) && gwasSNPs.contains(snp)) {
                    if (snpsWithMultipleEffects == null || snpsWithMultipleEffects.contains(snp)) {
                        if (isTF) {
                            nrTransEQTLWTFCis++;
                            uniqueSNPsWTFCis.add(snp);
                            tfsThatCauseTransFX.add(ensembl);
                        } else {
                            nrTransEQTLWCis++;
                            uniqueSNPsWoFCis.add(snp);
                            nonTFsThatCauseTransFX.add(ensembl);
                        }
                    }
                } else {
                    if (useAllCisEQTLSNPsAsBackground) {
                        if (!gwasSNPs.contains(snp)) {
                            if (isTF) {
                                nrGWASWTFCis++;
                                uniqueGWASSNPsWFCis.add(snp);
                                tfsThatDoNotCauseTransFX.add(ensembl);
                            } else {
                                nrGWASWCis++;
                                uniqueGWASSNPsWoFCis.add(snp);
                                nonTFsThatDoNotCauseTransFX.add(ensembl);
                            }
                        }
                    } else if (gwasSNPs.contains(snp)) {
                        if (isTF) {
                            nrGWASWTFCis++;
                            uniqueGWASSNPsWFCis.add(snp);
                            tfsThatDoNotCauseTransFX.add(ensembl);
                        } else {
                            nrGWASWCis++;
                            uniqueGWASSNPsWoFCis.add(snp);
                            nonTFsThatDoNotCauseTransFX.add(ensembl);
                        }
                    }

                }

            }
            visitedSNPs.add(snp);
//            }



            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        System.out.println("Nr TF total: " + uniqueEnsemblGenesThatareTF.size());
        System.out.println("Nr of unique NonTF genes: " + uniqueEnsemblGenesThatareNotTF.size());

        System.out.println("Probes: ");
        System.out.println("Trans w/ TF CIS: " + tfsThatCauseTransFX.size());
        System.out.println("Trans w/o TF CIS: " + nonTFsThatCauseTransFX.size());
        if (useAllCisEQTLSNPsAsBackground) {
            System.out.println("Cis-EQTLs w/ TF CIS: " + tfsThatDoNotCauseTransFX.size());
            System.out.println("Cis-EQTLs w/o TF CIS: " + nonTFsThatDoNotCauseTransFX.size());
        } else {
            System.out.println("GWAS w/ TF CIS: " + tfsThatDoNotCauseTransFX.size());
            System.out.println("GWAS w/o TF CIS: " + nonTFsThatDoNotCauseTransFX.size());
        }


        double x = ChiSquare.getX(tfsThatCauseTransFX.size(), nonTFsThatCauseTransFX.size(), tfsThatDoNotCauseTransFX.size(), nonTFsThatDoNotCauseTransFX.size());
        double p = ChiSquare.getP(1, x);

        FisherExactTest fet = new FisherExactTest();
        double fetp = fet.getFisherPValue(tfsThatCauseTransFX.size(), nonTFsThatCauseTransFX.size(), tfsThatDoNotCauseTransFX.size(), nonTFsThatDoNotCauseTransFX.size());
        System.out.println(fetp);

        System.out.println("SNPs: ");
        System.out.println("Unique SNPs w/ TF CIS: " + uniqueSNPsWTFCis.size());
        System.out.println("Unique SNPs w/o TF CIS: " + uniqueSNPsWoFCis.size());
        System.out.println("Unique GWAS SNPs w/ TF CIS " + uniqueGWASSNPsWFCis.size());
        System.out.println("Unique GWAS SNPs w/o TF CIS " + uniqueGWASSNPsWoFCis.size());

        System.out.println("X: " + x);
        System.out.println("P: " + p);

        if (printAllTFs) {
            System.out.println("---");
            for (String probe : tfsThatCauseTransFX) {
                System.out.println(probe);
            }
            System.out.println("---");
            for (String probe : nonTFsThatCauseTransFX) {
                System.out.println(probe);
            }


            System.out.println("\n\n\n\n\n\n\n\n\n\n\n");
            System.out.println("---");
            System.out.println("\n\n\n\n\n\n\n\n\n\n\n");
            for (String probe : tfsThatDoNotCauseTransFX) {
                System.out.println(probe);
            }
            for (String probe : nonTFsThatDoNotCauseTransFX) {
                System.out.println(probe);
            }
        }


    }

    private void run2(String funcpredSNPList, String transeqtls, String gwasnps) throws IOException {
        TextFile tf = new TextFile(funcpredSNPList, TextFile.R);
        HashSet<String> snpswithfunc = new HashSet<String>();
        snpswithfunc.addAll(tf.readAsArrayList());
        tf.close();

        System.out.println(snpswithfunc.size() + " func snps");

        TextFile tftrans = new TextFile(transeqtls, TextFile.R);
        HashSet<String> transeqtlsnps = new HashSet<String>();
        String[] elems = tftrans.readLineElems(TextFile.tab);
        elems = tftrans.readLineElems(TextFile.tab);
        while (elems != null) {
            transeqtlsnps.add(elems[1]);
            elems = tftrans.readLineElems(TextFile.tab);
        }
        tftrans.close();

        System.out.println(transeqtlsnps.size() + " transeqtl snps");

        HashSet<String> gwasSNPsWOTrans = new HashSet<String>();
        TextFile tfgwas = new TextFile(gwasnps, TextFile.R);
        String[] listofgwassnps = tfgwas.readAsArray();
        for (String s : listofgwassnps) {
            if (!transeqtlsnps.contains(s)) {
                gwasSNPsWOTrans.add(s);
            }
        }
        tfgwas.close();

        System.out.println(listofgwassnps.length + "\t" + gwasSNPsWOTrans.size() + " GWAS SNPs / GWAS SNPs w/o trans");



        int nrTransWOFunc = 0;
        int nrTransWFunc = 0;
        int nrGWASWOFunc = 0;
        int nrGWASWFunc = 0;

        for (String s : transeqtlsnps) {
            if (snpswithfunc.contains(s)) {
                nrTransWFunc++;
            } else {
                nrTransWOFunc++;
            }
        }

        for (String s : gwasSNPsWOTrans) {
            if (snpswithfunc.contains(s)) {
                nrGWASWFunc++;
            } else {
                nrGWASWOFunc++;
            }
        }


        System.out.println("Trans W: " + nrTransWFunc);
        System.out.println("Trans W/o: " + nrTransWOFunc);
        System.out.println("GWAS W: " + nrGWASWFunc);
        System.out.println("GWAS W/o: " + nrGWASWOFunc);

        FisherExactTest fet = new FisherExactTest();
        double p = fet.getFisherPValue(nrTransWFunc, nrTransWOFunc, nrGWASWFunc, nrGWASWOFunc);
        System.out.println("FET: " + p);
    }
}
