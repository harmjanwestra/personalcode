/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class eQTLFileDetermineValuesForTable1 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
//        String cisFile = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR0.05.txt-UpdatedEnsemblAnnotation.txt";
//        String cisFileProbeLevelFDR = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR0.05-ProbeLevel.txt-UpdatedEnsemblAnnotation.txt";
//        String transFile = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt-UpdatedEnsemblAnnotation.txt";
//        String transSNPfile = "/Volumes/iSnackHD/SkyDrive/AllTransSNPs.txt";
//        String outdir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/";
//        
        String cisFile = "d:\\SkyDrive\\latesteQTLs\\cisFDR0.05.txt.gz";
        String cisFileProbeLevelFDR = "d:\\SkyDrive\\latesteQTLs\\cisProbeLevelFDR0.05.txt.gz";
        String transFile = "d:\\SkyDrive\\latesteQTLs\\transFDR0.5.txt.gz";
        String transSNPfile = "d:\\SkyDrive\\AllTransSNPs.txt";
        String outdir = "d:\\SkyDrive\\latesteQTLs\\MetaAnalysisResults\\";


        try {
            eQTLFileDetermineValuesForTable1.populate(cisFile, cisFileProbeLevelFDR, transFile, transSNPfile, null, outdir);
//            System.out.println("");
//            System.out.println("----");
//            System.out.println("FDR 0.5");
//            System.out.println("");
//            cisFile = "d:\\SkyDrive\\latesteQTLs\\cisFDR0.5.txt";
//            cisFileProbeLevelFDR = "d:\\SkyDrive\\latesteQTLs\\cisProbeLevelFDR0.5.txt";
//            transFile = "d:\\SkyDrive\\latesteQTLs\\transFDR0.5.txt";
//            transSNPfile = "d:\\SkyDrive\\AllTransSNPs.txt";
//            outdir = "d:\\SkyDrive\\latesteQTLs\\MetaAnalysisResultsFDR0.5\\";
//
//            eQTLFileDetermineValuesForTable1.populate(cisFile, cisFileProbeLevelFDR, transFile, transSNPfile, null, outdir);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void populate(String cislocClassicFDR, String cislocProbeLevelFDR, String transloc, String listOfTransSNPs, String ensemblAnnotation, String outdir) throws IOException {

        Gpio.createDir(outdir);
//        HashSet<String> probesMappingToExon = new HashSet<String>();
//        HashSet<String> probesMappingToGene = new HashSet<String>();
//        HashMap<String, String> probeToGene = new HashMap<String, String>();
//        
//        TextFile ensemblFile = new TextFile(ensemblAnnotation, TextFile.R);
//        String[] ensemblElems = ensemblFile.readLineElems(TextFile.tab);
//        while(ensemblElems!=null){
//            
//            ensemblElems = ensemblFile.readLineElems(TextFile.tab);
//        }
//        ensemblFile.close();
//        

        HashMap<Pair<String, String>, Boolean> pairEffectDirections = new HashMap<Pair<String, String>, Boolean>();

        HashSet<String> allTransTestedSNPs = new HashSet<String>();
        TextFile tf = new TextFile(listOfTransSNPs, TextFile.R);
        allTransTestedSNPs.addAll(tf.readAsArrayList());
        tf.close();
        System.out.println("All tested trans SNPs: " + allTransTestedSNPs.size());

        HashMap<String, String> probeToGeneMap = new HashMap<String, String>();
        TextFile cisProbeLevel = new TextFile(cislocProbeLevelFDR, TextFile.R);
        String[] data = cisProbeLevel.readLineElems(TextFile.tab);
        HashSet<String> cisGenesProbeLevel = new HashSet<String>();
        HashSet<String> cisProbesProbeLevel = new HashSet<String>();
        HashSet<String> cisSNPSProbeLevel = new HashSet<String>();
        HashSet<Pair<String, String>> cisSNPProbePairsProbeLevel = new HashSet<Pair<String, String>>();
        HashSet<String> cisProbesNotMappingToGenesProbeLevel = new HashSet<String>();

        data = cisProbeLevel.readLineElems(TextFile.tab);
        double probeLevelFDRPValueThreshold = 0;
        double metaZ = 0;
        while (data != null) {

            String probe = data[4];
            double pval = Double.parseDouble(data[0]);
            if (pval > probeLevelFDRPValueThreshold) {
                probeLevelFDRPValueThreshold = pval;
                metaZ = Double.parseDouble(data[eQTLTextFile.METAZ]);
            }
            String snp = data[1];
            String gene = data[eQTLTextFile.HUGO];
            if (!gene.equals("-")) {
                String[] geneElems = gene.split(",");
                cisGenesProbeLevel.addAll(Arrays.asList(geneElems));
                probeToGeneMap.put(probe, gene);
                Pair<String, String> pair = new Pair<String, String>(snp, gene);
                double effect = Double.parseDouble(data[eQTLTextFile.METAZ]);
                boolean isPositive = true;
                if (effect < 0) {
                    isPositive = false;
                }
                pairEffectDirections.put(pair, isPositive);
            } else {
                cisProbesNotMappingToGenesProbeLevel.add(probe);
            }
            cisSNPSProbeLevel.add(snp);
            cisProbesProbeLevel.add(probe);
            Pair<String, String> pair = new Pair<String, String>(snp, probe);


            cisSNPProbePairsProbeLevel.add(pair);

            data = cisProbeLevel.readLineElems(TextFile.tab);
        }
        cisProbeLevel.close();

        System.out.println("FDR Threshold: " + probeLevelFDRPValueThreshold);
        System.out.println("MetaZ Threshold: " + metaZ);
        System.out.println("");

        TextFile cis = new TextFile(cislocClassicFDR, TextFile.R);
        data = cis.readLineElems(TextFile.tab);
        HashSet<String> cisGenes = new HashSet<String>();
        HashSet<String> cisProbes = new HashSet<String>();
        HashSet<String> cisProbesNotMappingToGenes = new HashSet<String>();
        HashSet<String> cisSNPS = new HashSet<String>();
        HashSet<Pair<String, String>> cisSNPProbePairs = new HashSet<Pair<String, String>>();

        HashSet<String> cisGenesAtFDR = new HashSet<String>();
        HashSet<String> cisProbesAtFDR = new HashSet<String>();
        HashSet<String> cisProbesNotMappingToGenesAtFDR = new HashSet<String>();
        HashSet<String> cisSNPSAtFDR = new HashSet<String>();
        HashSet<Pair<String, String>> cisSNPProbePairsAtFDR = new HashSet<Pair<String, String>>();

        data = cis.readLineElems(TextFile.tab);
        while (data != null) {
            double pval = Double.parseDouble(data[0]);
            String probe = data[4];
            String snp = data[1];
            String gene = data[eQTLTextFile.HUGO];

            if (pval <= probeLevelFDRPValueThreshold) {
                cisSNPSAtFDR.add(snp);
                cisProbesAtFDR.add(probe);
                cisSNPProbePairsAtFDR.add(new Pair<String, String>(snp, probe));
                if (!gene.equals("-")) {
                    String[] geneElems = gene.split(",");
                    cisGenesAtFDR.addAll(Arrays.asList(geneElems));
                    probeToGeneMap.put(probe, gene);
                    Pair<String, String> pair = new Pair<String, String>(snp, gene);
                    double effect = Double.parseDouble(data[eQTLTextFile.METAZ]);
                    boolean isPositive = true;
                    if (effect < 0) {
                        isPositive = false;
                    }
                    pairEffectDirections.put(pair, isPositive);
                } else {
                    cisProbesNotMappingToGenesAtFDR.add(probe);
                }
            }

            if (!gene.equals("-")) {
                String[] geneElems = gene.split(",");
                cisGenes.addAll(Arrays.asList(geneElems));
            } else {
                cisProbesNotMappingToGenes.add(probe);
            }
            cisSNPS.add(snp);
            cisProbes.add(probe);
            Pair<String, String> pair = new Pair<String, String>(snp, probe);

            cisSNPProbePairs.add(pair);
            data = cis.readLineElems(TextFile.tab);
        }
        cis.close();



        TextFile trans = new TextFile(transloc, TextFile.R);
        data = trans.readLineElems(TextFile.tab);
        HashSet<String> transGenes = new HashSet<String>();
        HashSet<String> transProbes = new HashSet<String>();
        HashSet<String> transSNPS = new HashSet<String>();
        HashSet<Pair<String, String>> transSNPProbePairs = new HashSet<Pair<String, String>>();
        HashSet<String> transProbesNotMappingToGenes = new HashSet<String>();

        data = trans.readLineElems(TextFile.tab);
        double transMetaZ = 0;
        double transPVal = 0;
        while (data != null) {

            String probe = data[4];

            transMetaZ = Double.parseDouble(data[eQTLTextFile.METAZ]);
            transPVal = Double.parseDouble(data[0]);

            String snp = data[1];
            String gene = data[eQTLTextFile.HUGO];
            if (!gene.equals("-")) {
                String[] geneElems = gene.split(",");
                transGenes.addAll(Arrays.asList(geneElems));
                probeToGeneMap.put(probe, gene);
                Pair<String, String> pair = new Pair<String, String>(snp, gene);
                double effect = Double.parseDouble(data[eQTLTextFile.METAZ]);
                boolean isPositive = true;
                if (effect < 0) {
                    isPositive = false;
                }
                pairEffectDirections.put(pair, isPositive);
            } else {
                transProbesNotMappingToGenes.add(probe);
            }
            transSNPS.add(snp);
            transProbes.add(probe);
            Pair<String, String> pair = new Pair<String, String>(snp, probe);
            transSNPProbePairs.add(pair);
            data = trans.readLineElems(TextFile.tab);
        }
        trans.close();

        System.out.println("Trans MetaZ: " + transMetaZ);
        System.out.println("Trans PVal: " + transPVal);
        System.out.println("");
        System.out.println("cis snpprobe pairs:\t" + cisSNPProbePairs.size());
        System.out.println("cis snps:\t" + cisSNPS.size());
        System.out.println("cis probes:\t" + cisProbes.size());
        System.out.println("cis probes w/o genes:\t" + cisProbesNotMappingToGenes.size());
        System.out.println("cis genes:\t" + cisGenes.size());
        System.out.println("");
        System.out.println("cis snpprobe pairs at ProbeLevel FDR:\t" + cisSNPProbePairsAtFDR.size());
        System.out.println("cis snps at ProbeLevel FDR:\t" + cisSNPSAtFDR.size());
        System.out.println("cis probes at ProbeLevel FDR:\t" + cisProbesAtFDR.size());
        System.out.println("cis probes w/o genes at ProbeLevel FDR:\t" + cisProbesNotMappingToGenesAtFDR.size());
        System.out.println("cis genes at ProbeLevel FDR:\t" + cisGenesAtFDR.size());
        System.out.println("");
        System.out.println("cis snpprobe pairs:\t" + cisSNPProbePairsProbeLevel.size());
        System.out.println("cis snps probe level FDR:\t" + cisSNPSProbeLevel.size());
        System.out.println("cis probes probe level FDR:\t" + cisProbesProbeLevel.size());
        System.out.println("cis probes probe level FDR w/o genes:\t" + cisProbesNotMappingToGenesProbeLevel.size());
        System.out.println("cis genes probe level FDR:\t" + cisGenesProbeLevel.size());
        System.out.println("");
        System.out.println("trans snpprobe pairs:\t" + transSNPProbePairs.size());
        System.out.println("trans snps:\t" + transSNPS.size());
        System.out.println("trans probes:\t" + transProbes.size());
        System.out.println("trans probes probe level FDR w/o genes:\t" + transProbesNotMappingToGenes.size());
        System.out.println("trans genes:\t" + transGenes.size());
        System.out.println("");

        // now determine the overlap for trans-eQTL SNPs

        int significantTransSNPsWithCisEffect = 0;
        HashSet<Pair<String, String>> cisTransPairs = new HashSet<Pair<String, String>>();
        TextFile transOutFile = new TextFile(outdir + "CisAndTransEffectsPerSNP.txt", TextFile.W);
        TextFile cisTransGenePairsOut = new TextFile(outdir + "cisTransGenePairs.txt", TextFile.W);
        transOutFile.writeln("SNP\tnrCisProbes\tnrCisGenes\tCisGenes\tnrTransProbes\tnrTransGenes\tTransGenes");
        for (String transSNP : transSNPS) {
            if (cisSNPSAtFDR.contains(transSNP)) {
                significantTransSNPsWithCisEffect++;
                // write trans pairs
                String transOut = transSNP;
                int nrCis = 0;
                int nrTrans = 0;

                ArrayList<String> probes = new ArrayList<String>();
                ArrayList<String> genes = new ArrayList<String>();

                for (Pair<String, String> cisPair : cisSNPProbePairsAtFDR) {
                    if (cisPair.getLeft().equals(transSNP)) {
                        String geneName = probeToGeneMap.get(cisPair.getRight());
                        if (geneName == null) {
                            geneName = "-";
                        }
                        probes.add(cisPair.getRight());
                        genes.add(geneName);

//                        cisTransPairs.add(new Pair<String, String>())
                        nrCis++;
                    }
                }
                HashSet<String> uniqueCisGenes = new HashSet<String>();
                for (String gene : genes) {
                    if (!gene.equals("-")) {
                        uniqueCisGenes.add(gene);
                    }
                }

                transOut += "\t" + probes.size() + "\t" + uniqueCisGenes.size() + "\t" + Strings.concat(probes, Strings.semicolon) + "\t" + Strings.concat(genes, Strings.semicolon);
                probes = new ArrayList<String>();
                genes = new ArrayList<String>();

                for (Pair<String, String> transPair : transSNPProbePairs) {
                    if (transPair.getLeft().equals(transSNP)) {
                        String geneName = probeToGeneMap.get(transPair.getRight());
                        if (geneName == null) {
                            geneName = "-";
                        }
                        probes.add(transPair.getRight());
                        genes.add(geneName);

                        nrTrans++;
                    }
                }
                HashSet<String> uniqueTransGenes = new HashSet<String>();
                for (String gene : genes) {
                    if (!gene.equals("-")) {
                        uniqueTransGenes.add(gene);
                    }
                }

                for (String cisGene : uniqueCisGenes) {
                    for (String transGene : uniqueTransGenes) {
                        Pair<String, String> pair = new Pair<String, String>(cisGene, transGene, "\t");
                        if (!cisTransPairs.contains(pair)) {
                            boolean cisDirection = pairEffectDirections.get(new Pair<String, String>(transSNP, cisGene));
                            boolean transDirection = pairEffectDirections.get(new Pair<String, String>(transSNP, transGene));
                            boolean opposite = false;
                            if (cisDirection != transDirection) {
                                opposite = true;
                            }

                            cisTransGenePairsOut.writeln(transSNP + "\t" + cisGene + "\t" + transGene + "\t" + cisDirection + "\t" + transDirection + "\t" + opposite);
                            cisTransPairs.add(pair);
                        }


                    }
                }

                transOut += "\t" + probes.size() + "\t" + uniqueCisGenes.size() + "\t" + Strings.concat(probes, Strings.semicolon) + "\t" + Strings.concat(genes, Strings.semicolon);
                transOutFile.writeln(transOut);
            }
        }
        transOutFile.close();


//        for (Pair<String, String> p : cisTransPairs) {
//            cisTransGenePairsOut.writeln(p + "");
//        }
        cisTransGenePairsOut.close();

        int nrTransTestedSNPsWithCisEffect = 0;

        HashSet<String> uniqueCisGenesForTransSNPs = new HashSet<String>();
        HashSet<String> uniqueNonGeneProbesForTransSNPs = new HashSet<String>();
        HashSet<String> uniqueCisProbesForTransSNPs = new HashSet<String>();
        HashSet<Pair<String, String>> uniqueCisSNPProbePairsForTransSNPs = new HashSet<Pair<String, String>>();
        for (String transSNP : allTransTestedSNPs) {
            if (cisSNPSAtFDR.contains(transSNP)) {
                nrTransTestedSNPsWithCisEffect++;

                for (Pair<String, String> cisPair : cisSNPProbePairsAtFDR) {

                    if (cisPair.getLeft().equals(transSNP)) {
                        uniqueCisSNPProbePairsForTransSNPs.add(cisPair);
                        String geneName = probeToGeneMap.get(cisPair.getRight());
                        if (geneName == null) {
                            geneName = "-";
                            uniqueNonGeneProbesForTransSNPs.add(cisPair.getRight());
                        } else {
                            uniqueCisGenesForTransSNPs.add(geneName);
                        }
                        uniqueCisProbesForTransSNPs.add(cisPair.getRight());

                    }
                }
            }
        }



        System.out.println("Nr cis SNP-Probe pairs for trans SNP: " + uniqueCisSNPProbePairsForTransSNPs.size());
        System.out.println("Nr Trans SNPs with Cis: " + nrTransTestedSNPsWithCisEffect);
        System.out.println("Nr Probes for those SNPs " + uniqueCisProbesForTransSNPs.size());
        System.out.println("Nr Probes w/o genes for those SNPs " + uniqueNonGeneProbesForTransSNPs.size());
        System.out.println("Nr genes for those SNPs " + uniqueCisGenesForTransSNPs.size());


        int nrCisTransSNPs = 0;
        int nrCisTransProbes = 0;
        for (String transSNP : allTransTestedSNPs) {
            if (cisSNPSAtFDR.contains(transSNP) && transSNPS.contains(transSNP)) {
                nrCisTransSNPs++;
            }
        }
        HashSet<String> uniqueCisTransGenes = new HashSet<String>();
        HashSet<String> uniqueCisTransProbesWithoutGenename = new HashSet<String>();
        for (String cisProbe : cisProbesAtFDR) {
            if (transProbes.contains(cisProbe)) {
                nrCisTransProbes++;
                String geneName = probeToGeneMap.get(cisProbe);
                if (geneName == null) {
                    uniqueCisTransProbesWithoutGenename.add(cisProbe);
                } else {
                    uniqueCisTransGenes.add(geneName);
                }
            }
        }
        System.out.println("");
        System.out.println("Nr cistrans snps: " + nrCisTransSNPs);
        System.out.println("Nr cistrans probes: " + nrCisTransProbes);
        System.out.println("Nr cistrans genes: " + uniqueCisTransGenes.size());
        System.out.println("Nr cistrans probes w/o genename: " + uniqueCisTransProbesWithoutGenename.size());

        int nrSNPsTransOnly = 0;
        int nrSNPsCisOnly = 0;
        for (String transSNP : allTransTestedSNPs) {
            if (!cisSNPSAtFDR.contains(transSNP) && transSNPS.contains(transSNP)) {
                nrSNPsTransOnly++;
            }
            if (cisSNPSAtFDR.contains(transSNP) && !transSNPS.contains(transSNP)) {
                nrSNPsCisOnly++;
            }
        }
        int nrProbesCisOnly = 0;
        int nrProbesTransOnly = 0;
        int nrProbesWoGenesCisOnly = 0;
        HashSet<String> uniqueGenesTransOnly = new HashSet<String>();

        int nrProbesWoGenesTransOnly = 0;
        HashSet<String> uniqueGenesCisOnly = new HashSet<String>();
        for (String transProbe : transProbes) {
            if (!uniqueCisProbesForTransSNPs.contains(transProbe)) {
                nrProbesTransOnly++;
                String geneName = probeToGeneMap.get(transProbe);
                if (geneName == null) {
                    nrProbesWoGenesTransOnly++;
                } else {
                    uniqueGenesTransOnly.add(geneName);
                }

            }
        }
        for (String transProbe : uniqueCisProbesForTransSNPs) {
            if (!transProbes.contains(transProbe)) {
                nrProbesCisOnly++;
                String geneName = probeToGeneMap.get(transProbe);
                if (geneName == null) {
                    nrProbesWoGenesCisOnly++;
                } else {
                    uniqueGenesCisOnly.add(geneName);
                }
            }
        }

        System.out.println("");
        System.out.println("Nr cisonly snps: " + nrSNPsCisOnly);
        System.out.println("Nr cisonly probes: " + nrProbesCisOnly);
        System.out.println("Nr cisonly genes: " + uniqueGenesCisOnly.size());
        System.out.println("Nr cisonly probes w/o genename: " + nrProbesWoGenesCisOnly);

        System.out.println("");
        System.out.println("Nr transonly snps: " + nrSNPsTransOnly);
        System.out.println("Nr transonly probes: " + nrProbesTransOnly);
        System.out.println("Nr transonly genes: " + uniqueGenesTransOnly.size());
        System.out.println("Nr transonly probes w/o genename: " + nrProbesWoGenesTransOnly);

    }
}
