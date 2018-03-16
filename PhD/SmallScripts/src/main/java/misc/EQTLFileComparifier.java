/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import umcg.genetica.containers.Chromosome;
import umcg.genetica.containers.Exon;
import umcg.genetica.containers.Gene;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Transcript;
import umcg.genetica.containers.Triple;
import umcg.genetica.ensembl.Features;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.ProbeAnnotation;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class EQTLFileComparifier {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here
            String affyFile = "/Volumes/iSnackHD/AeroFS/TransEQTLMegaAnalysis/RickJansen/2014-02-10-CisFX/eQTLsFDR0.05-ProbeLevel_rs.txt";
            String illuminaFile = "/Volumes/iSnackHD/AeroFS/TransEQTLMegaAnalysis/RickJansen/eQTLMetaAnalysisResultst-HT12v3.txt";
            //"/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt-FilteredForProbeLevelFDR.txt.gz";

//        String f1 = "/Volumes/iSnackHD/Data/Projects/RickJansen/2014-02-11-TransFX/Trans-eQTLsFDR0.05-ProbeLevel.txt_rs.txt";
//        String f2 = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt";
            String outfileLoc = "/Volumes/iSnackHD/AeroFS/TransEQTLMegaAnalysis/RickJansen/2014-02-18-Cis-FDR0.05-GeneLevelComparison.txt";

            String illuminaAnnotation = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-18-HT12v3.txt";
            String affymetrixAnnotation = "/Volumes/iSnackHD/AeroFS/TransEQTLMegaAnalysis/RickJansen/E_Annotation.txt";
            String ensemblFeatures = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl70_HG19/structures/structures_b70.txt";

            EQTLFileComparifier comparificator = new EQTLFileComparifier();

            comparificator.yeNeweCode(illuminaFile, affyFile, outfileLoc, illuminaAnnotation, affymetrixAnnotation, ensemblFeatures);
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    public void yeOldeCode(String probefile, String f1, String f2, String outfileLoc) {
        try {

            HashMap<String, String> probeToGene = new HashMap<String, String>();
            TextFile probeannotation = new TextFile(probefile, TextFile.R);
            probeannotation.readLine();
            String[] pbelems = probeannotation.readLineElems(TextFile.tab);
            while (pbelems != null) {
                String probe = pbelems[6];
                String gene = pbelems[2];
                probeToGene.put(probe, gene);
                pbelems = probeannotation.readLineElems(TextFile.tab);
            }
            probeannotation.close();

            HashSet<String> overlapSNPs = new HashSet<String>();
            HashSet<String> overlapGenes = new HashSet<String>();

            HashSet<Pair<String, String>> overlapEQTLs = new HashSet<Pair<String, String>>();

            HashSet<String> snpsF1 = new HashSet<String>();
            HashSet<Pair<String, String>> eqtlsF1 = new HashSet<Pair<String, String>>();
            HashSet<String> genesF1 = new HashSet<String>();

            HashSet<String> snpsF2 = new HashSet<String>();
            HashSet<String> genesF2 = new HashSet<String>();
            HashSet<Pair<String, String>> eqtlsF2 = new HashSet<Pair<String, String>>();

            TextFile tf = new TextFile(f1, TextFile.R);
            String[] header = tf.readLineElems(TextFile.tab);
            int geneCol1 = -1;
            int alleleCol1 = -1;
            int allAlleleCol1 = -1;
            int zscoreCol1 = -1;

            for (int i = 0; i < header.length; i++) {
                if (header[i].equals("HGNCName")) {
                    geneCol1 = i;
                }
                if (header[i].equals("AlleleAssessed")) {
                    alleleCol1 = i;
                }
                if (header[i].equals("SNPType")) {
                    allAlleleCol1 = i;
                }
                if (header[i].equals("OverallZScore")) {
                    zscoreCol1 = i;
                }
            }
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                snpsF1.add(elems[1]);
                genesF1.addAll(Arrays.asList(elems[geneCol1].split(";")));
                eqtlsF1.add(new Pair<String, String>(elems[1], elems[geneCol1]));
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            System.out.println(snpsF1.size() + " snps in f1");
            System.out.println(genesF1.size() + " genes in f1");
            System.out.println(eqtlsF1.size() + " eqtls in f1");

            TextFile tf2 = new TextFile(f2, TextFile.R);
            String[] header2 = tf2.readLineElems(TextFile.tab);
            int geneCol2 = -1;
            int alleleCol2 = -1;
            int allAlleleCol2 = -1;
            int zscoreCol2 = -1;
            for (int i = 0; i < header2.length; i++) {
                if (header2[i].equals("HGNCName")) {
                    geneCol2 = i;
                }

                if (header2[i].equals("AlleleAssessed")) {
                    alleleCol2 = i;
                }
                if (header2[i].equals("SNPType")) {
                    allAlleleCol2 = i;
                }
                if (header2[i].equals("OverallZScore")) {
                    zscoreCol2 = i;
                }
            }

            String[] elems2 = tf2.readLineElems(TextFile.tab);
            while (elems2 != null) {

                String probe = elems2[4];
                String genes = probeToGene.get(probe);
                String[] allGenes = genes.split(";");
                genesF2.addAll(Arrays.asList(allGenes));
                snpsF2.add(elems2[1]);

                for (String gene : allGenes) {
                    if (gene.trim().length() > 1) {
                        Pair<String, String> eqtl = new Pair<String, String>(elems2[1], gene);
                        if (eqtlsF1.contains(eqtl)) {
                            overlapEQTLs.add(eqtl);
                        }
                        if (genesF1.contains(gene)) {
                            overlapGenes.add(gene);
                        }
                    }
                }

                if (snpsF1.contains(elems2[1])) {
                    overlapSNPs.add(elems2[1]);
                }

                elems2 = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();
            System.out.println(snpsF2.size() + " SNPs in f2");
            System.out.println(genesF2.size() + " genes in f2");

            System.out.println(overlapSNPs.size() + " overlap in SNPs");
            System.out.println(overlapGenes.size() + " overlap in Genes");
            System.out.println(overlapEQTLs.size() + " overlap in eQTLs");

            HashMap<Pair<String, String>, Triple<String, String, Double>> sharedEffects = new HashMap<Pair<String, String>, Triple<String, String, Double>>();
            tf.open();
            tf.readLine();
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                Pair<String, String> eqtl = new Pair<String, String>(elems[1], elems[geneCol1]);
                String assessedAllele = elems[alleleCol1];
                String allAlleles = elems[allAlleleCol1];

                String z = elems[zscoreCol1];
                Double zD = Double.parseDouble(z);

                sharedEffects.put(eqtl, new Triple(allAlleles, assessedAllele, zD));

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            tf2.open();
            tf2.readLine();
            elems2 = tf2.readLineElems(TextFile.tab);
            TextFile outfile = new TextFile(outfileLoc, TextFile.W);
            outfile.writeln("SNP\tGene1\tAlleles1\tAlleleAssessed1\tAlleles2\tAllelesAssessed2\tFlipFx\tZ1\tZ2\tZ2Flipped\tSameDirection\tDistanceBetweenProbes");
            while (elems2 != null) {
//                Pair<String, String> eqtl = new Pair<String, String>(elems[1], elems[geneCol1]);
//                String assessedAllele = elems[alleleCol1];
//                String allAlleles = elems[allAlleleCol1];
//
//                String z = elems[zscoreCol1];
//                Double zD = Double.parseDouble(z);

//                sharedEffects.put(eqtl, new Triple(allAlleles, assessedAllele, zD));
                String probe = elems2[4];
                String genes = probeToGene.get(probe);
                String[] allGenes = genes.split(";");
                for (String gene : allGenes) {
                    if (gene.trim().length() > 1) {
                        Pair<String, String> eqtl = new Pair<String, String>(elems2[1], gene);

                        Triple<String, String, Double> eqtl1 = sharedEffects.get(eqtl);

                        if (eqtl1 != null) {

                            Boolean flipAllele = BaseAnnot.flipalleles(eqtl1.getLeft(), eqtl1.getMiddle(), elems2[allAlleleCol2], elems2[alleleCol2]);
                            Double flipZ = null;
                            Double z = Double.parseDouble(elems2[zscoreCol2]);

                            if (flipAllele == null) {
                                flipZ = null;
                            } else {
                                if (flipAllele == false) {
                                    flipZ = z;
                                } else {
                                    flipZ = -z;
                                }
                            }

                            boolean sameDirection = false;
                            if (flipZ != null) {
                                if (flipZ >= 0 && eqtl1.getRight() >= 0) {
                                    sameDirection = true;
                                }
                                if (flipZ <= 0 && eqtl1.getRight() <= 0) {
                                    sameDirection = true;
                                }
                            }

                            String outln = elems2[1] + "\t"
                                    + gene + "\t"
                                    + eqtl1.getLeft() + "\t"
                                    + eqtl1.getMiddle() + "\t"
                                    + elems2[allAlleleCol2] + "\t"
                                    + elems2[alleleCol2] + "\t"
                                    + flipAllele + "\t"
                                    + eqtl1.getRight() + "\t"
                                    + elems2[zscoreCol2] + "\t"
                                    + flipZ + "\t"
                                    + sameDirection;

                            outfile.writeln(outln);

                            break;
                        }

                    }

                }

                elems2 = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();
            outfile.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void yeNeweCode(String illuminaFile, String affyFile, String outfileLoc, String illuminaAnnotationFile, String affymetrixAnnotationFile, String ensemblFeatures) throws IOException {
        Features feat = new Features();
        feat.loadAnnotation(ensemblFeatures);

        ProbeAnnotation illuminaAnnotation = new ProbeAnnotation(illuminaAnnotationFile);
        ProbeAnnotation affyAnnotation = new ProbeAnnotation(affymetrixAnnotationFile);

        HashSet<String> snps = getOverlappingSNPs(illuminaFile, affyFile);

        System.out.println(snps.size() + " overlap between SNPs");
        Pair<HashMap<String, HashSet<Pair<String, Double>>>, HashMap<String, Pair<String, String>>> illuminaData = getEQTLsPerSNP(illuminaFile, snps);
        Pair<HashMap<String, HashSet<Pair<String, Double>>>, HashMap<String, Pair<String, String>>> affyData = getEQTLsPerSNP(affyFile, snps);

        HashMap<String, HashSet<Pair<String, Double>>> illuminaEQTLs = illuminaData.getLeft();
        HashMap<String, Pair<String, String>> illuminaAlleles = illuminaData.getRight();
        HashMap<String, HashSet<Pair<String, Double>>> affyEQTLs = affyData.getLeft();
        HashMap<String, Pair<String, String>> affyAlleles = affyData.getRight();

// get probe to exon annotation :)
        System.out.println("Loading annotation for Illumina");
        Triple<HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>> annotationForIllumina = loadAnnotationForProbes(illuminaEQTLs, illuminaAnnotation, snps, feat);
        System.out.println("Loading annotation for Affy");
        Triple<HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>> annotationForAffy = loadAnnotationForProbes(affyEQTLs, affyAnnotation, snps, feat);

        // integrate results.
        System.out.println("SNP\tAffy\tIllu");

        int nrEffectsSameSNP = 0;
        int nrEffectsSameSNPSameDirection = 0;
        int nrEffectsSameGene = 0;
        int nrEffectsSameGeneSameDirection = 0;
        int nrEffectsSameTranscript = 0;
        int nrEffestsSameTranscriptSameDirection = 0;
        int nrEffectsSameExon = 0;
        int nrEffectsSameExonSameDirection = 0;

        int[] distanceBins = new int[500]; // bins of 1kb
        int[] distanceBinsSameDirection = new int[500];

        int maxDistance = distanceBins.length * 1000;

        TextFile outf = new TextFile(outfileLoc, TextFile.W);
        String header = (""
                + "SNP\t"
                + "IlluminaAlleles-Assessed\t"
                + "AffyAlleles-Assessed\t"
                + "IlluminaProbe\t"
                + "IlluminaProbeHasAnnotation\t"
                + "IlluminaProbeGenes\t"
                + "IlluminaProbeTranscripts\t"
                + "IlluminaProbeExons\t"
                + "IlluminaZScore\t"
                + "IlluminaProbePosition\t"
                + "AffyProbe\t"
                + "AffyProbeHasAnnotation\t"
                + "AffyProbeGenes\t"
                + "AffyProbeTranscripts\t"
                + "AffyProbeExons\t"
                + "AffyZScore\t"
                + "AffyProbePosition\t"
                + "DifferenceBetweenProbePositions\t"
                + "GenesShared\t"
                + "TranscriptShared\t"
                + "ExonShared\t"
                + "ZFlipped\t"
                + "SameDirection");

        outf.writeln(header);
        TextFile tfGenesOut = new TextFile(outfileLoc + "-Genes.txt", TextFile.W);
        TextFile tfExonsOut = new TextFile(outfileLoc + "-Exons.txt", TextFile.W);
        TextFile tfTranscriptsOut = new TextFile(outfileLoc + "-Transcripts.txt", TextFile.W);
        tfGenesOut.writeln(header);
        tfExonsOut.writeln(header);
        tfTranscriptsOut.writeln(header);

        for (String s : snps) {

//            System.out.println(s + "\t" + affyEQTLs.get(s).size() + "\t" + illuminaEQTLs.get(s).size());
            HashSet< Pair<String, Double>> illuEQTLForSNP = illuminaEQTLs.get(s);
            HashSet< Pair<String, Double>> affyEQTLForSNP = affyEQTLs.get(s);

            Pair<String, String> illuAllele = illuminaAlleles.get(s);
            Pair<String, String> affyAllele = affyAlleles.get(s);

            Boolean flipAllele = BaseAnnot.flipalleles(illuAllele.getLeft(), illuAllele.getRight(), affyAllele.getLeft(), affyAllele.getRight());

            String snpStr = s + "\t" + illuAllele.getLeft() + "-" + illuAllele.getRight() + "\t" + affyAllele.getLeft() + "-" + affyAllele.getRight();

            for (Pair<String, Double> p1 : illuEQTLForSNP) {
                String illuProbe = p1.getLeft();
                Integer illuProbeID = illuminaAnnotation.getProbeToProbeId().get(illuProbe);

                if (illuProbeID == null) {
                    System.err.println("Illu probe not present in Annotation shizzle: " + illuProbe);
                } else {

                    int posIlluStart = illuminaAnnotation.getChrStart()[illuProbeID];
                    int posIlluStop = illuminaAnnotation.getChrEnd()[illuProbeID];

                    double illuZ = p1.getRight();

                    int midpointDistanceIllu = posIlluStart + (posIlluStart + posIlluStop) / 2;

                    HashSet<String> genesForIlluProbe = annotationForIllumina.getLeft().get(illuProbe);
                    HashSet<String> transcriptsForIlluProbe = annotationForIllumina.getMiddle().get(illuProbe);
                    HashSet<String> exonsForIlluProbe = annotationForIllumina.getRight().get(illuProbe);

                    boolean illuProbeHasAnnotation = true;
                    if (genesForIlluProbe.isEmpty() && transcriptsForIlluProbe.isEmpty() && exonsForIlluProbe.isEmpty()) {
                        illuProbeHasAnnotation = false;
                    }

                    String illuStr = snpStr + "\t"
                            + p1.getLeft() + "\t"
                            + illuProbeHasAnnotation + "\t"
                            + convertHashToString(genesForIlluProbe) + "\t"
                            + convertHashToString(transcriptsForIlluProbe) + "\t"
                            + convertHashToString(exonsForIlluProbe) + "\t"
                            + p1.getRight() + "\t"
                            + midpointDistanceIllu;

                    for (Pair<String, Double> p2 : affyEQTLForSNP) {

                        String affyProbe = p2.getLeft();
                        Integer affyProbeID = affyAnnotation.getProbeToProbeId().get(affyProbe);
                        if (affyProbeID == null) {
                            System.out.println("Affy probe not in filez0r: " + affyProbe);
                        } else {

                            HashSet<String> genesForAffyProbe = annotationForAffy.getLeft().get(affyProbe);
                            HashSet<String> transcriptsForAffyProbe = annotationForAffy.getMiddle().get(affyProbe);
                            HashSet<String> exonsForAffyProbe = annotationForAffy.getRight().get(affyProbe);

                            boolean affyProbeHasAnnotation = true;
                            if (genesForAffyProbe.isEmpty() && transcriptsForAffyProbe.isEmpty() && exonsForAffyProbe.isEmpty()) {
                                affyProbeHasAnnotation = false;
                            }

                            HashSet<String> sharedGenes = compareHashes(genesForIlluProbe, genesForAffyProbe);
                            HashSet<String> sharedTranscripts = compareHashes(transcriptsForIlluProbe, transcriptsForAffyProbe);
                            HashSet<String> sharedExons = compareHashes(exonsForIlluProbe, exonsForAffyProbe);

                            double z = p2.getRight();

                            int posAffyStart = affyAnnotation.getChrStart()[affyProbeID];
                            int posAffyStop = affyAnnotation.getChrEnd()[affyProbeID];

                            int midpointDistanceAffy = posAffyStart + (posAffyStart + posAffyStop) / 2;

                            Double flipZ = z;
                            if (flipAllele == null) {
                                flipZ = null;
                            } else {
                                if (flipAllele) {
                                    flipZ = z * -1;
                                }

                                int sameDirection = 0;
                                if (illuZ * flipZ >= 0) {
                                    sameDirection = 1;
                                }

                                int distanceBetweenProbes = Math.abs(midpointDistanceIllu - midpointDistanceAffy);

                                if (distanceBetweenProbes > maxDistance) {
                                    distanceBetweenProbes = maxDistance;
                                }

                                int distanceBin = (int) Math.floor(((double) distanceBetweenProbes / maxDistance) * distanceBins.length);
                                if (distanceBin >= distanceBins.length) {
                                    distanceBin = distanceBins.length - 1;
                                }

                                String outline = illuStr + "\t"
                                        + affyProbe + "\t"
                                        + affyProbeHasAnnotation + "\t"
                                        + convertHashToString(genesForAffyProbe) + "\t"
                                        + convertHashToString(transcriptsForAffyProbe) + "\t"
                                        + convertHashToString(exonsForAffyProbe) + "\t"
                                        + z + "\t"
                                        + midpointDistanceAffy + "\t"
                                        + distanceBetweenProbes + "\t"
                                        + convertHashToString(sharedGenes) + "\t"
                                        + convertHashToString(sharedTranscripts) + "\t"
                                        + convertHashToString(sharedExons) + "\t"
                                        + flipZ + "\t"
                                        + sameDirection;

                                nrEffectsSameSNP++;
                                distanceBins[distanceBin]++;
                                if (sameDirection > 0) {
                                    nrEffectsSameSNPSameDirection++;
                                    distanceBinsSameDirection[distanceBin]++;
                                }

                                if (!sharedExons.isEmpty()) {
                                    nrEffectsSameExon++;
                                    tfExonsOut.writeln(outline);

                                    if (sameDirection > 0) {
                                        nrEffectsSameExonSameDirection++;
                                    }
                                }
                                if (!sharedTranscripts.isEmpty()) {
                                    nrEffectsSameTranscript++;
                                    tfTranscriptsOut.writeln(outline);
                                    if (sameDirection > 0) {
                                        nrEffestsSameTranscriptSameDirection++;
                                    }
                                }
                                if (!sharedGenes.isEmpty()) {
                                    nrEffectsSameGene++;
                                    tfGenesOut.writeln(outline);
                                    if (sameDirection > 0) {
                                        nrEffectsSameGeneSameDirection++;
                                    }
                                }
                                outf.writeln(outline);
                            }

                        }

                    }
                }
            }

        }
        outf.close();

        TextFile summary = new TextFile(outfileLoc + "-Summary.txt", TextFile.W);

        summary.writeln("Comparison\tTotal\tNrSameDirection\tPercentageSameDirection");
        summary.writeln("nrEffectsSameSNP\t" + nrEffectsSameSNP + "\t" + nrEffectsSameSNPSameDirection + "\t" + ((double) nrEffectsSameSNPSameDirection / nrEffectsSameSNP));
        summary.writeln("nrEffectsSameGene\t" + nrEffectsSameGene + "\t" + nrEffectsSameGeneSameDirection + "\t" + ((double) nrEffectsSameGeneSameDirection / nrEffectsSameGene));
        summary.writeln("nrEffectsSameTranscript\t" + nrEffectsSameTranscript + "\t" + nrEffestsSameTranscriptSameDirection + "\t" + ((double) nrEffestsSameTranscriptSameDirection / nrEffectsSameTranscript));
        summary.writeln("nrEffectsSameExon\t" + nrEffectsSameExon + "\t" + nrEffectsSameExonSameDirection + "\t" + ((double) nrEffectsSameExonSameDirection / nrEffectsSameExon));

        summary.writeln("");
        summary.writeln("---");
        summary.writeln("Bin\tDistance\tCt\tSameDirection\tPerc");
        for (int i = 0; i < distanceBins.length; i++) {
            summary.writeln(i + "\t" + (i * 1000) + "\t" + distanceBins[i] + "\t" + distanceBinsSameDirection[i] + "\t" + ((double) distanceBinsSameDirection[i] / distanceBins[i]));
        }

        summary.close();

    }

    private Pair<HashMap<String, HashSet<Pair<String, Double>>>, HashMap<String, Pair<String, String>>> getEQTLsPerSNP(String eqtlfile, Set<String> snps) throws IOException {
        // need probe, alleles, zscore
        HashMap<String, HashSet<Pair<String, Double>>> eqtls = new HashMap<String, HashSet<Pair<String, Double>>>();
        HashMap<String, Pair<String, String>> eqtlsAlleles = new HashMap<String, Pair<String, String>>();

        TextFile tf = new TextFile(eqtlfile, TextFile.R);

        String[] header = tf.readLineElems(TextFile.tab);
        int geneCol1 = -1;
        int alleleCol1 = -1;
        int allAlleleCol1 = -1;
        int zscoreCol1 = -1;

        for (int i = 0; i < header.length; i++) {
            if (header[i].equals("HGNCName")) {
                geneCol1 = i;
            }
            if (header[i].equals("AlleleAssessed")) {
                alleleCol1 = i;
            }
            if (header[i].equals("SNPType")) {
                allAlleleCol1 = i;
            }
            if (header[i].equals("OverallZScore")) {
                zscoreCol1 = i;
            }
        }
        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {

            String snp = elems[1];
            if (snps.contains(snp)) {

                Double z = Double.parseDouble(elems[zscoreCol1]);

                String probe = elems[4];

                String assessedAllele = elems[alleleCol1];
                String allAlleles = elems[allAlleleCol1];
                eqtlsAlleles.put(snp, new Pair<String, String>(allAlleles, assessedAllele));

                HashSet<Pair<String, Double>> pairSet = eqtls.get(snp);
                if (pairSet == null) {
                    pairSet = new HashSet<Pair<String, Double>>();
                }
                pairSet.add(new Pair<String, Double>(probe, z));
                eqtls.put(snp, pairSet);
                ctr++;
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        System.out.println(ctr + " eQTLs loaded from file: " + eqtlfile);

        return new Pair<HashMap<String, HashSet<Pair<String, Double>>>, HashMap<String, Pair<String, String>>>(eqtls, eqtlsAlleles);

    }

    private HashSet<String> getOverlappingSNPs(String illuminaFile, String affyFile) throws IOException {
        TextFile tf = new TextFile(illuminaFile, TextFile.R);
        tf.readLine();
        Set<String> snps1 = tf.readAsSet(1, TextFile.tab);
        tf.close();

        TextFile tf2 = new TextFile(affyFile, TextFile.R);
        tf2.readLine();
        Set<String> snps2 = tf2.readAsSet(1, TextFile.tab);
        tf2.close();

        HashSet<String> overlap = new HashSet<String>();
        for (String s : snps1) {
            if (snps2.contains(s)) {
                overlap.add(s);
            }
        }
        return overlap;

    }

    private Triple<HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>>
            loadAnnotationForProbes(HashMap<String, HashSet<Pair<String, Double>>> alleQTLs, ProbeAnnotation annotation, HashSet<String> snps, Features feat) {

        HashMap<String, HashSet<String>> exonMap = new HashMap<String, HashSet<String>>();
        HashMap<String, HashSet<String>> transcriptMap = new HashMap<String, HashSet<String>>();
        HashMap<String, HashSet<String>> geneMap = new HashMap<String, HashSet<String>>();

        HashSet<String> probes = new HashSet<String>();
        for (String snp : snps) {
            HashSet<Pair<String, Double>> eqtls = alleQTLs.get(snp);
            for (Pair<String, Double> p : eqtls) {

                String probe = p.getLeft();
                probes.add(probe);
            }
        }
        System.out.println(probes.size() + " probes");
        int ctr = 0;
        for (String probe : probes) {
            Integer id = annotation.getProbeToProbeId().get(probe);

            ctr++;
            if (id != null) {

                short chr = annotation.getChr()[id];
                HashSet<String> exonsForProbe = exonMap.get(probe);
                HashSet<String> genesForProbe = geneMap.get(probe);
                HashSet<String> transForProbe = transcriptMap.get(probe);
                if (exonsForProbe == null) {
                    exonsForProbe = new HashSet<String>();
                }
                if (genesForProbe == null) {
                    genesForProbe = new HashSet<String>();
                }
                if (transForProbe == null) {
                    transForProbe = new HashSet<String>();
                }

                if (chr > 0) {

                    Chromosome c = feat.getChromosomeHash().get("" + chr);
                    if (c == null) {
                        System.err.println("Chromosome not found: " + chr + "\t for eQTL: " + ctr);
                    }
                    ArrayList<Gene> genes = c.getGenes();
                    int pStart = annotation.getChrStart()[id];
                    int pEnd = annotation.getChrEnd()[id];

                    for (Gene g : genes) {

                        int gStart = g.getStart();
                        int gEnd = g.getEnd();
                        boolean overlap = determineFeatureOverlap(gStart, gEnd, pStart, pEnd);

                        if (overlap) {
                            genesForProbe.add(g.getName());
                            // map transcripts
                            HashMap<String, Transcript> transcripts = g.getTranscripts();
                            Set<Map.Entry<String, Transcript>> transcriptSet = transcripts.entrySet();
                            for (Map.Entry<String, Transcript> entry : transcriptSet) {
                                Transcript t = entry.getValue();
                                gStart = t.getStart();
                                gEnd = t.getEnd();
                                overlap = determineFeatureOverlap(gStart, gEnd, pStart, pEnd);
                                if (overlap) {
                                    transForProbe.add(t.getName());
                                    // map exons
                                    Exon[] exons = t.getExonsRanked();
                                    for (Exon e : exons) {
                                        gStart = e.getStart();
                                        gEnd = e.getEnd();
                                        overlap = determineFeatureOverlap(gStart, gEnd, pStart, pEnd);
                                        if (overlap) {
                                            exonsForProbe.add(e.getName());
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // and that is how you annotate probe sequences.
                // science, biatch!
                exonMap.put(probe, exonsForProbe);
                transcriptMap.put(probe, transForProbe);
                geneMap.put(probe, genesForProbe);
            }

            if (ctr % 1000 == 0) {
                System.out.println(ctr);
            }
        }

        return new Triple<HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>>(geneMap, transcriptMap, exonMap);

    }

    private boolean determineFeatureOverlap(int fStart, int fEnd, int pStart, int pEnd) {
        boolean hit = false;
        if (pStart <= fStart && pEnd > fStart) {
            //           |########|
            //         |probe|
            hit = true;
        } else if (pStart >= fStart && pStart < fEnd && pEnd < fEnd) {
            //           |########|
            //             |probe|
            hit = true;
        } else if (pStart >= fStart && pStart < fEnd && pEnd > fEnd) {
            //           |########|
            //             |probe|
            hit = true;
        }
        return hit;
    }

    private HashSet<String> compareHashes(HashSet<String> hash1, HashSet<String> hash2) {
        HashSet<String> overlap = new HashSet<String>();
        for (String s : hash1) {
            if (hash2.contains(s)) {
                overlap.add(s);
            }
        }
        return overlap;
    }

    private String convertHashToString(HashSet<String> hash) {
        String[] data = hash.toArray(new String[0]);
        return Strings.concat(data, Strings.comma);
    }

}
