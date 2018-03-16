/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.snps;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.containers.Chromosome;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class SelectSNPsWithinACertainWindow {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here/*

        /*
         * // use LD pruned input (!1!)
         SelectSimilarSNPs q = new SelectSimilarSNPs();
         String reference = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
         //            String reference = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
         String querySNPFile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-GWASSNPTransEQTLEnrichment/GWASCatalogClumped.clumped";
         String cisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz";
         // cisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/QC/TruePositives.txt";
         String outdir = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-GWASSNPTransEQTLEnrichment/";
         String geneAnnotation = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl54_HG18/structures/structures_b54.txt";
         int nrOfSets = 100;
         double r2 = 0.2;
         */
        SelectSNPsWithinACertainWindow s = new SelectSNPsWithinACertainWindow();
        String query = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt.gz";
        String refer = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
        // refer = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
        String outfi = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-TransFineMappingAroundGWASSNPs/SNPsAroundGWASSNPs-TestRun.txt";
        String probet = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt";
        String source = "Probe";
        String cisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz";

        String dest = "";

        int windowsi = 250000;

        for (int i = 0; i < 1; i++) {
            try {
                String dir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-GWASSNPTransEQTLEnrichment/eQTLs/QuerySet/eQTLsFDR0.05.txt";
                String outfis = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-GWASSNPTransEQTLEnrichment/querySNPs.txt";
                s.run(dir, refer, outfis, windowsi, cisEQTLFile);
            } catch (IOException ex) {
                Logger.getLogger(SelectSNPsWithinACertainWindow.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

//        try {
//            s.run(query, refer, outfi, windowsi, cisEQTLFile);
//            dest = "HumanHT-12_V3_0_R2_11283641_A.txt";
//            s.replaceProbeNames(probet, source, dest, outfi + "-combos.txt", outfi + "-combos-HT12v3.txt");
//
//            dest = "HumanHT-12_V4_0_R1_15002873_B.txt";
//            s.replaceProbeNames(probet, source, dest, outfi + "-combos.txt", outfi + "-combos-HT12v4.txt");
//
//            dest = "H8v2ConvToHT12";
//            s.replaceProbeNames(probet, source, dest, outfi + "-combos.txt", outfi + "-combos-H8v2.txt");
//        } catch (IOException ex) {
//            Logger.getLogger(SelectSNPsWithinACertainWindow.class.getName()).log(Level.SEVERE, null, ex);
//        }
    }
    private HashSet<String> eQTLSNPs;

    public void run(String queryList, String reference, String outfile, int windowsize, String cisEQTLFile) throws IOException {
        // use the cis-eQTLs to determine the distance to the gene for each SNP.
        TextFile efile = new TextFile(cisEQTLFile, TextFile.R);
        efile.readLine();
        String[] elems = efile.readLineElemsReturnObjects(TextFile.tab);

        int lnctr = 0;

        if (eQTLSNPs == null) {
            eQTLSNPs = new HashSet<String>();
            while (elems != null) {
                String snp = elems[1];

                eQTLSNPs.add(snp);

                lnctr++;
                if (lnctr % 1000000 == 0) {
                    System.out.println(lnctr + " eQTLs parsed");
                }
                elems = efile.readLineElemsReturnObjects(TextFile.tab);
            }
        }
        efile.close();

        System.out.println("Loaded SNPs: " + eQTLSNPs.size());


        TriTyperGenotypeData genotypes = new TriTyperGenotypeData(reference);
        SNPLoader loader = genotypes.createSNPLoader();

        TextFile qfile = new TextFile(queryList, TextFile.R);
        HashSet<String> query = new HashSet<String>();
        HashMap<String, HashSet<String>> probesForQSNPs = new HashMap<String, HashSet<String>>();
        elems = qfile.readLineElems(TextFile.tab);
        HashSet<String> uniqueProbes = new HashSet<String>();
        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];
            HashSet<String> probes = probesForQSNPs.get(snp);
            if (probes == null) {
                probes = new HashSet<String>();
            }
            probes.add(probe);
            uniqueProbes.add(probe);
            query.add(snp);
            probesForQSNPs.put(snp, probes);
            elems = qfile.readLineElems(TextFile.tab);
        }
        qfile.close();

        System.out.println(uniqueProbes.size() + " probes in query");
        System.out.println(query.size() + " snps in query. Will select SNPs in window of: " + windowsize);
        DetermineLD ldcalc = new DetermineLD();
        TextFile tf = new TextFile(outfile, TextFile.W);
        TextFile tfCombos = new TextFile(outfile + "-combos.txt", TextFile.W);
        String[] availableSNPs = genotypes.getSNPs();

        int combocounter = 0;

        HashSet<String> printedSNPProbeCombos = new HashSet<String>();
        for (String querySNP : query) {
            Integer querySNPid = genotypes.getSnpToSNPId().get(querySNP);

            if (querySNPid != null) {
                byte chr = genotypes.getChr(querySNPid);
                int snppos = genotypes.getChrPos(querySNPid);

                SNP snpObj1 = genotypes.getSNPObject(querySNPid);
                loader.loadGenotypes(snpObj1);

                HashSet<String> probesForsnp = probesForQSNPs.get(querySNP);
                int ctr = 0;
                int ctr2 = 0;
                for (int snp = 0; snp < availableSNPs.length; snp++) {
                    if (snp != querySNPid) {
                        byte chr2 = genotypes.getChr(snp);
                        if (chr == chr2) {
                            String snp2 = availableSNPs[snp];
                            if (eQTLSNPs.contains(snp2)) {
                                int snppos2 = genotypes.getChrPos(snp);
                                int distance = Math.abs(snppos - snppos2);
                                if (distance < windowsize) {
                                    for (String probe : probesForsnp) {
                                        String combo = availableSNPs[snp] + "\t" + probe;
                                        if (!printedSNPProbeCombos.contains(combo)) {
                                            tfCombos.writeln(combo);
                                            printedSNPProbeCombos.add(combo);
                                            combocounter++;
                                        }
                                        ctr2++;
                                    }
                                    SNP snpObj2 = genotypes.getSNPObject(snp);
                                    loader.loadGenotypes(snpObj2);

                                    double r2 = ldcalc.getRSquared(snpObj1, snpObj2, genotypes, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);

                                    snpObj2.clearGenotypes();

                                    tf.writeln(querySNP + "\t" + availableSNPs[snp] + "\t" + distance + "\t" + r2);
                                    ctr++;
                                }
                            }

                        }

                    }
                }
                System.out.println(querySNP + "\t" + probesForsnp.size() + "\t" + ctr + "\t" + ctr2 + "\t" + combocounter);
                snpObj1.clearGenotypes();
            } else {
                System.out.println("ERROR: " + querySNP + " not in reference!?");
            }
        }
        tfCombos.close();
        tf.close();


    }

    public void replaceProbeNames(String pb, String source, String dest, String in, String out) throws IOException {
        ProbeTranslation pt = new ProbeTranslation();
        HashMap<String, String> probeToProbe = pt.getProbeTranslation(pb, source, dest);

        TextFile infile = new TextFile(in, TextFile.R);
        TextFile outfile = new TextFile(out, TextFile.W);
        String[] elems = infile.readLineElems(TextFile.tab);

        while (elems != null) {
            String snp = elems[0];
            String pro = elems[1];

            String newpro = probeToProbe.get(pro);
            if (newpro == null || newpro.equals("-")) {
            } else {
                outfile.writeln(snp + "\t" + newpro);
            }

            elems = infile.readLineElems(TextFile.tab);
        }
        infile.close();
        outfile.close();




    }
}
