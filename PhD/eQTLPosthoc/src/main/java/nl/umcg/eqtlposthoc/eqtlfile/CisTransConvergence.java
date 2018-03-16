/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CisTransConvergence {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//        determinenumberoftestedsnppairspertrait(args);
//        determineOverlapBetweenPossiblePairsAndRealPairs("/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2012-03-ConvergenceAnalysis/PossibleSNPPairsPerTrait.txt",
//                "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2012-03-ConvergenceAnalysis/ConvergenceTransFDR0.5.txt",
//                "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2012-03-ConvergenceAnalysis/SignificantPairsPerTraitFDR0.5.txt");
        run(args);
    }

    public static void determinenumberoftestedsnppairspertrait(String[] args) {
        try {
            String gwasCatalogLoc = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt";
            String thousandGenomes = "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/Merged/";
            thousandGenomes = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
            double LDThreshold = 0.05;
            String probeTranslationFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";

            String outfileName = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2012-03-ConvergenceAnalysis/PossibleSNPPairsPerTrait.txt";

            String cisFileName = "";
            String transFileName = "";
            double pvalCutoffCisFDR = 1;

            String snpsTested = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-TestedGWASSNPs.txt";
            HashSet<String> limitSNPs = new HashSet<String>();
            TextFile tf = new TextFile(snpsTested, TextFile.R);
            limitSNPs.addAll(tf.readAsArrayList());
            tf.close();

            // now iterate through all GWAS traits..
            HashMap<GWASTrait, ArrayList<String>> snpsPerTrait = new HashMap<GWASTrait, ArrayList<String>>();

            ProbeTranslation pt = new ProbeTranslation();
            HashMap<String, String> probeToGene = pt.getProbeTranslation(probeTranslationFile, "Probe", "HUGO");

            HashMap<String, String> probeToHT12 = pt.getProbeTranslation(probeTranslationFile, "Probe", "HumanHT-12_V3_0_R2_11283641_A.txt");

            GWASCatalog catalog = new GWASCatalog();
            catalog.read(gwasCatalogLoc);

            HashMap<String, HashSet<String>> transFXPerSNP = new HashMap<String, HashSet<String>>();
            HashMap<String, HashSet<String>> cisFXPerSNP = new HashMap<String, HashSet<String>>();


            TriTyperGenotypeData gDs = new TriTyperGenotypeData();
            gDs.load(thousandGenomes);
            SNPLoader loader = gDs.createSNPLoader();

            DetermineLD ldcalc = new DetermineLD();

            GWASTrait[] traits = catalog.getTraits();

            TextFile out = new TextFile(outfileName, TextFile.W);
            TextFile outCtr = new TextFile(outfileName + "NumberOfPairsPerTrait.txt", TextFile.W);

            for (GWASTrait trait : traits) {
                GWASSNP[] traitsnps = trait.getSNPs();
                ArrayList<String> snps = new ArrayList<String>();
                for (GWASSNP snp : traitsnps) {
                    if (limitSNPs.contains(snp.getName())) {
                        snps.add(snp.getName());
                    }
                }

                int ctr = 0;
                if (snps.size() > 1) {
                    String[] snparr = snps.toArray(new String[0]);
                    for (int i = 0; i < snparr.length; i++) {
                        Integer snpObj1Id = gDs.getSnpToSNPId().get(snparr[i]);
                        if (snpObj1Id != null) {
                            SNP snpObj1 = gDs.getSNPObject(snpObj1Id);
                            loader.loadGenotypes(snpObj1);
                            for (int j = i + 1; j < snparr.length; j++) {
                                Integer snpObj2Id = gDs.getSnpToSNPId().get(snparr[j]);
                                if (snpObj2Id != null) {
                                    SNP snpObj2 = gDs.getSNPObject(snpObj2Id);

                                    // check whether both are on the same chr or within 5Mb
                                    if (snpObj1.getChr() != snpObj2.getChr()) {
                                        out.writeln(trait.getName() + "\t" + snpObj1.getName() + "\t" + snpObj2.getName());
                                        ctr++;
                                    } else {
                                        if (snpObj1.getChr() == 6 && snpObj2.getChr() == 6) {
                                            if ((snpObj1.getChrPos() > 20000000 && snpObj1.getChrPos() < 40000000) && (snpObj2.getChrPos() > 20000000 && snpObj2.getChrPos() < 40000000)) {
                                                // don't do anything
                                            }
                                        } else if (Math.abs(snpObj1.getChrPos() - snpObj2.getChrPos()) < 5000000) {
                                            // determine LD
                                            loader.loadGenotypes(snpObj2);

                                            double r2 = ldcalc.getRSquared(snpObj1, snpObj2, gDs, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                            if (!Double.isNaN(r2)) {
                                                if (r2 < LDThreshold) {
                                                    out.writeln(trait.getName() + "\t" + snpObj1.getName() + "\t" + snpObj2.getName());
                                                    ctr++;
                                                }
                                            }

                                            snpObj2.clearGenotypes();
                                        } else {
                                            out.writeln(trait.getName() + "\t" + snpObj1.getName() + "\t" + snpObj2.getName());
                                            ctr++;
                                        }
                                    }
                                }
                            }
                            snpObj1.clearGenotypes();
                        }
                    }

                }
                if (ctr > 0) {
                    outCtr.writeln(trait.getName() + "\t" + ctr);
                }
            }
            outCtr.close();
            out.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void run(String[] args) {
        // TODO code application logic here

        try {
//            String cisFileName = "/Volumes/iSnackHD/SkyDrive/latesteQTLs/cisFDR0.5.txt";
//            String transFileName = "/Volumes/iSnackHD/SkyDrive/latesteQTLs/transFDR0.05.txt.gz";



            HashSet<String> snpsAllowed = new HashSet<String>();

            TextFile tfgwas = new TextFile("/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLude.txt", TextFile.R);
            snpsAllowed.addAll(tfgwas.readAsArrayList());
            tfgwas.close();


            double[] fdrthresholds = new double[]{0.05, 0.5};

            for (double fdr : fdrthresholds) {
                System.out.println("");
                System.out.println("FDR: " + fdr);
                for (int perm = 0; perm < 21; perm++) {
                    System.out.println("Perm: " + perm);
                    System.out.println("");
                    System.out.println("");
                    System.out.println("");
                    System.out.println("");
                    System.out.println("");

                    String gwasCatalogLoc = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt";
                    String thousandGenomes = "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/Merged/";
                    thousandGenomes = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
                    double LDThreshold = 0.05;
                    String probeTranslationFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";

                    String cisFileName = "";
                    String transFileName = "";
                    String outfileName = "";
                    double pvalCutoffCisFDR = 1;

                    if (perm == 0) {
                        cisFileName = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR.txt.gz";
                        transFileName = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR" + fdr + ".txt.gz";
                        outfileName = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/convergence/ConvergenceTransFDR" + fdr + ".txt";
                        pvalCutoffCisFDR = 1.3141052405861965E-4;//0.00352058461224396;  // 0.00352058461224396

                    } else {
                        cisFileName = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR.txt.gz";
                        transFileName = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/TopFXInPermutedFiles/PermutationRound" + perm + "-FDR" + fdr + ".txt";
                        outfileName = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/convergence/ConvergenceTransRandomSet" + perm + "-FDR" + fdr + ".txt";
                        pvalCutoffCisFDR = 1.3141052405861965E-4;//0.00352058461224396;  // 0.00352058461224396
                    }




                    ProbeTranslation pt = new ProbeTranslation();
                    HashMap<String, String> probeToGene = pt.getProbeTranslation(probeTranslationFile, "Probe", "HUGO");

                    HashMap<String, String> probeToHT12 = pt.getProbeTranslation(probeTranslationFile, "Probe", "HumanHT-12_V3_0_R2_11283641_A.txt");

                    GWASCatalog catalog = new GWASCatalog();
                    catalog.read(gwasCatalogLoc);

                    HashMap<String, HashSet<String>> transFXPerSNP = new HashMap<String, HashSet<String>>();
                    HashMap<String, HashSet<String>> cisFXPerSNP = new HashMap<String, HashSet<String>>();

                    TextFile transFile = new TextFile(transFileName, TextFile.R);
                    String[] elems = transFile.readLineElems(TextFile.tab);
                    elems = transFile.readLineElems(TextFile.tab);
                    int transctr = 0;
                    while (elems != null) {
                        String snp = elems[1];
                        String probe = elems[4];

                        if (snpsAllowed == null || snpsAllowed.contains(snp)) {

                            HashSet<String> trans = transFXPerSNP.get(snp);
                            if (trans == null) {
                                trans = new HashSet<String>();
                            }
                            trans.add(probe);
                            transFXPerSNP.put(snp, trans);
                            transctr++;
                        }
                        elems = transFile.readLineElems(TextFile.tab);
                    }
                    transFile.close();
                    System.out.println(transctr + " trans effects loaded");

                    TextFile cisFile = new TextFile(cisFileName, TextFile.R);
                    elems = cisFile.readLineElems(TextFile.tab);
                    elems = cisFile.readLineElems(TextFile.tab);
                    int cisCtr = 0;
                    while (elems != null) {

                        double pval = Double.parseDouble(elems[0]);

                        if (pval <= pvalCutoffCisFDR) {
                            String snp = elems[1];
                            String probe = elems[4];
                            if (snpsAllowed == null || snpsAllowed.contains(snp)) {

                                HashSet<String> cis = cisFXPerSNP.get(snp);
                                if (cis == null) {
                                    cis = new HashSet<String>();
                                }
                                cis.add(probe);
                                cisFXPerSNP.put(snp, cis);
                                cisCtr++;
                            }
                        }
                        elems = cisFile.readLineElems(TextFile.tab);

                    }
                    transFile.close();
                    System.out.println(cisCtr + " cis effects loaded");

                    GWASSNP test = catalog.getSnpToObj().get("rs1738074");
                    System.out.println(Strings.concat(test.getAssociatedTraitsArray(), Strings.tab));

                    TriTyperGenotypeData genotypeData = new TriTyperGenotypeData();
                    genotypeData.load(thousandGenomes);
                    SNPLoader loader = genotypeData.createSNPLoader();

                    DetermineLD ldcalc = new DetermineLD();

                    TextFile out = new TextFile(outfileName, TextFile.W);
                    TextFile outtraits = new TextFile(outfileName + "NumberOfConvergentEffectsPerTrait.txt", TextFile.W);
//            out.writeln("Trait\tSNP 1\tSNP 2\tLD\tProbe\tGene\tCisOrTrans\tQC");
                    out.writeln("Trait\tSNP 1\tSNP 1 Chromosome\tSNP 1 Chromosome position\tSNP 2\tSNP 2 Chromosome\tSNP 2 Chromosome position\tIllumina HT12v3 Array Address ID\tGene\tCisOrTrans\tQC");
                    GWASTrait[] traits = catalog.getTraits();
                    ProgressBar pb = new ProgressBar(traits.length);
//            for (GWASTrait t : traits) {
//            GWASSNP[] snps = catalog.getSnpsArray();

                    // convergence can be trans only as well....
                    // we need a pair-wise transcomparison for this
                    HashSet<String> pairsOutputed = new HashSet<String>();
                    int notIn1KG = 0;

                    HashSet<String> convergentTraits = new HashSet<String>();
                    for (GWASTrait t : traits) {
                        GWASSNP[] snps = t.getSNPs();
                        HashSet<String> convergentSNPs = new HashSet<String>();
                        HashSet<String> convertentProbes = new HashSet<String>();

                        for (int s = 0; s < snps.length; s++) {
                            String snp1Name = snps[s].getName();
                            HashSet<String> transFx = transFXPerSNP.get(snp1Name);
                            HashSet<String> cisFx = cisFXPerSNP.get(snp1Name);
                            GWASTrait[] traits1 = snps[s].getAssociatedTraitsArray();

                            for (int s2 = s + 1; s2 < snps.length; s2++) {
                                String snp2Name = snps[s2].getName();
                                GWASTrait[] traits2 = snps[s2].getAssociatedTraitsArray();
                                HashSet<String> transFx2 = transFXPerSNP.get(snp2Name);
                                HashSet<String> cisFx2 = cisFXPerSNP.get(snp2Name);

                                String QCString = "";
                                Integer snpId1 = genotypeData.getSnpToSNPId().get(snp1Name);
                                Integer snpId2 = genotypeData.getSnpToSNPId().get(snp2Name);
                                boolean canDetermineLD = true;
                                boolean failsQC = false;

                                if (snpId1 == null) {
                                    canDetermineLD = false;
                                    QCString += snp1Name + " not present in 1KG;";
                                    failsQC = true;
                                }
                                if (snpId2 == null) {
                                    canDetermineLD = false;
                                    QCString += snp2Name + " not present in 1KG";
                                    failsQC = true;
                                }

                                if (snpId1 != null && snpId2 != null) {
                                    SNP snpObj1 = genotypeData.getSNPObject(snpId1);
                                    SNP snpObj2 = genotypeData.getSNPObject(snpId2);
                                    if (canDetermineLD) {

                                        if (snpObj1.getChr() != snpObj2.getChr()) {
                                            QCString += "SNP1 and SNP2 are on different chromosomes (" + ChrAnnotation.parseByte(snpObj1.getChr()) + " / " + ChrAnnotation.parseByte(snpObj2.getChr()) + ")";
                                            failsQC = false;
                                        } else {
                                            if (snpObj1.getChr() == 6 && snpObj2.getChr() == 6) {
                                                if ((snpObj1.getChrPos() > 20000000 && snpObj1.getChrPos() < 40000000) && (snpObj2.getChrPos() > 20000000 && snpObj2.getChrPos() < 40000000)) {
                                                    failsQC = true;
                                                }
                                            }
                                            if (Math.abs(snpObj1.getChrPos() - snpObj2.getChrPos()) < 5000000) {
                                                loader.loadGenotypes(snpObj1);
                                                loader.loadGenotypes(snpObj2);

//                                                if (snpObj1.getMAF() >= 0.05 && snpObj2.getMAF() > 0.05) {
                                                double r2 = ldcalc.getRSquared(snpObj1, snpObj2, genotypeData, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                                if (Double.isNaN(r2)) {
                                                    failsQC = true;
                                                } else {
                                                    QCString += "" + r2;
                                                    if (r2 > LDThreshold) {
                                                        failsQC = true;
                                                    }
                                                }
//                                                } else {
//                                                    System.out.println("MAF is not sufficient for SNP 1 or SNP2: " + snpObj1.getName() + " - " + snpObj1.getMAF() + "\t" + snpObj2.getName() + " - " + snpObj2.getMAF());
//                                                    failsQC = true;
//                                                }


                                                snpObj1.clearGenotypes();
                                                snpObj2.clearGenotypes();
                                            } else {
                                                QCString += "distance between SNP1 and SNP2 > 5Mb";
                                                failsQC = false;
                                            }
                                        }
                                    } else {
//                            // check from the GWAS catalog where the SNPs actually are...
//                            if (snps[s].getChr() != snps[s2].getChr()) {
//                                QCString += "SNP1 and SNP2 are on different chromosomes (" + ChrAnnotation.parseByte(snps[s].getChr()) + " / " + ChrAnnotation.parseByte(snps[s2].getChr()) + ")";
//                                failsQC = false;
//                            }

                                        failsQC = true;
                                    }

                                    if (!failsQC) {




                                        if (transFx != null && transFx2 != null) {
                                            for (String probe : transFx) {
                                                if (transFx2.contains(probe)) {
                                                    String output = t.getName() + "\t" + snp1Name + "\t" + ChrAnnotation.parseByte(snpObj1.getChr()) + "\t" + snpObj1.getChrPos() + "\t" + snp2Name + "\t" + ChrAnnotation.parseByte(snpObj2.getChr()) + "\t" + snpObj2.getChrPos() + "\t" + probeToHT12.get(probe) + "\t" + probeToGene.get(probe) + "\tTrans-Trans\t" + QCString;
                                                    out.writeln(output);
//                                        System.out.println(output);
                                                    convergentTraits.add(t.getName());
                                                    convergentSNPs.add(snp1Name);
                                                    convergentSNPs.add(snp2Name);
                                                    convertentProbes.add(probe);
                                                }
                                            }
                                        }

                                        if (cisFx != null && transFx2 != null) {
                                            for (String probe : cisFx) {
                                                if (transFx2.contains(probe)) {
                                                    String output = t.getName() + "\t" + snp1Name + "\t" + ChrAnnotation.parseByte(snpObj1.getChr()) + "\t" + snpObj1.getChrPos() + "\t" + snp2Name + "\t" + ChrAnnotation.parseByte(snpObj2.getChr()) + "\t" + snpObj2.getChrPos() + "\t" + probeToHT12.get(probe) + "\t" + probeToGene.get(probe) + "\tCis-Trans\t" + QCString;
                                                    out.writeln(output);
//                                        System.out.println(output);
                                                    convergentTraits.add(t.getName());
                                                    convergentSNPs.add(snp1Name);
                                                    convergentSNPs.add(snp2Name);
                                                    convertentProbes.add(probe);
                                                }
                                            }
                                        }

                                        if (cisFx2 != null && transFx != null) {
                                            for (String probe : cisFx2) {
                                                if (transFx.contains(probe)) {
                                                    String output = t.getName() + "\t" + snp1Name + "\t" + ChrAnnotation.parseByte(snpObj1.getChr()) + "\t" + snpObj1.getChrPos() + "\t" + snp2Name + "\t" + ChrAnnotation.parseByte(snpObj2.getChr()) + "\t" + snpObj2.getChrPos() + "\t" + probeToHT12.get(probe) + "\t" + probeToGene.get(probe) + "\tTrans-Cis\t" + QCString;
                                                    out.writeln(output);
//                                        System.out.println(output);
                                                    convergentTraits.add(t.getName());
                                                    convergentSNPs.add(snp1Name);
                                                    convergentSNPs.add(snp2Name);
                                                    convertentProbes.add(probe);
                                                }
                                            }
                                        }

                                    }

                                }



                            }


                        }
                        if (convergentSNPs.size() > 0) {
                            System.out.println(t.getName() + "\tConvergent SNPs: " + convergentSNPs.size() + " / " + t.getSNPs().length + "\tConvergent Probes: " + convertentProbes.size());
                            outtraits.writeln(t.getName() + "\tConvergent SNPs: " + convergentSNPs.size() + " / " + t.getSNPs().length + "\tConvergent Probes: " + convertentProbes.size());
                        }
                    }

                    System.out.println("Convergent traits: " + convergentTraits.size());

//            for (int s = 0; s < snps.length; s++) {
//                String snp1Name = snps[s].getName();
//                HashSet<String> transFx = transFXPerSNP.get(snp1Name);
//                HashSet<String> cisFx = cisFXPerSNP.get(snp1Name);
//                GWASTrait[] traits1 = snps[s].getAssociatedTraitsArray();
//
//                for (int s2 = s + 1; s2 < snps.length; s2++) {
//                    String snp2Name = snps[s2].getName();
//                    GWASTrait[] traits2 = snps[s2].getAssociatedTraitsArray();
//                    HashSet<String> transFx2 = transFXPerSNP.get(snp2Name);
//                    HashSet<String> cisFx2 = cisFXPerSNP.get(snp2Name);
//
//                    String QCString = "";
//                    Integer snpId1 = genotypeData.getSnpToSNPId().get(snp1Name);
//                    Integer snpId2 = genotypeData.getSnpToSNPId().get(snp2Name);
//                    boolean canDetermineLD = true;
//                    boolean failsLD = false;
//
//                    if (snpId1 == null) {
//                        canDetermineLD = false;
//                        QCString += snp1Name + " not present in 1KG;";
//                    }
//                    if (snpId2 == null) {
//                        canDetermineLD = false;
//                        QCString += snp2Name + " not present in 1KG";
//                    }
//                    if (snpId1 != null && snpId2 != null) {
//                        SNP snpObj1 = genotypeData.getSNPObject(snpId1);
//                        SNP snpObj2 = genotypeData.getSNPObject(snpId2);
//                        if (snpObj1.getChr() != snpObj2.getChr()) {
//                            QCString += "SNP1 and SNP2 are on different chromosomes";
//                        } else {
//                            if (Math.abs(snpObj1.getChrPos() - snpObj2.getChrPos()) < 5000000) {
//                                loader.loadGenotypes(snpObj1);
//                                loader.loadGenotypes(snpObj2);
//                                double r2 = ldcalc.getRSquared(snpObj1, snpObj2, genotypeData, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
//                                QCString += "" + r2;
//                                snpObj1.clearGenotypes();
//                                snpObj2.clearGenotypes();
//                            } else {
//                                QCString += "distance between SNP1 and SNP2 > 5Mb";
//                            }
//                        }
//                    }
//
//
//                    if (transFx != null && transFx2 != null) {
//                        for (String probe : transFx) {
//                            if (transFx2.contains(probe)) {
//                                String output = snp1Name + "\t" + Strings.concat(traits1, Strings.semicolon) + "\t" + snp2Name + "\t" + Strings.concat(traits2, Strings.semicolon) + "\t" + probe + "\t" + probeToGene.get(probe) + "\tTrans-Trans\t" + QCString;
//                                out.writeln(output);
//                                System.out.println(output);
//                            }
//                        }
//                    }
//
//                    if (cisFx != null && transFx2 != null) {
//                        for (String probe : cisFx) {
//                            if (transFx2.contains(probe)) {
//                                String output = snp1Name + "\t" + Strings.concat(traits1, Strings.semicolon) + "\t" + snp2Name + "\t" + Strings.concat(traits2, Strings.semicolon) + "\t" + probe + "\t" + probeToGene.get(probe) + "\tCis-Trans\t" + QCString;
//                                out.writeln(output);
//                                System.out.println(output);
//                            }
//                        }
//                    }
//
//                    if (cisFx2 != null && transFx != null) {
//                        for (String probe : cisFx2) {
//                            if (transFx.contains(probe)) {
//                                String output = snp1Name + "\t" + Strings.concat(traits1, Strings.semicolon) + "\t" + snp2Name + "\t" + Strings.concat(traits2, Strings.semicolon) + "\t" + probe + "\t" + probeToGene.get(probe) + "\tTrans-Cis" + QCString;
//                                out.writeln(output);
//                                System.out.println(output);
//                            }
//                        }
//                    }
//
//                }
//
//
//            }
//            pb.close();
                    outtraits.close();
                    out.close();
                    loader.close();

                }
            }


        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static String toString(HashSet<String> cisFx, HashMap<String, String> probeToGene) {
        if (cisFx != null) {
            String[] probes = cisFx.toArray(new String[0]);
            String[] fx = new String[cisFx.size()];
            for (int i = 0; i < fx.length; i++) {
                fx[i] = probeToGene.get(probes[i]);
            }
            return Strings.concat(fx, Strings.semicolon);
        } else {
            return null;
        }
    }

    private static void determineOverlapBetweenPossiblePairsAndRealPairs(String allPairsPerTrait, String significantPairsPerTrait, String outputfilename) {
        try {

            TextFile tf = new TextFile(significantPairsPerTrait, TextFile.R);
            tf.readLine(); // skip header
            String[] elems = tf.readLineElems(TextFile.tab);
            HashSet<Pair<String, String>> significantPairs = new HashSet<Pair<String, String>>();
            while (elems != null) {
                String snp1 = elems[1];
                String snp2 = elems[4];

                significantPairs.add(new Pair<String, String>(snp1, snp2));
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            System.out.println(significantPairs.size() + " significant SNP pairs");

            TextFile tf2 = new TextFile(allPairsPerTrait, TextFile.R);

            elems = tf2.readLineElems(TextFile.tab);
            HashMap<String, HashSet<Pair<String, String>>> pairsPerTrait = new HashMap<String, HashSet<Pair<String, String>>>();
            HashSet<String> traits = new HashSet<String>();
            while (elems != null) {

                String trait = elems[0];
                String snp1 = elems[1];
                String snp2 = elems[2];


                HashSet<Pair<String, String>> pairs = pairsPerTrait.get(trait);
                if (pairs == null) {
                    pairs = new HashSet<Pair<String, String>>();
                }

                pairs.add(new Pair<String, String>(snp1, snp2));
                pairsPerTrait.put(trait, pairs);

                traits.add(trait);
                elems = tf2.readLineElems(TextFile.tab);
            }

            tf2.close();


            TextFile out = new TextFile(outputfilename, TextFile.W);
            for (String trait : traits) {
                int ctr = 0;
                HashSet<Pair<String, String>> pairs = pairsPerTrait.get(trait);
                for (Pair<String, String> p : pairs) {
                    String snp11 = p.getLeft();
                    String snp21 = p.getRight();
                    boolean match = false;
                    for (Pair<String, String> p2 : significantPairs) {
                        String snp12 = p2.getLeft();
                        String snp22 = p2.getRight();

                        if (snp12.equals(snp11) && snp22.equals(snp21)) {
                            // got one
                            match = true;
                        } else if (snp22.equals(snp11) && snp12.equals(snp21)) {
                            // got one
                            match = true;
                        } else {
                            // no match
                        }
                    }

                    if (match) {
                        ctr++;
                    }
                }

                if (ctr > 0) {
                    out.writeln(trait + "\t" + pairs.size() + "\t" + ctr + "\t" + ((double) ctr / pairs.size()));
                }

            }

            out.close();




        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
