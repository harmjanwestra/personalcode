package nl.harmjanwestra.playground.trityper;

import gnu.trove.map.hash.TObjectIntHashMap;
import org.apache.commons.cli.*;
import org.checkerframework.checker.units.qual.A;
import umcg.genetica.console.MultiThreadProgressBar;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.text.Strings;

import java.io.Console;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class LDCalculator {

    public static void main(String[] args) {
        LDCalculator ld = new LDCalculator();


        Options options = new Options();
        options.addRequiredOption("m", "mode", true, "mode: [compareall|traitprune|traitprune2]; compareall: no grouping per trait; traitprune: prune per trait and ld; traitprune2: ld and prune per trait");
        options.addRequiredOption("g", "gwas", true, "GWAS snp set");
        options.addRequiredOption("e", "eqtl", true, "eQTL snp set");
        options.addOption("ef", "eqtlfile", true, "eQTL file");

        options.addOption("gl", "gwaslist", true, "GWAS list file");
        options.addOption("ga", "gwasassoc", true, "GWAS assoc file");

        options.addRequiredOption("t", "trityper", true, "TriTyper genotype reference");
        options.addOption("i", "individuals", true, "Individuals to select from TriTyper reference");
        options.addOption("l", "ldthreshold", true, "LD threshold");
        options.addOption("lp", "ldprunethreshold", true, "LD prune threshold");
        options.addOption("lpd", "ldtprunedistance", true, "LD distance threshold");
        options.addOption("s6", "skipchr6", false, "Skip chr 6");

        options.addRequiredOption("o", "output", true, "Output dir/file");

        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine cmd = parser.parse(options, args);

            String mode = cmd.getOptionValue("m");
            switch (mode) {
                case "compareall": {
                    String setfile1 = cmd.getOptionValue("e");
                    String setfile2 = cmd.getOptionValue("g");
                    String trityper = cmd.getOptionValue("t");
                    String output = cmd.getOptionValue("o");

                    double ldthreshold = 0.8;
                    String ldThreshStr = cmd.getOptionValue("l");
                    if (ldThreshStr != null) {
                        ldthreshold = Double.parseDouble(ldThreshStr);
                    }
                    String individualsubset = cmd.getOptionValue("i");
                    boolean skipchr6 = cmd.hasOption("s6");
                    ld.run(setfile1, setfile2, trityper, individualsubset, ldthreshold, output, skipchr6);
                }
                break;
                case "traitprune": {
                    String gwasset = cmd.getOptionValue("g");
                    String eqtlset = cmd.getOptionValue("e");
                    String gwasassocfile = cmd.getOptionValue("ga");
                    String gwaslistfile = cmd.getOptionValue("gl");
                    String trityperdir = cmd.getOptionValue("t");
                    String output = cmd.getOptionValue("o");

                    double ldthreshold = 0.8;
                    String ldThreshStr = cmd.getOptionValue("l");
                    if (ldThreshStr != null) {
                        ldthreshold = Double.parseDouble(ldThreshStr);
                    }

                    double prunethreshold = 0.2;
                    String prunethresholdStr = cmd.getOptionValue("lp");
                    if (prunethresholdStr != null) {
                        prunethreshold = Double.parseDouble(prunethresholdStr);
                    }

                    int prunedistance = 1000000;
                    String prunedistanceStr = cmd.getOptionValue("lpd");
                    if (prunedistanceStr != null) {
                        prunedistance = Integer.parseInt(prunedistanceStr);
                    }

                    String individualsubset = cmd.getOptionValue("i");
                    boolean skipchr6 = cmd.hasOption("s6");
                    ld.runPerGWAS(eqtlset, gwasset, gwasassocfile, gwaslistfile, trityperdir, individualsubset, prunedistance, prunethreshold, ldthreshold, output, skipchr6);
                }
                break;
                case "traitprune2": {
                    String gwasset = cmd.getOptionValue("g");
                    String eqtlset = cmd.getOptionValue("e");
                    String eqtlfile = cmd.getOptionValue("ef");
                    String gwasassocfile = cmd.getOptionValue("ga");
                    String gwaslistfile = cmd.getOptionValue("gl");
                    String trityperdir = cmd.getOptionValue("t");
                    String output = cmd.getOptionValue("o");

                    double ldthreshold = 0.8;
                    String ldThreshStr = cmd.getOptionValue("l");
                    if (ldThreshStr != null) {
                        ldthreshold = Double.parseDouble(ldThreshStr);
                    }

                    double prunethreshold = 0.2;
                    String prunethresholdStr = cmd.getOptionValue("lp");
                    if (prunethresholdStr != null) {
                        prunethreshold = Double.parseDouble(prunethresholdStr);
                    }

                    int prunedistance = 1000000;
                    String prunedistanceStr = cmd.getOptionValue("lpd");
                    if (prunedistanceStr != null) {
                        prunedistance = Integer.parseInt(prunedistanceStr);
                    }

                    String individualsubset = cmd.getOptionValue("i");
                    boolean skipchr6 = cmd.hasOption("s6");
                    ld.runPerGWAS2(eqtlset, eqtlfile, gwasset, gwasassocfile, gwaslistfile, trityperdir, individualsubset, prunedistance, prunethreshold, ldthreshold, output, skipchr6);
                }
                break;

                default:
                    System.out.println("Invalid mode:");
                    HelpFormatter formatter = new HelpFormatter();
                    formatter.printHelp("ant", options);
                    break;
            }


        } catch (ParseException e) {
            System.out.println(e.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("ant", options);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void waitForEnter(String message, Object... args) {
        Console c = System.console();
        if (c != null) {
            // printf-like arguments
            if (message != null)
                c.format(message, args);
            c.format("\nPress ENTER to proceed.\n");
            c.readLine();
        }
    }


    public void runPerGWAS2(String eqtlset, String eqtlfile, String gwasset, String gwasAssocFile, String gwasListFile,
                            String trityperdir, String individualSubset,
                            int prunedistance, double prunethreshold,
                            double ldthreshold, String outputfile, boolean skipchr6) throws IOException {

        System.out.println("Skipping chr6: " + skipchr6);
        System.out.println("RSquared threshold: " + ldthreshold);
        System.out.println("Prune threshold: " + prunethreshold);
        System.out.println("Prune distance: " + prunedistance);

        HashSet<String> gwasSet = loadSet(gwasset);
        HashSet<String> eqtlSet = loadSet(eqtlset);

        System.out.println("gwasSet: " + gwasset + " has " + gwasSet.size() + " snps");
        System.out.println("eqtlSet: " + eqtlset + " has " + eqtlSet.size() + " snps");

        // load eQTL info
        HashMap<String, ArrayList<eQTLResult>> eqtlsPerSNP = loadEQTLInfo(eqtlfile);

        // load GWAS info
        Pair<HashMap<String, String>, HashMap<String, HashMap<String, GWASResult>>> gwasinfo = loadGWASInfo(gwasAssocFile, gwasListFile);
        HashMap<String, String> gwasIdToTrait = gwasinfo.getLeft();
        HashMap<String, HashMap<String, GWASResult>> gwasSNPs = gwasinfo.getRight();

        System.out.println("Loading data from: " + trityperdir);

        HashSet<String> allSNPsSet = new HashSet<>();
        allSNPsSet.addAll(gwasSet);
        allSNPsSet.addAll(eqtlSet);
        ArrayList<String> allSNPs = new ArrayList<>();
        allSNPs.addAll(allSNPsSet);
        System.out.println(allSNPs.size() + " SNPs total");
        Collections.sort(allSNPs);
        // preload genotypes

        TriTyperGenotypeData data = new TriTyperGenotypeData(trityperdir);
        if (individualSubset != null) {
            HashSet<String> indsToInclude = loadSet(individualSubset);
            System.out.println("Filtering individuals using: " + individualSubset + " having " + indsToInclude.size() + " individuals");
            String[] individuals = data.getIndividuals();
            Boolean[] include = new Boolean[individuals.length];
            for (int i = 0; i < individuals.length; i++) {
                if (indsToInclude.contains(individuals[i])) {
                    include[i] = true;
                } else {
                    include[i] = false;
                }
            }
            data.setIsIncluded(include);
        }


        HashMap<String, SNP> genotypes = loadSNPs(allSNPs, data, skipchr6);
        System.out.println(genotypes.size() + " snps found in total.");

        TextFile tfq = new TextFile(outputfile + "-eQTLSNPsPresent.txt", TextFile.W);
        for (String s : eqtlSet) {
            if (genotypes.containsKey(s)) {
                tfq.writeln(s);
            }
        }
        tfq.close();

        tfq = new TextFile(outputfile + "-GWASSNPsPresent.txt", TextFile.W);
        for (String s : gwasSet) {
            if (genotypes.containsKey(s)) {
                tfq.writeln(s);
            }
        }
        tfq.close();

        // get a list of pruned eQTL SNPs
        ArrayList<String> alleqtl = new ArrayList<>();
        alleqtl.addAll(eqtlSet);

        TextFile outPruneEQTL = new TextFile(outputfile + "-PrunedSNPsPerEQTL.txt", TextFile.W);
        outPruneEQTL.writeln("SNP\tAliases");
        HashSet<String> visitedEQTLSNPs = new HashSet<>();
        for (int i = 0; i < alleqtl.size(); i++) {
            String snp1 = alleqtl.get(i);
            SNP snpobj1 = genotypes.get(snp1);
            if (snpobj1 != null) {
                if (!visitedEQTLSNPs.contains(snp1)) {

                    DetermineLD ld = new DetermineLD();
                    ArrayList<String> aliases = new ArrayList<>();
                    visitedEQTLSNPs.add(snp1);
                    for (int j = i + 1; j < alleqtl.size(); j++) {
                        String snp2 = alleqtl.get(j);
                        SNP snpobj2 = genotypes.get(snp2);
                        if (snpobj2 != null) {
                            if (!visitedEQTLSNPs.contains(snp2)) {
                                if (snpobj1.getChr() == snpobj2.getChr() && Math.abs(snpobj1.getChrPos() - snpobj2.getChrPos()) < prunedistance) {
                                    // check LD
                                    Pair<Double, Double> ldinfo = ld.getLD(snpobj1, snpobj2, data, 1, false);
                                    double rsq = ldinfo.getRight();
                                    if (rsq > ldthreshold) {
                                        aliases.add(snp2);
                                        visitedEQTLSNPs.add(snp2);
                                    }
                                }
                            }
                        }
                    }
                    if (aliases.isEmpty()) {
                        outPruneEQTL.writeln(snp1 + "\t" + "-");
                    } else {
                        outPruneEQTL.writeln(snp1 + "\t" + Strings.concat(aliases, Strings.semicolon));
                    }

                }
            }

        }
        outPruneEQTL.close();

        // for each GWAS, determine loaded SNPs, and their associated p-values
        ArrayList<String> allgwas = new ArrayList<>();
        allgwas.addAll(gwasSNPs.keySet());

        TextFile outAll = new TextFile(outputfile + "-all.txt", TextFile.W);
        String header = "GWASId\tTrait\tTraitP\tTraitSNP1\tMAF\tHWEP\tSNP2\tMAF\tHWEP\tRSQ\tDPR\tENSG\tGene";
        outAll.writeln(header);

        TextFile outSum = new TextFile(outputfile + "-summary.txt", TextFile.W);
        String header2 = "GWASID\tTrait\tNrGWASSNPsTotal\tNrPrunedGWASSNPsLinked\tNrPrunedGWASSNPs\tPercGWASSNPsLinked\tnrEQTLSNPsLinked\tTotalEQTLSNPs\tPercEQTLSNPsLinked";
        outSum.writeln(header2);

        TextFile outPrune = new TextFile(outputfile + "-PrunedSNPsPerTrait.txt", TextFile.W);
        String pruneheader = "GWASID\tTrait\tIndexVariant\tIndexVariantP\tIndexVariantAlleles\tIndexVariantEffect\tLinkedEQTLSNP\tLD(rsq)\tLinkedEQTLGenes\tLinkedEQTLHgncIDs\tLinkedEQTLAlleles\tLinkedEQTLZScores\tLinkedEQTLP\tGWASClusterSize\tSNPsInGWASCluster";
        outPrune.writeln(pruneheader);


        // count eQTL SNPs with genotypes
        int nrEQTLSNPs = 0;
        for (String snp : eqtlSet) {
            if (genotypes.containsKey(snp)) {
                nrEQTLSNPs++;
            }
        }

        System.out.println(allgwas.size() + " GWASes to test");
        int finalNrEQTLSNPs = nrEQTLSNPs;

//        IntStream.range(0, allgwas.size()).parallel().forEach(g -> {
        ProgressBar pb = new ProgressBar(allgwas.size(), "Calculating LD.");
        IntStream.range(0, allgwas.size()).forEach(g -> {
            String gwasID = allgwas.get(g);
            HashMap<String, GWASResult> gwasSNPsPValues = gwasSNPs.get(gwasID);

            // find gwasSNPsPValues that are within 1mb of each other
            ArrayList<String> allGWASSNPs = new ArrayList<>();
            allGWASSNPs.addAll(gwasSNPsPValues.keySet());
            Collections.sort(allGWASSNPs);
//            HashSet<String> visited = new HashSet<>();

            // determine linked variants
            HashMap<String, HashSet<String>> allLinkedSNPs = new HashMap<>();
            HashMap<String, HashMap<String, Double>> ldforGWASSNPs = new HashMap<>();
            IntStream.range(0, allGWASSNPs.size()).forEach(i -> {
                DetermineLD ld = new DetermineLD();
                String gwasSNPId = allGWASSNPs.get(i);
                SNP gwasSNPObj = genotypes.get(gwasSNPId);

                HashSet<String> linkedSNPs = new HashSet<>();
                if (gwasSNPObj != null) {
                    byte chr = gwasSNPObj.getChr();
                    int pos = gwasSNPObj.getChrPos();
                    for (String eqtlSNP : eqtlSet) {
                        SNP eQTLSNPObj = genotypes.get(eqtlSNP);
                        if (eQTLSNPObj != null) {
                            byte chr2 = eQTLSNPObj.getChr();
                            int pos2 = eQTLSNPObj.getChrPos();
                            if (chr == chr2 && Math.abs(pos - pos2) < 5000000) {
                                Pair<Double, Double> ldvals = ld.getLD(gwasSNPObj, eQTLSNPObj, data, 1, false);
                                double rsq = ldvals.getRight();
                                double dpr = ldvals.getLeft();

                                synchronized (ldforGWASSNPs) {
                                    HashMap<String, Double> ldsnps = ldforGWASSNPs.get(gwasSNPId);
                                    if (ldsnps == null) {
                                        ldsnps = new HashMap<>();
                                    }
                                    ldsnps.put(eqtlSNP, rsq);
                                    ldforGWASSNPs.put(gwasSNPId, ldsnps);
                                }


                                if (rsq > ldthreshold) {
                                    linkedSNPs.add(eqtlSNP);
                                    ArrayList<eQTLResult> genedata = eqtlsPerSNP.get(eqtlSNP);
                                    ArrayList<String> genes = new ArrayList<>();
                                    ArrayList<String> hugos = new ArrayList<>();
                                    String ensg = "-";
                                    String hugo = "-";
                                    if (genedata != null) {
                                        for (eQTLResult p : genedata) {
                                            genes.add(p.gene);
                                            hugos.add(p.hugo);
                                        }
                                        ensg = Strings.concat(genes, Strings.semicolon);
                                        hugo = Strings.concat(hugos, Strings.semicolon);
                                    }

                                    try {
                                        outAll.writelnsynced(gwasID + "\t" + gwasIdToTrait.get(gwasID) + "\t" + gwasSNPsPValues.get(gwasSNPObj.getName())
                                                + "\t" + gwasSNPObj.getName() + "\t" + gwasSNPObj.getMAF() + "\t" + gwasSNPObj.getHWEP()
                                                + "\t" + eqtlSNP + "\t" + eQTLSNPObj.getMAF() + "\t" + eQTLSNPObj.getHWEP()
                                                + "\t" + rsq + "\t" + dpr
                                                + "\t" + ensg + "\t" + hugo);
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }

                                }
                            }
                        }
                    }
                }
                synchronized (allLinkedSNPs) {
                    HashSet<String> linked = allLinkedSNPs.get(gwasSNPId);
                    if (linked == null) {
                        linked = new HashSet<>();
                    }
                    linked.addAll(linkedSNPs);
                    allLinkedSNPs.put(gwasSNPId, linked);
                }
            });


            // prune / clump GWAS set
            HashSet<String> visited = new HashSet<String>();
            int nrGWASClustersLinked = 0;
            int nrGWASClusters = 0;
            HashSet<String> eqtlSNPsLinkedToGWASSNPClusters = new HashSet<String>();

            for (int s = 0; s < allGWASSNPs.size(); s++) {
                String gwasSNP1 = allGWASSNPs.get(s);
                SNP gwasSNP1Obj = genotypes.get(gwasSNP1);

                String indexVariant = gwasSNP1;

                if (gwasSNP1Obj != null && !visited.contains(gwasSNP1)) {
                    byte chr = gwasSNP1Obj.getChr();
                    Integer pos = gwasSNP1Obj.getChrPos(); // Integer.parseInt(elems[1]);

                    // check whether there are any SNPs linked to the current one
                    // determine which gwasSNPsPValues are within prunedistance of each other that are in LD
                    // if so, gather all linked eQTLs
                    ArrayList<Pair<String, Double>> gwasCluster = new ArrayList<>();
                    HashSet<String> eqtlSNPsLinkedToCluster = new HashSet<String>();
                    GWASResult gwasResult = gwasSNPsPValues.get(gwasSNP1);
                    Pair<String, Double> pair = new Pair<String, Double>(gwasSNP1, gwasResult.p, Pair.SORTBY.RIGHT);
                    gwasCluster.add(pair);
                    visited.add(gwasSNP1);
                    for (int s2 = s + 1; s2 < allGWASSNPs.size(); s2++) {
                        String gwasSNP2 = allGWASSNPs.get(s2);
                        SNP gwasSNP2Obj = genotypes.get(gwasSNP2);
                        if (gwasSNP2Obj != null && !visited.contains(gwasSNP2)) {
                            byte chr2 = gwasSNP2Obj.getChr();
                            Integer pos2 = gwasSNP2Obj.getChrPos(); // Integer.parseInt(elems2[1]);
                            if (chr == chr2) {
                                int distance = Math.abs(pos - pos2);
                                if (distance < prunedistance) {
                                    // determine LD
                                    DetermineLD ld = new DetermineLD();
                                    SNP snpjobj = genotypes.get(gwasSNP2);
                                    Pair<Double, Double> ldvals = ld.getLD(gwasSNP1Obj, snpjobj, data, 1, false);
                                    double rsq = ldvals.getRight();
                                    if (rsq > prunethreshold) {
                                        GWASResult gwasResult2 = gwasSNPsPValues.get(gwasSNP2);
                                        Pair<String, Double> pair2 = new Pair<String, Double>(gwasSNP2, gwasResult2.p, Pair.SORTBY.RIGHT);
                                        gwasCluster.add(pair2);
                                        visited.add(gwasSNP2);

                                        // get all linked eQTLs
                                        HashSet<String> linked = allLinkedSNPs.get(gwasSNP2);
                                        // add to 'gwasCluster' set
                                        if (linked != null) {
                                            eqtlSNPsLinkedToCluster.addAll(linked);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // prune gwasCluster
                    ArrayList<SNP> remainingSNPs = new ArrayList<>();
                    if (gwasCluster.size() == 1) {
                        String snpi = gwasCluster.get(0).getLeft();
                        SNP snpiobj = genotypes.get(snpi);
                        remainingSNPs.add(snpiobj);
                        HashSet<String> linked = allLinkedSNPs.get(snpi);
                        if (linked != null) {
                            eqtlSNPsLinkedToCluster.addAll(linked);
                        }
                    } else if (gwasCluster.size() > 1) {
                        // sort gwasCluster on P
                        Collections.sort(gwasCluster);
                        int min = Math.min(5, gwasCluster.size());
                        for (int i = 0; i < min; i++) {
                            System.out.println(gwasCluster.get(i));
                        }
                        System.out.println();
                        indexVariant = gwasCluster.get(0).getLeft();
//                        System.out.println(gwasCluster.size() + " SNPs in gwasCluster " + clusterctr.get() + ", " + remainingSNPs.size() + " remain after pruning for trait: " + gwasID + "\t" + gwasIdToTrait.get(gwasID));
//                        waitForEnter("Check this out!");
                    }

                    // eqtlSNPsLinkedToCluster has all linked SNPs for this gwasCluster.
                    if (!eqtlSNPsLinkedToCluster.isEmpty()) {
                        nrGWASClustersLinked++;
                        eqtlSNPsLinkedToGWASSNPClusters.addAll(eqtlSNPsLinkedToCluster);
                    }
                    nrGWASClusters++;

                    String clusterStr = "";
                    for (Pair<String, Double> clusterpair : gwasCluster) {
                        String snp = clusterpair.getLeft();
                        if (!snp.equals(gwasSNP1)) {
                            if (clusterStr.length() == 0) {
                                clusterStr += snp;
                            } else {
                                clusterStr += ";" + snp;
                            }
                        }
                    }


                    for (String eqtl : eqtlSNPsLinkedToCluster) {
                        ArrayList<eQTLResult> eqtlinfo = eqtlsPerSNP.get(eqtl);
                        if (eqtlinfo != null) {
                            for (eQTLResult r : eqtlinfo) {
                                // write index variant, clustered SNPs, and linked eQTL SNPs
                                // get list of genes as well?
                                // 	GWASID	Trait	IndexVariant	IndexVariantP	IndexVariantAlleles	IndexVariantEffect	nrEQTLSNPsLinked
                                // 	LinkedEQTLSNPs	LinkedEQTLGenes	LinkedEQTLHgncIDs	LinkedEQTLAlleles	LinkedEQTLZScores	GWASClusterSize	SNPsInGWASCluster
                                HashMap<String, Double> ldsnps = ldforGWASSNPs.get(indexVariant);
                                double rsq = ldsnps.get(r.snp);

                                // check direction?

                                String[] allelesGWAS = gwasSNPsPValues.get(gwasSNP1).alleles.split("/");
                                String assessedGWAS = "";
                                String gwasAlleleStr = gwasSNPsPValues.get(gwasSNP1).alleles;
                                direction dir = direction.UNKNOWN;
                                if (allelesGWAS.length > 1) {
                                    assessedGWAS = allelesGWAS[1];
                                    Boolean alleleflip = BaseAnnot.flipalleles(gwasSNPsPValues.get(gwasSNP1).alleles, assessedGWAS, r.alleles, r.assessed);
                                    if (alleleflip != null) {
                                        try {
                                            double beta = Double.parseDouble(gwasSNPsPValues.get(gwasSNP1).beta);
                                            double z = Double.parseDouble(r.z);
                                            if (alleleflip) {
                                                z *= -1;
                                            }
                                            if ((beta >= 0 && z >= 0) || (beta < 0 && z <= 0)) {
                                                dir = direction.SAME;
                                            } else {
                                                dir = direction.OPPOSITE;
                                            }
                                        } catch (NumberFormatException e) {

                                        }
                                    }
                                } else {

                                }

                                String outln = gwasID + "\t" + gwasIdToTrait.get(gwasID)
                                        + "\t" + indexVariant + "\t" + gwasSNPsPValues.get(gwasSNP1).p + "\t" + gwasAlleleStr + "\t" + gwasSNPsPValues.get(gwasSNP1).beta
                                        + "\t" + r.snp + "\t" + rsq + "\t" + r.gene + "\t" + r.hugo + "\t" + r.alleles + "\t" + r.assessed + "\t" + r.z + "\t" + r.p
                                        + "\t" + gwasCluster.size() + "\t" + clusterStr;
                                try {
                                    outPrune.writelnsynced(outln);
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                            }

                        }
                    }


                } // end: if snpobj!=null and notvisited(snp)
            } // end loop: iterate GWAS SNPs


            // write a summary
            double perc1 = (double) nrGWASClustersLinked / nrGWASClusters;
            if (Double.isNaN(perc1)) {
                perc1 = 0;
            }
            double perc2 = (double) eqtlSNPsLinkedToGWASSNPClusters.size() / finalNrEQTLSNPs;
            if (Double.isNaN(perc2)) {
                perc2 = 0;
            }

            String outln = gwasID + "\t" + gwasIdToTrait.get(gwasID) + "\t" + gwasSNPsPValues.size()
                    + "\t" + nrGWASClustersLinked + "\t" + nrGWASClusters + "\t" + perc1
                    + "\t" + eqtlSNPsLinkedToGWASSNPClusters.size() + "\t" + finalNrEQTLSNPs + "\t" + perc2;

            try {
                outSum.writelnsynced(outln);
            } catch (IOException e) {
                e.printStackTrace();
            }
            pb.iterateSynched();
        });
        pb.close();
        outAll.close();
        outSum.close();
        outPrune.close();

    }

    enum direction {
        SAME,
        OPPOSITE,
        UNKNOWN
    }

    class eQTLResult {
        String p;
        String snp;
        String gene;
        String hugo;
        String z;
        String alleles;
        String assessed;

    }

    private HashMap<String, ArrayList<eQTLResult>> loadEQTLInfo(String eqtlfile) throws IOException {
        TextFile tf = new TextFile(eqtlfile, TextFile.R);
        HashMap<String, ArrayList<eQTLResult>> output = new HashMap<>();
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            String snp = elems[1];
            String gene = elems[4];
            String z = elems[10];
            String alleles = elems[8];
            String assessed = elems[9];
            String hugo = elems[16];

            eQTLResult r = new eQTLResult();
            r.p = elems[0];
            r.alleles = alleles;
            r.assessed = assessed;
            r.z = z;
            r.snp = snp;
            r.gene = gene;
            r.hugo = hugo;
            ArrayList<eQTLResult> pairs = output.get(snp);
            if (pairs == null) {
                pairs = new ArrayList<>();
            }
            pairs.add(r);
            output.put(snp, pairs);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return output;

    }

    public void runPerGWAS(String eqtlset, String gwasset, String gwasAssocFile, String gwasListFile,
                           String trityperdir, String individualSubset,
                           int prunedistance, double prunethreshold,
                           double ldthreshold, String outputfile, boolean skipchr6) throws IOException {

        System.out.println("RSquared threshold: " + ldthreshold);
        System.out.println("Prune threshold: " + prunethreshold);
        System.out.println("Prune distance: " + prunedistance);

        HashSet<String> gwasSet = loadSet(gwasset);
        HashSet<String> eqtlSet = loadSet(eqtlset);
        System.out.println("gwasSet: " + gwasset + " has " + gwasSet.size() + " snps");
        System.out.println("eqtlSet: " + eqtlset + " has " + eqtlSet.size() + " snps");

        // load GWAS info
        Pair<HashMap<String, String>, HashMap<String, HashMap<String, GWASResult>>> gwasinfo = loadGWASInfo(gwasAssocFile, gwasListFile);
        HashMap<String, String> gwasIdToTrait = gwasinfo.getLeft();
        HashMap<String, HashMap<String, GWASResult>> gwasSNPs = gwasinfo.getRight();

        System.out.println("Loading data from: " + trityperdir);

        HashSet<String> allSNPsSet = new HashSet<>();
        allSNPsSet.addAll(gwasSet);
        allSNPsSet.addAll(eqtlSet);
        ArrayList<String> allSNPs = new ArrayList<>();
        allSNPs.addAll(allSNPsSet);
        System.out.println(allSNPs.size() + " SNPs total");
        Collections.sort(allSNPs);
        // preload genotypes

        TriTyperGenotypeData data = new TriTyperGenotypeData(trityperdir);
        if (individualSubset != null) {
            HashSet<String> indsToInclude = loadSet(individualSubset);
            System.out.println("Filtering individuals using: " + individualSubset + " having " + indsToInclude.size() + " individuals");
            String[] individuals = data.getIndividuals();
            Boolean[] include = new Boolean[individuals.length];
            for (int i = 0; i < individuals.length; i++) {
                if (indsToInclude.contains(individuals[i])) {
                    include[i] = true;
                } else {
                    include[i] = false;
                }
            }
            data.setIsIncluded(include);
        }

        HashMap<String, SNP> genotypes = loadSNPs(allSNPs, data, skipchr6);
        System.out.println(genotypes.size() + " snps found in total.");

        // for each GWAS, determine loaded SNPs, and their associated p-values
        ArrayList<String> allgwas = new ArrayList<>();
        allgwas.addAll(gwasSNPs.keySet());

        TextFile outAll = new TextFile(outputfile + "-all.txt", TextFile.W);
        String header = "GWASId\tTrait\tTraitP\tTraitSNP1\tMAF\tHWEP\tSNP2\tMAF\tHWEP\tRSQ\tDPR";
        outAll.writeln(header);

        TextFile outSum = new TextFile(outputfile + "-summary.txt", TextFile.W);
        String header2 = "GWASID\tTrait\tNrGWASSNPsTotal\tNrPrunedGWASSNPsLinked\tNrPrunedGWASSNPs\tPercGWASSNPsLinked\tnrEQTLSNPsLinked\tTotalEQTLSNPs\tPercEQTLSNPsLinked";
        outSum.writeln(header2);

        TextFile outPrune = new TextFile(outputfile + "-PrunedSNPsPerTrait.txt", TextFile.W);
        outPrune.writeln("Trait\tNrTotalSNPs\tNrPrunedSNPs\tPrunedSNPs");

        // count eQTL SNPs with genotypes
        int nrEQTLSNPs = 0;
        for (String snp : eqtlSet) {
            if (genotypes.containsKey(snp)) {
                nrEQTLSNPs++;
            }
        }

        System.out.println(allgwas.size() + " GWASes to test");
        int finalNrEQTLSNPs = nrEQTLSNPs;

//        IntStream.range(0, allgwas.size()).parallel().forEach(g -> {
        ProgressBar pb = new ProgressBar(allgwas.size(), "Calculating LD.");
        IntStream.range(0, allgwas.size()).forEach(g -> {
            String gwasID = allgwas.get(g);
            HashMap<String, GWASResult> snps = gwasSNPs.get(gwasID);

            // find snps that are within 1mb of each other
            ArrayList<String> allGWASSNPs = new ArrayList<>();
            allGWASSNPs.addAll(snps.keySet());
            Collections.sort(allGWASSNPs);
            HashSet<String> visited = new HashSet<>();

            HashSet<SNP> remainingSNPsForGWAS = new HashSet<>();

            AtomicInteger clusterctr = new AtomicInteger();
            // System.out.println(snps.size() + " SNPs before pruning for trait: " + gwasID + "\t" + gwasIdToTrait.get(gwasID));
            // prune snps for GWAS set
            for (int s = 0; s < allGWASSNPs.size(); s++) {
                String snp1 = allGWASSNPs.get(s);
                SNP snp1obj = genotypes.get(snp1);
                if (snp1obj != null && !visited.contains(snp1)) {
                    String[] elems = snp1.split(":");
                    byte chr = snp1obj.getChr();
                    Integer pos = snp1obj.getChrPos();// Integer.parseInt(elems[1]);

                    // determine which snps are within 1mb of each other
                    ArrayList<Pair<String, Double>> cluster = new ArrayList<>();
                    GWASResult gwasResult = snps.get(snp1);
                    Pair<String, Double> pair = new Pair<String, Double>(snp1, gwasResult.p, Pair.SORTBY.RIGHT);
                    cluster.add(pair);
                    visited.add(snp1);
                    for (int s2 = s + 1; s2 < allGWASSNPs.size(); s2++) {
                        String snp2 = allGWASSNPs.get(s2);
                        SNP snp2obj = genotypes.get(snp2);
                        if (snp2obj != null && !visited.contains(snp2)) {
                            String[] elems2 = allGWASSNPs.get(s2).split(":");
                            byte chr2 = snp2obj.getChr();
                            Integer pos2 = snp2obj.getChrPos();// Integer.parseInt(elems[1]);
                            if (chr == chr2) {
                                int distance = Math.abs(pos - pos2);
                                if (distance < prunedistance) {
                                    GWASResult gwasResult2 = snps.get(snp2);
                                    Pair<String, Double> pair2 = new Pair<String, Double>(snp2, gwasResult2.p, Pair.SORTBY.RIGHT);
                                    cluster.add(pair2);
                                    visited.add(snp2);
                                }
                            }
                        }
                    }


//                    // debug sorting
//                    for (int i = 0; (i < cluster.size() || i < 10); i++) {
//                        System.out.println(cluster.get(i));
//                    }

                    // prune cluster
                    ArrayList<SNP> remainingSNPs = new ArrayList<>();


                    if (cluster.size() == 1) {
                        String snpi = cluster.get(0).getLeft();
                        SNP snpiobj = genotypes.get(snpi);
                        remainingSNPs.add(snpiobj);
                    } else if (cluster.size() > 1) {
                        // sort cluster on P

                        Collections.sort(cluster);

                        DetermineLD ld = new DetermineLD();

                        HashSet<String> visitedPrune = new HashSet<>();
                        for (int i = 0; i < cluster.size(); i++) {
                            String snpi = cluster.get(i).getLeft();
                            if (!visitedPrune.contains(snpi)) {
                                visitedPrune.add(snpi);
                                SNP snpiobj = genotypes.get(snpi);
                                remainingSNPs.add(snpiobj);
                                for (int j = i + 1; j < cluster.size(); j++) {
                                    String snpj = cluster.get(j).getLeft();
                                    if (!visitedPrune.contains(snpj)) {
                                        // check LD
                                        SNP snpjobj = genotypes.get(snpj);
                                        Pair<Double, Double> ldvals = ld.getLD(snpiobj, snpjobj, data, 1, false);
                                        double rsq = ldvals.getRight();
//                                        System.out.println(snpi + "\t" + snpj + "\t" + rsq + "\t" + (rsq > prunethreshold));
                                        if (rsq > prunethreshold) {
                                            visitedPrune.add(snpj); // don't consider this snp as an independent one
                                        }
                                    }
                                }
                            }

                        }
//                        System.out.println(cluster.size() + " SNPs in cluster " + clusterctr.get() + ", " + remainingSNPs.size() + " remain after pruning for trait: " + gwasID + "\t" + gwasIdToTrait.get(gwasID));
//                        waitForEnter("Check this out!");

                    }

                    remainingSNPsForGWAS.addAll(remainingSNPs);
                    clusterctr.getAndIncrement();
                } // end: if snpobj!=null and notvisited(snp)
            } // end loop: iterate GWAS SNPs

            // calculate LD with eQTL SNPs for remaining SNPs;
            ArrayList<SNP> allRemainingSNPs = new ArrayList<>();
            allRemainingSNPs.addAll(remainingSNPsForGWAS);
            // System.out.println(allRemainingSNPs.size() + " SNPs remain after pruning for trait: " + gwasID + "\t" + gwasIdToTrait.get(gwasID));

            if (allRemainingSNPs.size() > 0) {
                ArrayList<String> allRemainingSNPStr = new ArrayList<>();
                for (SNP s : remainingSNPsForGWAS) {
                    allRemainingSNPStr.add(s.getName());
                }
                Collections.sort(allRemainingSNPStr);

                try {
                    outPrune.writelnsynced(gwasID + "\t" + allGWASSNPs.size() + "\t" + remainingSNPsForGWAS.size() + "\t" + Strings.concat(allRemainingSNPStr, Strings.semicolon));
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            AtomicInteger nrGWASSNPsLinked = new AtomicInteger();
            AtomicInteger nrEQTLSNPsLinked = new AtomicInteger();


            IntStream.range(0, allRemainingSNPs.size()).forEach(i -> {
                DetermineLD ld = new DetermineLD();
                SNP snp1obj = allRemainingSNPs.get(i);
                boolean linked = false;
                if (snp1obj != null) {
                    byte chr = snp1obj.getChr();
                    int pos = snp1obj.getChrPos();
                    for (String snp2 : eqtlSet) {
                        SNP snp2obj = genotypes.get(snp2);
                        if (snp2obj != null) {
                            byte chr2 = snp2obj.getChr();
                            int pos2 = snp2obj.getChrPos();
                            if (chr == chr2 && Math.abs(pos - pos2) < 5000000) {
                                Pair<Double, Double> ldvals = ld.getLD(snp1obj, snp2obj, data, 1, false);
                                double rsq = ldvals.getRight();
                                double dpr = ldvals.getLeft();
                                try {
                                    if (rsq > ldthreshold) {
                                        outAll.writelnsynced(gwasID + "\t" + gwasIdToTrait.get(gwasID) + "\t" + snps.get(snp1obj.getName()) + "\t" + snp1obj.getName() + "\t" + snp1obj.getMAF() + "\t" + snp1obj.getHWEP() + "\t" +
                                                snp2 + "\t" + snp2obj.getMAF() + "\t" + snp2obj.getHWEP() + "\t" + rsq + "\t" + dpr
                                        );
                                        linked = true;
                                        nrEQTLSNPsLinked.getAndIncrement();
                                    }
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                            }
                        }
                    }
                }
                if (linked) {
                    nrGWASSNPsLinked.getAndIncrement();
                }

            });

            // write a summary
            double perc1 = (double) nrGWASSNPsLinked.get() / allRemainingSNPs.size();
            if (Double.isNaN(perc1)) {
                perc1 = 0;
            }
            double perc2 = (double) nrEQTLSNPsLinked.get() / finalNrEQTLSNPs;
            if (Double.isNaN(perc2)) {
                perc2 = 0;
            }

            String outln = gwasID + "\t" + gwasIdToTrait.get(gwasID) + "\t" + snps.size()
                    + "\t" + nrGWASSNPsLinked + "\t" + allRemainingSNPs.size() + "\t" + perc1
                    + "\t" + nrEQTLSNPsLinked + "\t" + finalNrEQTLSNPs + "\t" + perc2;
            try {
                outSum.writelnsynced(outln);
            } catch (IOException e) {
                e.printStackTrace();
            }
            pb.iterateSynched();
        });
        pb.close();
        outAll.close();
        outSum.close();
        outPrune.close();

    }

    class GWASResult {
        double p;
        String snp;
        String alleles;
        String beta;
        String se;
    }

    private Pair<HashMap<String, String>, HashMap<String, HashMap<String, GWASResult>>> loadGWASInfo(String gwasAssocFile, String gwasListFile) throws IOException {
        System.out.println("Loading GWAS list: " + gwasListFile);
        TextFile tf = new TextFile(gwasListFile, TextFile.R);
        tf.readLine();
        HashMap<String, String> idToTrait = new HashMap<>();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String trait = elems[3];
            String id = elems[0];
            idToTrait.put(id, trait);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        System.out.println("Loading GWAS assoc: " + gwasAssocFile);
        tf = new TextFile(gwasAssocFile, TextFile.R);
        tf.readLine();

        HashMap<String, HashMap<String, GWASResult>> gwassnps = new HashMap<>();

        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            // ID      RsID    OtherAllele     EffectAllele    EffectAlleleFreq        Beta    SE      Pvalue
            String id = elems[0];
            HashMap<String, GWASResult> snps = gwassnps.get(id);
            if (snps == null) {
                snps = new HashMap<>();
            }


            String allele = elems[2] + "/" + elems[3];
            String beta = elems[5];
            String se = elems[6];

            String snp = elems[1];
            double p = 1;
            try {
                p = Double.parseDouble(elems[7]);
            } catch (NumberFormatException e) {

            }

            GWASResult r = new GWASResult();
            r.alleles = allele;
            r.beta = beta;
            r.se = se;
            r.snp = snp;
            r.p = p;

            snps.put(snp, r);
            gwassnps.put(id, snps);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return new Pair<>(idToTrait, gwassnps);
    }

    public void run(String set1file, String set2file, String trityperdir, String individualSubset, double threshold, String outputfile, boolean skipchr6) throws IOException {

        System.out.println("RSquared threshold: " + threshold);
        HashSet<String> set1 = loadSet(set1file);
        HashSet<String> set2 = loadSet(set2file);
        System.out.println("set1: " + set1file + " has " + set1.size() + " snps");
        System.out.println("set2: " + set2file + " has " + set2.size() + " snps");

        System.out.println("Loading data from: " + trityperdir);


        HashSet<String> foundSet1 = new HashSet<String>();
        HashSet<String> foundSet2 = new HashSet<String>();
        int ldcalculations = 0;


        HashSet<String> allSNPsSet = new HashSet<>();
        allSNPsSet.addAll(set1);
        allSNPsSet.addAll(set2);
        ArrayList<String> allSNPs = new ArrayList<>();
        allSNPs.addAll(allSNPsSet);
        System.out.println(allSNPs.size() + " SNPs total");
        Collections.sort(allSNPs);
        // preload genotypes

        TriTyperGenotypeData data = new TriTyperGenotypeData(trityperdir);
        if (individualSubset != null) {
            HashSet<String> indsToInclude = loadSet(individualSubset);
            System.out.println("Filtering individuals using: " + individualSubset + " having " + indsToInclude.size() + " individuals");
            String[] individuals = data.getIndividuals();
            Boolean[] include = new Boolean[individuals.length];
            for (int i = 0; i < individuals.length; i++) {
                if (indsToInclude.contains(individuals[i])) {
                    include[i] = true;
                } else {
                    include[i] = false;
                }
            }
            data.setIsIncluded(include);
        }

        HashMap<String, SNP> snps = loadSNPs(allSNPs, data, skipchr6);
        System.out.println(snps.size() + " snps found in total.");


        int ctr = 0;


        ArrayList<String> allset1SNPs = new ArrayList<>();
        allset1SNPs.addAll(set1);


        ProgressBar pb = new ProgressBar(allset1SNPs.size(), "Performing LD calculations");
        TextFile output = new TextFile(outputfile, TextFile.W);
        output.writeln("SNP1\tMAF1\tHWE-P1\tSNP2\tMAF2\tHWE-P1\tRsquared\tDprime");
        AtomicInteger identicalSNPs = new AtomicInteger();
        int[] ldbins = new int[11];
        IntStream.range(0, allset1SNPs.size()).parallel().forEach(i -> {
            DetermineLD ld = new DetermineLD();
            String snp1 = allset1SNPs.get(i);
            SNP snp1obj = snps.get(snp1);
            if (snp1obj != null) {

                String[] elems = snp1.split(":");
                byte chr = snp1obj.getChr();
                Integer pos = snp1obj.getChrPos(); // Integer.parseInt(elems[1]);
                for (String snp2 : set2) {
                    SNP snp2obj = snps.get(snp2);
                    if (snp2obj != null) {

                        String[] elems2 = snp2.split(":");
                        byte chr2 = snp2obj.getChr();
                        Integer pos2 = snp2obj.getChrPos(); // Integer.parseInt(elems[1]);
                        if (chr == chr2 && Math.abs(pos - pos2) < 5000000) {
                            Pair<Double, Double> ldvals = ld.getLD(snp1obj, snp2obj, data, 1, false);
                            double rsq = ldvals.getRight();
                            double dpr = ldvals.getLeft();

                            if (snp1.equals(snp2)) {
//                                System.out.println(snp1 + "\t" + rsq);

                                int bin = (int) Math.floor(rsq * 10);
                                ldbins[bin]++;
                                identicalSNPs.getAndIncrement();
                            }


                            try {
                                if (rsq > threshold) {
                                    output.writelnsynced(snp1 + "\t" + snp1obj.getMAF() + "\t" + snp1obj.getHWEP() + "\t" +
                                            snp2 + "\t" + snp2obj.getMAF() + "\t" + snp2obj.getHWEP() + "\t" + rsq + "\t" + dpr
                                    );
                                }
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }
                    }
                }

            }
            pb.iterateSynched();
        });

        System.out.println(identicalSNPs + " identical SNPs tested");
        for (int i = 0; i < ldbins.length; i++) {
            System.out.println(i + "\t" + ldbins[i]);
        }
        output.close();
        pb.close();


        System.out.println(foundSet1.size() + " SNPs found out of " + set1.size() + " defined in " + set1file);
        System.out.println(foundSet2.size() + " SNPs found out of " + set2.size() + " defined in " + set2file);
        System.out.println(ldcalculations + " calculations performed");

        // write SNPs that were not found
        writeNotFound(snps.keySet(), set1, outputfile + "-set1-found.txt", outputfile + "-set1-notfound.txt");
        writeNotFound(snps.keySet(), set2, outputfile + "-set2-found.txt", outputfile + "-set2-notfound.txt");

    }

    private HashMap<String, SNP> loadSNPs(ArrayList<String> allSNPs, TriTyperGenotypeData data, boolean skipchr6) throws IOException {


        // sort SNPs according to dataset
        TObjectIntHashMap<String> snpids = data.getSnpToSNPId();
        ArrayList<Pair<String, Integer>> ids = new ArrayList<>();
        for (String snp : allSNPs) {
            int id = snpids.get(snp);
            if (id >= 0) {
                ids.add(new Pair<String, Integer>(snp, id, Pair.SORTBY.RIGHT));
            }
        }
        Collections.sort(ids);

//        int min = Math.min(10, ids.size());
//        for (int i = 0; i < min; i++) {
//            System.out.println(ids.get(i));
//        }


        SNPLoader loader = data.createSNPLoader(1);
        HashMap<String, SNP> output = new HashMap<>();
        ProgressBar pb = new ProgressBar(allSNPs.size(), "Loading genotypes");
        int ctr = 0;
        for (Pair<String, Integer> snppair : ids) {
            Integer id1 = snppair.getRight();
            if (id1 >= 0) {
                SNP snpobj1 = data.getSNPObject(id1);
                if (!skipchr6 || snpobj1.getChr() != 6) {
                    loader.loadGenotypes(snpobj1);
                    if (loader.hasDosageInformation()) {
                        loader.loadDosage(snpobj1);
                    }
                    output.put(snppair.getLeft(), snpobj1);
                }
            }
            ctr++;
            pb.iterate();
        }
        loader.close();
        pb.close();
        return output;
    }


    private void writeNotFound(Set<String> foundSet1, HashSet<String> set1, String sfound, String snotfound) throws IOException {
        TextFile tf = new TextFile(snotfound, TextFile.W);
        TextFile tf2 = new TextFile(sfound, TextFile.W);
        for (String snp : set1) {
            if (!foundSet1.contains(snp)) {
                tf.writeln(snp);
            } else {
                tf2.writeln(snp);
            }
        }
        tf.close();
        tf2.close();
    }


    private HashSet<String> loadSet(String set1file) throws IOException {
        TextFile tf = new TextFile(set1file, TextFile.R);
        String ln = tf.readLine();
        HashSet<String> set = new HashSet<>();
        while (ln != null) {

            set.add(ln);


            ln = tf.readLine();
        }
        tf.close();
        return set;
    }
}
