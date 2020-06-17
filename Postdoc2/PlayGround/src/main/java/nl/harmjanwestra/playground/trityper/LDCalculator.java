package nl.harmjanwestra.playground.trityper;

import org.checkerframework.checker.units.qual.A;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

import java.io.Console;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class LDCalculator {

    public static void main(String[] args) {
        LDCalculator ld = new LDCalculator();


        if (args.length < 4) {
            System.out.println("Usage: compareall gwasset.txt eqtlset.txt trityperdir output [individualsubset] [ldthreshold:0.8]\n" +
                    "Or\t" +
                    "Usage: traitprune gwasset.txt eqtlset.txt gwasassocfile.txt.gz gwaslistfile.txt.gz trityperdir output [individualsubset] [ldthreshold:0.8] [ldprunethreshold:0.2] [ldprunedistance:1E6]");
        } else {

            String mode = args[0];
            try {
                if (mode.equals("compareall")) {
                    String setfile1 = args[1];
                    String setfile2 = args[2];
                    String trityper = args[3];
                    String output = args[4];

                    double threshold = 0.8;

                    String individualsubset = null;
                    if (args.length > 5) {
                        individualsubset = args[6];
                    }
                    if (args.length > 6) {
                        threshold = Double.parseDouble(args[6]);
                    }
                    ld.run(setfile1, setfile2, trityper, individualsubset, threshold, output);
                } else if (mode.equals("traitprune")) {
                    String gwasset = args[1];
                    String eqtlset = args[2];
                    String gwasassocfile = args[3];
                    String gwaslistfile = args[4];
                    String trityperdir = args[5];
                    String output = args[6];

                    double ldthreshold = 0.8;
                    double prunethreshold = 0.2;
                    int prunedistance = 1000000;

                    String individualsubset = null;
                    if (args.length > 7) {
                        individualsubset = args[7];
                    }
                    if (args.length > 8) {
                        ldthreshold = Double.parseDouble(args[8]);
                    }
                    if (args.length > 9) {
                        prunethreshold = Double.parseDouble(args[9]);
                    }
                    if (args.length > 10) {
                        prunedistance = Integer.parseInt(args[10]);
                    }
                    ld.runPerGWAS(eqtlset, gwasset, gwasassocfile, gwaslistfile, trityperdir, individualsubset, prunedistance, prunethreshold, ldthreshold, output);
                } else {
                    System.out.println("Select traitprune or compareall");
                }
            } catch (IOException e) {
                e.printStackTrace();
            }


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

    public void runPerGWAS(String eqtlset, String gwasset, String gwasAssocFile, String gwasListFile,
                           String trityperdir, String individualSubset,
                           int prunedistance, double prunethreshold,
                           double ldthreshold, String outputfile) throws IOException {

        System.out.println("RSquared threshold: " + ldthreshold);
        System.out.println("Prune threshold: " + prunethreshold);
        System.out.println("Prune distance: " + prunedistance);

        HashSet<String> set1 = loadSet(gwasset);
        HashSet<String> set2 = loadSet(eqtlset);
        System.out.println("set1: " + gwasset + " has " + set1.size() + " snps");
        System.out.println("set2: " + eqtlset + " has " + set2.size() + " snps");

        // load GWAS info
        Pair<HashMap<String, String>, HashMap<String, HashMap<String, Double>>> gwasinfo = loadGWASInfo(gwasAssocFile, gwasListFile);
        HashMap<String, String> gwasIdToTrait = gwasinfo.getLeft();
        HashMap<String, HashMap<String, Double>> gwasSNPs = gwasinfo.getRight();

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

        HashMap<String, SNP> genotypes = loadSNPs(allSNPs, data);
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

        // count eQTL SNPs with genotypes
        int nrEQTLSNPs = 0;
        for (String snp : set2) {
            if (genotypes.containsKey(snp)) {
                nrEQTLSNPs++;
            }
        }

        System.out.println(allgwas.size() + " GWASes to test");
        int finalNrEQTLSNPs = nrEQTLSNPs;

//        IntStream.range(0, allgwas.size()).parallel().forEach(g -> {
        IntStream.range(0, allgwas.size()).forEach(g -> {
            String gwasID = allgwas.get(g);
            HashMap<String, Double> snps = gwasSNPs.get(gwasID);

            // find snps that are within 1mb of each other
            ArrayList<String> allGWASSNPs = new ArrayList<>();
            allGWASSNPs.addAll(snps.keySet());
            Collections.sort(allGWASSNPs);
            HashSet<String> visited = new HashSet<>();

            HashSet<SNP> remainingSNPsForGWAS = new HashSet<>();

            AtomicInteger clusterctr = new AtomicInteger();
            System.out.println(snps.size() + " SNPs before pruning for trait: " + gwasID + "\t" + gwasIdToTrait.get(gwasID));
            // prune snps for GWAS set
            for (int s = 0; s < allGWASSNPs.size(); s++) {
                String snp1 = allGWASSNPs.get(s);
                SNP snp1obj = genotypes.get(snp1);
                if (snp1obj != null && !visited.contains(snp1)) {
                    String[] elems = snp1.split(":");
                    String chr = elems[0];
                    Integer pos = Integer.parseInt(elems[1]);

                    // determine which snps are within 1mb of each other
                    ArrayList<Pair<String, Double>> cluster = new ArrayList<>();
                    Double p = snps.get(snp1);
                    Pair<String, Double> pair = new Pair<String, Double>(snp1, p, Pair.SORTBY.RIGHT);
                    cluster.add(pair);
                    visited.add(snp1);
                    for (int s2 = s + 1; s2 < allGWASSNPs.size(); s2++) {
                        String snp2 = allGWASSNPs.get(s2);
                        SNP snp2obj = genotypes.get(snp2);
                        if (snp2obj != null && !visited.contains(snp2)) {
                            String[] elems2 = allGWASSNPs.get(s2).split(":");
                            String chr2 = elems2[0];
                            Integer pos2 = Integer.parseInt(elems2[1]);
                            if (chr.equals(chr2)) {
                                int distance = Math.abs(pos - pos2);
                                if (distance < prunedistance) {
                                    Double p2 = snps.get(snp2);
                                    Pair<String, Double> pair2 = new Pair<String, Double>(snp2, p2, Pair.SORTBY.RIGHT);
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
            System.out.println(allRemainingSNPs.size() + " SNPs remain after pruning for trait: " + gwasID + "\t" + gwasIdToTrait.get(gwasID));

            AtomicInteger nrGWASSNPsLinked = new AtomicInteger();
            AtomicInteger nrEQTLSNPsLinked = new AtomicInteger();


            IntStream.range(0, allRemainingSNPs.size()).forEach(i -> {
                DetermineLD ld = new DetermineLD();
                SNP snp1obj = allRemainingSNPs.get(i);
                boolean linked = false;
                if (snp1obj != null) {
                    byte chr = snp1obj.getChr();
                    int pos = snp1obj.getChrPos();
                    for (String snp2 : set2) {
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

        });
        outAll.close();
        outSum.close();


    }


    private Pair<HashMap<String, String>, HashMap<String, HashMap<String, Double>>> loadGWASInfo(String gwasAssocFile, String gwasListFile) throws IOException {
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

        HashMap<String, HashMap<String, Double>> gwassnps = new HashMap<>();

        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String id = elems[0];
            HashMap<String, Double> snps = gwassnps.get(id);
            if (snps == null) {
                snps = new HashMap<>();
            }

            String snp = elems[1];
            double p = 1;
            try {
                p = Double.parseDouble(elems[7]);
            } catch (NumberFormatException e) {

            }
            snps.put(snp, p);
            gwassnps.put(id, snps);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return new Pair<>(idToTrait, gwassnps);
    }

    public void run(String set1file, String set2file, String trityperdir, String individualSubset, double threshold, String outputfile) throws IOException {

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

        HashMap<String, SNP> snps = loadSNPs(allSNPs, data);
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
                String chr = elems[0];
                Integer pos = Integer.parseInt(elems[1]);
                for (String snp2 : set2) {
                    SNP snp2obj = snps.get(snp2);
                    if (snp2obj != null) {

                        String[] elems2 = snp2.split(":");
                        String chr2 = elems2[0];
                        Integer pos2 = Integer.parseInt(elems2[1]);
                        if (chr.equals(chr2) && Math.abs(pos - pos2) < 5000000) {
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

    private HashMap<String, SNP> loadSNPs(ArrayList<String> allSNPs, TriTyperGenotypeData data) throws IOException {


        SNPLoader loader = data.createSNPLoader(1);
        HashMap<String, SNP> output = new HashMap<>();
        ProgressBar pb = new ProgressBar(allSNPs.size(), "Loading genotypes");
        int ctr = 0;
        for (String snp : allSNPs) {
            Integer id1 = data.getSnpToSNPId().get(snp);
            if (id1 >= 0) {

                SNP snpobj1 = data.getSNPObject(id1);
                loader.loadGenotypes(snpobj1);
                if (loader.hasDosageInformation()) {
                    loader.loadDosage(snpobj1);
                }
                output.put(snp, snpobj1);
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
