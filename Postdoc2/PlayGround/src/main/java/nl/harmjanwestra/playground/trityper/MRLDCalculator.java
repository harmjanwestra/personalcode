package nl.harmjanwestra.playground.trityper;

import nl.harmjanwestra.playground.transeqtl.SupplementaryTables;
import org.apache.commons.cli.*;
import org.checkerframework.checker.units.qual.A;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class MRLDCalculator {

    public static void main(String[] args) {
        MRLDCalculator ld = new MRLDCalculator();


        Options options = new Options();
        options.addRequiredOption("e", "eqtlfile", true, "eQTL snp gene combos");
        options.addRequiredOption("m", "mrfile", true, "mr snp gene combos");
        options.addRequiredOption("t", "trityper", true, "TriTyper genotype reference");
        options.addRequiredOption("i", "individuals", true, "Individuals to select from TriTyper reference");
        options.addRequiredOption("o", "output", true, "Output dir/file");

        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine cmd = parser.parse(options, args);
            String eqtlfile = cmd.getOptionValue("e");
            String mrfile = cmd.getOptionValue("m");
            String trityperdir = cmd.getOptionValue("t");
            String individualsubset = cmd.getOptionValue("i");
            String output = cmd.getOptionValue("o");

            ld.run(eqtlfile, mrfile, trityperdir, individualsubset, output);

        } catch (ParseException e) {
            System.out.println(e.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("ant", options);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void run(String eqtlfile, String mrFile, String trityperdir, String individualSubset, String output) throws IOException {

        TextFile tf = new TextFile(mrFile, TextFile.R);
        HashMap<String, String> snpset = new HashMap<>();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length > 1) {
                String snp = elems[0];
                String gene = elems[1];
                snpset.put(snp, null);
            } else {
                System.out.println(elems.length + " < 2 in " + mrFile);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(snpset.size() + " SNPs in " + mrFile);

        // get appropriate SNP ids
        System.out.println("Reading ids from " + trityperdir + "SNPs.txt.gz");
        tf = new TextFile(trityperdir + "SNPs.txt.gz", TextFile.R);
        String ln = tf.readLine();
        int update = 0;
        int ctr = 0;
        while (ln != null) {
            elems = ln.split(":");
            String rs = elems[2];
            if (snpset.containsKey(rs)) {
                snpset.put(rs, ln);
                update++;
            }
            ctr++;
            if (ctr % 10000 == 0) {
                System.out.print(ctr + " snps read, " + update + " updated\r");
            }
            ln = tf.readLine();
        }
        tf.close();
        System.out.println();
        System.out.println(update + " SNPs updated");

        // load eqtls present in mr file
        HashSet<String> allSNPsHash = new HashSet<>();
        HashMap<String, HashSet<String>> mrEQTLs = new HashMap<>();
        tf = new TextFile(mrFile, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length > 1) {
                String snp = elems[0];
                String gene = elems[1];
                HashSet<String> snps = mrEQTLs.get(gene);

                String alias = snpset.get(snp);
                if (snps == null) {
                    snps = new HashSet<>();
                }
                if (alias != null) {
                    snps.add(alias);
                    allSNPsHash.add(alias);
                } else {
                    snps.add(snp);
                    allSNPsHash.add(snp);
                }

                allSNPsHash.add(snp);
                mrEQTLs.put(gene, snps);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(mrEQTLs.size() + " eqtls in mr file");

        // load actual eqtls
        System.out.println("Loading eqtls: " + eqtlfile);
        HashMap<String, HashSet<String>> eqtls = new HashMap<>();
        tf = new TextFile(eqtlfile, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[0];
            String gene = elems[1];

            HashSet<String> snpsforgene = eqtls.get(gene);
            if (snpsforgene == null) {
                snpsforgene = new HashSet<String>();
            }
            snpsforgene.add(snp);
            eqtls.put(gene, snpsforgene);
            allSNPsHash.add(snp);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        // preload genotypes
        System.out.println("Loading genotypes");
        LDCalculator ldcalc = new LDCalculator();
        TriTyperGenotypeData data = new TriTyperGenotypeData(trityperdir);
        if (individualSubset != null) {
            HashSet<String> indsToInclude = ldcalc.loadSet(individualSubset);
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

        ArrayList<String> allSNPs = new ArrayList<String>();
        allSNPs.addAll(allSNPsHash);
        HashMap<String, SNP> genotypes = ldcalc.loadSNPs(allSNPs, data, false);
        System.out.println(genotypes.size() + " SNP genotypes read");

        // perform LD calculations
        System.out.println("Performing LD calculations.");
        TextFile tfo = new TextFile(output, TextFile.W);
        tfo.writeln("MRSNP\trsid\tGene\tnrEQTLSNPsForGene\tEQTLSNP\tLD\tException");
        for (String gene : mrEQTLs.keySet()) {
            HashSet<String> mrSNPs = mrEQTLs.get(gene);
            HashSet<String> eqtlsnps = eqtls.get(gene);
            if (eqtlsnps == null) {
                for (String mrSNP : mrSNPs) {
                    String[] mrsnpelems = mrSNP.split(":");
                    String rs = mrsnpelems[0];
                    if (mrsnpelems.length > 2) {
                        rs = mrsnpelems[2];
                    }
                    tfo.writeln(mrSNP + "\t" + rs + "\t" + gene + "\tN/A\t0\t0\teSNPNotFound");
                }
            } else {
                for (String eqtlsnp : eqtlsnps) {
                    SNP esnpObj = genotypes.get(eqtlsnp);
                    if (esnpObj == null) {
                        for (String mrSNP : mrSNPs) {
                            String[] mrsnpelems = mrSNP.split(":");
                            String rs = mrsnpelems[0];
                            if (mrsnpelems.length > 2) {
                                rs = mrsnpelems[2];
                            }
                            tfo.writeln(mrSNP + "\t" + rs + "\t" + gene + "\t" + eqtlsnps.size() + "\t" + eqtlsnp + "\t0\teSNPNoGenotypes");
                        }
                    } else {
                        for (String mrSNP : mrSNPs) {
                            SNP mrsnpObj = genotypes.get(mrSNP);
                            String[] mrsnpelems = mrSNP.split(":");
                            String rs = mrsnpelems[0];
                            if (mrsnpelems.length > 2) {
                                rs = mrsnpelems[2];
                            }
                            if (mrsnpObj == null) {
                                tfo.writeln(mrSNP + "\t" + rs + "\t" + gene + "\t" + eqtlsnps.size() + "\t" + eqtlsnp + "\t0\tmrSNPNoGenotypes");
                            } else {
                                // calculate LD
                                DetermineLD ld = new DetermineLD();
                                Pair<Double, Double> ldvals = ld.getLD(esnpObj, mrsnpObj, data, 1, false);
                                double rsq = ldvals.getRight();
                                double dpr = ldvals.getLeft();
                                tfo.writeln(mrSNP + "\t" + rs + "\t" + gene + "\t" + eqtlsnps.size() + "\t" + eqtlsnp + "\t" + rsq + "\t-");
                            }
                        }
                    }
                }

            }
        }
        tfo.close();

    }

}
