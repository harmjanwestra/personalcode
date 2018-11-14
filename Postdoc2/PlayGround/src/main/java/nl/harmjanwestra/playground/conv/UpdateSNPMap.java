package nl.harmjanwestra.playground.conv;

import nl.harmjanwestra.utilities.enums.Chromosome;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class UpdateSNPMap {


    public static void main(String[] args) {
        UpdateSNPMap s = new UpdateSNPMap();

        if (args.length < 3) {
            System.out.println("useage: snpmap indir keepmissing");
        } else {
            try {
//			s.runSNPMap(args[0], args[1]);
                s.runSNPsSameBuildSingleDir(args[0], args[1], Boolean.parseBoolean(args[2]));
//			s.bimToBed(args[0], args[1], args[2]);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }


    public void runSNPsSameBuildSingleDir(String snpmap, String indir, boolean keepmissing) throws IOException {


        HashMap<Chromosome, HashMap<Integer, String>> posToRs = new HashMap<>();
        {
            TextFile tf = new TextFile(indir + "/SNPMappings.txt.gz", TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            HashMap<Chromosome, HashSet<Integer>> positionsToLoad = new HashMap<Chromosome, HashSet<Integer>>();
            int ctr = 0;
            while (elems != null) {

                Chromosome chr = Chromosome.parseChr(elems[0]);
                Integer pos = Integer.parseInt(elems[1]);

                HashSet<Integer> set = positionsToLoad.get(chr);
                if (set == null) {
                    set = new HashSet<>();
                }
                set.add(pos);

                positionsToLoad.put(chr, set);

                elems = tf.readLineElems(TextFile.tab);
                ctr++;
                if (ctr % 10000 == 0) {
                    System.out.print(ctr + " postions parsed\r");
                }

            }
            tf.close();

            for (Chromosome key : positionsToLoad.keySet()) {
                System.out.println(key + "\t" + positionsToLoad.get(key).size() + " positions loaded");
            }

            TextFile tf2 = new TextFile(snpmap, TextFile.R);
            System.out.println("Parsing snpmap " + snpmap);
            elems = tf2.readLineElems(TextFile.tab);
            ctr = 0;
            while (elems != null) {
                Chromosome chr = Chromosome.parseChr(elems[0]);
                Integer pos = Integer.parseInt(elems[1]);
                String rs = elems[2];

                HashSet<Integer> set = positionsToLoad.get(chr);
                if (set != null && set.contains(pos)) {
                    HashMap<Integer, String> set2 = posToRs.get(chr);
                    if (set2 == null) {
                        set2 = new HashMap<>();

                    }
                    set2.put(pos, rs);
                    posToRs.put(chr, set2);
                }
                elems = tf2.readLineElems(TextFile.tab);
                ctr++;
                if (ctr % 10000 == 0) {
                    System.out.print(ctr + " rsids parsed\r");
                }
            }
            tf2.close();

            for (Chromosome key : posToRs.keySet()) {
                System.out.println(key + "\t" + posToRs.get(key).size() + " rsids loaded");
            }
        }


        System.out.println("parsing: " + indir + "/SNPs.txt.gz");
        TextFile tf = new TextFile(indir + "/SNPs.txt.gz", TextFile.R);
        TextFile tfosnps = new TextFile(indir + "/SNPs-update.txt.gz", TextFile.W);
        TextFile tfosnpmap = new TextFile(indir + "/SNPMappings-update.txt.gz", TextFile.W);

        String ln = tf.readLine();
        int ctr = 0;
        while (ln != null) {

            String[] elems = ln.split(":");
            Chromosome chr = Chromosome.parseChr(elems[0]);
            String[] poselems = elems[1].split("_");
            Integer pos = Integer.parseInt(poselems[0]);

            HashMap<Integer, String> set = posToRs.get(chr);
            String rs = set.get(pos);
            if (rs == null) {

                tfosnps.writeln(ln);
                tfosnpmap.writeln(chr.getNumber() + "\t" + pos + "\t" + ln);
            } else {
                if (poselems.length > 1) {
                    tfosnps.writeln(rs + "_" + poselems[1]);
                    tfosnpmap.writeln(chr.getNumber() + "\t" + pos + "\t" + rs + "_" + poselems[1]);
                } else {
                    tfosnps.writeln(rs);
                    tfosnpmap.writeln(chr.getNumber() + "\t" + pos + "\t" + rs);
                }
            }

            ln = tf.readLine();
            ctr++;
            if (ctr % 10000 == 0) {
                System.out.print(ctr + " rsids parsed\r");
            }

        }

        tf.close();
        tfosnpmap.close();
        tfosnps.close();
    }

    public void runSNPs(String snpmap, String indir, boolean keepmissing) throws IOException {

        HashMap<String, String> posToRS = new HashMap<>();
        HashMap<String, String> rsToPos = new HashMap<>();

        System.out.println("Parsing: " + snpmap);
        TextFile tf = new TextFile(snpmap, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            Chromosome chr = Chromosome.parseChr(elems[0]);
            String pos = elems[1];
            String rsid = elems[2];
            posToRS.put(chr.getNumber() + ":" + pos, rsid);
            rsToPos.put(rsid, chr.getNumber() + ":" + pos);

            if (ctr % 10000 == 0) {
                System.out.print(ctr + " snps loaded.\r");
            }
            ctr++;
            elems = tf.readLineElems(TextFile.tab);
        }
        System.out.println();
        tf.close();

        System.out.println(posToRS.size() + " positions loaded");
        System.out.println(rsToPos.size() + " rsids loaded");

        for (int i = 1; i < 23; i++) {

            TextFile tfs = new TextFile(indir + "/chr" + i + "/SNPs.txt.gz", TextFile.R);
            System.out.println("Processing: " + tfs.getFileName());
            TextFile tfso = new TextFile(indir + "/chr" + i + "/SNPs-update.txt.gz", TextFile.W);
            TextFile tfsm = new TextFile(indir + "/chr" + i + "/SNPMappings-update.txt.gz", TextFile.W);
            int updated = 0;
            int removed = 0;
            int total = 0;
            String ln = tfs.readLine();
            while (ln != null) {
                if (ln.startsWith("rs")) {
                    if (!rsToPos.containsKey(ln)) {
                        tfso.writeln(ln + "-retired");
                        tfsm.writeln("0\t0\t" + ln + "-retired");
                        removed++;
                    } else {
                        String pos = rsToPos.get(ln);
                        String[] poselems = pos.split(":");
                        tfso.writeln(ln);
                        tfsm.writeln(poselems[0] + "\t" + poselems[1] + "\t" + ln);
                        updated++;
                    }
                } else {
                    String[] snpelems = ln.split(":");
                    if (snpelems.length == 2) {
                        Chromosome chr = Chromosome.parseChr(snpelems[0]);
                        String rs = posToRS.get(chr.getNumber() + ":" + snpelems[1]);

                        if (rs == null) {
                            if (keepmissing) {
                                tfso.writeln(chr.getNumber() + ":" + snpelems[1]);
                                tfsm.writeln(chr.getNumber() + "\t" + snpelems[1] + "\t" + chr.getNumber() + ":" + snpelems[1]);
                            } else {
                                tfso.writeln(ln + "-retired");
                                tfsm.writeln("0\t0\t" + ln + "-retired");
                            }
                            removed++;
                        } else {
                            tfso.writeln(rs);
                            tfsm.writeln(chr.getNumber() + "\t" + snpelems[1] + "\t" + rs);
                            updated++;
                        }
                    } else {
                        tfso.writeln(ln + "-retired");
                        tfsm.writeln("0\t0\t" + ln + "-retired");
                        removed++;
                    }

                }
                total++;
                ln = tfs.readLine();
            }
            tfso.close();
            tfsm.close();
            tfs.close();

            System.out.println(total + " written. " + removed + " removed. " + updated + " updated.");

        }
    }

    public void runSNPMap(String snpmap, String indir) throws IOException {

        HashMap<String, String> posToRS = new HashMap<>();
        HashMap<String, String> rsToPos = new HashMap<>();

        System.out.println("Parsing: " + snpmap);
        TextFile tf = new TextFile(snpmap, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            Chromosome chr = Chromosome.parseChr(elems[0]);
            String pos = elems[1];
            String rsid = elems[2];
            posToRS.put(chr.getNumber() + ":" + pos, rsid);
            rsToPos.put(rsid, chr.getNumber() + ":" + pos);

            if (ctr % 10000 == 0) {
                System.out.print(ctr + " snps loaded.\r");
            }
            ctr++;
            elems = tf.readLineElems(TextFile.tab);
        }
        System.out.println();
        tf.close();

        System.out.println(posToRS.size() + " positions loaded");
        System.out.println(rsToPos.size() + " rsids loaded");

        for (int i = 1; i < 23; i++) {

            TextFile tfs = new TextFile(indir + "/chr" + i + "/SNPMappings.txt.gz", TextFile.R);
            System.out.println("Processing: " + tfs.getFileName());
            TextFile tfso = new TextFile(indir + "/chr" + i + "/SNPs-update.txt.gz", TextFile.W);
            TextFile tfsm = new TextFile(indir + "/chr" + i + "/SNPMappings-update.txt.gz", TextFile.W);
            int updated = 0;
            int removed = 0;
            int total = 0;
            String ln = tfs.readLine();
            while (ln != null) {
                String[] lnelems = ln.split("\t");
                ln = lnelems[2];
                if (ln.startsWith("rs")) {
                    if (!rsToPos.containsKey(ln)) {
                        tfso.writeln(ln + "-retired");
                        tfsm.writeln("0\t0\t" + ln + "-retired");
                        removed++;
                    }
                } else {
                    String[] snpelems = ln.split(":");
                    if (snpelems.length == 2) {
                        Chromosome chr = Chromosome.parseChr(snpelems[0]);
                        String rs = posToRS.get(chr.getNumber() + ":" + snpelems[1]);
                        if (rs == null) {
                            tfso.writeln(ln + "-retired");
                            tfsm.writeln("0\t0\t" + ln + "-retired");
                            removed++;
                        } else {
                            tfso.writeln(rs);
                            tfsm.writeln(chr.getNumber() + "\t" + snpelems[1] + "\t" + rs);
                            updated++;
                        }
                    } else {
                        tfso.writeln(ln + "-retired");
                        tfsm.writeln("0\t0\t" + ln + "-retired");
                        removed++;
                    }

                }
                total++;
                ln = tfs.readLine();
            }
            tfso.close();
            tfsm.close();
            tfs.close();

            System.out.println(total + " written. " + removed + " removed. " + updated + " updated.");

        }
    }

    public void convertRsToPos(String snpfile, String annotfile, String output) {

    }

    public void bimToBed(String snpfile, String bim, String bed) throws IOException {
        TextFile tf = new TextFile(snpfile, TextFile.R);
        ArrayList<String> snps = tf.readAsArrayList();
        tf.close();

        TextFile tf2 = new TextFile(bim, TextFile.R);
        String[] elems = tf2.readLineElems(TextFile.tab);
        HashMap<String, String> rsToPos = new HashMap<>();

        while (elems != null) {
            Chromosome chr = Chromosome.parseChr(elems[0]);
            String rs = elems[1];
            String pos = elems[3];
            rsToPos.put(rs, "chr" + chr.getNumber() + "\t" + pos + "\t" + pos);
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();

        int found = 0;
        int notfound = 0;
        TextFile bedout = new TextFile(bed, TextFile.W);
        for (String snp : snps) {
            String[] split = snp.split(":");
            if (split.length == 2) {
                Chromosome chr = Chromosome.parseChr(split[0]);
                bedout.writeln("chr" + chr.getNumber() + "\t" + split[1] + "\t" + split[1] + "\t" + snp);
            } else {
                String outstr = rsToPos.get(snp);
                if (outstr != null) {
                    bedout.writeln(outstr + "\t" + snp);
                    found++;
                } else {
                    notfound++;
                    bedout.writeln("chr0\t0\t0\t" + split);
                }
            }
        }
        bedout.writeln();
        System.out.println(found + " snps found in bim. " + notfound + " not found");
    }

}
