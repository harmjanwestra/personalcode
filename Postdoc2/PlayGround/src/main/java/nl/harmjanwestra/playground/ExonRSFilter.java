package nl.harmjanwestra.playground;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.sets.FeatureSets;

import java.io.IOException;
import java.util.*;

public class ExonRSFilter {

    public static void main(String[] args) {

        String file = "D:\\Sync\\SyncThing\\Data\\Ref\\dbsnp\\common_all_20180418.vcf.gz";
        String outf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ExonRSIds\\gencode.v24.chr_patch_hapl_scaff.annotation.exon.interval_list_rsids.txt.gz";

        String interval = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ExonRSIds\\gencode.v24.chr_patch_hapl_scaff.annotation.exon.interval_list";

        try {
            TextFile tf = new TextFile(interval, TextFile.R);
            HashMap<Chromosome, TreeSet<Feature>> intervals = new HashMap<>();
            String ln = tf.readLine();


            while (ln != null) {

                if (!ln.startsWith("@")) {
                    String[] elems = ln.split("\t");
                    Chromosome chr = Chromosome.parseChr(elems[0]);
                    if (chr.isAutosome()) {
                        int start = Integer.parseInt(elems[1]);
                        int stop = Integer.parseInt(elems[2]);
                        Feature f = new Feature(chr, start, stop);

                        TreeSet<Feature> set = intervals.get(chr);
                        if (set == null) {
                            set = new TreeSet<>(new FeatureComparator(true));
                        }

                        set.add(f);
                        intervals.put(chr, set);
                    }
                }

                ln = tf.readLine();
            }
            tf.close();

            for (Chromosome chr : Chromosome.values()) {
                TreeSet<Feature> set = intervals.get(chr);
                if (set != null) {
                    System.out.println(chr + "\t" + set.size());

                }
            }

            TextFile out = new TextFile(outf, TextFile.W);
            TextFile vcfin = new TextFile(file, TextFile.R);

            ln = vcfin.readLine();
            TreeSet<Feature> setchr = null;
            Chromosome prevch = null;
            int lnctr = 0;
            while (ln != null) {

                if (!ln.startsWith("#")) {
                    String[] elems = ln.split("\t");
                    Chromosome chr = Chromosome.parseChr(elems[0]);
                    if (chr.isAutosome()) {
                        int start = Integer.parseInt(elems[1]);
                        int stop = start + 1;
                        Feature fstart = new Feature(chr, start - 100000, start);
                        Feature f = new Feature(chr, start, stop);
                        Feature fstop = new Feature(chr, start, start + 100000);

                        if (setchr == null) {
                            setchr = intervals.get(chr);
                            prevch = chr;
                            System.out.println("Parsing: " + chr);
                        } else if (!prevch.equals(chr)) {

                            setchr = intervals.get(chr);
                            prevch = chr;
                            System.out.println("Parsing: " + chr);
                        }
                        if (setchr != null) {
                            boolean overlap = false;

                            NavigableSet<Feature> subset = setchr.subSet(fstart, true, fstop, true);

                            for (Feature f2 : subset) {
                                if (f.overlaps(f2)) {
                                    overlap = true;
                                    break;
                                }
                            }
                            if (overlap) {
                                out.writeln(Strings.concat(elems, Strings.tab, 0, 4));
                            }
                        }

                    }
                }
                ln = vcfin.readLine();
                lnctr++;
                if (lnctr % 10000 == 0) {
                    if (prevch != null) {
                        System.out.print("\r" + lnctr + " lines, \tcurrent chr: " + prevch);
                    } else {
                        System.out.print("\r" + lnctr + " lines");
                    }
                }
            }
            vcfin.close();
            out.close();


        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
