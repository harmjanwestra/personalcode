package nl.harmjanwestra.playground.biogen.freeze2dot1.countsamples;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.*;

public class DNASamplesPerTissuePerDataset {


    public static void main(String[] args) {

        // note to self: this thing is not working properly
        String individuals = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-08-LinkFiles\\2020-02-08-individualfiles.txt";
        String popfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-08-LinkFiles\\2020-02-08-GenotypePopulationAssignments.txt";
        DNASamplesPerTissuePerDataset s = new DNASamplesPerTissuePerDataset();
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-17-SamplesBeforeQC\\DNACounts\\";

//        String sampleLinks = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiless\\2020-02-17-output\\links2-ExpressionSamplesLinked.txt";
        String sampleLinks = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-linkFiles2\\2020-02-17-output\\links2-ExpressionSamplesLinked.txt";
//        sampleLinks = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-02-17-output-exclusionsFixed\\links2-ExpressionSamplesLinked.txt";
        sampleLinks = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\links2-ExpressionSamplesLinked.txt";
        try {
//            s.countAvailableGenotypesPerDatasetAndPopulation(individuals, popfile, output);
            boolean countUniqueGenotypes = true;
//            s.countLinkedSamples(sampleLinks, popfile, countUniqueGenotypes);

            String tissueFolder = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\splitpertissue\\";

            s.countLinkedSamplesPerTissue(tissueFolder);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void countLinkedSamplesPerTissue(String tissueFolder) throws IOException {
        String[] tissues = new String[]{
                "amygdala",
                "basalganglia",
                "cerebellum",
                "cortex",
                "hippocampus",
                "hypothalamus",
                "spinalcord"
        };
        String[] populations = new String[]{
                "EUR", "AFR", "EAS", "SAS", "AMR"
        };

        String[] datasets = new String[]{
                "AMPAD-MAYO-V2",
                "AMPAD-MSBB-V2",
                "AMPAD-ROSMAP-V2",
                "Bipseq_1M",
                "Bipseq_h650",
                "BrainGVEX-V2",
                "Braineac",
                "CMC",
                "CMC_HBCC_set1",
                "CMC_HBCC_set2",
                "CMC_HBCC_set3",
                "GTEx",
                "GVEX",
                "LIBD_1M",
                "LIBD_5M",
                "LIBD_h650",
                "NABEC-H550",
                "NABEC-H610",
                "TargetALS",
                "UCLA_ASD",
                "ENA"
        };

        int[][][] cts = new int[datasets.length][tissues.length][3];
        for (int d = 0; d < datasets.length; d++) {
            for (int t = 0; t < tissues.length; t++) {
                for (int p = 0; p < populations.length; p++) {
                    int actualp = p;
                    if (p > 1) {
                        actualp = 2;
                    }
                    // amygdala.txt-dedup-gte.txt-EUR.txt-AMPAD-MAYO-V2.txt
                    // amygdala.txt-dedup-gte.txt-EAS.txt-AMPAD-MAYO-V2.txt

                    String file = tissueFolder + tissues[t] + ".txt-dedup-gte.txt-" + populations[p] + ".txt-" + datasets[d] + ".txt";
                    if (file.contains("cortex") && file.contains("EAS")) {

                    }
                    if (Gpio.exists(file)) {
                        TextFile tf = new TextFile(file, TextFile.R);
                        int ln = tf.countLines();
                        tf.close();
                        cts[d][t][actualp] += ln;
                    } else {
                        System.out.println("Cannot find: " + file);
                    }
                }
            }
        }

        for (int d = 0; d < datasets.length; d++) {

            String[] dsElems = datasets[d].split("\\\\");

            String ln = dsElems[dsElems.length - 1];
            for (int t = 0; t < tissues.length; t++) {
                for (int p = 0; p < cts[d][t].length; p++) {
                    ln += "\t" + cts[d][t][p];
                }
            }
            System.out.println(ln);
        }


    }

    private HashMap<String, String> loadPopulationInfo(String popfile) throws IOException {
        HashMap<String, String> popPerSample = new HashMap<>();
        TextFile tf = new TextFile(popfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            popPerSample.put(elems[0], elems[1]);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        return popPerSample;
    }

    public void countLinkedSamples(String sampleLinks, String popfile, boolean countUniqueGenotypes) throws IOException {
        HashMap<String, String> popPerSample = loadPopulationInfo(popfile);
        String[] populations = new String[]{
                "EUR", "AFR", "EAS", "SAS", "AMR"
        };

        HashMap<String, HashMap<String, Integer>> dsPopCtr = new HashMap<>();

        TextFile tf = new TextFile(sampleLinks, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        int linkCtr = 0;
        int linkCtrWithoutPopulation = 0;
        HashMap<String, HashSet<String>> seenGenotypesPerDataset = new HashMap<>();
        while (elems != null) {
            String gt = elems[1];
            String ds = elems[4];

            HashSet<String> seenIds = seenGenotypesPerDataset.get(ds);
            if (seenIds == null) {
                seenIds = new HashSet<>();
                seenGenotypesPerDataset.put(ds, seenIds);
            }

            if (!countUniqueGenotypes || !seenIds.contains(gt)) {
                seenIds.add(gt);
                String pop = popPerSample.get(gt);

                if (pop == null) {
                    System.out.println("Error: no population for sample: " + gt + " from " + ds);
                    linkCtrWithoutPopulation++;
                }
                HashMap<String, Integer> popctr = dsPopCtr.get(ds);
                if (popctr == null) {
                    popctr = new HashMap<>();
                }
                Integer ct = popctr.get(pop);
                if (ct == null) {
                    ct = 0;
                }
                ct++;
                popctr.put(pop, ct);
                dsPopCtr.put(ds, popctr);
            }
            linkCtr++;
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(linkCtr + " total lines.");
        System.out.println(linkCtrWithoutPopulation + " total lines without population.");
        ArrayList<String> keys = new ArrayList<String>();
        keys.addAll(dsPopCtr.keySet());
        Collections.sort(keys);
        for (String ds : keys) {
            String[] dsElems = ds.split("\\\\");
            String ln = dsElems[dsElems.length - 1];
            HashMap<String, Integer> popctr = dsPopCtr.get(ds);
            for (int p = 0; p < populations.length; p++) {
                Integer ct = popctr.get(populations[p]);
                if (ct == null) {
                    ct = 0;
                }
                ln += "\t" + ct;
            }
            System.out.println(ln);
        }

    }

    public void countAvailableGenotypesPerDatasetAndPopulation(String genotypedIndividualsFile, String popfile, String output) throws IOException {
        HashMap<String, String> popPerSample = loadPopulationInfo(popfile);

        String[] populations = new String[]{
                "EUR", "AFR", "EAS", "SAS", "AMR"
        };


        // load all genotyped individuals. Make sure to exclude those that have pihat > 0.125, but keep track of which samples were excluded
        // and also which genotyped individuals are present in multiple files
        HashMap<String, HashSet<String>> availableGenotypeIds = new HashMap<>();
        HashMap<String, HashSet<String>> allExclusions = new HashMap<>();
        TextFile tf2 = new TextFile(genotypedIndividualsFile, TextFile.R);
        String[] elems2 = tf2.readLineElems(TextFile.tab);
        TextFile gtOut = new TextFile(output + "AvailableGenotypeIDs.txt", TextFile.W);
        gtOut.writeln("GenotypeID\tDatasetFile");
        while (elems2 != null) {
            String sampleFile = elems2[0];
            Set<String> exclusions = null;
            if (elems2.length >= 2) {
                String exclusionfile = elems2[1];
                TextFile tf3 = new TextFile(exclusionfile, TextFile.R);
                exclusions = tf3.readAsSet(1, TextFile.tab);
                for (String s : exclusions) {
                    HashSet<String> exclusionset = allExclusions.get(s);
                    if (exclusionset == null) {
                        exclusionset = new HashSet<>();
                    }
                    exclusionset.add(exclusionfile);
                    allExclusions.put(s, exclusionset);
                }
                tf3.close();
            }

            TextFile tf3 = new TextFile(sampleFile, TextFile.R);
            String ln = tf3.readLine();
            HashSet<String> uniqueIds = new HashSet<String>();
            HashMap<String, Integer> popctr = new HashMap<String, Integer>();
            HashSet<String> uniqueIdsNotExcluded = new HashSet<String>();
            HashMap<String, Integer> popctrNotExcluded = new HashMap<String, Integer>();
            while (ln != null) {
                uniqueIds.add(ln);
                String pop = popPerSample.get(ln);

                if (pop == null) {
                    System.out.println("Error no population for sample: " + ln);
                }

                Integer ct = popctr.get(pop);
                if (ct == null) {
                    ct = 0;
                }
                ct++;
                popctr.put(pop, ct);
                if (exclusions == null || !exclusions.contains(ln)) {
                    HashSet<String> fileset = availableGenotypeIds.get(ln);
                    uniqueIdsNotExcluded.add(ln);

                    Integer ct2 = popctrNotExcluded.get(pop);
                    if (ct2 == null) {
                        ct2 = 0;
                    }
                    ct2++;
                    popctrNotExcluded.put(pop, ct2);
                    if (fileset == null) {
                        fileset = new HashSet<>();
                    }
                    fileset.add(sampleFile);
                    gtOut.writeln(ln + "\t" + sampleFile);
                    availableGenotypeIds.put(ln, fileset);
                }
                ln = tf3.readLine();
            }
            tf3.close();

            String popStr = "";
            String popStrNotExcluded = "";
            for (int p = 0; p < populations.length; p++) {
                String pop = populations[p];
                Integer ct = popctr.get(pop);
                if (ct == null) {
                    ct = 0;
                }
                popStr += "\t" + ct;

                Integer ctnotexcluded = popctrNotExcluded.get(pop);
                if (ctnotexcluded == null) {
                    ctnotexcluded = 0;
                }
                popStrNotExcluded += "\t" + ctnotexcluded;
            }

            System.out.println(sampleFile + "\t" + uniqueIds.size() + popStr + "\t" + uniqueIdsNotExcluded.size() + popStrNotExcluded);
            elems2 = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();
        gtOut.close();
    }
}
