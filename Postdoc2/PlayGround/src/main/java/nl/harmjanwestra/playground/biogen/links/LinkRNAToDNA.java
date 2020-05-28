package nl.harmjanwestra.playground.biogen.links;


import org.apache.poi.ss.formula.functions.T;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class LinkRNAToDNA {

    public static void main(String[] args) {
        LinkRNAToDNA l = new LinkRNAToDNA();

        String origlink = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\2019-04-13-RNASeqIdToGenotypeID.txt";
        String freeze1links = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\Freeze1Links.txt";
        String psychencodelinks = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\PsychEncodeIndividualToGT.txt";

        String[] origlinkPhenotypeFile = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-linkFiles\\2020-02-17-OriginalGenotypeExpressionLinks-withSamplesInExpressionData.txt",

        };
        String individualfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-08-LinkFiles\\2020-02-08-individualfiles.txt";


        individualfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-08-LinkFiles\\2020-02-08-individualfiles-woCMCHBCC.txt";

        String tissuemapfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-Tissues\\plots\\2020-02-06-eigenvectors-tissueClassification_residualSampleAssignmentAll.txt";
        String popfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-08-LinkFiles\\2020-02-08-GenotypePopulationAssignments.txt";
        String rnaseqoutliers = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-08-LinkFiles\\2020-02-11-OutlierSamplesAndThoseWithoutCovariates.txt";

        String cmcLinks = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\CMC_HumansampleIDkey_GTE.txt";

        String outdir = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\";
        String output = outdir + "links-";
        String output2 = outdir + "links2-";

        boolean excludepihat = true;

        try {
            Gpio.createDir(outdir);
            String cmcIndToDNA = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\CMC\\CMC_individualIdToGenotype.txt";
            String cmcIndToRNA = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\CMC\\CMC_individualIdToRNASeq.txt";
            String cmcDNAToRNA = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\CMC\\CMC_GenotypeIdToRNASeq-throughIndIds.txt";
            l.linkCMCSamples(cmcIndToDNA, cmcIndToRNA, cmcDNAToRNA);
            String newCMCLinkFile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\CMC_HumansampleIDkey_GTE-actualIDs.txt";
            l.cmcLinkCheck(cmcDNAToRNA, individualfile, newCMCLinkFile);
//            System.exit(0);
            String listOfExpressionSamples = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-tissues\\2020-02-17-listOfRNASeqSamples.txt";
            l.run2(listOfExpressionSamples, individualfile, origlinkPhenotypeFile, freeze1links, psychencodelinks, newCMCLinkFile, excludepihat, output2);

//			String newList = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-LinkFiles\\2020-05-25-output-exclusionsFixed-popcheck\\links2-ExpressionSamplesLinked.txt";
//			String oldList = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-linkFiles2\\2020-02-17-output\\links2-ExpressionSamplesLinked.txt";
//			compareLists(newList, oldList);
//			System.exit(-1);
        } catch (IOException e) {
            e.printStackTrace();
        }
//		System.exit(0);

//		try {
//			l.run(individualfile, origlink, freeze1links, psychencodelinks, rnaseqoutliers, output);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}

        boolean preferresequencedsamples = true;
        String file = outdir + "links2-ExpressionSamplesLinked.txt";
        String defupfileout = outdir + "links2-ExpressionSamplesLinked-dedup.txt";
        try {
            // remove duplicate genotype samples
            l.dedupRNA(file, defupfileout);
            String pathstrippeddedupfile = outdir + "links2-ExpressionSamplesLinked-dedup-pathstripped.txt";
            l.stripPaths(defupfileout, pathstrippeddedupfile);

//			System.exit(-1);


            String splittissuefolder = outdir + "splitpertissue\\";

            Gpio.createDir(splittissuefolder);


            String[] lsof = Gpio.getListOfFiles(splittissuefolder);
            for (String s : lsof) {
                Gpio.delete(s);
            }
            l.splitBasedOnTissueType(pathstrippeddedupfile, tissuemapfile, splittissuefolder);

            String[] tissues = new String[]{
                    "amygdala",
                    "hippocampus",
                    "basalganglia",
                    "hypothalamus",
                    "cerebellum",
                    "cortex",
                    "spinalcord"
            };
            String[] populations = new String[]{
                    "EUR", "AFR", "EAS", "SAS", "AMR"
            };
            for (String tissue : tissues) {
                String[] files = new String[]{outdir + "splitpertissue\\" + tissue + ".txt"};

                for (String linkfile : files) {
                    l.dedupDNA(linkfile, preferresequencedsamples);
                }

                files = new String[]{outdir + "splitpertissue\\" + tissue + ".txt-dedup-gte.txt"};
                for (String linkfile : files) {
                    l.splitPerPopulation(linkfile, popfile);
                }

                for (String pop : populations) {
                    String popgtefile = outdir + "splitpertissue\\" + tissue + ".txt-dedup-gte.txt-" + pop + ".txt";
                    if (Gpio.exists(popgtefile)) {
                        files = new String[]{popgtefile};
                    }
                    for (String linkfile : files) {
                        l.splitPerDataset(linkfile, pathstrippeddedupfile);
                    }
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
//		String[] inds = new String[]{
//
//		}

//		try {
//			l.run(inds, origlink, output);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
    }

    private void linkCMCSamples(String cmcIndToDNA, String cmcIndToRNA, String cmcDNAToRNA) throws IOException {
        HashMap<String, HashSet<String>> indToDNA = new HashMap<>();
        TextFile tf = new TextFile(cmcIndToDNA, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String ind = elems[0];
            String dna = elems[1];
            if (!dna.equals("NA")) {
                HashSet<String> set = indToDNA.get(ind);
                if (set == null) {
                    set = new HashSet<>();
                }
                set.add(dna);
                indToDNA.put(ind, set);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        HashMap<String, HashSet<String>> indToRNA = new HashMap<>();
        tf = new TextFile(cmcIndToRNA, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String ind = elems[0];
            String rna = elems[1];
            if (!rna.equals("NA")) {
                HashSet<String> set = indToDNA.get(ind);
                if (set == null) {
                    set = new HashSet<>();
                }
                set.add(rna);
                indToRNA.put(ind, set);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        // link the IDs
        TextFile output = new TextFile(cmcDNAToRNA, TextFile.W);
        TextFile output2 = new TextFile(cmcDNAToRNA + "DNAsWithoutRNA.txt", TextFile.W);
        for (String key : indToDNA.keySet()) {
            HashSet<String> genotypes = indToDNA.get(key);
            HashSet<String> rnas = indToRNA.get(key);
            if (rnas != null) {
                for (String geno : genotypes) {
                    for (String rna : rnas) {
                        output.writeln(geno + "\t" + rna + "\t" + key);
                    }
                }
            } else {
                output2.writeln(indToDNA.get(key) + "\t" + key);
            }
        }
        output.close();
        output2.close();
        output2 = new TextFile(cmcDNAToRNA + "RNAsWithoutDNA.txt", TextFile.W);
        for (String key : indToRNA.keySet()) {
            HashSet<String> genotypes = indToDNA.get(key);
            if (genotypes == null) {
                output2.writeln(indToRNA.get(key) + "\t" + key);
            }
        }
        output2.close();


    }

    private void cmcLinkCheck(String cmcLinks, String individualfile, String s) throws IOException {

        HashMap<String, String> genotypeIds = new HashMap<String, String>();
        HashMap<String, String> aliases = new HashMap<String, String>();
        TextFile tf = new TextFile(individualfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        TextFile out2 = new TextFile(s + "-allIds.txt", TextFile.W);
        while (elems != null) {
            String ind = elems[0];
            if (ind.contains("CMC")) {
                TextFile tf2 = new TextFile(ind, TextFile.R);
                String id = tf2.readLine();
                while (id != null) {
                    genotypeIds.put(id, id);
                    out2.writeln(id);
                    String[] idelems = id.split("_");
                    if (idelems.length > 2) {
                        String newId = idelems[1] + "_" + idelems[2];
                        genotypeIds.put(newId, id);
                        aliases.put(newId, id);
                        out2.writeln(newId);
                    }
                    id = tf2.readLine();
                }
                tf2.close();
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        out2.close();
        tf.close();
        System.out.println(genotypeIds.size() + " genotypes loaded.");
        TextFile out = new TextFile(s, TextFile.W);
        TextFile in = new TextFile(cmcLinks, TextFile.R);
        elems = in.readLineElems(TextFile.tab);
        HashSet<String> idsFound = new HashSet<>();
        while (elems != null) {
            String gt = elems[0];
            if (genotypeIds.containsKey(gt)) {
                idsFound.add(genotypeIds.get(gt));
                out.writeln(genotypeIds.get(gt) + "\t" + elems[1]);
            }
            elems = in.readLineElems(TextFile.tab);
        }
        in.close();
        out.close();


    }

    private static void compareLists(String newList, String oldList) throws IOException {
        HashSet<String> list1 = new HashSet<>();
        TextFile tf = new TextFile(oldList, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String q = elems[0] + "\t" + elems[1] + "\t" + elems[4];
            list1.add(q);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        HashSet<String> list2 = new HashSet<>();
        tf = new TextFile(newList, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String q = elems[0] + "\t" + elems[1] + "\t" + elems[4];

            if (list2.contains(q)) {
                System.out.println("Duplicate: " + Strings.concat(elems, Strings.tab));
            }
            list2.add(q);

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        System.out.println(list1.size());
        System.out.println(list2.size());

        System.out.println();
        System.out.println("Compare old to new:");
        for (String s : list1) {
            if (!s.startsWith("CMC_HBCC") && !s.startsWith("R")) {
                if (!list2.contains(s)) {
                    System.out.println(s + "\tnot present in new list");
                }
            }
        }


        System.out.println("Compare new to old:");
        for (String s : list2) {
            if (!s.startsWith("CMC_HBCC") && !s.startsWith("R")) {
                if (!list1.contains(s)) {
                    System.out.println(s + "\tnot present in old list");
                }
            }
        }
    }

    private void run2(String listOfExpressionSamples, String genotypedIndividualsFile,
                      String[] originalLinks,
                      String freeze1Links,
                      String psychEncodeLinks,
                      String cmcLinks, boolean excludepihat, String output) throws IOException {

        // load supposed links between samples
        HashMap<String, String> rnaToGenotype = new HashMap<String, String>();
        HashMap<String, String> rnaToDataset = new HashMap<String, String>();

        for (String linkfile : originalLinks) {
            TextFile tf = new TextFile(linkfile, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                rnaToGenotype.put(elems[0], elems[1]);
                if (elems.length > 3) {
                    rnaToDataset.put(elems[0], elems[3] + "\t" + elems[2]);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
        }


        HashMap<String, String> freeze1RnaToGenotype = new HashMap<String, String>();
        TextFile tf = new TextFile(freeze1Links, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            freeze1RnaToGenotype.put(elems[0], elems[1]);
//            freeze1RnaToGenotype.put(elems[0], elems[2] + " / " + elems[3]);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        HashMap<String, HashSet<String>> cmcRnaToGenotype = new HashMap<>();
        tf = new TextFile(cmcLinks, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems[1].contains("CMC_HBCC_RNA_PFC_3207")) {
                System.out.println("?");
            }
            HashSet<String> samples = cmcRnaToGenotype.get(elems[1]);
            if (samples == null) {
                samples = new HashSet<>();
            }
            samples.add(elems[0]);

//            freeze1RnaToGenotype.put(elems[0], elems[2] + " / " + elems[3]);
            cmcRnaToGenotype.put(elems[1], samples);

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        // pysch encode is coded as dna -> exp
        HashMap<String, HashSet<String>> psychEncodeRnaToGenotype = new HashMap<>();
        tf = new TextFile(psychEncodeLinks, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            HashSet<String> samples = psychEncodeRnaToGenotype.get(elems[0]);
            if (samples == null) {
                samples = new HashSet<>();
            }
            samples.add(elems[1]);

//            freeze1RnaToGenotype.put(elems[0], elems[2] + " / " + elems[3]);
            psychEncodeRnaToGenotype.put(elems[0], samples);

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        // load all genotyped individuals. Make sure to exclude those that have pihat > 0.125, but keep track of which samples were excluded
        // and also which genotyped individuals are present in multiple files
        HashMap<String, HashSet<String>> availableGenotypeIds = new HashMap<>();
        HashMap<String, HashSet<String>> allExclusions = new HashMap<>();
        TextFile tf2 = new TextFile(genotypedIndividualsFile, TextFile.R);
        String[] elems2 = tf2.readLineElems(TextFile.tab);
        TextFile gtOut = new TextFile(output + "AvailableGenotypeIDs.txt", TextFile.W);
        gtOut.writeln("GenotypeID\tDatasetFile");
        while (elems2 != null) {
            Set<String> exclusions = null;
            if (excludepihat) {
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
                }
            }

            String sampleFile = elems2[0];
            TextFile tf3 = new TextFile(sampleFile, TextFile.R);
            String ln = tf3.readLine();
            while (ln != null) {

                if (exclusions == null || !exclusions.contains(ln)) {
                    HashSet<String> fileset = availableGenotypeIds.get(ln);
                    if (ln.contains("4256126122")) {
                        System.out.println("Got it");
                    }
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

            elems2 = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();
        gtOut.close();
        // now iterate the expression samples and see which samples we can and can't link
        TextFile tf3 = new TextFile(listOfExpressionSamples, TextFile.R);
        TextFile expSamplesNotLinked = new TextFile(output + "ExpressionSamplesNotLinked.txt", TextFile.W);
        expSamplesNotLinked.writeln("RnaID\tRNADataset\tGenotypeQuery\tReason");
        TextFile expSamplesLinked = new TextFile(output + "ExpressionSamplesLinked.txt", TextFile.W);
        expSamplesLinked.writeln("RnaID\tGenotypeID\tMetaCohort\tRNADataset\tGenotypeIndividualsDataset");

        String rnaSample = tf3.readLine();

        while (rnaSample != null) {
            String rnaDataset = rnaToDataset.get(rnaSample);
            if (rnaSample.contains("CMC_HBCC_RNA_PFC_3207")) {
                System.out.println("Watch this");
            }
            String genotypeImLookingFor = rnaToGenotype.get(rnaSample);
            if (genotypeImLookingFor == null) {
                System.err.println("Error no original genotype assignment for: " + rnaSample + "\t" + rnaDataset);
                System.exit(0);
            }

            HashSet<String> genotypeDatasets = availableGenotypeIds.get(genotypeImLookingFor);
            if (genotypeDatasets == null) {
                // genotype individual not found

                // check whether we can recover it using Freeze 1 links
                String genotypeImActuallyLookingFor = freeze1RnaToGenotype.get(rnaSample);
                boolean found = false;
                if (!found) {
                    if (genotypeImActuallyLookingFor != null) {
                        genotypeDatasets = availableGenotypeIds.get(genotypeImActuallyLookingFor);
                        if (genotypeDatasets != null) {
                            for (String genotypeDataset : genotypeDatasets) {
                                expSamplesLinked.writeln(rnaSample + "\t" + genotypeImActuallyLookingFor + "\t" + rnaDataset + "\t" + genotypeDataset);
                            }
                            found = true;
                        }
                    }
                }
                if (!found) {
                    // check old psychencode links
                    HashSet<String> genotypesImActuallyLookingFor = psychEncodeRnaToGenotype.get(rnaSample);
                    if (genotypesImActuallyLookingFor != null) {
                        for (String genotype : genotypesImActuallyLookingFor) {
                            genotypeDatasets = availableGenotypeIds.get(genotype);
                            if (genotypeDatasets != null) {
                                for (String genotypeDataset : genotypeDatasets) {
                                    expSamplesLinked.writeln(rnaSample + "\t" + genotypeImActuallyLookingFor + "\t" + rnaDataset + "\t" + genotypeDataset);
                                }
                                found = true;
                            }
                        }
                    }
                }

                // do some dataset specific fixes...
                if (!found) {
                    // braineac
                    genotypeImActuallyLookingFor = "0_" + genotypeImLookingFor + "_1";
                    genotypeDatasets = availableGenotypeIds.get(genotypeImActuallyLookingFor);
                    if (genotypeDatasets != null) {
                        for (String genotypeDataset : genotypeDatasets) {
                            expSamplesLinked.writeln(rnaSample + "\t" + genotypeImActuallyLookingFor + "\t" + rnaDataset + "\t" + genotypeDataset);
                        }
                        found = true;
                    }
                }

                if (!found) {
                    // NABEC
                    genotypeImActuallyLookingFor = "0_" + genotypeImLookingFor;
                    genotypeDatasets = availableGenotypeIds.get(genotypeImActuallyLookingFor);
                    if (genotypeDatasets != null) {
                        for (String genotypeDataset : genotypeDatasets) {
                            expSamplesLinked.writeln(rnaSample + "\t" + genotypeImActuallyLookingFor + "\t" + rnaDataset + "\t" + genotypeDataset);
                        }
                        found = true;
                    }
                    if (!found) {
                        // NABEC-h610
                        genotypeImActuallyLookingFor = "0_" + genotypeImLookingFor + "_" + genotypeImLookingFor;
                        genotypeDatasets = availableGenotypeIds.get(genotypeImActuallyLookingFor);
                        if (genotypeDatasets != null) {
                            for (String genotypeDataset : genotypeDatasets) {
                                expSamplesLinked.writeln(rnaSample + "\t" + genotypeImActuallyLookingFor + "\t" + rnaDataset + "\t" + genotypeDataset);
                            }
                            found = true;
                        }
                    }

                    if (!found) {
                        genotypeImActuallyLookingFor = "0_" + genotypeImLookingFor.replaceAll("-", "_");
                        genotypeDatasets = availableGenotypeIds.get(genotypeImActuallyLookingFor);
                        if (genotypeDatasets != null) {
                            for (String genotypeDataset : genotypeDatasets) {
                                expSamplesLinked.writeln(rnaSample + "\t" + genotypeImActuallyLookingFor + "\t" + rnaDataset + "\t" + genotypeDataset);
                            }
                            found = true;
                        }
                    }
                }

                // now find the last CMC_HBCC samples
                if (!found) {
                    if (rnaSample.startsWith("CMC_HBCC_RNA_PFC_3207")) {
                        System.out.println("Got it");
                    }

                    // check old psychencode links ; some of them list the actual genotype ID... :(
                    HashSet<String> genotypesImActuallyLookingFor = psychEncodeRnaToGenotype.get(genotypeImLookingFor);
                    if (genotypesImActuallyLookingFor == null) {
                        genotypesImActuallyLookingFor = new HashSet<>();
                    }
                    // fix: add original CMC links as well
                    HashSet<String> cmc = cmcRnaToGenotype.get(rnaSample);
                    if (cmc != null) {
                        genotypesImActuallyLookingFor.addAll(cmc);
                    }
                    if (genotypesImActuallyLookingFor != null) {
                        for (String genotype : genotypesImActuallyLookingFor) {
                            // check if we can find it immediately
                            genotypeDatasets = availableGenotypeIds.get(genotype);
                            if (genotypeDatasets != null) {
                                for (String genotypeDataset : genotypeDatasets) {
                                    expSamplesLinked.writeln(rnaSample + "\t" + genotype + "\t" + rnaDataset + "\t" + genotypeDataset);
                                }
                                found = true;
                            }

                            // else, do some magic:
                            // actual genotype has a family ID stuck in front.
                            if (!found) {
                                // iterate all individuals
                                for (String ind : availableGenotypeIds.keySet()) {
                                    String[] indElems = ind.split("_");
                                    if (indElems.length == 3) {
                                        String indWithoutFamId = indElems[1] + "_" + indElems[2];
                                        if (indWithoutFamId.equals(genotype)) {
                                            // found it!
                                            if (indWithoutFamId.contains("4463344373")) {

                                                // should be 4040296008_A, is
                                                System.out.println("found the gt as well");
                                            }
                                            String genotypeImActuallyActuallyLookingFor = ind;
                                            genotypeDatasets = availableGenotypeIds.get(genotypeImActuallyActuallyLookingFor);
                                            if (genotypeDatasets != null) {
                                                for (String genotypeDataset : genotypeDatasets) {
                                                    expSamplesLinked.writeln(rnaSample + "\t" + genotypeImActuallyActuallyLookingFor + "\t" + rnaDataset + "\t" + genotypeDataset);
                                                }
                                                found = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (!found) {
                    // TargetALS
                    genotypeImActuallyLookingFor = genotypeImLookingFor + "-b38";
                    genotypeDatasets = availableGenotypeIds.get(genotypeImActuallyLookingFor);
                    if (genotypeDatasets != null) {
                        for (String genotypeDataset : genotypeDatasets) {
                            expSamplesLinked.writeln(rnaSample + "\t" + genotypeImActuallyLookingFor + "\t" + rnaDataset + "\t" + genotypeDataset);
                        }
                        found = true;
                    }
                }


                if (!found) {
                    // check whether it is in the exclusions list
                    HashSet<String> exclusionFiles = allExclusions.get(genotypeImLookingFor);
                    if (exclusionFiles == null) {
                        expSamplesNotLinked.writeln(rnaSample + "\t" + rnaDataset + "\t" + genotypeImLookingFor + "\tGenotypeID not found");
                    } else {
                        ArrayList<String> exclusionlist = new ArrayList<>();
                        exclusionlist.addAll(exclusionFiles);
                        expSamplesNotLinked.writeln(rnaSample + "\t" + rnaDataset + "\t" + genotypeImLookingFor + "\tGenotypeQC: " + Strings.concat(exclusionlist, Strings.semicolon));
                    }
                }


            } else {
                for (String genotypeDataset : genotypeDatasets) {
                    expSamplesLinked.writeln(rnaSample + "\t" + genotypeImLookingFor + "\t" + rnaDataset + "\t" + genotypeDataset);
                }
            }
            rnaSample = tf3.readLine();
        }
        expSamplesNotLinked.close();
        expSamplesLinked.close();
        tf3.close();

    }

    private void removeDuplicateAndRelatedDNA(String defupfileout, String duplicateDNAIds, String pathstrippeddedupfilewithoutrelated) throws IOException {
        HashSet<String> idsToExclude = new HashSet<>();
        TextFile tf = new TextFile(duplicateDNAIds, TextFile.R);
        idsToExclude.addAll(tf.readAsArrayList());
        tf.close();
        TextFile tf2 = new TextFile(defupfileout, TextFile.R);
        TextFile outf = new TextFile(pathstrippeddedupfilewithoutrelated, TextFile.W);

        String[] elems = tf2.readLineElems(TextFile.tab);
        while (elems != null) {

            String dna = elems[1];
            if (!idsToExclude.contains(dna)) {
                outf.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = tf2.readLineElems(TextFile.tab);
        }

        outf.close();
        tf2.close();

    }

    private void stripPaths(String defupfileout, String pathstrippeddedupfile) throws IOException {

        TextFile in = new TextFile(defupfileout, TextFile.R);
        TextFile out = new TextFile(pathstrippeddedupfile, TextFile.W);

        String[] elems = in.readLineElems(TextFile.tab);
        while (elems != null) {

            String path = elems[4];
            String[] pathelems = path.split("\\\\");
            String name = pathelems[pathelems.length - 1];
            name = name.replace("-Individuals.txt", "");
            elems[4] = name;

            out.writeln(Strings.concat(elems, Strings.tab));
            elems = in.readLineElems(TextFile.tab);
        }
        in.close();
        out.close();
    }

    private void splitPerDataset(String linkfile, String defupfileout) throws IOException {
        HashMap<String, String> samplePerDataset = new HashMap<>();
        HashMap<String, TextFile> dsout = new HashMap<>();

        TextFile tf = new TextFile(defupfileout, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);


        while (elems != null) {

            String sample = elems[0];
            String dataset = elems[4];
            TextFile tfo = dsout.get(dataset);
            if (tfo == null) {
                dsout.put(dataset, new TextFile(linkfile + "-" + dataset + ".txt", TextFile.W));
            }
            samplePerDataset.put(sample, dataset);

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(linkfile, TextFile.R);

        elems = tf2.readLineElems(TextFile.tab);
        TextFile sampleToDs = new TextFile(linkfile + "-SampleToDataset.txt", TextFile.W);
        while (elems != null) {

            String sample = elems[1];
            String ds = samplePerDataset.get(sample);
            sampleToDs.writeln(sample + "\t" + ds);
            if (ds == null) {
                System.out.println("Could not find ds for " + sample);

            } else {
                dsout.get(ds).writeln(elems[0] + "\t" + elems[1]);
            }

            elems = tf2.readLineElems(TextFile.tab);
        }
        sampleToDs.close();
        tf2.close();
        for (String key : dsout.keySet()) {
            dsout.get(key).close();
        }

    }

    private void splitPerPopulation(String file, String popfile) throws IOException {
        HashMap<String, String> popPerSample = new HashMap<>();
        HashMap<String, TextFile> popout = new HashMap<>();
        TextFile tf = new TextFile(popfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            popPerSample.put(elems[0], elems[1]);
            TextFile tfout = popout.get(elems[1]);
            if (tfout == null) {
                popout.put(elems[1], new TextFile(file + "-" + elems[1] + ".txt", TextFile.W));
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(file, TextFile.R);

        elems = tf2.readLineElems(TextFile.tab);
        while (elems != null) {

            String gt = elems[0];
            String pop = popPerSample.get(gt);
//            if (gt.startsWith("SAMN") || gt.startsWith("SAMEA")) {
//                pop = "EUR";
//            }

            if (pop == null) {
                System.out.println("Could not find population for sample: " + gt);
            } else {
                popout.get(pop).writeln(elems[0] + "\t" + elems[1]);
            }

            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();

        for (String key : popout.keySet()) {
            popout.get(key).close();
        }
    }

    private void dedupDNA(String file, boolean preferreseqenced) throws IOException {

        if (preferreseqenced) {


            HashMap<String, HashSet<String>> links = new HashMap<>();
            TextFile tf = new TextFile(file, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                HashSet<String> rnas = links.get(elems[1]);
                if (rnas == null) {
                    rnas = new HashSet<>();
                }
                rnas.add(elems[0]);
                links.put(elems[1], rnas);

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile tfo = new TextFile(file + "-dedup-gte.txt", TextFile.W);
            for (String key : links.keySet()) {
                HashSet<String> rnas = links.get(key);
                String reseqsample = null;
                String othersample = null;
                for (String s : rnas) {
                    if (s.contains("reseq")) {
                        reseqsample = s;
                    } else {
                        othersample = s;
                    }
                }

                if (reseqsample == null) {
                    // pick a random sample
                    tfo.writeln(key + "\t" + othersample);
                } else {
                    tfo.writeln(key + "\t" + reseqsample);
                }

            }
            tfo.close();

        } else {
            HashSet<String> visitedInd = new HashSet<>();
            TextFile tf = new TextFile(file, TextFile.R);
            TextFile tfo = new TextFile(file + "-dedup-gte.txt", TextFile.W);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                if (!visitedInd.contains(elems[1])) {
                    tfo.writeln(elems[1] + "\t" + elems[0]);
                    visitedInd.add(elems[1]);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            tfo.close();
        }


    }


    private void splitBasedOnTissueType(String defupfileout, String tissuemapfile, String splittissuefolder) throws IOException {
        Gpio.createDir(splittissuefolder);
        HashMap<String, String> tissuemap = new HashMap<>();
        HashMap<String, TextFile> tissueout = new HashMap<>();
        TextFile tf = new TextFile(tissuemapfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            String sample = elems[0];
            String tissue = elems[1];

            tissuemap.put(sample, tissue);
            TextFile tissuetf = tissueout.get(tissue);
            if (tissuetf == null) {
                tissueout.put(tissue, new TextFile(splittissuefolder + tissue + ".txt", TextFile.W));
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(defupfileout, TextFile.R);
        elems = tf2.readLineElems(TextFile.tab);
        while (elems != null) {

            String sample = elems[0];
            String tissue = tissuemap.get(sample);
            if (tissue == null) {
                System.err.println("Tissue assignment missing for sample " + sample);
            } else {
                tissueout.get(tissue).writeln(sample + "\t" + elems[1]);
            }
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();

        for (String key : tissueout.keySet()) {
            tissueout.get(key).close();
        }

    }

    private void dedupRNA(String file, String fileout) throws IOException {
        HashSet<String> visitedRNA = new HashSet<>();
        TextFile tf = new TextFile(file, TextFile.R);
        TextFile tfo = new TextFile(fileout, TextFile.W);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            if (!visitedRNA.contains(elems[0])) {
                tfo.writeln(Strings.concat(elems, Strings.tab));
                visitedRNA.add(elems[0]);
            }

            elems = tf.readLineElems(TextFile.tab);
        }
        tfo.close();
        tf.close();
    }

    public class Quad {
        String rna;
        String dna;
        String ds;
        String meta;
    }

    public HashMap<String, String> readAsMap(String file) throws IOException {
        HashMap<String, String> map = new HashMap<>();
        TextFile tfq = new TextFile(file, TextFile.R);
        String[] elems = tfq.readLineElems(TextFile.tab);
        while (elems != null) {
            if (!elems[1].equals("NotMapped")) {
                map.put(elems[0], elems[1]);
            }
            elems = tfq.readLineElems(TextFile.tab);

        }
        tfq.close();
        return map;
    }

    public void run(String indfile, String origlink, String freeze1links, String psychencodelinks, String rnseqoutlierfile, String output) throws IOException {

        HashSet<String> rnaseqoutliers = null;
        if (rnseqoutlierfile != null) {
            rnaseqoutliers = new HashSet<>();
            TextFile tf = new TextFile(rnseqoutlierfile, TextFile.R);
            rnaseqoutliers.addAll(tf.readAsArrayList());
            tf.close();
        }

        ArrayList<Quad> data = new ArrayList<>();

        TextFile tf1 = new TextFile(origlink, TextFile.R);
        String[] elems1 = tf1.readLineElems(TextFile.tab);

        while (elems1 != null) {
            Quad q = new Quad();
            q.dna = elems1[1];
            q.rna = elems1[0];
            q.meta = elems1[2];
            q.ds = elems1[3];
            if (rnaseqoutliers == null || !rnaseqoutliers.contains(q.rna)) {
                data.add(q);
            }

            elems1 = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();

        HashMap<String, String> freeze1map = readAsMap(freeze1links);
        HashMap<String, String> psychmap = readAsMap(psychencodelinks);

        TextFile tf = new TextFile(indfile, TextFile.R);

        String[] lnelems = tf.readLineElems(TextFile.tab);
        HashSet<String> foundinds = new HashSet<>();
        TextFile allindOut = new TextFile(output + "AllInds.txt", TextFile.W);
        TextFile foundindOut = new TextFile(output + "FoundInds.txt", TextFile.W);

        while (lnelems != null) {

            System.out.println(lnelems[0]);
            HashSet<String> excludeInds = new HashSet<String>();
            if (lnelems.length > 1) {
                System.out.println(lnelems[1]);
                TextFile tf2 = new TextFile(lnelems[1], TextFile.R);
                String[] dupelems = tf2.readLineElems(TextFile.tab);
                while (dupelems != null) {
                    excludeInds.add(dupelems[3]);
                    dupelems = tf2.readLineElems(TextFile.tab);
                }
                tf2.close();
            }


            //
            HashSet<String> inds = new HashSet<>();
            String individualFile = lnelems[0];
            TextFile tf3 = new TextFile(individualFile, TextFile.R);
            String ln2 = tf3.readLine();
            int included = 0;
            int excluded = 0;
            while (ln2 != null) {
                if (excludeInds.contains(ln2)) {
                    excluded++;
                } else {
                    inds.add(ln2);
                    included++;
                    allindOut.writeln(ln2 + "\t" + lnelems[0]);
                }

                ln2 = tf3.readLine();
            }
            tf3.close();

            System.out.println("File has: " + included + " included, " + excluded + " excluded.");

            boolean found = false;
            for (Quad q : data) {
                found = false;
                String rna = q.rna;
                String f1 = freeze1map.get(rna);

                String dna = q.dna;
                if (inds.contains(dna)) {
                    foundinds.add(f1);
                    foundinds.add(dna);
                    foundindOut.writeln(rna + "\t" + dna + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
                    found = true;
                }
                if (!found && f1 != null) {

                    // try freeze 1 mapping
                    if (inds.contains(f1)) {
                        foundinds.add(f1);
                        foundinds.add(dna);
                        foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
                        found = true;
                    }
                }

                // targetals
                if (!found) {
                    f1 = dna + "-b38";
                    if (inds.contains(f1)) {
                        foundinds.add(f1);
                        foundinds.add(dna);
                        foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
                        found = true;
                    }
                }

                // NABEC
                if (!found) {
                    f1 = "0_" + dna;
                    if (inds.contains(f1)) {
                        foundinds.add(f1);
                        foundinds.add(dna);
                        foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
                        found = true;
                    }
                }

                if (!found) {
                    if (dna.contains("UMARY")) {
//						System.out.println("Try this!");
                    }
                    f1 = "0_" + dna.replaceAll("-", "_");
                    if (inds.contains(f1)) {
                        foundinds.add(f1);
                        foundinds.add(dna);
                        foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
                        found = true;
                    }
                }

                // HBCC
                if (!found) {

                    for (String ind : inds) {
                        String[] indelems = ind.split("_");
                        if (indelems.length == 3) {
                            String actual = indelems[1] + "_" + indelems[2];
                            if (actual.equals(dna)) {

                                foundinds.add(f1);
                                foundinds.add(dna);
                                foundindOut.writeln(rna + "\t" + ind + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
                                found = true;
                                break;
                            }
                        }
                    }
                }
                // braineac
                if (!found) {
                    f1 = "0_" + dna + "_1";
                    if (inds.contains(f1)) {
                        foundinds.add(dna);
                        foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
                        found = true;
                    }
                }

                // NABEC-H610
                if (!found) {
                    // SH-07-73 --> 0_SH-07-73_SH-07-73
                    f1 = "0_" + dna + "_" + dna;
                    if (inds.contains(f1)) {
                        foundinds.add(dna);
                        foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
                        found = true;
                    }
                }

                if (!found) {
                    if (rna.equals("Br1157_R2836") && individualFile.contains("CMC_HBCC_set1")) {
                        System.out.println("Found it");
                    }
                    String psychencodeind = psychmap.get(dna);
                    if (psychencodeind != null) {
                        for (String ind : inds) {
                            String[] indelems = ind.split("_");
                            if (indelems.length == 3) {
                                String actual = indelems[1] + "_" + indelems[2];
                                if (actual.equals(psychencodeind)) {
                                    foundinds.add(psychencodeind);
                                    foundinds.add(dna);
                                    foundindOut.writeln(rna + "\t" + ind + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
                                    found = true;
                                    break;
                                }
                            }
                        }

                    }
                }

            }

            System.out.println(foundinds.size() + " ids found sofar..");
            lnelems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        allindOut.close();
        foundindOut.close();

        TextFile notfoundOut = new TextFile(output + "NotFoundInds.txt", TextFile.W);
        int nrnoutfound = 0;
        for (Quad d : data) {
            if (!foundinds.contains(d.dna)) {
                notfoundOut.writeln(d.rna + "\t" + d.dna + "\t" + d.ds + "\t" + d.meta);
                nrnoutfound++;
            }
        }
        notfoundOut.close();
        System.out.println(nrnoutfound + " not found out of " + data.size());


    }

    public void run(String[] individualFiles, String origlink, String output) throws IOException {
        HashMap<String, String> rnaToGenotype = new HashMap<>();
        TextFile tf = new TextFile(origlink, TextFile.R);
        tf.readLine();

        HashSet<String> visitedRNA = new HashSet<>();

        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems.length > 2) {
                String id = elems[0] + "_" + elems[1];
                if (visitedRNA.contains(elems[1])) {
                    System.out.println("Dup RNA: " + elems[1]);
                }
                visitedRNA.add(elems[1]);
                String gt = elems[2];
                if (rnaToGenotype.containsKey(id)) {
                    System.out.println("Dup: " + id + "\t" + elems[0] + "\t" + elems[1]);
                } else {
                    rnaToGenotype.put(id, gt);
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


    }
}
