package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class Subset {

    public static void main(String[] args) {

//        String europeans = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-europeans.txt";
//        String dups = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\geneticsimilarity\\duplicates.txt";

        String europeans = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\AMPAD-Europeans.txt";
        String dups = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\geneticsimilarity\\duplicates.txt";
        Subset s = new Subset();


        try {
//			s.gtfToProbeAnnotationFile(europeans, dups, false);


            String europeanfilter = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\AMPAD-Europeans.txt";
            String duplicateDNAs = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\geneticsimilarity\\duplicates.txt";
            String[] mapfile = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\linksToRNA\\MayoCBE_final.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\linksToRNA\\MayoTCX_final.txt",
//                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\linksToRNA\\MSBB_final.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\linksToRNA\\ROSMAP_final.txt"
            };
            String[] rnaoutliers = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-AMPAD\\mayocbe-outliers.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-AMPAD\\mayotcx-outliers.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-AMPAD\\rosmap-outliers.txt",
            };

            for (int i = 0; i < mapfile.length; i++) {

//                s.filterdups(mapfile[i], europeanfilter, duplicateDNAs, rnaoutliers[i]);
            }

            // CMC
            europeanfilter = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\genotypepca\\CMC-Europeans-gtids.txt";
            duplicateDNAs = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\geneticsimilarity\\dnadup.txt";
            String rnastoremove = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\outliers.txt";
            mapfile = new String[]{
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\linkfiles\\genotype-rna-CMC.txt"
            };

            for (String map : mapfile) {
                s.filterdups(map, europeanfilter, duplicateDNAs, rnastoremove);
            }

            // TargetALS
            europeanfilter = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-europeans.txt";
            duplicateDNAs = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\geneticsimilarity\\duplicates.txt";
             rnastoremove = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\rnapca\\outliers-pc1lt0.041.txt";
            mapfile = new String[]{

                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\gtepertissue\\GTE-Cell_line.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\gtepertissue\\GTE-Cerebellum.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\gtepertissue\\GTE-Cortex.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\gtepertissue\\GTE-Group.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\gtepertissue\\GTE-Liver.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\gtepertissue\\GTE-MotorCortex.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\gtepertissue\\GTE-SpinalCord.txt"
            };
            for (String map : mapfile) {
//                s.filterdups(map, europeanfilter, duplicateDNAs, rnastoremove);
            }


        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private void filterdups(String mapfile, String dnafilter, String dnadupfilter, String rnafilter) throws IOException {


        HashSet<String> set = null;
        if (dnafilter != null) {
            TextFile tf2 = new TextFile(dnafilter, TextFile.R);
            ArrayList<String> list2 = tf2.readAsArrayList();
            tf2.close();
            set = new HashSet<String>();
            set.addAll(list2);
        }

        HashSet<String> setrna = null;
        if (rnafilter != null) {
            TextFile tf2 = new TextFile(rnafilter, TextFile.R);
            ArrayList<String> list2 = tf2.readAsArrayList();
            tf2.close();
            setrna = new HashSet<String>();
            setrna.addAll(list2);
        }
        HashSet<String> setdups = null;
        if (dnadupfilter != null) {
            TextFile tf2 = new TextFile(dnadupfilter, TextFile.R);
            ArrayList<String> list2 = tf2.readAsArrayList();
            tf2.close();
            setdups = new HashSet<String>();
            setdups.addAll(list2);
        }


        HashSet<String> uniquelyLinkedDNA = new HashSet<>();
        TextFile tf = new TextFile(mapfile, TextFile.R);
        TextFile tfo = new TextFile(mapfile + "-filtered.txt", TextFile.W);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            if (elems.length > 1 && elems[0].length() > 0 && elems[1].length() > 0) {
                String rna = elems[1];
                if (setrna == null || !setrna.contains(rna)) {
                    String dna = elems[0];
                    if (set == null || set.contains(dna)) {
                        if (setdups == null || !setdups.contains(dna)) {
                            uniquelyLinkedDNA.add(dna);
                            tfo.writeln(dna + "\t" + rna);
                        }
                    }
                }
            } else {
                System.out.println(Strings.concat(elems, Strings.semicolon));
            }

            elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();
        tfo.close();
        System.out.println(uniquelyLinkedDNA.size() + " unique DNAs in " + mapfile);

    }

    public void run(String file1, String subsetFile, boolean include) throws IOException {
        TextFile tf = new TextFile(file1, TextFile.R);
        ArrayList<String> list = tf.readAsArrayList();
        tf.close();

        TextFile tf2 = new TextFile(subsetFile, TextFile.R);
        ArrayList<String> list2 = tf2.readAsArrayList();
        tf2.close();

        HashSet<String> set = new HashSet<String>();
        set.addAll(list2);

        int overlap = 0;
        for (String s : list) {
            if (include) {
                if (set.contains(s)) {
                    overlap++;
                    System.out.println(s);
                }
            } else {
                if (!set.contains(s)) {
                    overlap++;
                    System.out.println(s);
                }
            }
        }

        System.out.println(overlap + " out of " + list.size());


    }
}
