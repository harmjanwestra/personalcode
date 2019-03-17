package nl.harmjanwestra.playground.biogen.datasets;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class NABEC {

    public static void main(String[] args) {


        try {
            String mapfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\links\\NABEC-rnaseq-to-sampleId.txt";
            HashMap<String, ArrayList<String>> fullmap = readMap(mapfile);

            String cageseqfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\rna\\NABEC-cageseq-samples.txt";
            fullmap = filterRNA(cageseqfile, false, fullmap);

            String rnaSamplesFailingQC = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-samples_to_filter.txt";
            fullmap = filterRNA(rnaSamplesFailingQC, false, fullmap);

            String rnaSamplesRNAPCAQC = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\rna\\NABEC-rnaseq-samples_pc1gt0.0420.txt";
            fullmap = filterRNA(rnaSamplesRNAPCAQC, true, fullmap);

            String hap610inds = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\dna\\hap610-Individuals.txt";
            HashMap<String, ArrayList<String>> map610 = filterDNA(hap610inds, true, fullmap);
            String hap610indsEur = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\dna\\hap610europeans.txt";
            map610 = filterDNA(hap610indsEur, true, map610);
            String hap610indsEurDups = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\dna\\hap610geneticdups.txt";
            map610 = filterDNA(hap610indsEurDups, false, map610);

            String hap610gte = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\links\\hap610-nonunique.txt";
            // write

            // random pick

            System.out.println();

            String hap550inds = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\dna\\hap550-Individuals.txt";
            HashMap<String, ArrayList<String>> map550 = filterDNA(hap550inds, true, fullmap);

            // filter for europeans
            String hap550indsEur = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\dna\\hap550europeans.txt";
            map550 = filterDNA(hap550indsEur, true, map550);

            // filter for duplicates
            String hap550indsEurDups = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\dna\\hap610geneticdups.txt";
            map550 = filterDNA(hap550indsEurDups, false, map550);

            String hap550gte = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-02-26-NABEC\\links\\hap550-nonunique.txt";


        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private static HashMap<String, ArrayList<String>> readMap(String mapfile) throws IOException {
        TextFile tf = new TextFile(mapfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);


        HashMap<String, ArrayList<String>> map = new HashMap<String, ArrayList<String>>();
        int nrsamples = 0;
        while (elems != null) {
            if (elems.length > 1) {
                String dnaname = elems[1].split("_")[0];
                String rna = elems[0];


                ArrayList<String> samples = map.get(dnaname);
                if (samples == null) {
                    samples = new ArrayList<>();
                }
                samples.add(rna);
                nrsamples++;

                map.put(dnaname, samples);

            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println("Map: " + mapfile);
        System.out.println("DNAs: " + map.size() + "\tRNAs: " + nrsamples);
        return map;
    }

    private static HashMap<String, ArrayList<String>> filterDNA(String list, boolean include, HashMap<String, ArrayList<String>> map) throws IOException {
        HashMap<String, ArrayList<String>> tmpmap = new HashMap<>();
        HashSet<String> set = getSet(list);
        int nrsamples = 0;
        for (String key : map.keySet()) {
            if ((include && set.contains(key)) || (!include && !set.contains(key))) {
                tmpmap.put(key, map.get(key));
                nrsamples += map.get(key).size();
            }
        }
        System.out.println(tmpmap.size() + " DNAs remain, and " + nrsamples + " RNAs");
        return tmpmap;
    }

    private static HashMap<String, ArrayList<String>> filterRNA(String list, boolean include, HashMap<String, ArrayList<String>> map) throws IOException {
        HashMap<String, ArrayList<String>> tmpmap = new HashMap<>();
        HashSet<String> set = getSet(list);

        int nrsamples = 0;
        for (String key : map.keySet()) {
            ArrayList<String> tmpList = new ArrayList<>();
            ArrayList<String> prevlist = map.get(key);
            for (String s : prevlist) {
                if ((include && set.contains(s)) || (!include && !set.contains(s))) {
                    tmpList.add(s);
                    nrsamples++;
                }
            }
            if (!tmpList.isEmpty()) {
                tmpmap.put(key, tmpList);
            }
        }
        System.out.println(tmpmap.size() + " DNAs remain, and " + nrsamples + " RNAs");
        return tmpmap;
    }

    private static HashSet<String> getSet(String list) throws IOException {
        HashSet<String> out = new HashSet<>();
        TextFile tf = new TextFile(list, TextFile.R);
        out.addAll(tf.readAsArrayList());
        tf.close();
        System.out.println(out.size() + " read from " + list);
        return out;
    }
}
