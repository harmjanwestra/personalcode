package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class ReplaceMDSIds {


    public static void main(String[] args) {

        ReplaceMDSIds r = new ReplaceMDSIds();

        try {
//			r.run("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-MDS\\Mayo\\mayomds.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\MAYO-DNAToRNA-CBE-samples.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-MDS\\Mayo\\mayomds-CBE.txt"
////			);
//			r.run("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-MDS\\Mayo\\mayomds.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\MAYO-DNAToRNA-TCX-samples.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-MDS\\Mayo\\mayomds-TCX.txt"
//			);
//			r.run("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-MDS\\Msbb\\mssbmds.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMSBB\\MSBB-DNAToRNA.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-MDS\\Msbb\\mssbmds-EXP.txt"
//			);
//			r.run("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-MDS\\Rosmap\\rosmapmds.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADRosmap\\ROSMAP-DNAToRNA.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-MDS\\Rosmap\\rosmapmds-EXP.txt"
//			);

            // cmc
            r.run("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\mds\\mds.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\linkfiles\\genotype-rna-CMC.txt-filtered.txt",
                    "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-CMC\\mds\\mds-rewrite.txt");

            String[] tissues = new String[]{
                    "Any",
                    "BACCO",
                    "BAMYG",
                    "BBGCA",
                    "BBGNA",
                    "BBGPU",
                    "BCERE",
                    "BCHEM",
                    "BCORT",
                    "BFRCO",
                    "BHIPP",
                    "BHYPO",
                    "BSPCO",
                    "BSUNI"
            };
            for (String tissue : tissues) {
                r.run("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\mds\\gtex.mds.txt",
                        "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\initialGTEs\\GTEx-" + tissue + ".txt",
                        "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\mds\\gtex.mds-rewrite-" + tissue + ".txt");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void run(String mds, String map, String out) throws IOException {
        TextFile tf = new TextFile(map, TextFile.R);
        Map<String, String> gte = tf.readAsHashMap(0, 1);
        tf.close();

        TextFile m = new TextFile(mds, TextFile.R);
        TextFile mo = new TextFile(out, TextFile.W);
        mo.writeln(m.readLine());

        String[] elems = m.readLineElems(TextFile.tab);
        while (elems != null) {
            String sample = elems[0];
            if (gte.containsKey(sample)) {
                String conv = gte.get(sample);
                elems[0] = conv;
                mo.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = m.readLineElems(TextFile.tab);
        }
        mo.close();
        m.close();
    }


    public void rewriteWithDuplicates(String mds, String mapfile, boolean concatIds, String concatStr, String out) throws IOException {
        HashMap<String, ArrayList<String>> map = new HashMap<String, ArrayList<String>>();
        TextFile tf = new TextFile(mapfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            if (elems.length > 1) {
                String dna = elems[0];
                ArrayList<String> list = map.get(dna);
                if (list == null) {
                    list = new ArrayList<>();
                }

                list.add(elems[1]);
                map.put(dna, list);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(mds, TextFile.R);
        TextFile tfo = new TextFile(out, TextFile.W);
//        String ln
        String ln = tf2.readLine();
        int ctr = 0;
        while (ln != null) {
            while (ln.contains("  ")) {
                ln = ln.replaceAll("  ", " ");
            }
            if (ln.startsWith(" ")) {
                ln = ln.substring(1);
            }
            elems = ln.split(" ");
            if (ctr == 0) {
                String outStr = "Sample\tMDS1\tMDS2\tMDS3\tMDS4";
                tfo.writeln(outStr);
            } else {
                String id = elems[1];
                if (concatIds) {
                    id = elems[0] + concatStr + elems[1];
                }

                ArrayList<String> set = map.get(id);
                if (set != null) {
                    for (String sample : set) {
                        String outstr = sample + "\t" + Strings.concat(elems, Strings.tab, 2, 5);
                        tfo.writeln(outstr);
                    }
                }
            }

            ln = tf2.readLine();
            ctr++;
        }
        tf2.close();
    }
}
