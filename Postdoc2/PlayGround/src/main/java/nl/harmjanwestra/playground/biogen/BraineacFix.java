package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class BraineacFix {

    public static void main(String[] args) {
        BraineacFix f = new BraineacFix();
        try {
            f.fixlinks();
            f.fixmds();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void fixmds() throws IOException {
        String[] tissues = new String[]{
                "HIPP", "PUTM", "SNIG"
        };

        for (String t : tissues) {
            TextFile link = new TextFile("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\initiallinks.txt-EUR-" + t + ".txt", TextFile.R);
            HashMap<String, String> map = new HashMap<>();
            String[] elems = link.readLineElems(TextFile.tab);
            while (elems != null) {
                map.put(elems[0], elems[1]);
                elems = link.readLineElems(TextFile.tab);
            }

            link.close();
            TextFile tf = new TextFile("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\mds\\brainfix.mds", TextFile.R);
            TextFile tfo = new TextFile("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\mds\\brainfix-mds-" + t + ".txt", TextFile.W);
            tfo.writeln(tf.readLine());
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String id = elems[0];
                String rna = map.get(id);
                if (rna != null) {
                    tfo.writeln(rna + "\t" + Strings.concat(elems, Strings.tab, 3, elems.length));
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            tfo.close();
        }
    }

    public void fixlinks() {
        String in = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\SampleInfoBiogen.csv";
        String rnasamples = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\rnaseqsamples.txt";
        String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\initiallinks.txt";
        try {

            ArrayList<String> samples = new ArrayList<String>();
            TextFile tf1 = new TextFile(rnasamples, TextFile.R);
            String ln = tf1.readLine();
            while (ln != null) {
                samples.add(ln);
                ln = tf1.readLine();
            }
            tf1.close();


            TextFile tf = new TextFile(in, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.comma);
            TextFile tfo = new TextFile(out, TextFile.W);
            while (elems != null) {
                String rna = elems[0];
                String dna = elems[2];
                for (String s : samples) {
                    if (s.contains(rna)) {
                        System.out.println(rna + "\t" + s + "\t" + dna);
                        tfo.writeln(dna + "\t" + s);
                    }
                }

                elems = tf.readLineElems(TextFile.comma);
            }

            tfo.close();
            tf.close();

            String eursamples = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\dnapca\\Braineac-Europeans.txt";
            TextFile tf5 = new TextFile(eursamples, TextFile.R);
            HashSet<String> eurs = new HashSet<String>();
            String ln2 = tf5.readLine();
            while (ln2 != null) {
                eurs.add(ln2);
                ln2 = tf5.readLine();
            }

            tf5.close();

            ArrayList<String> eurrnasamples = new ArrayList<>();
            {
                TextFile tf3 = new TextFile(out, TextFile.R);
                TextFile tf4 = new TextFile(out + "-EUR.txt", TextFile.W);
                String[] elems3 = tf3.readLineElems(TextFile.tab);

                while (elems3 != null) {
                    if (eurs.contains(elems3[0])) {
                        tf4.writeln(Strings.concat(elems3, Strings.tab));
                    } else {
                        String dna = elems3[0].replaceAll("/", "_");
                        if (eurs.contains(dna)) {
                            elems3[0] = dna;
                            tf4.writeln(Strings.concat(elems3, Strings.tab));
                            eurrnasamples.add(elems3[1]);
                        }
                    }
                    elems3 = tf3.readLineElems(TextFile.tab);
                }

                tf4.close();
                tf3.close();
            }

            HashMap<String, String> rnaToTissue = new HashMap<String, String>();
            HashSet<String> tissues = new HashSet<>();
            TextFile tf6 = new TextFile("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-19-Braineac\\brainregions.txt", TextFile.R);
            String[] elems6 = tf6.readLineElems(TextFile.tab);
            while (elems6 != null) {
                tissues.add(elems6[1]);
                String rna = elems6[0];
                for (String s : eurrnasamples) {
                    if (s.contains(rna)) {
                        rnaToTissue.put(s, elems6[1]);
                    }
                }
                elems6 = tf6.readLineElems(TextFile.tab);
            }
            tf.close();


            for (String tissue : tissues) {
                TextFile tf3 = new TextFile(out + "-EUR.txt", TextFile.R);
                TextFile tf4 = new TextFile(out + "-EUR-" + tissue + ".txt", TextFile.W);
                String[] elems3 = tf3.readLineElems(TextFile.tab);
                HashSet<String> visited = new HashSet<String>();
                HashSet<String> visitedrna = new HashSet<String>();
                while (elems3 != null) {
                    String tissuerna = rnaToTissue.get(elems3[1]);
                    if (tissuerna != null && tissue.equals(tissuerna)) {
                        if (!visited.contains(elems3[0]) && !visitedrna.contains(elems3[1])) {
                            tf4.writeln(Strings.concat(elems3, Strings.tab));
                            visited.add(elems3[0]);
                            visitedrna.add(elems3[1]);
                        }
                    }
                    elems3 = tf3.readLineElems(TextFile.tab);
                }
                tf3.close();
                tf4.close();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
