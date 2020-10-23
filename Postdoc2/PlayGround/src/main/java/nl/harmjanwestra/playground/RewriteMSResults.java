package nl.harmjanwestra.playground;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;

public class RewriteMSResults {

    public static void main(String[] args) {

        try {
            TextFile tf = new TextFile("U:\\IEUGWAS\\2020-08-MSGWAS\\discovery_metav3.0.meta.gz", TextFile.R);
            String[] elems = tf.readLineElems(Strings.whitespace);
            HashMap<String, String> replacements = new HashMap<>();
            while (elems != null) {
                String chr = elems[0];
                String pos = elems[1];
                String a1 = elems[3];
                String a2 = elems[4];
                String id = chr + ":" + pos + ":" + a1 + ":" + a2;
                replacements.put(id, null);
                elems = tf.readLineElems(Strings.whitespace);
            }
            tf.close();

            System.out.println(replacements.size() + " keys");

            TextFile tf2 = new TextFile("D:\\Sync\\SyncThing\\Data\\Ref\\dbsnp\\All_b37_b151_20180418.vcf.gz", TextFile.R);
            elems = tf2.readLineElems(TextFile.tab);
            int nrRepl = 0;
            int ctr = 0;
            while (elems != null) {
                if (elems.length > 4) {
                    String chr = elems[0];
                    String pos = elems[1];
                    String rs = elems[2];
                    String a1 = elems[3];
                    String a2 = elems[4];
                    if (a1.length() < 2 && a2.length() < 2) {
                        String id = chr + ":" + pos + ":" + a1 + ":" + a2;
                        if (replacements.containsKey(id)) {
                            replacements.put(id, rs);
                            nrRepl++;
                        } else {
                            id = chr + ":" + pos + ":" + a2 + ":" + a1;
                            if (replacements.containsKey(id)) {
                                replacements.put(id, rs);
                                nrRepl++;
                            }
                        }
                    }
                    ctr++;
                    if (ctr % 1000000 == 0) {
                        System.out.println(ctr + " lines, " + nrRepl + " replacements");
                    }
                }

                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            tf = new TextFile("U:\\IEUGWAS\\2020-08-MSGWAS\\discovery_metav3.0.meta.gz", TextFile.R);
            TextFile tfo = new TextFile("U:\\IEUGWAS\\2020-08-MSGWAS\\discovery_metav3.0.meta-wRSId.gz", TextFile.W);
            elems = tf.readLineElems(Strings.whitespace);
            tfo.writeln(Strings.concat(elems, Strings.tab) + "\tRSId");
            while (elems != null) {
                String chr = elems[0];
                String pos = elems[1];
                String a1 = elems[3];
                String a2 = elems[4];
                String id = chr + ":" + pos + ":" + a1 + ":" + a2;
                String repl = replacements.get(id);
                if (repl == null) {
                    repl = "-";
                }
                tfo.writeln(Strings.concat(elems, Strings.tab) + "\t" + repl);
                elems = tf.readLineElems(Strings.whitespace);
            }
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }


}
