package nl.harmjanwestra.playground.cis;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class FilterTopFX {
    public static void main(String[] args) {

        String in = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38.txt.gz";
        String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38-topfxpergene.txt.gz";

        try {
            TextFile tf = new TextFile(in, TextFile.R);
            TextFile tfo = new TextFile(out, TextFile.W);
            tfo.writeln(tf.readLine());
            HashSet<String> geneseen = new HashSet<String>();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String gene = elems[4];
                if (!geneseen.contains(gene)) {
                    tfo.writeln(Strings.concat(elems, Strings.tab));
                    geneseen.add(gene);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            tfo.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
