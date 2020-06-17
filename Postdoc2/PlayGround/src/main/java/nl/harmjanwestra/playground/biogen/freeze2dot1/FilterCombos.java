package nl.harmjanwestra.playground.biogen.freeze2dot1;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class FilterCombos {


    public static void main(String[] args) {

        String allowedGenes = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.proteincoding.txt.gz";
        String allowedSNPs = null;

        // bios eqtl combos
        String input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisBIOS\\010517-BIOS-independent-eQTLs-MetaBrain2dot1IDs.txt";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisBIOS\\010517-BIOS-independent-eQTLs-MetaBrain2dot1IDs-proteincoding.txt";
        FilterCombos f = new FilterCombos();

        try {
//            f.filter(input, 0, 1, allowedSNPs, allowedGenes, output);
            // filter gene list
            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisBIOS\\010517-BIOS-testedGenes-MetaBrain2dot1IDs.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisBIOS\\010517-BIOS-testedGenes-MetaBrain2dot1IDs-intersectMetaBrain2dot1.txt";
            f.filterList(input, allowedGenes, output);
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    private void filterList(String input, String allowedSet, String output) throws IOException {
        HashSet<String> geneset = loadSet(allowedSet);
        TextFile tf = new TextFile(input, TextFile.R);
        TextFile tfo = new TextFile(output, TextFile.W);

        String ln = tf.readLine();
        HashSet<String> visited = new HashSet<>();
        while (ln != null) {
            if (geneset.contains(ln)) {
                if (!visited.contains(ln)) {
                    tfo.writeln(ln);
                    visited.add(ln);
                }

            }
            ln = tf.readLine();
        }

        tf.close();
        tfo.close();

    }

    public void filter(String input, int snpcol, int genecol, String snpfilter, String genefilter, String output) throws IOException {

        HashSet<String> geneset = loadSet(genefilter);
        HashSet<String> snpset = loadSet(snpfilter);

        TextFile tfo = new TextFile(output, TextFile.W);
        TextFile tf = new TextFile(input, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String gene = elems[genecol];
            String snp = elems[snpcol];
            boolean snpok = (snpset == null || snpset.contains(snp));
            boolean geneok = (geneset == null || geneset.contains(gene));
            if (snpok && geneok) {
                tfo.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tfo.close();
        tf.close();

    }

    private HashSet<String> loadSet(String filterfile) throws IOException {
        if (filterfile == null) {
            return null;
        } else {
            HashSet<String> output = new HashSet<>();
            TextFile tf = new TextFile(filterfile, TextFile.R);
            String ln = tf.readLine();
            while (ln != null) {
                output.add(ln);
                ln = tf.readLine();
            }
            tf.close();
            return output;
        }
    }

}
