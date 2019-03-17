package nl.harmjanwestra.playground.cis;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

public class SplitCisEQTLsByGene {

    public static void main(String[] args) {
        SplitCisEQTLsByGene c = new SplitCisEQTLsByGene();

        String efile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.txt.gz";
        String efilesort = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-sorted.txt";
        String outdir = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\splitpergene\\";
        String genefile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\sortedGeneSNPCombos.txt.gz-genes.txt";

        try {
            c.run(efile, efilesort, genefile, outdir);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void run(String eqtlfile, String eqtlfilesorted, String genefile, String outdir) throws IOException {


        TextFile tf1 = new TextFile(eqtlfile, TextFile.R);
        String header = tf1.readLine();
        tf1.close();
        TextFile tf = new TextFile(eqtlfilesorted, TextFile.R);
        String ln = tf.readLine();
        String prevgene = null;
        TextFile out = null;
        int ctr = 0;
        while (ln != null) {
            String[] subelems = Strings.subsplit(ln, Strings.tab, 4, 5);
            String gene = subelems[0];

            if (prevgene == null || !gene.equals(prevgene)) {
                String outfile = outdir + gene + ".txt.gz";
                if (out != null) {
                    out.close();
                }
                out = new TextFile(outfile, TextFile.W);
                out.writeln(header);
                out.writeln(ln);
            } else {
                out.writeln(ln);

            }
            prevgene = gene;
            ctr++;
            if (ctr % 10000 == 0) {
                System.out.print("\r" + gene + "\t" + ctr);
            }
            ln = tf.readLine();
        }
        if (out != null) {
            out.close();
        }
        tf.close();

    }

    private void runGene(String eqtlfile, String outdir, String v, AtomicInteger g, ProgressBar pb) throws IOException {
        TextFile tf = new TextFile(eqtlfile, TextFile.R);
        String header = tf.readLine();

        String outfile = outdir + v + ".txt.gz";
        TextFile tfo = new TextFile(outfile, TextFile.W);
        tfo.writeln(header);
        String ln = tf.readLine();
        int lnctr = 0;
        DecimalFormat formatter = (DecimalFormat) NumberFormat.getInstance(Locale.US);
        while (ln != null) {

            String[] subelems = Strings.subsplit(ln, Strings.tab, 4, 5);
            String gene = subelems[0];
            if (gene.equals(v)) {
                tfo.writeln(ln);
            }

            lnctr++;
            if (lnctr % 1000000 == 0) {

                System.out.print("\r" + v + "\t" + ((double) lnctr / 150000000));
            }

            ln = tf.readLine();
        }
        System.out.println();
        System.out.println(v + "\t" + lnctr + "\tdone");
        tfo.close();
        tf.close();
        g.getAndIncrement();
        int q = g.get();
        pb.set(q);
    }

}
