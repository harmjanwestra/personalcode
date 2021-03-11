package nl.harmjanwestra.playground.biogen.cellcounts;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;

public class MergeCellCountCorrelationsWithTransEQTLs {

    public static void main(String[] args) {

        String[] snp = {"7:12244161:rs1990622:A_G", "7:12223752:rs11974335:G_T", "7:12225245:rs10950398:G_A"};
        boolean[] flip = {false, true, false};
        String eqtl = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\trans\\eQTLs-crossMappingEQTLsRemoved-FDR0.05.txt";
        String correlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\cellcounts\\Cortex-EUR\\correlations\\cortex_eur_trans_coefficients.txt.gz";
        String celltype = "CellMapNNLS_Neuron";
        MergeCellCountCorrelationsWithTransEQTLs c = new MergeCellCountCorrelationsWithTransEQTLs();
        for (int i = 0; i < snp.length; i++) {
            String snpf = snp[i].replaceAll(":", "_");
            String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\cellcounts\\Cortex-EUR\\correlations\\" + snpf + ".txt";
            try {
                c.run(eqtl, correlfile, snp[i], celltype, flip[i], out);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void run(String eqtl, String correlfile, String snpquery, String celltype, boolean flip, String outf) throws IOException {

        TextFile tf = new TextFile(correlfile, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        String[] ctelems = tf.readLineElems(TextFile.tab);
        HashMap<String, Double> ctToCorrel = new HashMap<>();
        while (ctelems != null) {
            String ct = ctelems[0];
            if (ct.equals(celltype)) {
                for (int i = 1; i < ctelems.length; i++) {
                    ctToCorrel.put(header[i], Double.parseDouble(ctelems[i]));
                }
            }
            ctelems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        tf = new TextFile(eqtl, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        TextFile out = new TextFile(outf, TextFile.W);
        out.writeln("Gene\tZ\tCellTypeCorrel");
        while (elems != null) {
            String snp = elems[1];
            String gene = elems[4];
            if (snp.equals(snpquery)) {
                double z = Double.parseDouble(elems[10]);
                if (flip) {
                    z = -z;
                }
                double other = ctToCorrel.get(gene);
                out.writeln(gene + "\t" + z + "\t" + other);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        out.close();

    }
}
