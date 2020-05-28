package nl.harmjanwestra.playground.biogen.freeze2dot1;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.playground.biogen.locusdebug.FilterCorrelationBetweenGene;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Freeze2Comparison {

    public static void main(String[] args) {

        String freeze2cis = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\Freeze2CisEQTLs-MetaBrain2dot1IDs.txt.gz";
        String freeze2trans = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\Freeze2TransEQTLs-MetaBrain2dot1IDs.txt.gz";
        String freeze2dot1cis = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\freeze2dot1-cis\\eQTLsFDR-ProbeLevel.txt.gz";
        String freeze2dot1trans = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\freeze2dot1-trans\\eQTLsFDR.txt.gz";
        String patchgenes = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-PatchSequenceIssue\\2019-09-26-patch_genes.txt";
        String outputcis = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\PatchGenesVsNoPatchGenes-cis.png";
        String outputcistxt = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\PatchGenesVsNoPatchGenes-cis.txt";
        String outputtrans = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\PatchGenesVsNoPatchGenes-trans.png";
        String outputtranstxt = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\PatchGenesVsNoPatchGenes-trans.txt";


        Freeze2Comparison c = new Freeze2Comparison();
        try {
            c.compareToFreeze2(patchgenes, freeze2cis, freeze2dot1cis, outputcis, outputcistxt);
            c.compareToFreeze2(patchgenes, freeze2trans, freeze2dot1trans, outputtrans, outputtranstxt);
        } catch (IOException | DocumentException e) {
            e.printStackTrace();
        }


    }

    public class Gene {
        String name;
        HashSet<String> ensgids = new HashSet<String>();
        boolean hasautosomalcopies = false;
        boolean haspatchcopies = false;
    }

    class EQTL {
        String gene;
        String snp;
        String alleles;
        String alleleAssessed;
        Double z1;
        Double FDR;
    }

    public void compareToFreeze2(String patchgenes, String freeze2, String freeze2dot1, String output, String outputtxt) throws IOException, DocumentException {

        TextFile tf2 = new TextFile(patchgenes, TextFile.R);
        String[] header = tf2.readLineElems(TextFile.tab);

        HashMap<String, Gene> ensgToGene = new HashMap<String, Gene>();
        ArrayList<Gene> genes = new ArrayList<>();
        String[] elems = tf2.readLineElems(TextFile.tab);
        while (elems != null) {

            Gene g = new Gene();
            g.name = elems[0];

            for (int i = 1; i < elems.length; i++) {
                if (!elems[i].equals("-")) {
                    String[] ensgids = elems[i].split(";");

                    for (String s : ensgids) {
                        s = s.split("\\.")[0];
                        g.ensgids.add(s);
                        ensgToGene.put(s, g);
                    }
                    if (i < 23 && g.ensgids.size() > 1) {
                        g.hasautosomalcopies = true;
                    } else if (i >= 23) {
                        g.haspatchcopies = true;
                    }
                }
            }
            genes.add(g);

            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();


        HashMap<String, EQTL> strToEQTL = new HashMap<String, EQTL>();

        TextFile tf = new TextFile(freeze2, TextFile.R);
        tf.readLine();
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            Double z = Double.parseDouble(elems[10]);
            String snp = elems[1];
            String gene = elems[4];

            String combo = snp + "-" + gene;
            EQTL e = new EQTL();
            e.alleleAssessed = elems[9];
            e.alleles = elems[8];
            e.snp = snp;
            e.gene = gene;
            e.z1 = z;
            e.FDR = Double.parseDouble(elems[elems.length - 1]);
            strToEQTL.put(combo, e);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        tf = new TextFile(freeze2dot1, TextFile.R);
        tf.readLine();
        elems = tf.readLineElems(TextFile.tab);

        ArrayList<Double> x1 = new ArrayList<Double>();
        ArrayList<Double> x2 = new ArrayList<Double>();
        ArrayList<Double> y1 = new ArrayList<Double>();
        ArrayList<Double> y2 = new ArrayList<Double>();


        TextFile tfo = new TextFile(outputtxt, TextFile.W);

        while (elems != null) {
            Double z = Double.parseDouble(elems[10]);
            String snp = elems[1];
            String gene = elems[4];

            String combo = snp + "-" + gene;

            EQTL e = strToEQTL.get(combo);
            if (e != null) {
                String alleles = elems[8];
                String assessed = elems[9];
                Boolean flip = BaseAnnot.flipalleles(e.alleles, e.alleleAssessed, alleles, assessed);
                if (flip == null) {
                    // meh
                } else {
                    if (flip) {
                        z *= -1;
                    }

                    Double fdr = Double.parseDouble(elems[elems.length - 1]);
                    String hugo = elems[elems.length - 6];
                    String dotlessgene = gene.split("\\.")[0];
                    Gene g = ensgToGene.get(dotlessgene);
                    if (g != null) {
                        if (g.haspatchcopies) {
                            x2.add(e.z1);
                            y2.add(z);
                            tfo.writeln(snp + "\t" + dotlessgene + "\t" + e.z1 + "\t" + z + "\t" + e.FDR + "\t" + fdr + "\t" + hugo + "\tPATCH");
                        } else if (g.hasautosomalcopies) {
                            x1.add(e.z1);
                            y1.add(z);
                            tfo.writeln(snp + "\t" + dotlessgene + "\t" + e.z1 + "\t" + z + "\t" + e.FDR + "\t" + fdr + "\t" + hugo + "\tNoPatch");
                        } else {
                            x1.add(e.z1);
                            y1.add(z);
                            tfo.writeln(snp + "\t" + dotlessgene + "\t" + e.z1 + "\t" + z + "\t" + e.FDR + "\t" + fdr + "\t" + hugo + "\tMultiAutosomalCopies");
                        }
                    } else {
                        System.out.println("?" + dotlessgene);
                        x1.add(e.z1);
                        y1.add(z);
                        tfo.writeln(snp + "\t" + dotlessgene + "\t" + e.z1 + "\t" + z + "\t" + e.FDR + "\t" + fdr + "\t" + hugo + "\tNoPatch");
                    }
                }

            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        tfo.close();
        double[][] x = new double[2][];
        double[][] y = new double[2][];
        x[0] = Primitives.toPrimitiveArr(x1);
        x[1] = Primitives.toPrimitiveArr(x2);
        y[0] = Primitives.toPrimitiveArr(y1);
        y[1] = Primitives.toPrimitiveArr(y2);


        Grid grid = new Grid(500, 500, 1, 1, 100, 100);
        ScatterplotPanel scp = new ScatterplotPanel(1, 1);
        scp.setData(x, y);
        scp.setPlotElems(true, true);
        scp.setLabels("Freeze2", "Freeze2dot1");
        scp.setDatasetLabels(new String[]{"NonPatchGenes", "PatchGenes"});
        grid.addPanel(scp);
        grid.draw(output);


    }
}
