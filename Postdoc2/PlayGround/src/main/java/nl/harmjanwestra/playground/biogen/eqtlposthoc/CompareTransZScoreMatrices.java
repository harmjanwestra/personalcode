package nl.harmjanwestra.playground.biogen.eqtlposthoc;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class CompareTransZScoreMatrices {


    public static void main(String[] args) {

        String z1 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-07-eqtlgenTransZScorecompare\\eqtlgen\\ZScoreMatrix.txt.gz";
        String n1 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-07-eqtlgenTransZScorecompare\\eqtlgen\\ZScoreMatrixNrSamples.txt.gz";
        String z2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-07-eqtlgenTransZScorecompare\\metabrain\\ZScoreMatrix.txt.gz";
        String n2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-07-eqtlgenTransZScorecompare\\metabrain\\ZScoreMatrixNrSamples.txt.gz";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-07-eqtlgenTransZScorecompare\\snpCompare.txt";
        String output2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-07-eqtlgenTransZScorecompare\\geneCompare.txt";
        String eqtlsD1 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\cis\\eQTLsFDR0.05-ProbeLevel.txt.gz";
        String eqtlsD2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cis\\eQTLDump-sortZ-FDR.txt.gz-Significant-0.05.txt.gz";
        String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz";
        String gwasAssoc = "U:\\IEUGWAS\\2020-06-25-2020-05-03-allTopAssociations-wgwascatalog-wALS-wAlzheimer-wMetaBrain.txt.gz";
        String gwasList = "U:\\IEUGWAS\\2020-06-01-2020-05-03-gwaslist-wgwascatalog-wALS-wAlzheimer-wMetaBrain.txt.gz";

        String snpsToOutput = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-07-eqtlgenTransZScorecompare\\snpsToSelect.txt";

        CompareTransZScoreMatrices comp = new CompareTransZScoreMatrices();
        try {
            comp.runSNPs(z1, n1, z2, n2, output, gwasAssoc, gwasList, eqtlsD1, eqtlsD2, snpsToOutput, gtf, "eQTLGen", "MetaBrain-Cortex-EUR");
//            comp.runGenes(z1, z2, output2, gtf);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    class ZMatLine {
        String snp;
        String alleles;
        String assessed;
        double[] z;
    }

    class Zmat {
        HashMap<String, Integer> geneMap;
        HashMap<String, Integer> snpMap;
        ArrayList<ZMatLine> zMat;
        ArrayList<String> genes;
        ArrayList<String> snps;

        public Zmat(String file) throws IOException {
            TextFile tf = new TextFile(file, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            genes = new ArrayList<String>();
            geneMap = new HashMap<String, Integer>();
            for (int i = 3; i < elems.length; i++) {
                String gene = elems[i].split("\\.")[0];
                gene = gene.split("_")[0];
                genes.add(gene);
                geneMap.put(gene, i - 3);
            }
            elems = tf.readLineElems(TextFile.tab);
            System.out.println(file + "\tgenes: " + genes.size());
            zMat = new ArrayList<>();
            snpMap = new HashMap<>();
            snps = new ArrayList<>();
            int ctr = 0;

            while (elems != null) {
                String[] snpelems = elems[0].split(":");
                String snp = elems[0];
                if (snpelems.length > 1) {
                    snp = snpelems[2];
                }

                String alleles = elems[1];
                String assessed = elems[2];
                double[] z = new double[elems.length - 3];
                for (int i = 3; i < elems.length; i++) {
                    z[i - 3] = Double.parseDouble(elems[i]);
                }
                ZMatLine zml = new ZMatLine();
                zml.snp = snp;
                zml.alleles = alleles;
                zml.assessed = assessed;
                zml.z = z;
                zMat.add(zml);
                snpMap.put(snp, ctr);
                snps.add(snp);
                ctr++;

                if (ctr % 100 == 0) {
                    System.out.print("Loading: " + ctr + "\r");
                }

                elems = tf.readLineElems(TextFile.tab);
            }
            System.out.println();
            tf.close();
            System.out.println(file + "\tsnps: " + snps.size());
        }


    }

    public void runGenes(String zmat1file, String zmat2file, String output, String gtffile) throws IOException {
        GTFAnnotation gtf = new GTFAnnotation(gtffile);

        HashMap<String, String> ensToHGNC = new HashMap<>();
        for (Gene g : gtf.getGenes()) {
            String name = g.getName().split("\\.")[0];
            ensToHGNC.put(name, g.getGeneSymbol());
        }

        Zmat zmat1 = new Zmat(zmat1file);
        Zmat zmat2 = new Zmat(zmat2file);

        TextFile outf = new TextFile(output, TextFile.W);
        outf.writeln("Gene\tsymbol\tn\tr\trsq");
        ProgressBar pb = new ProgressBar(zmat1.genes.size());
        for (String gene : zmat1.genes) {
            Integer gidx1 = zmat1.geneMap.get(gene);
            Integer gidx2 = zmat2.geneMap.get(gene);

            if (gidx1 != null && gidx2 != null) {
                ArrayList<Double> x = new ArrayList<>();
                ArrayList<Double> y = new ArrayList<>();

                for (String snp : zmat1.snps) {
                    Integer idx2 = zmat2.snpMap.get(snp);
                    Integer idx1 = zmat1.snpMap.get(snp);
                    if (idx1 != null && idx2 != null) {
                        ZMatLine zml1 = zmat1.zMat.get(idx1);
                        ZMatLine zml2 = zmat2.zMat.get(idx2);
                        Boolean flip = BaseAnnot.flipalleles(zml1.alleles, zml1.assessed, zml2.alleles, zml2.assessed);
                        if (flip != null) {
                            double[] z1 = zml1.z;
                            double[] z2 = zml2.z;
                            double zv1 = z1[gidx1];
                            double zv2 = z2[gidx2];
                            if (!Double.isNaN(zv1) && !Double.isNaN(zv2)) {
                                x.add(zv1);
                                if (flip) {
                                    y.add(zv2 * -1);
                                } else {
                                    y.add(zv2);
                                }
                            }
                        }
                    }
                }

                double[] px = Primitives.toPrimitiveArr(x);
                double[] py = Primitives.toPrimitiveArr(y);
                if (px.length > 0) {
                    SpearmansCorrelation sp = new SpearmansCorrelation();
                    double r = sp.correlation(px, py);

                    String ln = gene + "\t" + ensToHGNC.get(gene) + "\t" + px.length + "\t" + r + "\t" + (r * r);
                    outf.writeln(ln);
                }
            }
            pb.iterate();
        }
        outf.close();
        pb.close();
    }

    public void runSNPs(String zmat1file, String nmat1file, String zmat2file, String nmat2file, String output, String gwasAssoc,
                        String gwasList, String eqtlfileD1, String eqtlfileD2, String snpfilterfile, String gtffile, String d1name, String d2name) throws IOException {

        GTFAnnotation gtf = new GTFAnnotation(gtffile);

        HashMap<String, String> ensToHGNC = new HashMap<>();
        for (Gene g : gtf.getGenes()) {
            String name = g.getName().split("\\.")[0];
            ensToHGNC.put(name, g.getGeneSymbol());
        }


        HashSet<String> snpsToOutput = null;
        if (snpfilterfile != null) {
            TextFile t1 = new TextFile(snpfilterfile, TextFile.R);
            ArrayList<String> al1 = t1.readAsArrayList();
            snpsToOutput = new HashSet<>();
            snpsToOutput.addAll(al1);
            t1.close();
        }

        HashMap<String, String> eqtlsD1 = loadEQTLs(eqtlfileD1, ensToHGNC);
        HashMap<String, String> eqtlsD2 = loadEQTLs(eqtlfileD2, ensToHGNC);

        HashMap<String, String> gwas = loadGWASAssoc(gwasAssoc, gwasList);


        // load matrices
        Zmat zmat1 = new Zmat(zmat1file);
        Zmat nmat1 = new Zmat(nmat1file);

        Zmat zmat2 = new Zmat(zmat2file);
        Zmat nmat2 = new Zmat(nmat2file);

        // compare
        TextFile outf = new TextFile(output, TextFile.W);
        outf.writeln("SNP\tD1CisEQTLs\tD2CisEQTLs\tn\tSpearman-r\tSpearman-rsq\tZ\tP\tTraitsAssociatedWithSNP");
        ProgressBar pb = new ProgressBar(zmat1.snps.size());
        int maxN = zmat1.genes.size();
        if (zmat2.genes.size() > maxN) {
            maxN = zmat2.genes.size();
        }
        Correlation.correlationToZScore(maxN);
        for (String snp : zmat1.snps) {

            Integer idx1 = zmat1.snpMap.get(snp);
            Integer idxn1 = nmat1.snpMap.get(snp);

            Integer idx2 = zmat2.snpMap.get(snp);
            Integer idxn2 = nmat2.snpMap.get(snp);

            if (idx1 != null && idx2 != null) {
                ZMatLine zml1 = zmat1.zMat.get(idx1);
                ZMatLine zml2 = zmat2.zMat.get(idx2);
                ZMatLine nml1 = nmat1.zMat.get(idxn1);
                ZMatLine nml2 = nmat2.zMat.get(idxn2);
                Boolean flip = BaseAnnot.flipalleles(zml1.alleles, zml1.assessed, zml2.alleles, zml2.assessed);
                if (flip != null) {
                    double[] z1 = zml1.z;
                    double[] z2 = zml2.z;

                    double[] n1 = nml1.z;
                    double[] n2 = nml2.z;

                    ArrayList<Double> x = new ArrayList<>();
                    ArrayList<Double> xn = new ArrayList<>();
                    ArrayList<Double> y = new ArrayList<>();
                    ArrayList<Double> yn = new ArrayList<>();
                    ArrayList<String> genes = new ArrayList<>();
                    for (int g = 0; g < zmat1.genes.size(); g++) {
                        String gene = zmat1.genes.get(g);
                        Integer idxg2 = zmat2.geneMap.get(gene);

                        if (idxg2 != null) {
                            double zv1 = z1[g];
                            double zv2 = z2[idxg2];
                            if (!Double.isNaN(zv1) && !Double.isNaN(zv2)) {

                                // also add sample size
                                Integer idxng1 = nmat1.geneMap.get(gene);
                                Integer idxng2 = nmat2.geneMap.get(gene);

                                xn.add(n1[idxng1]);
                                yn.add(n2[idxng2]);

                                genes.add(gene);
                                x.add(zv1);
                                if (flip) {
                                    y.add(zv2 * -1);
                                } else {
                                    y.add(zv2);
                                }
                            }
                        }
                    }

                    double[] px = Primitives.toPrimitiveArr(x);
                    double[] py = Primitives.toPrimitiveArr(y);
                    if (snpsToOutput != null && snpsToOutput.contains(snp)) {


                        TextFile outf2 = new TextFile(output + "-" + snp + ".txt", TextFile.W);
                        outf2.writeln("D1: " + d1name + "\tD2: " + d2name);
                        outf2.writeln(snp + "\tAllelesD1: " + zml1.alleles + ", Assessed: " + zml1.assessed + "\tAllelelesD2: " + zml2.alleles + ", Assessed: " + zml2.assessed + "\tFlipped: " + flip);

                        String header = "Gene\tHGNC\tD1Z\tD1N\tD2Z\tD2N";
                        outf2.writeln(header);
                        for (int i = 0; i < genes.size(); i++) {
                            String gene = genes.get(i);
                            String ln = gene + "\t" + ensToHGNC.get(gene) + "\t" + px[i] + "\t" + xn.get(i) + "\t" + py[i] + "\t" + yn.get(i);
                            outf2.writeln(ln);
                        }
                        outf2.close();
                    }
                    if (px.length > 0) {
                        SpearmansCorrelation sp = new SpearmansCorrelation();

                        double r = sp.correlation(px, py);
                        double z = Correlation.convertCorrelationToZScore(px.length, r);
                        double p = ZScores.zToP(z);
                        String traits = gwas.get(snp);
                        String ln = snp + "\t" + eqtlsD1.get(snp) + "\t" + eqtlsD2.get(snp) + "\t" + px.length + "\t" + r + "\t" + (r * r) + "\t" + z + "\t" + p + "\t" + traits;
                        outf.writeln(ln);
                    }
                }
            }
            pb.iterate();
        }
        outf.close();
        pb.close();
    }

    private HashMap<String, String> loadEQTLs(String eqtlfile, HashMap<String, String> ensgToHGNC) throws IOException {
        TextFile tf = new TextFile(eqtlfile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        HashMap<String, ArrayList<String>> genes = new HashMap<>();
        while (elems != null) {

            String snp = "";
            String[] snpelems = elems[1].split(":");
            if (snpelems.length < 3) {
                snp = snpelems[0];
            } else {
                snp = snpelems[2];
            }

            String gene = ensgToHGNC.get(elems[4].split("\\.")[0]);
            ArrayList<String> g = genes.get(snp);
            if (g == null) {
                g = new ArrayList<>();
            }
            g.add(gene);
            genes.put(snp, g);
            elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();
        HashMap<String, String> output = new HashMap<>();
        for (String key : genes.keySet()) {
            ArrayList<String> g = genes.get(key);
            output.put(key, Strings.concat(g, Strings.semicolon));
        }
        return output;
    }

    private HashMap<String, String> loadGWASAssoc(String gwasAssoc, String gwasList) throws IOException {
        TextFile tf = new TextFile(gwasList, TextFile.R);
        HashMap<String, String> idToTrait = new HashMap<String, String>();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String id = elems[0];
            String trait = elems[3];
            idToTrait.put(id, trait);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(gwasAssoc, TextFile.R);
        tf2.readLine();
        elems = tf2.readLineElems(TextFile.tab);
        HashMap<String, HashSet<String>> output = new HashMap<>();
        while (elems != null) {
            String id = elems[0];
            String[] snpelems = elems[1].split(":");
            String snp = elems[1];
            if (snpelems.length > 1) {
                snp = snpelems[2];
            }

            HashSet<String> traits = output.get(snp);
            if (traits == null) {
                traits = new HashSet<>();
            }
            traits.add(idToTrait.get(id));
            output.put(snp, traits);
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();

        HashMap<String, String> output2 = new HashMap<>();
        for (String key : output.keySet()) {
            HashSet<String> vals = output.get(key);
            ArrayList<String> alltraits = new ArrayList<>();
            alltraits.addAll(vals);
            output2.put(key, Strings.concat(alltraits, Strings.semicolon));
        }
        return output2;

    }
}
