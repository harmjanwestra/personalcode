package nl.harmjanwestra.playground.biogen.freeze2dot1.eqtlgen;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.playground.legacy.vcf.VCFTabix;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import org.broad.tribble.readers.TabixReader;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class CompareMAF {


    public static void main(String[] args) {
        String eqtlgen = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz";
        String metabrain = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-21-SNPSummaryStats\\2020-04-13-eqtls-rsidfix-cortex-cis-EUR-SNPQCLog-MAF.txt.gz";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-21-SNPSummaryStats\\2020-04-13-eqtls-rsidfix-cortex-cis-EUR-SNPQCLog-MAF-compareToEQTLgen.txt";

        String kgprefix = "";
        String kgsamples = "";

        CompareMAF n = new CompareMAF();
        try {
            n.run(eqtlgen, metabrain, output);
//            n.compare1kg(metabrain, kgprefix, kgsamples, output);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void compare1kg(String metabrain, String kgprefix, String kgsamples, String output) throws IOException, DocumentException {

        VCFTabix tabix = new VCFTabix(kgprefix);
        boolean[] filter = tabix.getSampleFilter(kgsamples);


        TextFile tf = new TextFile(metabrain, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        String[] elems = tf.readLineElems(TextFile.tab);

        ArrayList<Double> x = new ArrayList<Double>();
        ArrayList<Double> y = new ArrayList<Double>();

        TextFile out = new TextFile(output, TextFile.W);
        out.writeln("ID\tMAFMetaBrain\tMAF1kg\tFlipped\tequalid");
        while (elems != null) {

            String id = elems[0];

            String[] idelems = id.split(":");
            String metabrainrs = idelems[2];
            String[] metabrainalleles = idelems[idelems.length - 1].split("_");
            double metabrainmaf = Double.parseDouble(elems[elems.length - 1]);
            String minor = elems[elems.length - 2];
            int pos = Integer.parseInt(idelems[1]);
            SNPFeature feature = new SNPFeature(Chromosome.parseChr(idelems[0]), pos, pos + 1);

            TabixReader.Iterator variants = tabix.query(feature);
            String variant = variants.next();
            while (variant != null) {
                VCFVariant snp = new VCFVariant(variant, VCFVariant.PARSE.ALL, filter);
                String[] alleles = snp.getAlleles();
                double maf = snp.getMAF();

                int ct = compareAlleles(alleles, metabrainalleles);
                if (ct == 2) {
                    String minorallele = snp.getMinorAllele();
                    boolean flipped = false;
                    String rsid = snp.getId();
                    boolean equalid = metabrainrs.equals(rsid);
                    if (!minorallele.equals(minor)) {
                        // alleles are flipped
                        metabrainmaf = 1 - metabrainmaf;
                        flipped = true;
                    }
                    x.add(metabrainmaf);
                    y.add(maf);
                    out.writeln(id + "\t" + metabrainmaf + "\t" + maf + "\t" + flipped + "\t" + equalid);
                    break;
                }
                variant = variants.next();
            }

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        out.close();

        Grid g = new Grid(500, 500, 1, 1, 50, 50);
        ScatterplotPanel p = new ScatterplotPanel(1, 1);
        p.setData(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
        p.setAlpha(0.1f);
        p.setPlotElems(true, false);
        p.setDataRange(new Range(0, 0, 1.0, 1.0));
        g.addPanel(p);
        g.draw(output + ".png");
    }

    private int compareAlleles(String[] alleles, String[] metabrainalleles) {
        HashSet<String> allelesset = new HashSet<String>();
        for (String s : alleles) {
            allelesset.add(s);
        }
        int ctr = 0;
        for (String s : metabrainalleles) {
            if (allelesset.contains(s)) {
                ctr++;
            }
        }
        return ctr;

    }


    HashMap<String, double[]> snps = new HashMap<String, double[]>();

    public void run(String eqtlgen, String metabrain, String output) throws IOException {


        TextFile tf = new TextFile(eqtlgen, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String rsid = elems[0];
            Double maf = Double.NaN;
            try {
                maf = Double.parseDouble(elems[elems.length - 1]);
            } catch (NumberFormatException e) {

            }
            double[] d = new double[2];
            if (!Double.isNaN(maf)) {
                if (maf > 0.5) {
                    maf = 1 - maf;
                }
            }
            d[0] = maf;
            d[1] = Double.NaN;
            snps.put(rsid, d);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(metabrain, TextFile.R);
        tf2.readLine();
        elems = tf2.readLineElems(TextFile.tab);

        while (elems != null) {
            String id = elems[0];
            Double maf = Double.parseDouble(elems[elems.length - 1]);
            if (maf >= 0.01) {
                String[] idelems = id.split(":");
                String rsid = idelems[2];
                if (rsid.equals("nors")) {
                    rsid = id;
                }
                double[] d = snps.get(rsid);
                if (d == null) {
                    d = new double[2];
                    d[0] = Double.NaN;
                }
                d[1] = maf;
                snps.put(rsid, d);
            }
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();

        TextFile out = new TextFile(output, TextFile.W);
        TextFile outshared = new TextFile(output + "-shared.txt", TextFile.W);
        out.writeln("ID\tMAFeqtlgen\tMAFmetabrain");
        outshared.writeln("ID\tMAFeqtlgen\tMAFmetabrain");
        int shared = 0;
        int sharedMAFLt005 = 0;
        int sharedmaflt001 = 0;

        int missing2 = 0;
        int missing2MAFLt005 = 0;
        int missing2maflt001 = 0;

        int missing1 = 0;
        int missing1MAFLt005 = 0;
        int missing1maflt001 = 0;

        int[] missingbins = new int[100];
        int[] missingbins2 = new int[100];

        for (String key : snps.keySet()) {
            double[] d = snps.get(key);

            double maf1 = d[0];
            double maf2 = d[1];
            if (!Double.isNaN(maf1) && !Double.isNaN(maf2)) {
                shared++;
                if (maf1 < 0.05) {
                    sharedMAFLt005++;
                }
                if (maf1 < 0.01) {
                    sharedmaflt001++;
                }
                outshared.writeln(key + "\t" + d[0] + "\t" + d[1]);
            } else if (!Double.isNaN(maf1) && Double.isNaN(maf2)) {
                missing2++;

                int bin = (int) Math.floor(maf1 * 2 * 100);
                if (bin == missingbins.length) {
                    bin = bin - 1;
                }
                missingbins[bin]++;


                if (maf1 < 0.05) {
                    missing2MAFLt005++;
                }
                if (maf1 < 0.01) {
                    missing2maflt001++;
                }
            } else if (Double.isNaN(maf1) && !Double.isNaN(maf2)) {
                int bin = (int) Math.floor(maf2 * 2 * 100);
                if (bin == missingbins.length) {
                    bin = bin - 1;
                }
                missingbins2[bin]++;

                missing1++;
                if (maf2 < 0.05) {
                    missing1MAFLt005++;
                }
                if (maf2 < 0.01) {
                    missing1maflt001++;
                }
            }

            out.writeln(key + "\t" + d[0] + "\t" + d[1]);
        }
        out.close();
        outshared.close();
        System.out.println(shared + "\t" + sharedmaflt001 + "\t" + sharedMAFLt005);
        System.out.println(missing1 + "\t" + missing1maflt001 + "\t" + missing1MAFLt005);
        System.out.println(missing2 + "\t" + missing2maflt001 + "\t" + missing2MAFLt005);

        System.out.println();

        DecimalFormat format = new DecimalFormat("#.###");
        for (int i = 0; i < missingbins.length; i++) {
            double perc = (double) i / 100 / 2;
            System.out.println(format.format(perc) + "\t" + missingbins[i] + "\t" + missingbins2[i]);
        }

    }

}
