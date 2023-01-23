package nl.harmjanwestra.playground.biogen.downstreamer;

import com.google.common.util.concurrent.AtomicDouble;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class OverlapCoregScores {


    public static void main(String[] args) {
        String zmat = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\2020-05-26-Cortex-EUR-AFR-noENA-noPCA\\ZScoreMatrix-strip.txt.gz";
        String scoref = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MSFigure\\Downstreamer_MS_metabrain-coregulationscores.txt";
        scoref = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MSFigure\\Downstreamer_MS_metabrain-cortex-coregulationscores.txt";

        String snpfil = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MSFigure\\metabrain-MSGWAS2019SNPs.txt";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MSFigure\\2020-05-26-Cortex-EUR-AFR-noENA-noPCA\\MSGWASSNPs-CorrelWDownstreamer-cortexCoregscores.txt";

        String[] snpfiles = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\metabrain-MSGWAS2019SNPs.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\2020-ALSGWAS-SNPS.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\2020-ALSGWAS-SNPS.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\Ripke2014-SCZ-SNPs.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\Ripke2014-SCZ-SNPs.txt"
        };

        String[] scorefs = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\Downstreamer_MS_metabrain_snRNASeq-coregulationscores.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\Downstreamer-ALS_metabrain_snRNASeq-coregulationscores.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\Downstreamer-ALS_metabrain-coregulationscores.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\Downstreamer-SCZ_metabrain_snRNASeq-coregulationscores.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\Downstreamer-SCZ_metabrain-coregulationscores.txt"
        };

        String[] outputs = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\2020-05-26-Cortex-EUR-AFR-noENA-noPCA\\MSGWASSNPs-CorrelWDownstreamer-snRNAseq-Coregscores.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\2020-05-26-Cortex-EUR-AFR-noENA-noPCA\\ALSGWASSNPs-CorrelWDownstreamer-snRNAseq-Coregscores.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\2020-05-26-Cortex-EUR-AFR-noENA-noPCA\\ALSGWASSNPs-CorrelWDownstreamer-Coregscores.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\2020-05-26-Cortex-EUR-AFR-noENA-noPCA\\SCZGWASSNPs-CorrelWDownstreamer-snRNAseq-Coregscores.txt",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-12-01-MS-SCZ-ALS\\2020-05-26-Cortex-EUR-AFR-noENA-noPCA\\SCZGWASSNPs-CorrelWDownstreamer-Coregscores.txt"
        };

        int snpcol = 1;

        OverlapCoregScores r = new OverlapCoregScores();
        try {
            for (int i = 0; i < snpfiles.length; i++) {
                r.run(zmat, scorefs[i], snpfiles[i], snpcol, outputs[i]);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void run(String zmatrix, String scorefile, String snpfilter, int snpcol, String output) throws Exception {


        TextFile tfa = new TextFile(snpfilter, TextFile.R);
        HashSet<String> allowedSNPs = (HashSet<String>) tfa.readAsSet(snpcol, TextFile.tab);

        System.out.println(allowedSNPs + " snp ids");

        HashMap<String, Double> scores = new HashMap<String, Double>();
        TextFile tf = new TextFile(scorefile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            scores.put(elems[0], Double.parseDouble(elems[1]));
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(scores.size() + " scores");

        // correlate rows

        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(zmatrix);
        System.out.println("Done");

        TextFile tfo1 = new TextFile(output + "-perSNP.txt", TextFile.W);
        tfo1.writeln("SNP\tVals\tCorrelZ\tCorrelZSquare\tCorrelZSquareZSquare");
        ProgressBar pb = new ProgressBar(ds.rows(), "Iterating rows");
        ProgressBar finalPb = pb;
        IntStream.range(0, ds.rows()).parallel().forEach(r -> {
            String name = ds.getRowObjects().get(r);
            if (allowedSNPs.contains(name)) {
                ArrayList<Double> x = new ArrayList<>();
                ArrayList<Double> x2 = new ArrayList<>();
                ArrayList<Double> y = new ArrayList<>();
                ArrayList<Double> y2 = new ArrayList<>();

                for (int c = 0; c < ds.columns(); c++) {
                    String gene = ds.getColObjects().get(c).split("\\.")[0];
                    Double score = scores.get(gene);
                    if (score != null && !Double.isNaN(score)) {
                        Double val = ds.getElement(r, c);
                        if (!Double.isNaN(val)) {
                            y.add(val);
                            y2.add(val * val);
                            x.add(score);
                            x2.add(score * score);
                        }
                    }
                }
                SpearmansCorrelation sp = new SpearmansCorrelation();
                double corr = sp.correlation(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
                double corr2 = sp.correlation(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y2));
                double corr3 = sp.correlation(Primitives.toPrimitiveArr(x2), Primitives.toPrimitiveArr(y2));
                try {
                    tfo1.writelnsynced(name + "\t" + x.size() + "\t" + corr + "\t" + corr2 + "\t" + corr3);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            finalPb.iterateSynched();
        });
        pb.close();


        // calculate chi-squareds
        TextFile tfo2 = new TextFile(output + "-perGene.txt", TextFile.W);
        tfo2.writeln("Gene\tScore\tNrZScores\tSumchiSq\tAvgSumchiSq");
        ArrayList<Double> x = new ArrayList<>();
        ArrayList<Double> x2 = new ArrayList<>();
        ArrayList<Double> y = new ArrayList<>();
        ArrayList<Double> y2 = new ArrayList<>();
        pb = new ProgressBar(ds.columns(), "Iterating columns");
        for (int c = 0; c < ds.columns(); c++) {
            String gene = ds.getColObjects().get(c).split("\\.")[0];
            Double score = scores.get(gene);
            AtomicDouble sumxi = new AtomicDouble();
            AtomicInteger nrvals = new AtomicInteger();
            if (score != null && !Double.isNaN(score)) {
                int finalC = c;
                IntStream.range(0, ds.rows()).parallel().forEach(r -> {
                    String name = ds.getRowObjects().get(r);
                    if (allowedSNPs.contains(name)) {
                        Double val = ds.getElement(r, finalC);
                        if (!Double.isNaN(val)) {
                            val = (val * val);
                            sumxi.getAndAdd(val);
                            nrvals.getAndIncrement();
                        }
                    }
                });
            }

            if (score != null) {
                x.add(score);
                x2.add(score * score);
                y.add(sumxi.get());
                y2.add(sumxi.get() / nrvals.get());
                tfo2.writeln(gene + "\t" + score + "\t" + nrvals + "\t" + sumxi + "\t" + (sumxi.get() / nrvals.get()));
            }
            pb.set(c);

        }
        pb.close();

        SpearmansCorrelation sp = new SpearmansCorrelation();
        double corr = sp.correlation(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
        tfo1.writeln("SumChiSq\t" + x.size() + "\t" + corr + "\t-\t-");

        corr = sp.correlation(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y2));
        tfo1.writeln("SumChiSqAverage\t" + x.size() + "\t" + corr + "\t-\t-");
        tfo1.close();
        tfo2.close();

    }
}