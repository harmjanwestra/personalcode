package nl.harmjanwestra.playground.cis.perm;


import JSci.maths.ArrayMath;
import JSci.maths.statistics.BetaDistribution;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.facebook.presto.hive.$internal.jodd.io.StreamGobbler;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.Gene;
import org.w3c.dom.Text;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Binning;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;
import umcg.genetica.util.RunTimer;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.IntStream;

public class FDRFromMaxAbsZ {

    public static void main(String[] args) {
        String perm = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-03-13-FDRFromMaxAbsPermutedZ\\2019-03-12-MaxAbsPermutedZScorePerGene.txt.gz";
        String real = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-03-13-FDRFromMaxAbsPermutedZ\\2019-03-12-MaxAbsZScoresRealData.txt.gz";
        String maxPerPerm = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-03-13-FDRFromMaxAbsPermutedZ\\2019-03-12-MaxAbsZScoresPermData.txt.gz";
        String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-03-13-FDRFromMaxAbsPermutedZ\\beta\\pvals";

//        String perm = "/mnt/d/Sync/SyncThing/Postdoc2/2019-eQTLMeta/data/2019-03-13-FDRFromMaxAbsPermutedZ/2019-03-12-MaxAbsPermutedZScorePerGene.txt.gz";
//        String real = "/mnt/d/Sync/SyncThing/Postdoc2/2019-eQTLMeta/data/2019-03-13-FDRFromMaxAbsPermutedZ/2019-03-12-MaxAbsZScoresRealData.txt.gz";
//        String maxPerPerm = "/mnt/d/Sync/SyncThing/Postdoc2/2019-eQTLMeta/data/2019-03-13-FDRFromMaxAbsPermutedZ/2019-03-12-MaxAbsZScoresPermData.txt.gz";
//        String out = "/mnt/d/Sync/SyncThing/Postdoc2/2019-eQTLMeta/data/2019-03-13-FDRFromMaxAbsPermutedZ/2019-03-12-MaxAbsZScoresRealDataAndPermMerged.txt";
//        String scriptloc = "/mnt/d/Sync/SyncThing/Postdoc2/2019-eQTLMeta/data/2019-03-13-FDRFromMaxAbsPermutedZ/betadist.py";

        FDRFromMaxAbsZ p = new FDRFromMaxAbsZ();
        try {
//            p.run(real, perm, maxPerPerm, scriptloc, out);
//            p.betaDist(real, perm, maxPerPerm, out);

            String annot = "D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
            String outtss = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-03-13-FDRFromMaxAbsPermutedZ\\2019-03-12-MaxAbsZScoresRealDataWithTSS.txt";
//            p.calcTSS(real, annot, outtss);

            String outbetap = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-03-13-FDRFromMaxAbsPermutedZ\\2019-03-12-MaxAbsZScoresRealDataWithTSS-betaP.txt";
            String betapvalloc = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2019-03-13-FDRFromMaxAbsPermutedZ\\beta\\";
            p.mergeBetaPVals(outtss, perm, betapvalloc, outbetap);


        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public void calcTSS(String real, String annot, String out) throws IOException {
        GTFAnnotation gtf = new GTFAnnotation(annot);

        TextFile in = new TextFile(real, TextFile.R);
        TextFile of = new TextFile(out, TextFile.W);
        String header = in.readLine() + "\tTSSDist";
        of.writeln(header);
        String[] elems = in.readLineElems(TextFile.tab);

        while (elems != null) {
            String gene = elems[0];

            Gene geneobj = gtf.getStrToGene().get(gene);
            if (geneobj != null) {
                Chromosome chr = Chromosome.parseChr(elems[2]);
                Integer pos = Integer.parseInt(elems[3]);

                int dist = 0;
                if (geneobj.getStrand().equals(Strand.NEG)) {

                    dist = geneobj.getStop() - pos;

                } else {
                    dist = geneobj.getStart() - pos;

                }
                of.writeln(Strings.concat(elems, Strings.tab) + "\t" + dist);
            } else {
                of.writeln(Strings.concat(elems, Strings.tab) + "\tNaN");
            }
            elems = in.readLineElems(TextFile.tab);
        }
        in.close();
        of.close();
    }

    public void mergeBetaPVals(String real, String perm, String betapvalloc, String out) throws IOException {

        DoubleMatrixDataset<String, String> dsPerm = DoubleMatrixDataset.loadDoubleData(perm);

        ArrayList<HashMap<String, String>> geneToBetaP = new ArrayList();
        for (int i = -1; i < 10; i++) {
            String file = betapvalloc + "betapvals-real.txt";
            if (i >= 0) {
                file = betapvalloc + "betapvals-perm" + i + ".txt";
            }
            HashMap<String, String> pvals = new HashMap<>();
            TextFile f1 = new TextFile(file, TextFile.R);
            String[] elems = f1.readLineElems(TextFile.tab);
            while (elems != null) {
                pvals.put(elems[0], elems[2]);
                elems = f1.readLineElems(TextFile.tab);
            }
            f1.close();

            geneToBetaP.add(pvals);
        }


        TextFile tfout = new TextFile(out, TextFile.W);
        TextFile tf = new TextFile(real, TextFile.R);

        String bonusheader = "";
        for (int p = -1; p < 10; p++) {
            if (p == -1) {
                bonusheader += "\tPBetaDist-real";
            } else {
                bonusheader += "\tPBetaDist-perm" + p;
            }
        }
        String header = tf.readLine()+ "\tPzGivenPerm\tNZBiggerInPermutations\tNPermutedZTotal" + bonusheader;
        tfout.writeln(header);

        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String gene = elems[0];
            String[] pvalsBeta = new String[geneToBetaP.size()];
            for (int i = 0; i < geneToBetaP.size(); i++) {
                HashMap<String, String> zs = geneToBetaP.get(i);
                String p = zs.get(gene);
                pvalsBeta[i] = p;
            }


            DoubleMatrix1D row = dsPerm.getRow(gene);
            double[] vals = filterNaNAndSort(row);
            double[] pvals = converToP(vals);
            double z = Double.parseDouble(elems[4]);

            int n = 0;
            for (int d = 0; d < vals.length; d++) {
                if (vals[d] >= z) {
                    n++;
                }
            }
            double pderived = (double) n / vals.length;

            double realp = ZScores.zToP(Double.parseDouble(elems[4]));

            tfout.writeln(Strings.concat(elems, Strings.tab) + "\t" + pderived + "\t" + n + "\t" + vals.length + "\t" + realp + "\t" + Strings.concat(pvalsBeta, Strings.tab));
            elems = tf.readLineElems(TextFile.tab);
        }

        tfout.close();
        tf.close();

    }

    public void betaDist(String real, String perm, String maxZPerPermFile, String out) throws IOException {

        DoubleMatrixDataset<String, String> dsPerm = DoubleMatrixDataset.loadDoubleData(perm);
        DoubleMatrixDataset<String, String> dsMaxPerm = DoubleMatrixDataset.loadDoubleData(maxZPerPermFile);

        for (int f = 0; f < 2; f++) {

            if (f == 0) {
                TextFile tf = new TextFile(real, TextFile.R);
                TextFile tfo = new TextFile(out + "-real.txt", TextFile.W);
                String header = "Gene\tP\t" + Strings.concat(dsPerm.getColObjects(), Strings.tab);
//                + "\t" + permheader
//                + "\t" + Strings.concat(dsPerm.getColObjects(), Strings.tab);
                tfo.writeln(header);
                int ctr = 0;
                tf.readLine();
                String ln = tf.readLine();

                while (ln != null) {
                    String[] elems = Strings.tab.split(ln);
                    String gene = elems[0];
                    Double z = Double.parseDouble(elems[4]);
                    DoubleMatrix1D row = dsPerm.getRow(gene);
                    double[] vals = filterNaNAndSort(row);
                    double[] pvals = converToP(vals);

                    double obsp = ZScores.zToP(z);

                    tfo.writeln(gene + "\t" + obsp + "\t" + Strings.concat(pvals, Strings.tab));

                    ctr++;

                    ln = tf.readLine();
                }
//
                tfo.close();
                tf.close();
            } else {
                int nrperm = 10;
                for (int p = 0; p < nrperm; p++) {
                    TextFile tfo = new TextFile(out + "-perm" + p + ".txt", TextFile.W);
                    String header = "Gene\tP\t" + Strings.concat(dsPerm.getColObjects(), Strings.tab);
                    tfo.writeln(header);
                    for (int g = 0; g < dsMaxPerm.rows(); g++) {
                        DoubleMatrix1D genedata = dsMaxPerm.getRow(g);
                        String gene = dsMaxPerm.getRowObjects().get(g);
                        DoubleMatrix1D genepermdata = dsPerm.getRow(gene);

                        double z = genedata.getQuick(p);
                        int n = 0;
                        int ntotal = 0;
                        HashSet<Integer> disallowedCols = new HashSet<Integer>();
                        for (int q = p; q < genepermdata.cardinality(); q += nrperm) {
                            disallowedCols.add(q);
                        }

                        ArrayList<Double> pvalList = new ArrayList<>();
                        for (int q = p; q < genepermdata.cardinality(); q++) {
                            if (!disallowedCols.contains(q)) {
                                double otherz = genepermdata.get(q);
                                if (!Double.isNaN(otherz)) {
                                    pvalList.add(ZScores.zToP(otherz));
                                    ntotal++;
                                }
                            }
                        }

                        Collections.sort(pvalList);

                        double[] pvals = Primitives.toPrimitiveArr(pvalList.toArray(new Double[0]));
                        double obsp = ZScores.zToP(z);
                        tfo.writeln(gene + "\t" + obsp + "\t" + Strings.concat(pvals, Strings.tab));
                    }
                    tfo.close();
                }

            }

        }

    }

    public void run(String real, String perm, String maxZPerPermFile, String scriptloc, String out) throws IOException {


        DoubleMatrixDataset<String, String> dsPerm = DoubleMatrixDataset.loadDoubleData(perm);
        DoubleMatrixDataset<String, String> dsMaxPerm = DoubleMatrixDataset.loadDoubleData(maxZPerPermFile);


        // determine number of stronger z-scores, excluding certain permutation
        int nrperm = 10;

//        for (int p = 0; p < nrperm; p++) {
        if (1 == 1) {
            IntStream.range(-1, nrperm).parallel().forEach(p -> {
                        try {
                            if (p == -1) {
                                TextFile tf = new TextFile(real, TextFile.R);
                                TextFile tfo = new TextFile(out, TextFile.W);
                                String header = tf.readLine() + "\tP\tDerivedP\tNPermZ>=Z\tNTotal\tP(betacdf)";
//                + "\t" + permheader
//                + "\t" + Strings.concat(dsPerm.getColObjects(), Strings.tab);
                                tfo.writeln(header);
                                int ctr = 0;
                                String ln = tf.readLine();


                                while (ln != null) {
                                    String[] elems = Strings.tab.split(ln);
                                    String gene = elems[0];
                                    Double z = Double.parseDouble(elems[4]);

                                    DoubleMatrix1D row = dsPerm.getRow(gene);
                                    double[] vals = filterNaNAndSort(row);
                                    double[] pvals = converToP(vals);

                                    int n = 0;
                                    for (int d = 0; d < vals.length; d++) {
                                        if (vals[d] >= z) {
                                            n++;
                                        }
                                    }
                                    double pderived = (double) n / vals.length;
//                                    System.out.println(Strings.concat(pvals, Strings.comma));
//                                    System.exit(-1);

                                    double obsp = ZScores.zToP(z);
                                    double maxP = 1d;

                                    String pythoninput = Strings.concat(pvals, Strings.comma);
//                                    System.out.println("Writing");
//                                    writer.write(obsp + " " + maxP + " " + pythoninput + "\n");
//                                    System.out.println("done writing");


//                                    String errline = readererr.readLine();
//                                    while (errline != null) {
//                                        System.out.println(errline);
//                                        errline = readererr.readLine();
//                                    }

//                                    String line = reader.readLine();
//                                    System.out.println("Done reading.");
//                                    double pbetacdf = Double.parseDouble(line);
                                    // V1

                                    Process process = Runtime.getRuntime().exec("python " + scriptloc + " " + obsp + " " + maxP + " " + pythoninput);
                                    BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                                    BufferedReader readererr = new BufferedReader(new InputStreamReader(process.getErrorStream()));
                                    String line = reader.readLine();

//                                System.out.println("Errors:");
                                    String errline = readererr.readLine();
                                    while (errline != null) {
                                        System.out.println(errline);
                                        errline = readererr.readLine();
                                    }
                                    readererr.close();

//                System.out.println("Output");
                                    double pbetacdf = Double.parseDouble(line);
                                    reader.close();
                                    readererr.close();
//                                writer.close();
                                    process.destroy();

//                                    System.out.println(p + " - " + ctr + "/" + dsMaxPerm.rows() + "\t" + obsp + "\t" + pderived + "\t" + n + "\t" + vals.length + "\t" + pbetacdf);
                                    tfo.writeln(ln + "\t" + obsp + "\t" + pderived + "\t" + n + "\t" + vals.length + "\t" + pbetacdf);

                                    ctr++;

                                    if (ctr % 100 == 0) {
                                        System.out.println(p + "\t" + ctr / dsMaxPerm.rows());
                                    }
                                    ln = tf.readLine();
                                }
//
                                tfo.close();
                                tf.close();
                            } else {

                                TextFile outf = new TextFile(out + "MaxAbsZPvaluesPerm-" + p + ".txt", TextFile.W);
                                outf.writeln("Gene\tMaxAbsZ\tP\tP(maxabsz>abszotherperm)\tN(maxabsz>abszotherperm)\tNTotal\tP(betacdf)");

                                RunTimer t = new RunTimer();
                                t.start();


                                for (int g = 0; g < dsMaxPerm.rows(); g++) {
                                    DoubleMatrix1D genedata = dsMaxPerm.getRow(g);
                                    String gene = dsMaxPerm.getRowObjects().get(g);
                                    DoubleMatrix1D genepermdata = dsPerm.getRow(gene);

                                    double z = genedata.getQuick(p);
                                    int n = 0;
                                    int ntotal = 0;
                                    HashSet<Integer> disallowedCols = new HashSet<Integer>();
                                    for (int q = p; q < genepermdata.cardinality(); q += nrperm) {
                                        disallowedCols.add(q);
                                    }

                                    ArrayList<Double> pvalList = new ArrayList<>();
                                    for (int q = p; q < genepermdata.cardinality(); q++) {
                                        if (!disallowedCols.contains(q)) {
                                            double otherz = genepermdata.get(q);
                                            if (!Double.isNaN(otherz)) {
                                                if (otherz >= z) {
                                                    n++;
                                                }
                                                pvalList.add(ZScores.zToP(otherz));
                                                ntotal++;
                                            }
                                        }
                                    }

                                    Collections.sort(pvalList);


                                    double[] pvals = Primitives.toPrimitiveArr(pvalList.toArray(new Double[0]));
//                double maxP = ArrayMath.max(pvals) + 0.000001;
                                    double maxP = 1d;

                                    double obsp = ZScores.zToP(z);


                                    String pythoninput = Strings.concat(pvals, Strings.comma);
//                                    writer.write(obsp + " " + maxP + " " + pythoninput + "\n");
//
//
//                                    String errline = readererr.readLine();
//                                    while (errline != null) {
//                                        System.out.println(errline);
//                                        errline = readererr.readLine();
//                                    }
//
//                                    String line = reader.readLine();
//                                    double pbetacdf = Double.parseDouble(line);
//
//                                     V1
//                                    /*
                                    Process process = Runtime.getRuntime().exec("python " + scriptloc + " " + obsp + " " + maxP + " " + pythoninput);


                                    BufferedReader readererr = new BufferedReader(new InputStreamReader(process.getErrorStream()));
                                    BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                                    String line = reader.readLine();

//                                System.out.println("Errors:");
                                    String errline = readererr.readLine();
                                    while (errline != null) {
                                        System.out.println(errline);
                                        errline = readererr.readLine();
                                    }
                                    readererr.close();
                                    reader.close();
                                    process.destroy();
//                System.out.println("Output");
                                    double pbetacdf = Double.parseDouble(line);
                                    reader.close();
//*/

                                    String outln = gene + "\t" + z + "\t" + ZScores.zToP(z) + "\t" + ((double) n / ntotal) + "\t" + n + "\t" + ntotal + "\t" + pbetacdf;
//                                    System.out.println(p + " - " + g + "/" + dsMaxPerm.rows() + "\t" + outln);
                                    outf.writeln(outln);

                                    if (g % 100 == 0) {
                                        long diff = t.getTimeDiff();
                                        double tperitem = (double) diff / g;
                                        double remain = tperitem * (dsMaxPerm.rows() - g);

                                        System.out.println(p + "\t" + g + "\t" + dsMaxPerm.rows() + " \tT-: " + t.getTimeDesc((long) remain));
                                    }
//


                                }

                                outf.close();

                            }
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
            );
        }
        System.exit(-1);

    }

    private double[] converToP(double[] vals) {
        double[] p = new double[vals.length];
        for (int i = 0; i < p.length; i++) {
            p[i] = ZScores.zToP(vals[i]);
        }
        return p;
    }

    private double[] filterNaNAndSort(DoubleMatrix1D row) {
        ArrayList<Double> d = new ArrayList<Double>();
        for (int i = 0; i < row.cardinality(); i++) {
            double v = row.get(i);
            if (!Double.isNaN(v)) {
                d.add(v);
            }
        }
        Collections.sort(d);

        return Primitives.toPrimitiveArr(d);


    }
}
