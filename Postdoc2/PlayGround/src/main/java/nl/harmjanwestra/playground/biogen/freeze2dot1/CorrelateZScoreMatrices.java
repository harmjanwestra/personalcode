package nl.harmjanwestra.playground.biogen.freeze2dot1;

import com.google.common.util.concurrent.AtomicDouble;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class CorrelateZScoreMatrices {

    public static void main(String[] args) {

        String matrix1 = "";
        String matrix2 = "";
        String maffile = "";
        CorrelateZScoreMatrices z = new CorrelateZScoreMatrices();
        try {
//            z.correlate(matrix1, matrix2, maffile);
            String cisfile1 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\cis-cortex-eur-dump\\eQTLDump-sortZ-FDR-CohortInfoRemoved.txt.gz.txt.gz";
            String cisfile2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-MetaBrain2IDs.txt.gz";
            String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\cis-systematiccomparison-EUR\\FDRinBoth\\";
            String maffile1 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\SNPQCLog-fix.txt.gz";
            String maffile2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added-metabrain2dot1ids.txt.gz";
            z.mergeCisEQTLFiles(cisfile1, maffile1, cisfile2, maffile2, output, true, true, 0.05, true);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public void mergeCisEQTLFiles(String efile1, String maf1file, String efile2, String maf2file, String output, boolean requiresignificance1, boolean requiresignificance2, double fdrthreshold, boolean includenonoverlapping) throws IOException {

        Gpio.createDir(output);

        HashMap<String, HashMap<String, minimaleqtl>> map = new HashMap<>();
        TextFile tf = new TextFile(efile1, TextFile.R);
        System.out.println("Reading efile: " + efile1);
        String[] header = tf.readLineElems(TextFile.tab);
        int snpcol = -1;
        int genecol = -1;
        int ncol = -1;
        int allelecol = -1;
        int assesedcol = -1;
        int zcol = -1;
        int fdrcol = -1;
        for (int i = 0; i < header.length; i++) {
            if (header[i].equals("SNPName") || header[i].equals("")) {
                snpcol = i;
            } else if (header[i].equals("ProbeName")) {
                genecol = i;
            } else if (header[i].equals("SNPType")) {
                allelecol = i;
            } else if (header[i].equals("AlleleAssessed")) {
                assesedcol = i;
            } else if (header[i].equals("OverallZScore")) {
                zcol = i;
            } else if (header[i].equals("SumNumberOfSamples")) {
                ncol = i;
            } else if (header[i].equals("FDR")) {
                fdrcol = i;
            }
        }

        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        int inc = 0;
        while (elems != null) {
            String snp = Strings.cache(elems[snpcol]);

            String alleles = Strings.cache(elems[allelecol]);
            String assessed = Strings.cache(elems[assesedcol]);
            String gene = elems[genecol];
            float z = Float.parseFloat(elems[zcol]);
            float fdr = Float.parseFloat(elems[fdrcol]);
            minimaleqtl e = new minimaleqtl();

            e.alleles = alleles;
            e.assessed = assessed;
            e.fdr = fdr;
            e.n = Short.parseShort(elems[ncol]);
            e.z = z;
            boolean include = true;
            if (requiresignificance1 && fdr > fdrthreshold) {
                include = false;
                break;
            }
            if (include) {
                HashMap<String, minimaleqtl> lsofeqtl = map.get(gene);
                if (lsofeqtl == null) {
                    lsofeqtl = new HashMap<>();
                }
                lsofeqtl.put(snp, e);
                map.put(gene, lsofeqtl);
                inc++;
            }
            elems = tf.readLineElems(TextFile.tab);
            ctr++;
            if (ctr % 1000000 == 0) {
                System.out.print(ctr + " lines processed. " + map.size() + " genes so far, " + inc + " eqtls included\r");
            }
        }
        System.out.println();
        tf.close();
        System.out.print(ctr + " lines processed. " + map.size() + " genes so far, " + inc + " eqtls included\n");
        System.out.println();

        tf = new TextFile(efile2, TextFile.R);
        System.out.println("Reading efile: " + efile2);
        HashSet<String> overlappingGenes = new HashSet<String>();
        header = tf.readLineElems(TextFile.tab);
        snpcol = -1;
        genecol = -1;
        ncol = -1;
        allelecol = -1;
        assesedcol = -1;
        zcol = -1;
        fdrcol = -1;
        for (int i = 0; i < header.length; i++) {
            if (header[i].equals("SNPName") || header[i].equals("NotInDs-SNPName")) {
                snpcol = i;
            } else if (header[i].equals("ProbeName") || header[i].equals("NotInDs-ProbeName")) {
                genecol = i;
            } else if (header[i].equals("SNPType")) {
                allelecol = i;
            } else if (header[i].equals("AlleleAssessed")) {
                assesedcol = i;
            } else if (header[i].equals("OverallZScore")) {
                zcol = i;
            } else if (header[i].equals("SumNumberOfSamples")) {
                ncol = i;
            } else if (header[i].equals("FDR")) {
                fdrcol = i;
            }
        }

        elems = tf.readLineElems(TextFile.tab);
        ctr = 0;
        inc = 0;
        while (elems != null) {

            String snp = Strings.cache(elems[snpcol]);
            String gene = Strings.cache(elems[genecol]);
            HashMap<String, minimaleqtl> lsofeqtl = map.get(gene);

            if (lsofeqtl == null && includenonoverlapping) {
                lsofeqtl = new HashMap<>();
                String alleles = Strings.cache(elems[allelecol]);
                String assessed = Strings.cache(elems[assesedcol]);

                float z = Float.parseFloat(elems[zcol]);
                float fdr = Float.parseFloat(elems[fdrcol]);
                minimaleqtl e = new minimaleqtl();

                e.alleles = alleles;
                e.assessed = assessed;
                e.fdr2 = fdr;
                e.n2 = Short.parseShort(elems[ncol]);
                e.z2 = z;
                lsofeqtl.put(snp, e);
                map.put(gene, lsofeqtl);

            } else if (lsofeqtl != null) {
                minimaleqtl e = lsofeqtl.get(snp);
                if (e == null && includenonoverlapping) {
                    String alleles = Strings.cache(elems[allelecol]);
                    String assessed = Strings.cache(elems[assesedcol]);

                    float z = Float.parseFloat(elems[zcol]);
                    float fdr = Float.parseFloat(elems[fdrcol]);
                    e = new minimaleqtl();

                    e.alleles = alleles;
                    e.assessed = assessed;
                    e.fdr2 = fdr;
                    e.n2 = Short.parseShort(elems[ncol]);
                    e.z2 = z;
                    lsofeqtl.put(snp, e);
                    map.put(gene, lsofeqtl);
                } else if (e != null) {
                    float fdr = Float.parseFloat(elems[fdrcol]);
                    boolean include = true;
                    if (requiresignificance2 && fdr > fdrthreshold) {
                        include = false;
                        break;
                    }
                    if (include) {
                        String alleles = Strings.cache(elems[allelecol]);
                        String assessed = Strings.cache(elems[assesedcol]);
                        Boolean flip = BaseAnnot.flipalleles(e.alleles, e.assessed, alleles, assessed);
                        if (flip != null) {
                            float z = Float.parseFloat(elems[zcol]);
                            short n = Short.parseShort(elems[ncol]);
                            if (flip) {
                                z *= -1;
                            }

                            e.z2 = z;
                            e.fdr2 = fdr;
                            e.n2 = n;
                            e.match = true;
                            overlappingGenes.add(gene);
                            inc++;
                        }
                    }
                }
            }
            elems = tf.readLineElems(TextFile.tab);
            ctr++;
            if (ctr % 1000000 == 0) {
                System.out.print(ctr + " lines processed. " + map.size() + " genes so far, " + inc + " eqtls included.\r");
            }
        }
        System.out.print(ctr + " lines processed. " + map.size() + " genes so far, " + inc + " eqtls included.\r");
        System.out.println("Done");
        tf.close();

        for (String gene : map.keySet()) {
            if (!overlappingGenes.contains(gene)) {
                map.put(gene, null);
            }
        }

        HashMap<String, Float> mafmap1 = loadMafMAP(maf1file);
        HashMap<String, Float> mafmap2 = loadMafMAP(maf2file);

        int maxlen = 0;
        for (String gene : overlappingGenes) {
            HashMap<String, minimaleqtl> lsofeqtl = map.get(gene);
            if (lsofeqtl.size() > maxlen) {
                maxlen = lsofeqtl.size();
            }
        }

        Correlation.correlationToZScore(maxlen);

        // correlate overlapping genes
        Gpio.createDir(output + "/geneoutput/");
        TextFile outputcorrel = new TextFile(output + "Correlations.txt", TextFile.W);
        outputcorrel.writeln("gene" +
                "\tn" +
                "\tRb" +
                "\tRb-se" +
                "\tRb-z" +
                "\tRb-p" +
                "\trSpearman-beta" +
                "\tpSpearman-beta" +
                "\trPearson-beta" +
                "\tpPearson-beta" +
                "\trSpearman-zscores" +
                "\tpSpearman-zscores" +
                "\trPearson-zscores" +
                "\tpPearson-zscores" +
                "\tfdr1" +
                "\tfdr2");

        /*
          String lnout = gene + "\t" + b1t.length + "\t" + Strings.concat(rbv, Strings.tab) + "\t" + ZScores.pToZ(p) + "\t" + p
                                + "\t" + rsp + "\t" + psp
                                + "\t" + rpe + "\t" + ppe
                                + "\t" + rspz + "\t" + pspz
                                + "\t" + rpez + "\t" + ppez
                                + "\t" + minFDR1[0] + "\t" + minFDR2[0];
         */


        int gctr = 0;
        for (String gene : overlappingGenes) {
            HashMap<String, minimaleqtl> lsofeqtl = map.get(gene);
            // convert to beta + se
            double[] b1 = new double[lsofeqtl.size()];
            double[] b2 = new double[lsofeqtl.size()];
            double[] se1 = new double[lsofeqtl.size()];
            double[] se2 = new double[lsofeqtl.size()];
            double[] z1 = new double[lsofeqtl.size()];
            double[] z2 = new double[lsofeqtl.size()];

            final AtomicInteger q = new AtomicInteger();
            TextFile outputsum = new TextFile(output + "/geneoutput/" + gene + ".txt.gz", TextFile.W);
            outputsum.writeln("SNP\tAlleles\tAssessed" +
                    "\tz1\tn1\tfdr1\tmaf1\tb1\tse1" +
                    "\tz2\tn2\tfdr2\tmaf2\tb2\tse2");
            final double[] minFDR1 = {1};
            final double[] minFDR2 = {1};

            lsofeqtl.forEach((k, v) -> {
                if (!v.match) {
                    boolean include = true;
                    if (requiresignificance1 && v.fdr > fdrthreshold) {
                        include = false;
                    }
                    if (requiresignificance2 && v.fdr2 > fdrthreshold) {
                        include = false;
                    }
                    if (include) {
                        Float maf1 = mafmap1.get(k);
                        Float maf2 = mafmap2.get(k);
                        if (!Float.isNaN(v.z)) {
                            if (maf1 != null) {
                                double[] zToBeta = ZScores.zToBeta(v.z, maf1, v.n);
                                try {
                                    outputsum.writeln(k + "\t" + v.alleles + "\t" + v.assessed
                                            + "\t" + v.z + "\t" + v.n + "\t" + v.fdr + "\t" + maf1 + "\t" + zToBeta[0] + "\t" + zToBeta[1]
                                            + "\t" + v.z2 + "\t" + v.n2 + "\t" + v.fdr2 + "\t" + maf2 + "\t" + 0 + "\t" + 0);
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                            }
                        } else {
                            if (maf2 != null) {
                                double[] zToBeta = ZScores.zToBeta(v.z2, maf2, v.n2);
                                try {
                                    outputsum.writeln(k + "\t" + v.alleles + "\t" + v.assessed
                                            + "\t" + v.z + "\t" + v.n + "\t" + v.fdr + "\t" + maf1 + "\t" + 0 + "\t" + 0
                                            + "\t" + v.z2 + "\t" + v.n2 + "\t" + v.fdr2 + "\t" + maf2 + "\t" + zToBeta[0] + "\t" + zToBeta[1]);
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                            }
                        }
                    }
                } else {
                    boolean include = true;
                    if (requiresignificance1 && v.fdr > fdrthreshold) {
                        include = false;
                    }
                    if (requiresignificance2 && v.fdr2 > fdrthreshold) {
                        include = false;
                    }

                    if (include) {
                        Float maf1 = mafmap1.get(k);
                        Float maf2 = mafmap2.get(k);
                        if (maf1 != null && maf2 != null) {
                            double[] val1 = ZScores.zToBeta(v.z, maf1, v.n);
                            double[] val2 = ZScores.zToBeta(v.z2, maf2, v.n2);
                            int qlocal = q.get();
                            b1[qlocal] = val1[0];
                            se1[qlocal] = val1[1];
                            b2[qlocal] = val2[0];
                            se2[qlocal] = val2[1];
                            z1[qlocal] = v.z;
                            z2[qlocal] = v.z2;
                            q.getAndIncrement();

                            if (v.fdr < minFDR1[0]) {
                                minFDR1[0] = v.fdr;
                            }
                            if (v.fdr2 < minFDR2[0]) {
                                minFDR2[0] = v.fdr2;
                            }
                            try {
                                outputsum.writeln(k + "\t" + v.alleles + "\t" + v.assessed
                                        + "\t" + v.z + "\t" + v.n + "\t" + v.fdr + "\t" + maf1 + "\t" + val1[0] + "\t" + val1[1]
                                        + "\t" + v.z2 + "\t" + v.n2 + "\t" + v.fdr2 + "\t" + maf2 + "\t" + val2[0] + "\t" + val2[1]);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }
                    }
                }
            });
            outputsum.close();


            // trim
            int ql = q.get();
            if (ql > 20) {
                double[] b1t = trim(b1, ql);
                double[] b2t = trim(b2, ql);
                double[] se1t = trim(se1, ql);
                double[] se2t = trim(se2, ql);
                double[] z1t = trim(z1, ql);
                double[] z2t = trim(z2, ql);

                double[] theta = new double[se1t.length];
                if (b1t.length > 20) {
                    double[] rbv = Rb(b1t, se1t, b2t, se2t, theta);
                    // print correlation somewhere

                    SpearmansCorrelation correl = new SpearmansCorrelation();
                    double rsp = correl.correlation(b1t, b2t);
                    double psp = Correlation.convertCorrelationToZScore(b1t.length, rsp);

                    PearsonsCorrelation pc = new PearsonsCorrelation();
                    double rpe = pc.correlation(b1t, b2t);
                    double ppe = Correlation.convertCorrelationToZScore(b1t.length, rpe);


                    double rspz = correl.correlation(z1t, z2t);
                    double pspz = Correlation.convertCorrelationToZScore(b1t.length, rsp);


                    double rpez = pc.correlation(z1t, z2t);
                    double ppez = Correlation.convertCorrelationToZScore(b1t.length, rpe);


                    try {

                        double p = rbv[2];

                        String lnout = gene + "\t" + b1t.length + "\t" + Strings.concat(rbv, Strings.tab) + "\t" + ZScores.pToZ(p) + "\t" + p
                                + "\t" + rsp + "\t" + psp
                                + "\t" + rpe + "\t" + ppe
                                + "\t" + rspz + "\t" + pspz
                                + "\t" + rpez + "\t" + ppez
                                + "\t" + minFDR1[0] + "\t" + minFDR2[0];
                        System.out.println(gctr + "\t" + lnout);
                        outputcorrel.writeln(lnout);
                    } catch (Exception e) {
                        e.printStackTrace();
                        String lnout = gene + "\t" + b1t.length + "\t" + Strings.concat(rbv, Strings.tab) + "\t" + rsp + "\t" + psp + "\t" + minFDR1[0] + "\t" + minFDR2[0];
                        System.err.println(lnout);
                    }

                }
            }
            gctr++;
        }
        outputcorrel.close();
    }

    private double[] trim(double[] b1, int q) {
        double[] v = new double[q];
        System.arraycopy(b1, 0, v, 0, q);
        return v;
    }

    private HashMap<String, Float> loadMafMAP(String maf1) throws IOException {
        HashMap<String, Float> f = new HashMap<>();
        TextFile tf = new TextFile(maf1, TextFile.R);
        System.out.println("Loading MAF map: " + maf1);
        String[] header = tf.readLineElems(TextFile.tab);
        int snpcol = -1;
        int mafcol = -1;
        for (int i = 0; i < header.length; i++) {
            if (header[i].toLowerCase().equals("snp")) {
                snpcol = i;
            } else if (header[i].toLowerCase().equals("maf") || header[i].toLowerCase().equals("alleleb_all")) {
                mafcol = i;
            }
        }

        System.out.println(maf1 + "\tsnpcol: " + snpcol + ", mafcol: " + mafcol);

        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {

            boolean b = Strings.cacheContains(elems[snpcol]);
            if (b) {
                String snp = Strings.cache(elems[snpcol]);
                if (!elems[mafcol].equals("NA")) {
                    Float maf = Float.parseFloat(elems[mafcol]);
                    if (!Float.isNaN(maf) && !Float.isInfinite(maf)) {
                        if (maf > 0.5) {
                            maf = 1 - maf;
                        }
                        f.put(snp, maf);
                    }
                }
            }

            elems = tf.readLineElems(TextFile.tab);
            ctr++;
            if (ctr % 1000000 == 0) {
                System.out.print(ctr + " lines read.\r");
            }
        }
        tf.close();
        return f;
    }


    public class SNP {
        String id;
        String alleles;
        String assessed;
        double[] z;
    }

    public class Matrix {
        ArrayList<SNP> snps;
        HashMap<String, SNP> stringToSNP = new HashMap<>();
        ArrayList<String> genes;
        HashMap<String, Integer> geneToIdx = new HashMap<>();

        public Matrix(String matrix) throws IOException {
            snps = new ArrayList<>();
            genes = new ArrayList<>();

            TextFile tf = new TextFile(matrix, TextFile.R);
            String[] header = tf.readLineElems(TextFile.tab);

            int ctr = 0;
            for (int d = 3; d < header.length; d++) {
                genes.add(header[d]);
                geneToIdx.put(header[d], ctr);
                ctr++;
            }
            String[] elems = tf.readLineElems(TextFile.tab);

            while (elems != null) {
                String id = elems[1];
                String alleles = elems[2];
                String assessed = elems[3];
                SNP s = new SNP();
                s.id = id;
                s.alleles = alleles;
                s.assessed = assessed;
                s.z = new double[genes.size()];
                for (int d = 0; d < elems.length; d++) {
                    s.z[d] = Double.parseDouble(elems[d]);
                }
                snps.add(s);
                stringToSNP.put(id, s);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();


        }

    }

    public void correlate(String matrix1file, String matrix2file, String maffile) throws IOException {

        Matrix m1 = new Matrix(matrix1file);
        Matrix m2 = new Matrix(matrix2file);


    }

    public double[] Rb(double[] b1, double[] se1, double[] b2, double[] se2, double[] theta) {

        double[] squarese1 = square(se1);
        double[] squarese2 = square(se2);

        double var_b1 = Descriptives.variance(b1) - Descriptives.mean(squarese1);
        double var_b2 = Descriptives.variance(b2) - Descriptives.mean(squarese2);
        if (var_b1 < 0) {
            var_b1 = Descriptives.variance(b1);
        }
        if (var_b2 < 0) {
            var_b2 = Descriptives.variance(b2);
        }

        double cov_b1_b2 = JSci.maths.ArrayMath.covariance(b1, b2) - (Descriptives.mean(theta) * Math.sqrt(Descriptives.mean(squarese1) * Descriptives.mean(squarese2)));
        double r = cov_b1_b2 / Math.sqrt(var_b1 * var_b2);
        double n = b1.length;
        double[] r_jack = new double[b1.length];

        IntStream.range(0, b1.length).parallel().forEach(k -> {
            double[] b1_jack = jack(b1, k);
            double[] se1_jack = jack(se1, k);
            double[] b2_jack = jack(b2, k);
            double[] se2_jack = jack(se2, k);
            double[] theta_jack = jack(theta, k);

            double[] squarese1_jack = square(se1_jack);
            double[] squarese2_jack = square(se2_jack);
            double var_b1_jack = Descriptives.variance(b1_jack) - Descriptives.mean(squarese1_jack);
            double var_b2_jack = Descriptives.variance(b2_jack) - Descriptives.mean(squarese2_jack);

            if (var_b1_jack < 0) {
                var_b1_jack = Descriptives.variance(b1_jack);
            }
            if (var_b2_jack < 0) {
                var_b2_jack = Descriptives.variance(b2_jack);
            }

            double cov_e1_jack_e2_jack = Descriptives.mean(theta_jack) * Math.sqrt(Descriptives.mean(squarese1_jack) * Descriptives.mean(squarese2_jack));
            double cov_b1_b2_jack = JSci.maths.ArrayMath.covariance(b1_jack, b2_jack) - cov_e1_jack_e2_jack;
            double rtmp = cov_b1_b2_jack / Math.sqrt(var_b1_jack * var_b2_jack);
            r_jack[k] = rtmp;
        });


        // remove nans
        double meanR = Descriptives.mean(r_jack);
        double se_r = Math.sqrt((n - 1) / (n * difsqsum(r_jack, meanR)));
        double p = ChiSquare.getP(1, (r / se_r) * (r / se_r));


        return new double[]{r, se_r, p};
    }

    private double[] jack(double[] b1, int k) {
        double[] jack = new double[b1.length - 1];
        int ctr = 0;
        for (int d = 0; d < b1.length; d++) {
            if (d != k) {
                jack[ctr] = b1[d];
                ctr++;
            }

        }
        return jack;
    }

    private double difsqsum(double[] r_jack, double meanR) {
        double sq = 0;
        for (int d = 0; d < r_jack.length; d++) {
            double v = (r_jack[d] - meanR);
            if (!Double.isNaN(v) && !Double.isInfinite(v)) {
                sq += (v * v);
            }
        }
        return sq;
    }

    private double[] square(double[] se2) {
        double[] output = new double[se2.length];
        for (int d = 0; d < se2.length; d++) {
            output[d] = se2[d] * se2[d];
        }
        return output;
    }

}
