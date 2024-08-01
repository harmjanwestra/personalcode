package nl.harmjanwestra.playground.mbqtlvalidate;

import JSci.maths.statistics.FDistribution;
import cern.jet.random.tdouble.StudentT;
import nl.harmjanwestra.playground.legacy.vcf.VCFGenotypeData;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.util.ArrayList;
import java.util.HashMap;


public class MbQTLValidator {


    public static void main(String[] args) {

        if (args.length < 5) {
            System.out.println("Usage: vcf gte exp snpgene output");
            System.exit(-1);
        }
        String vcffile = args[0];
        String gtefile = args[1];
        String expfile = args[2];
        String snpgenefile = args[3];
        String output = args[4];
        MbQTLValidator v = new MbQTLValidator();
        try {
            v.run(vcffile, gtefile, expfile, snpgenefile, output);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }


    public void run(String vcffile, String gtefile, String expfile, String snpgenefile, String output) throws Exception {

        HashMap<String, String> snpGenePairs = new HashMap<>();

        TextFile tf = new TextFile(snpgenefile, TextFile.R);
        System.out.println("Parsing: " + snpgenefile);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[0];
            String gene = elems[1];

            snpGenePairs.put(snp, gene);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        System.out.println(snpGenePairs.size() + " snp/gene pairs");

        HashMap<String, String> genotypeRNAPairs = new HashMap<>();
        System.out.println("Parsing: " + gtefile);
        tf = new TextFile(gtefile, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String geno = elems[0];
            String rna = elems[1];
            genotypeRNAPairs.put(geno, rna);

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(genotypeRNAPairs.size() + " genotype/rna pairs");


        System.out.println("Loading exp data: " + expfile);
        DoubleMatrixDataset<String, String> expData = DoubleMatrixDataset.loadDoubleData(expfile);
        System.out.println("Expdata: " + expData.rows() + " x " + expData.columns());


        ArrayList<Integer> gtIndex = new ArrayList<>();
        ArrayList<Integer> rnaIndex = new ArrayList<>();

        TextFile outf = new TextFile(output, TextFile.W);
        String header = "SNP\tGene\tN" +
                "\tCorrelationJSci\tZ\tP" +
                "\tOLSSlope\tOLSP" +
                "\tOLSSlopeNoIntercept\tOLSPNoIntercept" +
                "\tOLSSlopeNorm\tOLSPNorm" +
                "\tOLSSlopeNormNoIntercept\tOLSPNormNoIntercept" +
                "\tfQTLcorr\tfQTLdf\tfQTLt\tfQTLp\tfQTLslope\tfQTLslope_se";
        outf.writeln(header);

        int lnctr = 0;
        int testctr = 0;

        System.out.println("Parsing: " + vcffile);
        TextFile vcf = new TextFile(vcffile, TextFile.R);
        String ln = vcf.readLine();

        boolean headerParsed = false;

        while (ln != null) {
            if (ln.startsWith("##")) {
                // ?
            } else if (ln.startsWith("#CHROM")) {
                System.out.println("Parsing header...");
                elems = ln.split("\t");
                for (int i = 9; i < elems.length; i++) {
                    String rna = genotypeRNAPairs.get(elems[i]);
//                    System.out.println(elems[i] + " --> " + rna);
                    if (rna != null) {

                        int rnaId = expData.getColIndex(rna);

                        if (rnaId >= 0) {
                            gtIndex.add(i);
                            rnaIndex.add(rnaId);
                        }
                    }
                }
                int nrGt = elems.length - 9;
                System.out.println(gtIndex.size() + " genotype samples overlapping RNA out of " + nrGt);

                if (gtIndex.isEmpty()) {
                    System.err.println("Error: no matching genotype samples loaded. Is header misformed?- " + lnctr);
                    System.exit(-1);
                }
            } else {
                if (gtIndex.isEmpty()) {
                    System.err.println("Error: no matching genotype samples loaded. Is header misformed? - " + lnctr);
                    System.exit(-1);
                }

                elems = ln.split("\t");
                if (elems.length > 9) {
                    String snpName = elems[2];
                    String geneName = snpGenePairs.get(snpName);
//                    System.out.println(snpName + "\t" + geneName);
                    String[] formatStr = elems[8].split(":");
                    int gtCol = 0;
                    for (int q = 0; q < formatStr.length; q++) {
                        if (formatStr[q].equals("GT")) {
                            gtCol = q;
                            break;
                        }
                    }

                    if (geneName != null) {
                        int geneId = expData.getRowIndex(geneName);
                        if (geneId >= 0) {
                            // test this eQTL: parse SNP first
                            double[] x = new double[gtIndex.size()];
                            double[] y = new double[gtIndex.size()];
                            double averageGt = 0;
                            int nrNotMissing = 0;
                            int nrMissing = 0;
                            for (int i = 0; i < gtIndex.size(); i++) {
                                int gtIdx = gtIndex.get(i);
                                int rnaIdx = rnaIndex.get(i);

                                y[i] = expData.getElementQuick(geneId, rnaIdx);

                                String gtStr = elems[gtIdx].split(":")[gtCol];
                                String[] gt = gtStr.split("/");
                                if (gt.length < 2) {
                                    gt = gtStr.split("\\|");
                                }
                                double dosage = -1;
                                if (!gt[0].equals(".")) {
                                    dosage = Double.parseDouble(gt[0]) + Double.parseDouble(gt[1]);
                                    nrNotMissing++;
                                    averageGt += dosage;
                                } else {
                                    nrMissing++;
                                }
                                x[i] = dosage;
                            }

                            if (nrMissing > 0) {
                                System.out.println("Warning: " + snpName + " does not have complete genotypes.");
                                averageGt /= nrNotMissing;
                                for (int i = 0; i < gtIndex.size(); i++) {
                                    if (x[i] == -1) {
                                        x[i] = averageGt;
                                    }
                                }
                            }


                            double[] xCopy = new double[x.length];
                            double[] yCopy = new double[x.length];
                            System.arraycopy(x, 0, xCopy, 0, x.length);
                            System.arraycopy(y, 0, yCopy, 0, y.length);

                            double[] xCopy2 = new double[x.length];
                            double[] yCopy2 = new double[y.length];
                            System.arraycopy(x, 0, xCopy2, 0, x.length);
                            System.arraycopy(y, 0, yCopy2, 0, y.length);

                            StudentT tDist = new cern.jet.random.tdouble.StudentT(x.length - 2, (new cern.jet.random.tdouble.engine.DRand()));
                            org.apache.commons.math3.stat.regression.SimpleRegression reg = new org.apache.commons.math3.stat.regression.SimpleRegression(true);
                            for (int q = 0; q < x.length; q++) {
                                reg.addData(x[q], y[q]);
                            }
                            double olsSlope = reg.getSlope();
                            double olsP = reg.getSignificance();

                            reg = new org.apache.commons.math3.stat.regression.SimpleRegression(false);
                            for (int q = 0; q < x.length; q++) {
                                reg.addData(x[q], y[q]);
                            }
                            double olsSlopeNoIntercept = reg.getSlope();
                            double olsPNoIntercept = reg.getSignificance();


                            normalize(xCopy);
                            normalize(yCopy);

                            reg = new org.apache.commons.math3.stat.regression.SimpleRegression(true);
                            for (int q = 0; q < x.length; q++) {
                                reg.addData(xCopy[q], yCopy[q]);
                            }
                            double olsSlopeNorm = reg.getSlope();
                            double olsPNorm = reg.getSignificance();

                            reg = new org.apache.commons.math3.stat.regression.SimpleRegression(false);
                            for (int q = 0; q < x.length; q++) {
                                reg.addData(xCopy[q], yCopy[q]);
                            }
                            double olsSlopeNormNoIntercept = reg.getSlope();
                            double olsPNormNoIntercept = reg.getSignificance();

                            double correlationJSci = JSci.maths.ArrayMath.correlation(x, y);
                            double z = ZScores.correlationToZ(correlationJSci, x.length, tDist);
                            double p = ZScores.zToP(z);

                            double[] correlationFastQTL = fastQTL(xCopy2, yCopy2);

                            if (correlationJSci != 0) {
                                String outln = snpName + "\t" + geneName + "\t" + x.length + "\t" + correlationJSci + "\t" + z + "\t" + p + "\t" + olsSlope + "\t" + olsP + "\t" + olsSlopeNoIntercept + "\t" + olsPNoIntercept + "\t" + olsSlopeNorm + "\t" + olsPNorm + "\t" + olsSlopeNormNoIntercept + "\t" + olsPNormNoIntercept + "\t" + Strings.concat(correlationFastQTL, Strings.tab);
                                outf.writeln(outln);
                            }
                            testctr++;
                        }


                    }
                }


            }
            lnctr++;
            ln = vcf.readLine();
            if (lnctr % 1000 == 0) {
                System.out.println(lnctr + " lines, " + testctr + " tests.");
            }
        }
        vcf.close();
        outf.close();

    }

    private void normalize(double[] arr) {
        double mean = JSci.maths.ArrayMath.mean(arr);
        double stdev = JSci.maths.ArrayMath.standardDeviation(arr);
//        System.out.println(mean + "\t" + stdev);
        for (int v = 0; v < arr.length; v++) {
            arr[v] = (arr[v] - mean) / stdev;
        }

    }


    private void normalizeFastQTL(double[] v) {
        int sample_count = v.length;
        double mean = 0.0;
        for (double value : v) {
            mean += value;
        }
        double sum = 0;
        mean /= sample_count;
        for (int i = 0; i < v.length; i++) {
            v[i] -= mean;
            sum += v[i] * v[i];
        }
        sum = Math.sqrt(sum);
        if (sum == 0) sum = 1;
        for (int s = 0; s < sample_count; s++) {
            v[s] /= sum;
        }
    }

    private double stdevFastQTL(double[] v) {
        double m_newS = 0;
        double m_oldS = 0;
        double m_newM = 0;
        double m_oldM = 0;
        double m_n = v.length;
        for (int i = 0; i < v.length; i++) {
            double x = v[i];
            m_newM = m_oldM + (x - m_oldM) / m_n;
            m_newS = m_oldS + (x - m_oldM) * (x - m_newM);
            m_oldM = m_newM;
            m_oldS = m_newS;
        }
        double var = m_newS / (m_n - 1);
        double stdev = Math.sqrt(var);
        return stdev;
    }

    private double correlFastQTL(double[] x, double[] y) {
        double corr = 0.0;
        for (int s = 0; s < x.length; s++) {
            corr += x[s] * y[s];
        }
        return corr;
    }

    private double getTStatFasgtQTL(double df, double corr) {
        double t = df * corr * corr / (1 - corr * corr);
        return t;
    }

    private double getPfromTFastQTL(double tstat2, double df) {
        //  double p = pf(tstat2, 1, df, 0, 0); // pf is from standard R lib
        FDistribution f = new FDistribution(1, df); // could be the other way around
        return f.cumulative(tstat2); // ?
    }

    private double getSlopteFastQTL(double nominal_correlation, double psd, double gsd) {
        if (gsd < 1e-16 || psd < 1e-16) return 0;
        else return nominal_correlation * psd / gsd;
    }


    private double[] fastQTL(double[] x, double[] y) {

        double stdevx = stdevFastQTL(x);
        double stdevy = stdevFastQTL(y);
        normalizeFastQTL(x);
        normalizeFastQTL(y);

        double corr = correlFastQTL(x, y);
        double df = x.length - 2;
        double t = getTStatFasgtQTL(df, corr);
        double p = getPfromTFastQTL(t, df);
        p = 1 - p; // other side...
        double slope = getSlopteFastQTL(corr, stdevy, stdevx);
        double slope_se = Math.abs(slope) / Math.sqrt(t);
        return new double[]{corr, df, t, p, slope, slope_se};
    }
}
