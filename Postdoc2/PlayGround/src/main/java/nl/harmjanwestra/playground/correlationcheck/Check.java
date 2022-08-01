package nl.harmjanwestra.playground.correlationcheck;

import nl.harmjanwestra.playground.qtltest.GeneExpressionData;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;

public class Check {
    public static void main(String[] args) {
        Check c = new Check();
//        if (args.length < 6) {
//            System.out.println("Usage: expbeforecorrection expaftercorrection covarfile sampletodataset genelimit outdir");
//        } else {
//            try {
//                c.run(args[0], args[1], args[2]
//                        , args[3], args[4], args[5]);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }

        if (args.length < 4) {
            System.out.println("Usage: expfile covarfile [sampletodataset] genelimit outdir");
        } else {
            try {

                String expfile = args[0];
                String covarfile = args[1];
                String sampletodataset = null;
                String genelimit = null;
                String outdir = null;
                if (args.length > 4) {
                    sampletodataset = args[2];
                    genelimit = args[3];
                    outdir = args[4];
                } else {
                    genelimit = args[2];
                    outdir = args[3];
                }


                c.correlateCovariates(expfile, covarfile, sampletodataset, genelimit, outdir);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void run(String expFileBeforeCorrection, String expFileAfterCorrection, String covarFile,
                    String sampleToDataset, String geneLimitFile, String outdir) throws IOException {

        System.out.println("Reading gene limit: " + geneLimitFile);
        TextFile tfg = new TextFile(geneLimitFile, TextFile.R);
        Set<String> geneselection = tfg.readAsSet(0, TextFile.tab);
        tfg.close();
        System.out.println(geneselection.size() + " genes to limit to");

        System.out.println("Reading sample to dataset: " + sampleToDataset);
        HashMap<String, ArrayList<String>> dsToSample = new HashMap<String, ArrayList<String>>();
        ArrayList<String> datasets = new ArrayList<>();
        TextFile tf = new TextFile(sampleToDataset, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        HashSet<String> sampleselection = new HashSet<>();
        while (elems != null) {
            String ds = elems[1];
            String rna = elems[0];
            ArrayList<String> rnas = dsToSample.get(ds);
            if (rnas == null) {
                rnas = new ArrayList<>();
                dsToSample.put(ds, rnas);
                datasets.add(ds);
            }
            sampleselection.add(rna);
            rnas.add(rna);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(datasets.size() + " datasets.");

        nl.harmjanwestra.playground.qtltest.GeneExpressionData ds1 = new GeneExpressionData(expFileBeforeCorrection, geneselection, sampleselection);
        nl.harmjanwestra.playground.qtltest.GeneExpressionData ds2 = new GeneExpressionData(expFileAfterCorrection, geneselection, sampleselection);
        nl.harmjanwestra.playground.qtltest.GeneExpressionData dsCov = new GeneExpressionData(covarFile, null, null); // samples are in genemap, rather than samplemap

        for (String ds : dsToSample.keySet()) {
            TextFile out = new TextFile(outdir + ds + ".txt", TextFile.W);
            out.writeln("Gene\tPre-MDS1\tPost-MDS1" +
                    "\tPre-MDS2\tPost-MDS2" +
                    "\tPre-MDS3\tPost-MDS3" +
                    "\tPre-MDS4\tPost-MDS4");
            ArrayList<String> samples = dsToSample.get(ds);
            System.out.println(ds);
            for (int i = 0; i < ds1.genes.length; i++) {
                String gene = ds1.genes[i];
                Integer otherGeneId = ds2.geneMap.get(gene);
                Integer cov1ID = dsCov.sampleMap.get("MDS1");
                Integer cov2ID = dsCov.sampleMap.get("MDS2");
                Integer cov3ID = dsCov.sampleMap.get("MDS3");
                Integer cov4ID = dsCov.sampleMap.get("MDS4");

                ArrayList<Double> covVals1 = new ArrayList<>();
                ArrayList<Double> covVals2 = new ArrayList<>();
                ArrayList<Double> covVals3 = new ArrayList<>();
                ArrayList<Double> covVals4 = new ArrayList<>();
                ArrayList<Double> exp1Vals = new ArrayList<>();
                ArrayList<Double> exp2Vals = new ArrayList<>();

                for (String s : samples) {
                    Integer id1 = ds1.sampleMap.get(s);
                    Integer id2 = ds2.sampleMap.get(s);
                    Integer idc = dsCov.geneMap.get(s);
                    if (id1 != null && id2 != null && idc != null) {
                        covVals1.add(dsCov.data[idc][cov1ID]);
                        covVals2.add(dsCov.data[idc][cov2ID]);
                        covVals3.add(dsCov.data[idc][cov3ID]);
                        covVals4.add(dsCov.data[idc][cov4ID]);
                        exp1Vals.add(ds1.data[i][id1]);
                        exp2Vals.add(ds2.data[otherGeneId][id2]);
                    }
                }

                // correlate
                double r11 = Correlation.correlate(Primitives.toPrimitiveArr(covVals1), Primitives.toPrimitiveArr(exp1Vals));
                double r12 = Correlation.correlate(Primitives.toPrimitiveArr(covVals1), Primitives.toPrimitiveArr(exp2Vals));
                double r21 = Correlation.correlate(Primitives.toPrimitiveArr(covVals2), Primitives.toPrimitiveArr(exp1Vals));
                double r22 = Correlation.correlate(Primitives.toPrimitiveArr(covVals2), Primitives.toPrimitiveArr(exp2Vals));
                double r31 = Correlation.correlate(Primitives.toPrimitiveArr(covVals3), Primitives.toPrimitiveArr(exp1Vals));
                double r32 = Correlation.correlate(Primitives.toPrimitiveArr(covVals3), Primitives.toPrimitiveArr(exp2Vals));
                double r41 = Correlation.correlate(Primitives.toPrimitiveArr(covVals4), Primitives.toPrimitiveArr(exp1Vals));
                double r42 = Correlation.correlate(Primitives.toPrimitiveArr(covVals4), Primitives.toPrimitiveArr(exp2Vals));

//                if (Math.abs(r12) > 0.2) {
//                    OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
//                    double[][] x = new double[covVals1.size()][1];
//                    double[] y = Primitives.toPrimitiveArr(exp1Vals);
//                    double[] y2 = Primitives.toPrimitiveArr(exp2Vals);
//                    for (int q = 0; q < covVals1.size(); q++) {
//                        x[q][0] = covVals1.get(q);
//                    }
//                    ols.newSampleData(y, x);
//                    double[] param = ols.estimateRegressionParameters();
//                    double[] residuals = ols.estimateResiduals();
//                    double rsq = ols.calculateRSquared();
//                    double rsq2 = ols.calculateAdjustedRSquared();
//
//
//                    double rr1 = Correlation.correlate(residuals, Primitives.toPrimitiveArr(covVals1));
//
//
//                    for (int q = 0; q < covVals1.size(); q++) {
//                        System.out.println(y[q] + "\t" + covVals1.get(q) + "\t" + residuals[q]);
//                    }
//
//                    System.out.println(ds + "\t" + gene + "\t" + r11 + "\t" + r12 + "\t" + Strings.concat(param, Strings.semicolon) + "\t" + rr1 + "\t" + rsq + "\t" + rsq2);
//                    ols.newSampleData(y2, x);
//
//                    param = ols.estimateRegressionParameters();
//                    residuals = ols.estimateResiduals();
//                    double rr2 = Correlation.correlate(residuals, Primitives.toPrimitiveArr(covVals1));
//                    rsq = ols.calculateRSquared();
//                    rsq2 = ols.calculateAdjustedRSquared();
////                    System.out.println(ds + "\t" + gene + "\t" + r11 + "\t" + r12 + "\t" + Strings.concat(param, Strings.semicolon) + "\t" + rr2 + "\t" + rsq + "\t" + rsq2);
////                    System.exit(-1);
//                }

                out.writeln(gene + "\t" + r11 + "\t" + r12
                        + "\t" + r21 + "\t" + r22
                        + "\t" + r31 + "\t" + r32
                        + "\t" + r41 + "\t" + r42);

                if (i % 1000 == 0) {
                    System.out.print(ds + "\t" + i + "\r");
                }
            }
            System.out.print(ds + "\tdone\n");
            out.close();
        }
    }

    public void correlateCovariates(String expfile,
                                    String covarFile,
                                    String sampleToDataset,
                                    String geneLimitFile,
                                    String outdir) throws IOException {

        System.out.println("Reading gene limit: " + geneLimitFile);
        TextFile tfg = new TextFile(geneLimitFile, TextFile.R);
        Set<String> geneselection = tfg.readAsSet(0, TextFile.tab);
        tfg.close();
        System.out.println(geneselection.size() + " genes to limit to");
        HashMap<String, ArrayList<String>> dsToSample = new HashMap<String, ArrayList<String>>();
        HashSet<String> sampleselection = null;
        if (sampleToDataset != null) {
            System.out.println("Reading sample to dataset: " + sampleToDataset);

            ArrayList<String> datasets = new ArrayList<>();
            TextFile tf = new TextFile(sampleToDataset, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            sampleselection = new HashSet<>();
            while (elems != null) {
                String ds = elems[1];
                String rna = elems[0];
                ArrayList<String> rnas = dsToSample.get(ds);
                if (rnas == null) {
                    rnas = new ArrayList<>();
                    dsToSample.put(ds, rnas);
                    datasets.add(ds);
                }
                sampleselection.add(rna);
                rnas.add(rna);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            System.out.println(datasets.size() + " datasets.");
        }

        nl.harmjanwestra.playground.qtltest.GeneExpressionData ds1 = new GeneExpressionData(expfile, geneselection, sampleselection);
        nl.harmjanwestra.playground.qtltest.GeneExpressionData dsCov = new GeneExpressionData(covarFile, null, null); // samples are in genemap, rather than samplemap

        if (sampleToDataset == null) {
            sampleselection = new HashSet<>();

            ArrayList<String> allSamples = new ArrayList<>();
            dsToSample.put("Overall", allSamples);
            for (String sample : ds1.samples) {
                sampleselection.add(sample);
                allSamples.add(sample);

            }

        }

        for (String ds : dsToSample.keySet()) {
            ArrayList<String> samples = dsToSample.get(ds);
            System.out.println(ds);
            double[][] output = new double[ds1.genes.length][dsCov.samples.length];
            for (int i = 0; i < ds1.genes.length; i++) {
                String gene = ds1.genes[i];
                for (int cov1ID = 0; cov1ID < dsCov.samples.length; cov1ID++) {
                    ArrayList<Double> covVals1 = new ArrayList<>();
                    ArrayList<Double> exp1Vals = new ArrayList<>();

                    for (String s : samples) {
                        Integer id1 = ds1.sampleMap.get(s);
                        Integer idc = dsCov.geneMap.get(s);
                        if (id1 != null && idc != null) {
                            covVals1.add(dsCov.data[idc][cov1ID]);
                            exp1Vals.add(ds1.data[i][id1]);

                        }
                    }

                    double r11 = Correlation.correlate(Primitives.toPrimitiveArr(covVals1), Primitives.toPrimitiveArr(exp1Vals));
                    output[i][cov1ID] = r11;
                }

                if (i % 1000 == 0) {
                    System.out.print(ds + "\t" + i + "\r");
                }
            }
            DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>();

            dsout.setMatrix(output);
            try {
                dsout.setRowObjects(Arrays.asList(ds1.genes));
                dsout.setColObjects(Arrays.asList(dsCov.samples));
            } catch (Exception e) {
                e.printStackTrace();
            }
            dsout.save(outdir + ds + ".txt");
            System.out.print(ds + "\tdone\n");
        }
    }

}
