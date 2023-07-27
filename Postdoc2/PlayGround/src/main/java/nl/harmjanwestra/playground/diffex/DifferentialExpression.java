package nl.harmjanwestra.playground.diffex;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.TTest;
import umcg.genetica.math.stats.WilcoxonMannWhitney;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

class DifferentialExpression {


    public static void main(String[] args) {

        if (args.length < 3) {
            System.out.println("Usage: phenotypefile phenotypename covariates");
//            System.out.println("Usage: expfile.txt.gz phenotypefile phenotypename outfile.txt.gz");
            System.exit(-1);
        }

        DifferentialExpression d = new DifferentialExpression();
        try {
//            d.run(args[0], args[1], args[2], args[3]);
            d.correlatePhenotypes(args[0],args[1],args[2]);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public void correlatePhenotypes(String phenotypefile, String phenotypename, String covariatefile) throws Exception {
        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(covariatefile);

        ArrayList<String> colSelection = new ArrayList<>();
        for (String s : ds.getColObjects()) {
            if (!s.contains("CMC_HBCC")) {
                colSelection.add(s);
            }
        }
        String[] colSelectionArr = colSelection.toArray(new String[0]);
        ds = ds.viewColSelection(colSelectionArr);

        TextFile tf = new TextFile(phenotypefile, TextFile.R);

        int samplecol = -1;
        int phenocol = -1;
        String[] header = tf.readLineElems(TextFile.tab);
        for (int i = 0; i < header.length; i++) {
            if (header[i].toLowerCase().equals("rnaseq_id")) {
                samplecol = i;
            }
            if (header[i].toLowerCase().equals(phenotypename.toLowerCase())) {
                phenocol = i;
            }
        }

        boolean ok = true;
        if (phenocol < 0) {
            System.out.println("Phenotype column with name " + phenotypename + " not found in " + phenotypefile);
            ok = false;
        }
        if (samplecol < 0) {
            System.out.println("Sample column with name " + samplecol + " not found in " + phenotypefile);
            ok = false;
        }
        if (!ok) {
            tf.close();
            System.exit(-1);
        }

        HashMap<String, String> phenotypes = new HashMap<>();
        HashMap<String, Integer> phenomap = new HashMap<>();
        int pctr = 0;
        LinkedHashMap<String, Integer> dsSamples = ds.getHashRows();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String sample = elems[samplecol];
            String pheno = elems[phenocol].toLowerCase();
            Integer index = dsSamples.get(sample);
            if (index != null) {
                phenotypes.put(sample, pheno);
                Integer phenoid = phenomap.get(pheno);
                if (phenoid == null) {
                    System.out.println("New phenotyope: " + pheno + " - " + pctr);
                    phenomap.put(pheno, pctr);
                    pctr++;
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        if (phenomap.size() > 2) {
            System.out.println("This code is dumb and can only handle a single phenotype");
            System.exit(-1);
        }


        double[] pheno = new double[phenotypes.size()];
        pctr = 0;
        String[] rowselection = new String[phenotypes.size()];
        for (int i = 0; i < ds.rows(); i++) {
            String sample = ds.getRowObjects().get(i);
            String p = phenotypes.get(sample);
            if (p != null) {
                int idx = phenomap.get(p);
                pheno[pctr] = idx;
                rowselection[pctr] = sample;
                pctr++;
            }
        }

        DoubleMatrixDataset<String, String> filtered = ds.viewRowSelection(rowselection);
        for (int c = 0; c < filtered.columns(); c++) {
            DoubleMatrix1D col = filtered.viewCol(c);
            double[] colarr = col.toArray();
            double correl = Correlation.correlate(pheno, colarr);
            System.out.println(filtered.getColObjects().get(c) + "\t" + correl + "\t" + colarr.length);
        }


        OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
        ols.newSampleData(pheno, filtered.getMatrixAs2dDoubleArray());
        double rsq = ols.calculateRSquared();
        System.out.println("multiple rsquared: " + rsq);


    }

    public void run(String expfile, String phenotypefile, String phenotypename, String outfile) throws Exception {
        // expect samples on columns
        System.out.println("Loading " + expfile);

        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(expfile);

        TextFile tf = new TextFile(phenotypefile, TextFile.R);

        int samplecol = -1;
        int phenocol = -1;
        String[] header = tf.readLineElems(TextFile.tab);
        for (int i = 0; i < header.length; i++) {
            if (header[i].toLowerCase().equals("rnaseq_id")) {
                samplecol = i;
            }
            if (header[i].toLowerCase().equals(phenotypename.toLowerCase())) {
                phenocol = i;
            }
        }

        boolean ok = true;
        if (phenocol < 0) {
            System.out.println("Phenotype column with name " + phenotypename + " not found in " + phenotypefile);
            ok = false;
        }
        if (samplecol < 0) {
            System.out.println("Sample column with name " + samplecol + " not found in " + phenotypefile);
            ok = false;
        }
        if (!ok) {
            tf.close();
            System.exit(-1);
        }

        HashMap<String, String> phenotypes = new HashMap<>();
        HashMap<String, Integer> phenomap = new HashMap<>();
        int pctr = 0;
        LinkedHashMap<String, Integer> dsSamples = ds.getHashCols();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String sample = elems[samplecol];
            String pheno = elems[phenocol].toLowerCase();
            Integer index = dsSamples.get(sample);
            if (index != null) {
                phenotypes.put(sample, pheno);
                Integer phenoid = phenomap.get(pheno);
                if (phenoid == null) {
                    System.out.println("New phenotyope: " + pheno + " - " + pctr);
                    phenomap.put(pheno, pctr);
                    pctr++;
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        System.out.println(phenomap.size() + " phenotypes");
        System.out.println(phenotypes.size() + " overlapping samples");

        if (phenomap.size() > 2) {
            System.out.println("This code is dumb and can only handle a single phenotype");
            System.exit(-1);
        }

        // split the samples into two groups
        ArrayList<Integer> grpA = new ArrayList<>();
        ArrayList<Integer> grpB = new ArrayList<>();

        ArrayList<String> dsSampleList = ds.getColObjects();
        for (int i = 0; i < dsSampleList.size(); i++) {
            String sample = dsSampleList.get(i);
            String pheno = phenotypes.get(sample);
            Integer phenoId = phenomap.get(pheno);
            if (phenoId != null) {
                if (phenoId == 0) {
                    grpA.add(i);
                } else {
                    grpB.add(i);
                }
            }
        }

        TextFile outf = new TextFile(outfile, TextFile.W);
        outf.writeln("Gene\tnA\tnB\tFC\tTTest-P\tTTest-log10P\tWilcoxon-P\tWilcoxon-log10P\tWilcoxon-AUC\tMWU-P");
        ProgressBar pb = new ProgressBar(ds.rows(), "Calculating differential expression.");
        for (int row = 0; row < ds.rows(); row++) {
            String gene = ds.getRowObjects().get(row);
            double[] a = new double[grpA.size()];
            double[] b = new double[grpB.size()];
            int ctr = 0;
            for (int i : grpA) {
                a[ctr] = ds.getElementQuick(row, i);
                ctr++;
            }
            ctr = 0;
            for (int i : grpB) {
                b[ctr] = ds.getElementQuick(row, i);
                ctr++;
            }
            double tp = TTest.test(a, b);
            WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();
            
            double wp = mwm.returnWilcoxonMannWhitneyPValue(a, b);
            double meanA = Descriptives.mean(a);
            double meanB = Descriptives.mean(b);



            double fc = meanB / meanA;
            String outln = gene
                    + "\t" + a.length
                    + "\t" + b.length
                    + "\t" + fc
                    + "\t" + tp
                    + "\t" + (-Math.log10(tp))
                    + "\t" + wp
                    + "\t" + (-Math.log10(wp))
                    + "\t" + mwm.getAUC();
            outf.writeln(outln);
            pb.set(row);
        }
        outf.close();
        pb.close();

    }

}