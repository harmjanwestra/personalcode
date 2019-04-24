package nl.harmjanwestra.playground.biogen.covariates;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;

public class CovariateProcessor {

    public static void main(String[] args) {
        String covmat = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\brain.phenotype_QC_covariates.txt";
        String outdir = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\pca\\";
        try {

            CovariateProcessor c = new CovariateProcessor();
//            c.correlateNonCategorical(covmat, outdir);
            String input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates.txt";
            String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-cohortvars.txt";
//            c.makeDummyVariables(input, 1, output);

            String quals = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-qualityscores.txt";
            String qualsout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-qualityscores-filter.txt";
//            c.filterNumberOfVars(quals, qualsout, 0.05);


            String allcovars = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates.txt";
            String allcovarsfiltered = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\2019-03-27-brain.phenotype_QC_covariates-noncategorical.txt";
            c.removeCategoricalColumns(allcovars,allcovarsfiltered);

        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void filterNumberOfVars(String quals, String qualsout, double maxpercmissing) throws Exception {

        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(quals);

        boolean[] includecol = new boolean[ds.columns()];
        int nrRemainingCols = 0;
        ArrayList<String> newcolobjs = new ArrayList<>();
        for (int col = 0; col < ds.columns(); col++) {
            double[] vals = ds.getMatrix().viewColumn(col).toArray();
            int nrmissing = 0;
            for (int r = 0; r < vals.length; r++) {
                if (Double.isNaN(vals[r])) {
                    nrmissing++;
                }
            }
            double percmissing = (double) nrmissing / vals.length;
            if (percmissing < maxpercmissing) {
                includecol[col] = true;
                newcolobjs.add(ds.getColObjects().get(col));
                nrRemainingCols++;
            } else {
                System.out.println("Excluding " + ds.getColObjects().get(col) + " - " + percmissing);
            }
        }

        DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<>(ds.rows(), nrRemainingCols);
        int colctr = 0;
        for (int col = 0; col < ds.columns(); col++) {
            if (includecol[col]) {
                dsout.getMatrix().viewColumn(colctr).assign(ds.getMatrix().viewColumn(col));
                colctr++;
            }
        }

        dsout.setRowObjects(ds.getRowObjects());
        dsout.setColObjects(newcolobjs);
        dsout.save(qualsout);

    }

    private void makeDummyVariables(String input, int i, String output) throws IOException {

        HashMap<String, Integer> cats = new HashMap<>();
        HashMap<String, Integer> catsctr = new HashMap<>();
        TextFile tf = new TextFile(input, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        String[] elems = tf.readLineElems(TextFile.tab);

        while (elems != null) {
            String cat = elems[i];
            if (cats.containsKey(cat)) {

                cats.put(cat, cats.size());
            }
            Integer ct = catsctr.get(cat);
            if (ct == null) {
                ct = 0;
            }
            ct++;
            catsctr.put(cat, ct);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        // make largest category the reference
        String maxcat = null;
        int maxcatct = 0;

        for (Map.Entry<String, Integer> catpair : catsctr.entrySet()) {
            if (catpair.getValue() > maxcatct) {
                maxcatct = catpair.getValue();
                maxcat = catpair.getKey();
            }
        }

        System.out.println(maxcat + "  is max category with  " + maxcatct + " samples.");

        // reindex categories
        HashMap<String, Integer> catsfinal = new HashMap<>();
        ArrayList<String> catlist = new ArrayList<>();
        int catctr = 0;
        for (Map.Entry<String, Integer> catpair : catsctr.entrySet()) {
            if (!catpair.getKey().equals(maxcat)) {
                catsfinal.put(catpair.getKey(), catctr);
                catlist.add(catpair.getKey());
                catctr++;
            }
        }

        String headerout = "Sample\t" + Strings.concat(catlist, Strings.tab);
        TextFile outf = new TextFile(output, TextFile.W);
        outf.writeln(headerout);
        tf.open();
        tf.readLine();
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {

            int[] membership = new int[catsfinal.size()];
            String cat = elems[i];
            Integer catid = catsfinal.get(cat);
            if (catid != null) {
                membership[catid] = 1;
            }
            String lnout = elems[0] + "\t" + Strings.concat(membership, Strings.tab);
            outf.writeln(lnout);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        outf.close();

    }

    private void correlateNonCategorical(String input, String outdir) throws Exception {


        removeCategoricalColumns(input, outdir + "NonCategorical.txt");

        removeNonVariantColumns(outdir + "NonCategorical.txt", 20, outdir + "NonCategoricalAndVariable.txt");

        // now we should have a double matrix dataset.. samples on the rows
        DoubleMatrixDataset<String, String> covariates = DoubleMatrixDataset.loadDoubleData(outdir + "NonCategoricalAndVariable.txt");

        // prepare correlation matrix over covariates

        double[][] cormat = new double[covariates.columns()][covariates.columns()];

        for (int i = 0; i < covariates.columns(); i++) {
            DoubleMatrix1D covA = covariates.getCol(i);
            for (int j = i + 1; j < covariates.columns(); j++) {
                DoubleMatrix1D covB = covariates.getCol(j);
                double r = removeNanAndCorrelate(covA, covB);
                cormat[i][j] = r;
                cormat[j][i] = r;

            }
            System.out.println(i + " / " + covariates.columns());

            cormat[i][i] = 1d;
        }

//        // remove cols with NaN
//        double[][] cormatNoNan = new double[covariates.columns() - nrWithNan][covariates.columns() - nrWithNan];
//        ArrayList<String> covariatesWithoutNan = new ArrayList<>();
//        int rctr = 0;
//        for (int i = 0; i < covariates.columns(); i++) {
//            if (!hashnan[i]) {
//                int colctr = 0;
//                for (int j = 0; j < covariates.columns(); j++) {
//                    if (!hashnan[j]) {
//                        cormatNoNan[rctr][colctr] = cormat[i][j];
//                        colctr++;
//                    }
//                }
//                covariatesWithoutNan.add(covariates.getColObjects().get(i));
//                rctr++;
//            }
//        }
//
//
//        DoubleMatrixDataset<String, String> correlationmatrix = new DoubleMatrixDataset<>(covariatesWithoutNan.size(), covariatesWithoutNan.size());
//        correlationmatrix.setRowObjects(covariatesWithoutNan);
//        correlationmatrix.save(outdir + "CorrelationMatrixWithoutNaN.txt");
//
//
//        HashSet<String> covariateHash = new HashSet<>();
//        covariateHash.addAll(covariatesWithoutNan);
//        covariates = DoubleMatrixDataset.loadSubsetOfTextDoubleData(outdir + "NonCategoricalAndVariable.txt", '\t', null, covariateHash);

//
//        double[][] finalCormat = new double[covariateNamesWithOutNaN.size()][covariateNamesWithOutNaN.size()];
//        int ictr = 0;
//        int jctr = 0;
//        for (int i = 0; i < covariateNames.size(); i++) {
//
//            if (!hasNan[i]) {
//                jctr = 0;
//
//                for (int j = i + 1; j < covariateNames.size(); j++) {
//                    if (!hasNan[j]) {
//                        finalCormat[ictr][jctr] = cormat[i][j];
//                        jctr++;
//                    }
//                }
//                ictr++;
//            }
//        }
//
//        DoubleMatrixDataset<String, String> corWithOutNan = new DoubleMatrixDataset<>();
//        try {
//            corWithOutNan.setMatrix(cormat);
//            corWithOutNan.setColObjects(covariateNamesWithOutNaN);
//            corWithOutNan.setRowObjects(covariateNamesWithOutNaN);
//            corWithOutNan.save(outdir + "PrunedCovariateMatrix-CorMatWithOutNaN.txt");
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//

//        // perform PCA on remainder of matrix
//        Jama.EigenvalueDecomposition eig = PCA.eigenValueDecomposition(cormat);
//        DoubleMatrixDataset<String, String> datasetEV = new DoubleMatrixDataset<>(corWithOutNan.getColObjects().size(), corWithOutNan.getColObjects().size());
//        try {
//            datasetEV.setColObjects(corWithOutNan.getColObjects());
//            datasetEV.setRowObjects(corWithOutNan.getColObjects());
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//        double[] eigenValues = eig.getRealEigenvalues();
//
//
//        int nrOfPCsToCalculate = ds.getColObjects().size();
//        System.out.println("PCA\tPCANr\tEigenValue\tExplainedVariance\tTotalExplainedVariance");
//
//        TextFile out = new TextFile(outdir + "PrunedCovariateMatrix-PCAOverCovariatesEigenvalues.txt.gz", TextFile.W);
//        double cumExpVarPCA = 0;
//
//        out.writeln("PCA\tPCANr\tEigenValue\tExplainedVariance\tTotalExplainedVariance");
//
//
//        for (int pca = 0; pca < nrOfPCsToCalculate; pca++) {
//            double expVarPCA = PCA.getEigenValueVar(eigenValues, pca);
//            double[] pca1ExpEigenVector = PCA.getEigenVector(eig, eigenValues, pca);
//            for (int s = 0; s < ds.getColObjects().size(); s++) {
//                datasetEV.getMatrix().setQuick(s, pca, pca1ExpEigenVector[s]);
//            }
//
//            int pcaNr = pca + 1;
//            cumExpVarPCA += expVarPCA;
//            out.write(pcaNr + "\t" + eigenValues[eigenValues.length - 1 - pca] + "\t" + expVarPCA + "\t" + cumExpVarPCA + "\n");
//            datasetEV.getColObjects().set(pca, "Comp" + String.valueOf(pcaNr));
//            System.out.println("PCA:\t" + pcaNr + "\t" + eigenValues[eigenValues.length - 1 - pca] + "\t" + expVarPCA + "\t" + cumExpVarPCA);
//        }
//        out.close();
//
//        datasetEV.save(outdir + "PrunedCovariateMatrix-PCAOverCovariatesEigenvectors.txt.gz");
//
//        System.out.println("Calculating PCs");
//        System.out.println("Initializing PCA matrix");
//        umcg.genetica.math.matrix.DoubleMatrixDataset<String, String> datasetPCAOverSamplesPCAs = new umcg.genetica.math.matrix.DoubleMatrixDataset<String, String>(ds.getRowObjects().size(), nrOfPCsToCalculate);
//        datasetPCAOverSamplesPCAs.rowObjects = ds.getRowObjects();
//        for (int s = 0; s < nrOfPCsToCalculate; s++) {
//            datasetPCAOverSamplesPCAs.colObjects.set(s, "Comp" + String.valueOf(s + 1));
//        }
//        for (int p = 0; p < ds.getRowObjects().size(); p++) {
//            for (int t = 0; t < nrOfPCsToCalculate; t++) {
//                datasetPCAOverSamplesPCAs.getRawData()[p][t] = 0;
//            }
//        }
//
//        ProgressBar pb = new ProgressBar(ds.getRowObjects().size(), "Calculating the PCA scores per sample: ");
//
//        int nrSamples = ds.getRowObjects().size();
//        int nrCovariates = ds.getColObjects().size();
//
//        double[][] evrawdata = datasetEV.getMatrix().toArray();
//        double[][] datasetPCAOverSamplesPCAsrawdata = datasetPCAOverSamplesPCAs.getRawData();
////        double[][] datasetrawdata = ds.getMatrix().toArray();
//
//        // multithread
//        Integer finalNrOfPCsToCalculate = nrOfPCsToCalculate;
//        IntStream.range(0, nrSamples).parallel().forEach(probe -> {
////			double[] probePCAs = datasetPCAOverSamplesPCAsrawdata[probe];
//            double[] probedata = datasetrawdata[probe];
//            for (int pc = 0; pc < finalNrOfPCsToCalculate; pc++) {
//                for (int sample = 0; sample < nrCovariates; sample++) {
//                    double probeCoefficient = evrawdata[sample][pc];
//                    datasetPCAOverSamplesPCAsrawdata[probe][pc] += probeCoefficient * probedata[sample];
//                }
//            }
//            pb.iterateSynched();
//        });


    }

    private void removeCategoricalColumns(String input, String outfile) throws IOException {

        Pair<boolean[], ArrayList<HashMap<String, Integer>>> catpair = getCategoricalColumns(input);
        boolean[] iscat = catpair.getLeft();
        int nrnoncat = 0;

        TextFile tf = new TextFile(input, TextFile.R);
        int nrSamples = tf.countLines() - 1;
        tf.close();
        tf.open();

        String[] header = tf.readLineElems(TextFile.tab);
        ArrayList<String> covariateNames = new ArrayList<>();
        for (int i = 0; i < header.length; i++) {
            if (!iscat[i]) {
                covariateNames.add(header[i]);
                nrnoncat++;
            }
        }


        double[][] vals = new double[nrSamples][nrnoncat];

        String[] elems = tf.readLineElems(TextFile.tab);
        ArrayList<String> samples = new ArrayList<>();
        HashSet<String> uniqueSamples = new HashSet<String>();
        int sampleCtr = 0;
        while (elems != null) {

            if (!uniqueSamples.contains(elems[0])) {

                double[] rowdata = new double[nrnoncat];
                int nrNan = 0;
                int colCtr = 0;
                for (int c = 0; c < elems.length; c++) {
                    if (!iscat[c]) {
                        double val = Double.NaN;
                        try {
                            val = Double.parseDouble(elems[c]);
                        } catch (NumberFormatException e) {
                        }

                        rowdata[colCtr] = val;
                        colCtr++;
                    }
                }

                // determine if row is variable or not

                vals[sampleCtr] = rowdata;
                samples.add(elems[0]);
                uniqueSamples.add(elems[0]);
                sampleCtr++;

            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        double[][] valspruned = new double[uniqueSamples.size()][];
        for (
                int c = 0; c < uniqueSamples.size(); c++) {
            valspruned[c] = vals[c];
        }


        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>();

        try {
            ds.setMatrix(valspruned); // samples on rows,
            ds.setRowObjects(samples);
            ds.setColObjects(covariateNames);
            ds.save(outfile);
        } catch (
                Exception e) {
            e.printStackTrace();
        }

    }

    private double removeNanAndCorrelate(DoubleMatrix1D colA, DoubleMatrix1D colB) {
        ArrayList<Double> x = new ArrayList<>();
        ArrayList<Double> y = new ArrayList<>();
        for (int c = 0; c < colA.cardinality(); c++) {
            double a = colA.getQuick(c);
            double b = colB.getQuick(c);

            if (!Double.isNaN(a) && !Double.isNaN(b)) {
                x.add(a);
                y.add(b);
            }
        }
        double[] xarr = Primitives.toPrimitiveArr(x);
        double[] yarr = Primitives.toPrimitiveArr(y);

        SpearmansCorrelation c = new SpearmansCorrelation();
        if (xarr.length > 500) {
            return c.correlation(xarr, yarr);
        } else {
            return Double.NaN;
        }


    }


    private Pair<boolean[], ArrayList<HashMap<String, Integer>>> getCategoricalColumns(String in) throws IOException {
        ArrayList<HashSet<String>> uniqueValuesPercolumn = getUniqueValsPerColumn(in);

        boolean[] isCategorical = new boolean[uniqueValuesPercolumn.size()];
        ArrayList<HashMap<String, Integer>> uniqueCatsPerCol = new ArrayList<>();

        int[] nrNonNan = new int[uniqueValuesPercolumn.size()];

        int nrCategorical = 0;
        for (int i = 0; i < isCategorical.length; i++) {
            uniqueCatsPerCol.add(new HashMap<>());
        }

        for (int i = 0; i < isCategorical.length; i++) {
            HashSet<String> uniqueVariables = uniqueValuesPercolumn.get(i);
            for (String s : uniqueVariables) {
                try {
                    double val = Double.parseDouble(s);
                    if (!Double.isNaN(val)) {
                        nrNonNan[i]++;
                    }

                } catch (NumberFormatException e) {
                    HashMap<String, Integer> uniqueCats = uniqueCatsPerCol.get(i);
                    if (uniqueCats == null) {
                        uniqueCats = new HashMap<>();
                    }
                    uniqueCats.put(s, uniqueCats.size());
                    uniqueCatsPerCol.set(i, uniqueCats);
                    isCategorical[i] = true;
                }
            }
        }

        System.out.println(nrCategorical + "  categorical variables");

        TextFile tf = new TextFile(in, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        tf.close();

        System.out.println();
        System.out.println("Categorical columns:");
        for (int i = 0; i < isCategorical.length; i++) {
            if (isCategorical[i]) {
                HashMap<String, Integer> uniqueCats = uniqueCatsPerCol.get(i);
                for (String key : uniqueCats.keySet()) {
                    System.out.println(header[i] + "\t" + key);
                }
                System.out.println();
            }
        }
        Pair<boolean[], ArrayList<HashMap<String, Integer>>> output = new Pair<>(isCategorical, uniqueCatsPerCol);
        return output;
    }

    private ArrayList<HashSet<String>> getUniqueValsPerColumn(String in) throws IOException {
        TextFile tf = new TextFile(in, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        ArrayList<HashSet<String>> uniqueValuesPercolumn = new ArrayList<HashSet<String>>();
        for (int i = 0; i < header.length; i++) {
            uniqueValuesPercolumn.add(new HashSet<>());
        }
        String[] elems = tf.readLineElemsReturnObjects(TextFile.tab);
        while (elems != null) {
            for (int i = 0; i < elems.length; i++) {
                if (!elems[i].equals("NA")) {
                    uniqueValuesPercolumn.get(i).add(elems[i]);
                }
            }
            elems = tf.readLineElemsReturnObjects(TextFile.tab);
        }
        tf.close();
        return uniqueValuesPercolumn;
    }


    public void removeNonVariantColumns(String in, int minNrOfUniqueValues, String out) throws IOException {
        ArrayList<HashSet<String>> uniqueValuesPercolumn = getUniqueValsPerColumn(in);

        TextFile tf = new TextFile(in, TextFile.R);
        String[] header = tf.readLineElems(TextFile.tab);
        int nrRemain = 0;
        boolean[] includecol = new boolean[header.length];
        for (int i = 0; i < includecol.length; i++) {
            HashSet<String> cols = uniqueValuesPercolumn.get(i);
            if (cols.size() > minNrOfUniqueValues) {
                includecol[i] = true;
                nrRemain++;
            }
        }

        System.out.println(nrRemain + " out of " + header.length + " columns remain");


        TextFile tfo = new TextFile(out, TextFile.W);
        String headerout = Strings.concat(header, includecol, Strings.tab);
        tfo.writeln(headerout);
        tf.open();
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);

        while (elems != null) {

            for (int i = 0; i < elems.length; i++) {
                if (elems[i].equals("NA")) {
                    elems[i] = "" + Double.NaN;
                }
            }

            String lnout = Strings.concat(elems, includecol, Strings.tab);
            tfo.writeln(lnout);

            elems = tf.readLineElems(TextFile.tab);
        }

        tfo.close();


    }
}
