package nl.harmjanwestra.methylation;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.VIF;
import umcg.genetica.text.Strings;

import java.util.*;

public class CovariateCorrect450K {

    public void correct(String input, String covariates, String output) throws Exception {
        DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(input);
        Pair<DoubleMatrixDataset<String, String>, ArrayList<String>> covarData = loadCovariateValues(covariates, it.getCols());
        ArrayList<String> sampleOrder = covarData.getRight();
        DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(sampleOrder, output);

        ArrayList<String> rowids = new ArrayList<>(it.getRows());
        ProgressBar pb = new ProgressBar(it.getNrRows(), "Removing covariates...");
        int rowid = 0;

        int[] colToCol = new int[it.getNrCols()];
        HashMap<String, Integer> colMap = new HashMap<String, Integer>();

        for (int i = 0; i < sampleOrder.size(); i++) {
            colMap.put(sampleOrder.get(i), i);
        }
        Set<String> curcols = it.getCols();
        int ctr = 0;
        for (String col : curcols) {
            Integer id = colMap.get(col);
            if (id == null) {
                colToCol[ctr] = -1;
            } else {
                colToCol[ctr] = id;
            }
            ctr++;
        }


        double[][] covariateDataMatrix = covarData.getLeft().viewDice().getMatrixAs2dDoubleArray();

        double[] data = new double[sampleOrder.size()];
        for (double[] row : it) {

            // select samples to include
            Arrays.fill(data, Double.NaN);
            select(row, colToCol, data);

            double[] y = data;
            OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
            ols.newSampleData(y, covariateDataMatrix);
            double[] yout = ols.estimateResiduals();

            writer.append(yout, rowids.get(rowid));
            pb.set(rowid);

            rowid++;
        }
        pb.close();

        it.close();
        writer.close();


    }

    private void select(double[] row, int[] colToCol, double[] data) {
        for (int c = 0; c < row.length; c++) {
            int id = colToCol[c];
            if (id != -1) {
                data[id] = row[c];
            }
        }
    }

    // NOTE: this new code switches around columns and rows for the covariate matrix
    private Pair<DoubleMatrixDataset<String, String>, ArrayList<String>> loadCovariateValues(String covariatesToRemove,
                                                                                             Set<String> datasetSamples) throws Exception {
        System.out.println("- Removing covariates as defined in: " + covariatesToRemove);
        TextFile covariates = new TextFile(covariatesToRemove, TextFile.R);
        int numRows = covariates.countLines() - 1; // minus the header :)
        int numCols = covariates.countCols(TextFile.tab) - 1; // minus the header's row identifier (if any)

        if (numRows == 0 || numCols == 0) {
            System.err.println("Covariate file is empty, but no covariates found in file! Is your file format correct?");
            System.err.println("The program is expecting the following: tab separated, one covariate per row, one sample per column, with sample identifiers identical to your --in file.");
            System.exit(0);
        } else {
            System.out.println("Covariate file has " + numRows + " rows and " + numCols + " columns");
        }


        // first hash up which samples are in the dataset
        HashMap<String, Integer> samplesInDatasetIndex = new HashMap<String, Integer>();
        String[] allSamplesInDataset = datasetSamples.toArray(new String[0]);
        for (int i = 0; i < allSamplesInDataset.length; i++) {
            samplesInDatasetIndex.put(allSamplesInDataset[i], i);
        }

        // read the column names from the covariate file
        // expect the samples on the columns
        String[] elems = covariates.readLineElemsReturnReference(TextFile.tab); // header

        int ctr = 0;
        boolean[] sampleInDatasetIncludedInCovariates = new boolean[datasetSamples.size()];
        ArrayList<String> columnNames = new ArrayList<String>();
        for (int i = 1; i < elems.length; i++) {
            Integer index = samplesInDatasetIndex.get(elems[i]);
            columnNames.add(elems[i]);
            if (index != null) {
                sampleInDatasetIncludedInCovariates[index] = true;
                ctr++;
            }
        }

        // read the covariate names, expect them to be on the rows
        ArrayList<String> rowNames = new ArrayList<String>();
        elems = covariates.readLineElemsReturnReference(TextFile.tab); // first line
        while (elems != null) {
            rowNames.add(elems[0]);
            elems = covariates.readLineElemsReturnReference(TextFile.tab);
        }
        covariates.close();

        boolean isTransposed = false;
        if (ctr == 0) {
            System.err.println("No matching samples detected between covariate file and dataset. Maybe your covariate file needs to be transposed? Will test that for you now:");
            for (String rowName : rowNames) {
                Integer index = samplesInDatasetIndex.get(rowName);
                if (index != null) {
                    sampleInDatasetIncludedInCovariates[index] = true;
                    ctr++;
                }
            }

            if (ctr == 0) {
                System.err.println("Transposing the data does not seem to resolve the issue. Please check your sample identifiers.");
                System.exit(0);
            } else {
                System.out.println("Transposing the covariate file reveals: " + ctr + " samples present.");
                isTransposed = true;

            }


        }

//        if (dataset.columns() != numSamples) {
//            System.out.println("Covariates loaded from: " + covariatesToRemove + ", but the number of samples does not correspond! " + numSamples + " in covariates file, " + dataset.columns() + " in dataset...");
//            System.out.println("Please note that missing samples will be removed from your eventual corrected --in file.");
//        }
        if (!isTransposed && ctr <= numRows + 2 || isTransposed && ctr <= numCols + 2) {
            System.err.println("Fewer samples present than minimally required for the normalization, (>covariates+3 samples needed).");
            System.exit(0);
        }
        if (ctr < datasetSamples.size()) {
            System.err.println("Covariates loaded from: " + covariatesToRemove + ", but not all samples present in covariates file! " + ctr + " present in covariates file, out of " + datasetSamples.size() + " in dataset...");
            System.out.println("Your dataset will be adjusted accordingly.");
        }
        int nrCovariates = numRows;
        if (isTransposed) {
            nrCovariates = numCols;
        }

        // make matrix with equal sample size
//        double[][] covariateValues = new double[nrCovariates][dataset.columns()];
        DoubleMatrixDataset<String, String> covariateValues = new DoubleMatrixDataset<>(nrCovariates, datasetSamples.size());
        covariateValues.getMatrix().assign(Double.NaN);

        int lineCtr = 0;
        covariates.open();
        String[] headerElems = covariates.readLineElemsReturnReference(TextFile.tab); // header
        elems = covariates.readLineElemsReturnReference(TextFile.tab);
        while (elems != null) {
            if (isTransposed) {
                String sampleName = elems[0];
                Integer sampleIdInDataset = samplesInDatasetIndex.get(sampleName);
                if (sampleIdInDataset != null) {
                    for (int i = 1; i < elems.length; i++) {
                        try {
                            covariateValues.setElementQuick(i - 1, sampleIdInDataset, Double.parseDouble(elems[i]));
                        } catch (NumberFormatException e) {
//                            System.out.println("WARNING: " + elems[i] + " is not a numeric value! in " + covariatesToRemove + " at line: " + (lineCtr + 1) + ".");
//                            covariateValues[i - 1][sampleIdInDataset] = Double.NaN;
//                            sampleInDatasetIncludedInCovariates[sampleIdInDataset] = false;
                        }
                    }
                }
            } else {
                for (int i = 1; i < elems.length; i++) {
                    String sampleName = headerElems[i];
                    Integer sampleIdInDataset = samplesInDatasetIndex.get(sampleName);
                    if (sampleIdInDataset != null) {
                        try {
                            covariateValues.setElementQuick(lineCtr, sampleIdInDataset, Double.parseDouble(elems[i]));
                        } catch (NumberFormatException e) {
//                            System.out.println("WARNING: " + elems[i] + " is not a numeric value at line: " + (lineCtr + 1) + "\tcolumn: " + i);
                        }
                    }
                }
            }
            elems = covariates.readLineElemsReturnReference(TextFile.tab);
            lineCtr++;
        }
        covariates.close();

        // investigate how many covariates there actually is data for.
        int covariateCtr = 0;

        boolean[] includeCovariate = new boolean[covariateValues.rows()];
        for (int row = 0; row < covariateValues.rows(); row++) {
            int nrColsFilled = 0;
            for (int col = 0; col < covariateValues.columns(); col++) {
                if (!Double.isNaN(covariateValues.getElementQuick(row, col))) {
                    nrColsFilled++;
                }
            }

            if (nrColsFilled == 0) {
                // there's no data for this covariate....
                includeCovariate[row] = false;
            } else {
                includeCovariate[row] = true;
                covariateCtr++;
            }
        }

        if (covariateCtr == 0) {
            System.err.println("ERROR: none of your covariates seem to have valid numerical values.. Please check your covariate file.");
            System.exit(0);
        }

        ArrayList<String> covariateNames = null;
        if (isTransposed) {
            covariateNames = columnNames;
        } else {
            covariateNames = rowNames;
        }

        if (covariateCtr != covariateValues.rows()) {
            // remove covariates with missing values
            System.out.println("Removing covariates that have no data at all.");
//            double[][] newCovariateData = new double[covariateCtr][dataset.columns()];
            DoubleMatrixDataset<String, String> newCovariateData = new DoubleMatrixDataset<>(covariateCtr, datasetSamples.size());
            ArrayList<String> newCovariateNames = new ArrayList<String>();
            int newCovariateCTR = 0;
            for (int row = 0; row < covariateValues.rows(); row++) {
                if (includeCovariate[row]) {
                    newCovariateNames.add(covariateNames.get(row));

                    for (int col = 0; col < covariateValues.columns(); col++) {
                        double val = covariateValues.getElementQuick(row, col);
                        newCovariateData.setElementQuick(newCovariateCTR, col, val);

                        // check whether we should include all samples, but don't remove yet: sync this with the expression/whatever dastaset
                        if (Double.isNaN(val)) {
                            sampleInDatasetIncludedInCovariates[col] = false;
                        }
                    }
                    newCovariateCTR++;
                } else {
                    System.out.println(covariateNames.get(row) + " removed.");
                }
            }


            nrCovariates = newCovariateCTR;
            covariateValues = newCovariateData;
            covariateNames = newCovariateNames;
        }

        // investigate how many samples there actually is data for.
        for (int row = 0; row < covariateValues.rows(); row++) {
            for (int col = 0; col < covariateValues.columns(); col++) {
                if (Double.isNaN(covariateValues.getElementQuick(row, col))) {
                    sampleInDatasetIncludedInCovariates[col] = false;
                }
            }
        }

        int sampleCtr = 0;
        for (int q = 0; q < sampleInDatasetIncludedInCovariates.length; q++) {
            if (sampleInDatasetIncludedInCovariates[q]) {
                sampleCtr++;
            }
        }

        // remove samples that have a missing value for at least one covariate
//        if (sampleCtr == sampleInDatasetIncludedInCovariates.length) {
//            System.out.println("There were no missing values or samples in your covariate file. Sample size will remain unchanged.");
//            DoubleMatrixDataset<String, String> covariateDataset = new DoubleMatrixDataset<String, String>(covariateValues, dataset.rowObjects, covariateNames);
//            return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(covariateDataset, dataset);
//        } else {

        System.out.println("Your covariate corrected dataset will have " + sampleCtr + " samples, after removing samples with missing covariate values.");


        // DoubleMatrixDataset<String, String> finalData = new DoubleMatrixDataset<>(dataset.rows(), sampleCtr);
        DoubleMatrixDataset<String, String> finalCovariates = new DoubleMatrixDataset<>(nrCovariates, sampleCtr);
        ArrayList<String> newColObjects = new ArrayList<String>();

        for (int col = 0; col < datasetSamples.size(); col++) {
            if (sampleInDatasetIncludedInCovariates[col]) {
                newColObjects.add(allSamplesInDataset[col]);
            }
        }

        for (int row = 0; row < covariateValues.rows(); row++) {
            int includedCovariateSampleCtr = 0;
            for (int col = 0; col < datasetSamples.size(); col++) {
                // replace covariate data...
                if (sampleInDatasetIncludedInCovariates[col]) {
                    finalCovariates.setElementQuick(row, includedCovariateSampleCtr, covariateValues.getElementQuick(row, col));
                    includedCovariateSampleCtr++;
                }
            }
        }

        finalCovariates.setRowObjects(covariateNames);
        finalCovariates.setColObjects(newColObjects);


        System.out.println("Checking variance of covariates");
        boolean inflated = true;
        int iter = 0;

        // remove zero variance covariates
        HashSet<Integer> skipRow = new HashSet<>();
        for (int row = 0; row < finalCovariates.rows(); row++) {
            double[] y = finalCovariates.getRow(row).toArray(); //[row];
            double var = JSci.maths.ArrayMath.variance(y);
            if (var == 0) {
                System.out.println("Variance for covariate " + finalCovariates.getRowObjects().get(row) + " == 0. Removing covariate.");
                skipRow.add(row);
            }
        }

        if (!skipRow.isEmpty()) {
            finalCovariates = excludeRows(finalCovariates, skipRow);
        }
        if (finalCovariates.rows() == 0) {
            System.err.println("ERROR: no covariates remain after removing zero variance covariates.");
            System.exit(-1);
        }

        // correct for variance inflation...
        VIF vif = new VIF();
        finalCovariates = vif.vifCorrect(finalCovariates.viewDice(), (1 - 1E-4)); // code here has covariates on the rows; move them to the columns instead.

        finalCovariates = finalCovariates.viewDice();
        System.out.println("");
        System.out.println("Remaining covariates: ");
        System.out.println(Strings.concat(finalCovariates.getRowObjects(), Strings.semicolon));
        System.out.println("");

        return new Pair<DoubleMatrixDataset<String, String>, ArrayList<String>>(finalCovariates, newColObjects);
//        }
    }

    private DoubleMatrixDataset<String, String> excludeRows(DoubleMatrixDataset<String, String> finalCovariates, HashSet<Integer> skipRow) throws Exception {
        // remove aliased row
        DoubleMatrixDataset<String, String> tmp = new DoubleMatrixDataset<>(finalCovariates.rows() - skipRow.size(), finalCovariates.columns());
        tmp.setColObjects(finalCovariates.getColObjects());
        ArrayList<String> rowObjs = new ArrayList<>();
        int rowctr = 0;
        for (int row = 0; row < finalCovariates.rows(); row++) {
            if (!skipRow.contains(row)) {
                for (int col = 0; col < finalCovariates.columns(); col++) {
                    tmp.setElementQuick(rowctr, col, finalCovariates.getElementQuick(row, col));
                }
                rowObjs.add(finalCovariates.getRowObjects().get(row));
                rowctr++;
            }
        }

        tmp.setRowObjects(rowObjs);
        return tmp;
    }
}
