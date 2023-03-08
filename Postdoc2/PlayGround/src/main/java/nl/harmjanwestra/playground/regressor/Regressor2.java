package nl.harmjanwestra.playground.regressor;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.VIF;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.stream.IntStream;

public class Regressor2 {

	double defaultVIFthreshold = (1 - 1E-4);

	public static void main(String[] args) {
		if (args.length < 3) {
			System.out.println("Usage: expfile.txt.gz covfile.txt.gz output.txt.gz");
		} else {
			try {
				Regressor2 r = new Regressor2();

				DoubleMatrixDataset<String, String> dataset = DoubleMatrixDataset.loadDoubleData(args[0]);
				r.adjustCovariates(dataset, args[2], args[1]);
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
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

	// NOTE: this new code switches around columns and rows for the covariate matrix
	private Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> loadCovariateValues(String covariatesToRemove,
																											   DoubleMatrixDataset<String, String> dataset) throws Exception {
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
		String[] allSamplesInDataset = dataset.getColObjects().toArray(new String[0]);
		for (int i = 0; i < allSamplesInDataset.length; i++) {
			samplesInDatasetIndex.put(allSamplesInDataset[i], i);
		}

		// read the column names from the covariate file
		// expect the samples on the columns
		String[] elems = covariates.readLineElemsReturnReference(TextFile.tab); // header

		int ctr = 0;
		boolean[] sampleInDatasetIncludedInCovariates = new boolean[dataset.columns()];
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
		if (ctr < dataset.columns()) {
			System.err.println("Covariates loaded from: " + covariatesToRemove + ", but not all samples present in covariates file! " + ctr + " present in covariates file, out of " + dataset.columns() + " in dataset...");
			System.out.println("Your dataset will be adjusted accordingly.");
		}
		int nrCovariates = numRows;
		if (isTransposed) {
			nrCovariates = numCols;
		}

		// make matrix with equal sample size
//        double[][] covariateValues = new double[nrCovariates][dataset.columns()];
		DoubleMatrixDataset<String, String> covariateValues = new DoubleMatrixDataset<>(nrCovariates, dataset.columns());
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
			DoubleMatrixDataset<String, String> newCovariateData = new DoubleMatrixDataset<>(covariateCtr, dataset.columns());
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


		DoubleMatrixDataset<String, String> finalData = new DoubleMatrixDataset<>(dataset.rows(), sampleCtr);
		DoubleMatrixDataset<String, String> finalCovariates = new DoubleMatrixDataset<>(nrCovariates, sampleCtr);
		ArrayList<String> newColObjects = new ArrayList<String>();

		for (int col = 0; col < dataset.columns(); col++) {
			if (sampleInDatasetIncludedInCovariates[col]) {
				newColObjects.add(dataset.getColObjects().get(col));
			}
		}

		for (int row = 0; row < dataset.rows(); row++) {
			int includedSampleCtr = 0;
			for (int col = 0; col < dataset.columns(); col++) {
				if (sampleInDatasetIncludedInCovariates[col]) {
					// include sample
					finalData.setElementQuick(row, includedSampleCtr, dataset.getElementQuick(row, col));
					includedSampleCtr++;
				}
			}
		}

		for (int row = 0; row < covariateValues.rows(); row++) {
			int includedCovariateSampleCtr = 0;
			for (int col = 0; col < dataset.columns(); col++) {
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
		vif.setDebug(true);
		finalCovariates = vif.vifCorrect(finalCovariates.viewDice(), defaultVIFthreshold); // code here has covariates on the rows; move them to the columns instead.

		finalCovariates = finalCovariates.viewDice(); // covariates on rows again
		System.out.println("");
		System.out.println("Remaining covariates: ");
		System.out.println(Strings.concat(finalCovariates.getRowObjects(), Strings.semicolon));
		System.out.println("");

//        finalCovariates.save(covariatesToRemove + "-asLoadedByNormalizer.txt");
		finalData.setRowObjects(dataset.getRowObjects());
		finalData.setColObjects(newColObjects);
//        finalData.save(dataset + "-SampleSizeCorrectedForCovariates.txt");
		return new Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>>(finalCovariates, finalData);
//        }
	}

	public String adjustCovariates(DoubleMatrixDataset<String, String> traitData,
								   String fileNamePrefix,
								   String covariatesToRemove
	) throws IOException, Exception {

		// load covariate data, and remove samples for which there is missing covariate data.
		Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> traitAndCovarData = loadCovariateValues(covariatesToRemove, traitData);
//		DoubleMatrixDataset<String, String> covariateDataset = covariateData.getLeft();

		//		DoubleMatrixDataset<String, String> covariateDatasetTranspose = covariateDataset.viewDice();


		traitData.setMatrix(traitAndCovarData.getRight().getMatrix());        // samples on columns, features on rows
		traitData.setColObjects(traitAndCovarData.getRight().getColObjects());
		traitData.setRowObjects(traitAndCovarData.getRight().getRowObjects());

		DoubleMatrixDataset<String, String> covariateDataset = traitAndCovarData.getLeft(); // covariates on rows, samples on columns
		DoubleMatrixDataset<String, String> covariateDatasetTranspose = covariateDataset.viewDice(); // covariates on columns, samples on rows
		double[][] covariateDataMatrix = covariateDatasetTranspose.getMatrixAs2dDoubleArray();

		// use Apache multivariate OLS in stead of PCA for covariate correction
		// move the samples to the rows, covariates on columns
//		double[][] covariateDataMatrix = new double[traitData.columns()][covariateDataset.rows()];
//
//		int[] sampleMap = new int[traitData.columns()];
//		for (int s = 0; s < covariateDataset.getColObjects().size(); s++) {
//			String sample = covariateDataset.getColObjects().get(s);
//			Integer id = traitData.getHashCols().get(sample);
//			sampleMap[s] = id;
//		}
//
//		for (int row = 0; row < covariateDataset.rows(); row++) {
//			for (int col = 0; col < covariateDataset.columns(); col++) {
//				Integer sampleid = sampleMap[col];
//				covariateDataMatrix[sampleid][row] = covariateDataset.getElementQuick(row, col);
//			}
//		}

		DoubleMatrixDataset<String, String> outputmat = new DoubleMatrixDataset<>(traitData.rows(), traitData.columns());
		outputmat.setRowObjects(traitData.getRowObjects());
		outputmat.setColObjects(traitData.getColObjects());
		DoubleMatrixDataset<String, String> finalTraitData = traitData;

		HashSet<Integer> dumpRows = new HashSet<Integer>();
		ProgressBar pb = new ProgressBar(traitData.rows(), "Calculating OLS residuals...");
		DoubleMatrixDataset<String, String> finalOutputmat = outputmat;
		IntStream.range(0, traitData.rows()).parallel().forEach(row -> {

			double[] y = finalTraitData.getRow(row).toArray();
			String rowId = traitData.getRowObjects().get(row);

			// y may contain missing values, same for covariateDataMatrix
			boolean[] isNonNan = new boolean[y.length];
			int nanctr = 0;
			ArrayList<Double> yNonNan = new ArrayList<>();
			HashSet<Integer> nanSamples = new HashSet<>();
			for (int i = 0; i < y.length; i++) {
				double v = y[i];
				if (Double.isNaN(v)) {
					nanctr++;
					nanSamples.add(i);
				} else {
					isNonNan[i] = true;
					yNonNan.add(v);
				}
			}

			if (nanctr > 0) {
				// row contains nan values, subset the covariate data
				try {

					DoubleMatrixDataset<String, String> covariatesNonNan = excludeRows(covariateDatasetTranspose, nanSamples);


					// subsetting may introduce collinearity, double check for nonvariant and collinear columns

//					System.out.println(yNonNan.size() + "\t" + covariatesNonNan.rows() + "\t" + covariatesNonNan.columns());

					if (covariatesNonNan.columns() >= (yNonNan.size() - 1)) {
						System.out.println("Warning: more covariates (" + covariatesNonNan.columns() + ") than measurements (" + yNonNan.size() + ") for " + rowId + " - removing from dataset");
						synchronized (dumpRows) {
							dumpRows.add(row);
						}
					} else {
						VIF v = new VIF();
						covariatesNonNan = v.vifCorrect(covariatesNonNan, defaultVIFthreshold);

						double[][] covariateDataNonNan = covariatesNonNan.getMatrixAs2dDoubleArray();
						double[] yNonNanArr = Primitives.toPrimitiveArr(yNonNan);

						OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
						ols.setNoIntercept(false);
						ols.newSampleData(yNonNanArr, covariateDataNonNan);
						double[] yout = ols.estimateResiduals();
						int cctr = 0;
						for (int c = 0; c < y.length; c++) {
							if (isNonNan[c]) {
								finalOutputmat.setElementQuick(row, c, yout[cctr]);
								cctr++;
							} else {
								finalOutputmat.setElementQuick(row, c, Double.NaN);
							}
						}
					}


				} catch (Exception e) {
					throw new RuntimeException(e);
				}

			} else {
				if (covariateDatasetTranspose.columns() >= (y.length - 1)) {
					System.out.println("Warning: more covariates (" + covariateDatasetTranspose.columns() + ") than measurements (" + y.length + ") for " + rowId + " - removing from dataset");
					synchronized (dumpRows) {
						dumpRows.add(row);
					}
				} else {
					OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
					ols.setNoIntercept(false);
					ols.newSampleData(y, covariateDataMatrix);
					double[] yout = ols.estimateResiduals();

					for (int c = 0; c < yout.length; c++) {
						finalOutputmat.setElementQuick(row, c, yout[c]);
					}
				}
			}
			pb.iterateSynched();
		});
		pb.close();

		// dump some rows from the outputmat
		if (!dumpRows.isEmpty()) {
			outputmat = excludeRows(finalOutputmat, dumpRows);
		}

		traitData.setMatrix(outputmat.getMatrix());
		fileNamePrefix += ".CovariatesRemovedOLS";


		System.out.println("Saving: " + fileNamePrefix);

		save(traitData, fileNamePrefix + ".txt.gz");

		System.out.println("Done");
		return null;
	}

	private void save(DoubleMatrixDataset<String, String> traitData, String file) throws IOException {
		TextFile out = new TextFile(file, true);
		String header = "-\t" + Strings.concat(traitData.getColObjects(), Strings.tab);

		out.writeln(header);
		ProgressBar pb = new ProgressBar(traitData.rows());
		for (int r = 0; r < traitData.rows(); r++) {

			double[] rowData = traitData.getRow(r).toArray();
			String outln = traitData.getRowObjects().get(r) + "\t" + Strings.concat(rowData, Strings.tab);
			out.writeln(outln);
			pb.set(r);
		}
		pb.close();
		out.close();
	}
}
