package nl.harmjanwestra.playground.regressor;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.stream.IntStream;

public class Regressor2Test {


	public static void main(String[] args) {

		if (args.length < 3) {
			System.out.println("Usage: exp.txt.gz cov.txt.gz out.txt.gz");
		} else {
			Regressor2Test test = new Regressor2Test();
			try {
				test.test(args[0], args[1], args[2]);
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
//		String i1 = "D:\\tmp\\Cortex-EUR-PSI-filtered-logit-outliersRemoved.txt.gz";
//		String i2 = "D:\\tmp\\Cortex-EUR-covarsAnd30PCs.txt.gz";
//		String i3 = "D:\\tmp\\Cortex-EUR-covarsAnd30PCs-correl.txt.gz";
//		Regressor2Test t = new Regressor2Test();
//		try {
//			t.test(i1, i2, i3);
//		} catch (Exception e) {
//			throw new RuntimeException(e);
//		}
	}

	public void test(String expfile, String covfile, String output) throws Exception {

		System.out.println("Loading header: " + expfile);
		HashSet<String> samplesInExp = new HashSet<>();
		TextFile tf = new TextFile(expfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		tf.close();
		System.out.println(elems.length + " elements");

		for (String sample : elems) {
			if (sample.strip().length() > 0) {
				samplesInExp.add(sample);
			}
		}

		System.out.println("Loading header: " + covfile);
		HashSet<String> sharedSamples = new HashSet<>();
		tf = new TextFile(covfile, TextFile.R);
		elems = tf.readLineElems(TextFile.tab);
		tf.close();
		System.out.println(elems.length + " elements");

		for (String sample : elems) {
			if (sample.strip().length() > 0) {
				if (samplesInExp.contains(sample)) {
					sharedSamples.add(sample);
				}
			}
		}

		System.out.println(sharedSamples.size() + " shared samples detected sofar.");

		boolean transpose = false;
		if (sharedSamples.isEmpty() || sharedSamples.size() < 3) {
			System.out.println("Error: no or very few shared columns");
			System.out.println("Are the samples on the rows for " + covfile + "?");
			System.out.println("Checking now...");
			sharedSamples = new HashSet<>();
			tf = new TextFile(covfile, TextFile.R);
			tf.readLine();
			String ln = tf.readLine();
			while (ln != null) {
				elems = ln.split("\t");
				if (elems.length > 0) {
					if (samplesInExp.contains(elems[0])) {
						sharedSamples.add(elems[0]);
					}
				}
				ln = tf.readLine();
			}
			tf.close();

			if (sharedSamples.isEmpty() || sharedSamples.size() < 3) {
				System.out.println("Nope, transposing did not help. Please look whether the exp file has samples on columns.");
				System.exit(0);
			} else {
				System.out.println(sharedSamples.size() + " shared samples after transposing.");
				transpose = true;
			}
		}

		// match covariates to exp data
		// assume samples on columns first
		DoubleMatrixDataset expds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(expfile, '\t', null, sharedSamples);
		System.out.println("Exp data: " + expds.rows() + " x " + expds.columns());
		DoubleMatrixDataset cov = null;
		if (!transpose) {
			cov = DoubleMatrixDataset.loadSubsetOfTextDoubleData(covfile, '\t', null, sharedSamples);
		} else {
			cov = DoubleMatrixDataset.loadSubsetOfTextDoubleData(covfile, '\t', sharedSamples, null);
			cov = cov.viewDice();
		}
		System.out.println("Cov data: " + cov.rows() + " x " + cov.columns());

		// reindex the covariates
		cov.reorderCols(expds.getHashCols());

		// correlate each covariate with each gene
		TextFile tfp = new TextFile(output, TextFile.W);
		String binheader = "";
		double v = 0;
		DecimalFormat format = new DecimalFormat("#.##");
		for (int i = 0; i < 12; i++) {
			v = i * 0.1;
			binheader += "\tbin-" + format.format(v);
		}
		tfp.writeln("Covariate\tnrNaNs" + binheader);
		DoubleMatrixDataset finalCov = cov;

		ProgressBar pb = new ProgressBar(cov.rows() * expds.rows(), "Testing covariates * genes correlations");

		IntStream.range(0, cov.rows()).parallel().forEach(covRow -> {
			int[] bins = new int[12];
			int nrnans = 0;
			double[] covarvals = finalCov.getRow(covRow).toArray();
			for (int expRow = 0; expRow < expds.rows(); expRow++) {
				double[] expvals = expds.getRow(expRow).toArray();
				// strip missing values

				Pair<Double, Integer> corvals = Correlation.correlateWithNaNValues(covarvals, expvals);
				double r = corvals.getLeft();
				int n = corvals.getRight();
				if (!Double.isNaN(r)) {
					int bin = (int) Math.floor(Math.abs(r) * 10);
					bins[bin]++;
				} else {
					nrnans++;
				}
				pb.iterateSynched();
			}
			String outln = finalCov.getRowObjects().get(covRow) + "\t" + nrnans + "\t" + Strings.concat(bins, Strings.tab);
			System.out.println(outln);
			try {
				tfp.writelnsynced(outln);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}

		});
		tfp.close();
		pb.close();
	}
}
