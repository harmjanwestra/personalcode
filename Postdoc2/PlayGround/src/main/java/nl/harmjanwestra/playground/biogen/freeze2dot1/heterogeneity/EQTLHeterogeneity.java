package nl.harmjanwestra.playground.biogen.freeze2dot1.heterogeneity;

import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ChiSquare;
import umcg.genetica.math.stats.Heterogeneity;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

public class EQTLHeterogeneity {

	public static void main(String[] args) {

//        String indir = "D:\\Sync\\TMP\\suptable2\\";
//        String[] in = new String[]{
//                indir + "basalganglia-EUR-30PCs-1000perm-IterationsMerged.txt",
//                indir + "cerebellum-EUR-60PCs-1000perm-IterationsMerged.txt",
//                indir + "cortex-AFR-40PCs-1000perm-IterationsMerged.txt",
//                indir + "cortex-EAS-30PCs-merged-withqval-significant.txt",
//                indir + "cortex-EUR-80PCs-1000perm-IterationsMerged.txt",
//                indir + "hippocampus-EUR-30PCs-1000perm-IterationsMerged.txt",
//                indir + "spinalcord-EUR-20PCs-1000perm-IterationsMerged.txt"
//        };
//
//        for (int d = 0; d < in.length; d++) {
//            EQTLHeterogeneity c = new EQTLHeterogeneity();
//            try {
//                String out = in[d] + "-Heterogeneity.txt";
//                c.runBetaQTLFile(in[d], out);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }

//        }

		String indir = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2021-12-transeqtls\\";
		String[] in = new String[]{
				indir + "Cortex-EUR-AFR-noENA-0PCs-final-eQTLs-crossMappingEQTLsRemoved-FDR0.05.txt.gz",
				indir + "Cortex-EUR-AFR-noENA-100PCs-final-eQTLs-crossMappingEQTLsRemoved-FDR0.05.txt.gz",
				indir + "Cortex-EUR-AFR-noENA-noAMPAD-0PCs-final-eQTLs-crossMappingEQTLsRemoved-FDR0.05.txt.gz",
				indir + "Cortex-EUR-AFR-noENA-noAMPAD-80PCs-final-eQTLs-crossMappingEQTLsRemoved-FDR0.05.txt.gz",
		};

		in = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2022-02-CellTypeGWAS\\eQTLsFDR0.05.txt"
		};

		in = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2022-02-Heterogeneity\\cortex-EURandAFR-noENA-100PCs-1000perm-IterationsMerged.txt.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2022-02-Heterogeneity\\cortex-EURandAFR-noENA-noAMPAD-80PCs-1000perm-IterationsMerged.txt.gz"
		};

		for (int d = 0; d < in.length; d++) {
			EQTLHeterogeneity c = new EQTLHeterogeneity();
			try {
				String out = in[d] + "-Heterogeneity.txt";
//				c.runBetaQTLFile(in[d], out);
				c.runBetaQTLFile(in[d], "d:\\test.txt");
//                c.runMetaQTLFile(in[d], out);
			} catch (IOException e) {
				e.printStackTrace();
			}

		}
	}

	public void runBetaQTLFile(String file, String output) throws IOException {
		Heterogeneity het = new Heterogeneity();
		// Heterogeneity.getISq(zscores, samplesizes);

		TextFile tf = new TextFile(file, TextFile.R);
		String header = tf.readLine();
		header += "\tISqR\tiSQZ\tiSQZfix";
		String[] elems = tf.readLineElems(TextFile.tab);
		TextFile tfo = new TextFile(output, TextFile.W);
		tfo.writeln(header);
		Locale currentLocale = Locale.US;
		DecimalFormatSymbols unusualSymbols = new DecimalFormatSymbols(currentLocale);
		DecimalFormat df = new DecimalFormat("#.###", unusualSymbols);
		while (elems != null) {
			String[] datasetCorStr = elems[17].split(";");
			String[] datasetNStr = elems[19].split(";");
			String[] datasetZStr = elems[18].split(";");

			if (datasetCorStr.length < 2) {
				System.out.println("Error: " + datasetCorStr.length);
				System.exit(0);
			}

			double[] corr = new double[datasetCorStr.length];
			double[] z = new double[datasetCorStr.length];
			int[] n = new int[datasetCorStr.length];
			for (int d = 0; d < corr.length; d++) {
				try {
					corr[d] = Double.parseDouble(datasetCorStr[d]);
					z[d] = Double.parseDouble(datasetZStr[d]);
					n[d] = Integer.parseInt(datasetNStr[d]);
				} catch (NumberFormatException e) {
					corr[d] = Double.NaN;
					z[d] = Double.NaN;
					n[d] = 0;
				}
			}
			// Triple<double[], double[], int[]> stripped = stripnan(corr, z, n);

			double[] heterogeneityR = Heterogeneity.isqFromCorrAndNFixedEffect(corr, n);
			Triple<Double, Double, Integer> heterogeneityZ = Heterogeneity.getISq(z, n);
			Triple<Double, Double, Integer> heterogeneityZ2 = getISqFix(z, n);
//            double q = heterogeneity[4];
			double i2R = heterogeneityR[7];
			double i2Z = heterogeneityZ.getLeft();
			double i2Z2 = heterogeneityZ2.getLeft();

			tfo.writeln(Strings.concat(elems, Strings.tab) + "\t" + df.format(i2R) + "\t" + df.format(i2Z) + "\t" + df.format(i2Z2));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfo.close();
	}

	public void runMetaQTLFile(String file, String output) throws IOException {
		Heterogeneity het = new Heterogeneity();
		// Heterogeneity.getISq(zscores, samplesizes);

		TextFile tf = new TextFile(file, TextFile.R);
		String header = tf.readLine();
		header += "\tISq";
		String[] elems = tf.readLineElems(TextFile.tab);
		TextFile tfo = new TextFile(output, TextFile.W);
		tfo.writeln(header);
		Locale currentLocale = Locale.US;
		DecimalFormatSymbols unusualSymbols = new DecimalFormatSymbols(currentLocale);
		DecimalFormat df = new DecimalFormat("#.###", unusualSymbols);
		while (elems != null) {

            /*
0 PValue
1 SNPName
2 SNPChr
3 SNPChrPos
4 ProbeName
5 ProbeChr
6 ProbeCenterChrPos
7 CisTrans
8 SNPType
9 AlleleAssessed
0 OverallZScore
11 DatasetsWhereSNPProbePairIsAvailableAndPassesQC
12 DatasetsZScores
13 DatasetsNrSamples
14 IncludedDatasetsMeanProbeExpression
15 IncludedDatasetsProbeExpressionVariance
16 HGNCName
17 IncludedDatasetsCorrelationCoefficient
18 Meta-Beta (SE)
19 Beta (SE)
20 FoldChange
21 FDR
             */

			String[] datasetCorStr = elems[17].split(";");
			String[] datasetNStr = elems[13].split(";");
			String[] datasetZStr = elems[12].split(";");

			double[] corr = new double[datasetCorStr.length];
			double[] z = new double[datasetCorStr.length];
			int[] n = new int[datasetCorStr.length];
			for (int d = 0; d < corr.length; d++) {
				try {
					corr[d] = Double.parseDouble(datasetCorStr[d]);
					z[d] = Double.parseDouble(datasetZStr[d]);
					n[d] = Integer.parseInt(datasetNStr[d]);
				} catch (NumberFormatException e) {
					corr[d] = Double.NaN;
					z[d] = Double.NaN;
					n[d] = 0;
				}
			}
			// Triple<double[], double[], int[]> stripped = stripnan(corr, z, n);

			double[] heterogeneity = Heterogeneity.isqFromCorrAndNFixedEffect(corr, n);
//            double q = heterogeneity[4];
			double i2 = heterogeneity[7];


			tfo.writeln(Strings.concat(elems, Strings.tab) + "\t" + df.format(i2));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfo.close();
	}


	/*
	 * This class was freely adapted from Abecasis' METAL package for meta-analysis
	 * (http://www.sph.umich.edu/csg/abecasis/Metal/)
	 * it should produce output that is highly correlated with the isqFromCorrAndNFixedEffect function below
	 * datasetZ: array of z-scores per dataset, NaN for missing
	 * datasetWeights: array of sample sizes (not squared)
	 */

	/*
	 * This class was freely adapted from Abecasis' METAL package for meta-analysis
	 * (http://www.sph.umich.edu/csg/abecasis/Metal/)
	 * it should produce output that is highly correlated with the isqFromCorrAndNFixedEffect function below
	 * datasetZ: array of z-scores per dataset, NaN for missing
	 * datasetWeights: array of sample sizes (not squared)
	 */
	public static Triple<Double, Double, Integer> getISqFix(double[] datasetZ, int[] datasetWeights) {

		double weightedZ = 0;
		int totalSample = 0;
		for (int d = 0; d < datasetZ.length; d++) {
			if (!Double.isNaN(datasetZ[d])) {
				weightedZ += Math.sqrt(datasetWeights[d]) * datasetZ[d];
				totalSample += datasetWeights[d];
			}
		}

		double hetSum = 0;
		int hetDf = 0;
		weightedZ/=Math.sqrt(totalSample);
		for (int d = 0; d < datasetZ.length; d++) {
			if (!Double.isNaN(datasetZ[d])) {
				double expectedZ = Math.sqrt(datasetWeights[d]) * weightedZ;

				hetSum += (datasetZ[d] - expectedZ) * (datasetZ[d] - expectedZ);
				hetDf++;
			}
		}

		double p = 1d;
		double i = 0d;

		double iExp = ((hetSum - hetDf + 1) / hetSum) * 100d;
		if (hetDf <= 1 || hetSum < 1E-7) {
			p = 1;
		} else {
			p = ChiSquare.getP(hetDf - 1, hetSum);
		}


		if (hetSum <= (hetDf - 1) || hetDf <= 1) {
			i = 0;
		} else {
			i = iExp;
		}

		if (i > 0) {
			i /= 100;
		}
		return new Triple<Double, Double, Integer>(i, p, totalSample);
	}
}
