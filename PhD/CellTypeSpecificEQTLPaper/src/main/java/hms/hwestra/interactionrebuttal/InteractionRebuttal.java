/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hms.hwestra.interactionrebuttal;

import eqtlmappingpipeline.interactionanalysis.InteractionAnalysisMultiThreaded;
import eqtlmappingpipeline.normalization.Normalizer;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.TTest;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author hwestra
 */
public class InteractionRebuttal {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {

		InteractionRebuttal r = new InteractionRebuttal();
		try {
			r.run();
			// TODO code application logic here
		} catch (IOException ex) {
			Logger.getLogger(InteractionRebuttal.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	public void run() throws IOException {

//        String file2 = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/EGCUTData/EGCUTValidSamplesNeutrosAndTop58Proxy.txt";
//        String file1 = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/2014-10-28-EGCUTPCs/ExpressionData.txt.gz.PCAOverSamplesEigenvectors.txt";
//        String out = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/1.2.PCsAndNeutros/PCsVersusNeutrophils";
//        correlate(file1, true, file2, false, false, null, out);
//        String annotation = "/Sync/AeroFS/cellTypeeQTL/DataFiles/2013-07-18-HT12v3.txt";
//        String in = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/1.3.Correlations58Probes/Top58ProbesAndCorrelation.txt";
//        String out = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/1.3.Correlations58Probes/Top58ProbesAndCorrelationWithGenes.txt";
//        annotate(in, annotation, out);
//        String fileIn = "/Data/tmp/2012-12-21-CisAssociationsProbeLevelFDR0.5.txt";
//        String fileOut = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/2.5.CisEQTLEffectSize/CisEQTLsWEffectSize.txt";
//        determineCorrelationFromZScore(fileIn, fileOut);
//        String file1 = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/2.5.CisEQTLEffectSize/CellTypeSpecificEQTLs.txt";
//        String file2 = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/2.5.CisEQTLEffectSize/CisEQTLsWEffectSize.txt";
//        String fileout = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/2.5.CisEQTLEffectSize/CellTypeSpecificEQTLs-WithR.txt";
//        intersectEQTLFiles(file1, file2, fileout);
//        Correlation.correlationToZScore(825);
//        double correlation = -0.16;
//        double z = Correlation.convertCorrelationToZScore(825, correlation);
//        double p = ZScores.zToP(z);
//        System.out.println(correlation + "\t" + z + "\t" + p);
//        String f="/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/1.1.Endophenotypes/EGCUTValidSamplesNeutrosMonocytesAgeGenderAndCellProxyWithNOD2WTop58ProbeProxy.txt";
//        String fout = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/1.1.Endophenotypes/plots/out";
//        plotBoxPlots(f, "Sex", fout);
		String normal1 = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/3.4.RobustSEInteraction/EGCUT/InteractionResults.txt";
		String robust1 = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/3.4.RobustSEInteraction/EGCUTRobust/InteractionResultsRobust.txt";
		String normal2 = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/3.4.RobustSEInteraction/Groningen/GRNGInteractionResults.txt";
		String robust2 = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/3.4.RobustSEInteraction/Groningen-Robust/GRNGInteractionResultsRobust.txt";
		String out = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/3.4.RobustSEInteraction/Merged.txt";
//        merge(normal1, robust1, normal2, robust2, out);
//



//		spearmanSE(normal1, robust1);
//		spearmanSE(normal2, robust2);
//        String expressionfile = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/EGCUTData/EGCUT-RawDataQNLog2.txt";
//        String gte = null; //"/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/pfff.txt";
//        String outputdir = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/1.5.IterativeProxyProbes-Top100Probes/";
//        String probeFile = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/2.2.ExpressionVsNeutros/Top100probes.txt";
//        int permutations = 1000;
//        String proxy = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/EGCUTData/EGCUTValidSamplesNeutrosOnly.txt";
//        Integer numprobesmax = 100;
//        iterativelyBuildProxiesAndCorrelate(expressionfile, gte, outputdir, probeFile, permutations, proxy, numprobesmax, true);

	}

	public void iterativelyBuildProxiesAndCorrelate(String expressionfile, String gte,
													String outputdir, String probefile, int permutations,
													String proxy, Integer numProbesMax, boolean justcorrelate) throws IOException {

		Set<String> samplesToUse = null;
		if (gte != null) {
			TextFile tfq = new TextFile(gte, TextFile.R);
			samplesToUse = tfq.readAsSet(1, TextFile.tab);
			tfq.close();
		}
		DoubleMatrixDataset<String, String> rawExpressionDataset = new DoubleMatrixDataset<String, String>(expressionfile, null, samplesToUse);
		int numTotalIter = 0;

		ArrayList<String> probes = null;
		if (probefile != null) {
			TextFile tf = new TextFile(probefile, TextFile.R);
			probes = tf.readAsArrayList();
			tf.close();
			numProbesMax = probes.size();
		} else {
			System.out.println("Selecting random probes");
			List<String> rows = rawExpressionDataset.rowObjects;
			probes = new ArrayList<String>(rows);
		}

		outputdir = Gpio.formatAsDirectory(outputdir);
		Gpio.createDir(outputdir);
		int iter = 5;
		System.out.println(probes.size() + " probes availables");

		int remainder = numProbesMax % iter;
		int numProbesMaxIter = numProbesMax + (iter - remainder);

		if (!justcorrelate) {

			for (int num = 0; num < numProbesMaxIter + 1; num += iter) {
				int probesToSelect = num;
				if (num == 0) {
					probesToSelect = 1;
				}
				if (num > numProbesMax) {
					probesToSelect = numProbesMax;
				}
				System.out.println("Selecting: " + probesToSelect + " probes");
				for (int permutation = 0; permutation < permutations; permutation++) {
					Collections.shuffle(probes);
					List<String> subsample = probes.subList(0, probesToSelect);

					// create output dir
					String outputdirPerm = outputdir + probesToSelect + "-Probes/Permutation-" + permutation + "/";
					outputdirPerm = Gpio.formatAsDirectory(outputdirPerm);
					Gpio.createDir(outputdirPerm);
					String subset = outputdirPerm + "probes.txt";
					TextFile probeout = new TextFile(subset, TextFile.W);
					probeout.writeList(subsample);
					probeout.close();

					// run normalizer
					prepareDataForCelltypeSpecificEQTLMapping(rawExpressionDataset, expressionfile, outputdirPerm, Double.NaN, subset, null, null, null, 4);
					// remove superfluous files
					// correlate with cell count
				}
				numTotalIter++;
			}
		}

		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(proxy); // samples on the rows...
		ds.transposeDataset(); // samples on the columns
		for (int row = 0; row < ds.nrRows; row++) {
			String pheno = ds.rowObjects.get(row);

			double[] x = ds.rawData[row];
			System.out.println("x length: " + x.length);

			TextFile statsout = new TextFile(outputdir + pheno + ".txt", TextFile.W);
			statsout.writeln("Num\tMeanPearson\tsdPearson\tMeanSpearman\tsdSpearman");
			SpearmansCorrelation sp = new SpearmansCorrelation();
			for (int num = 0; num < numProbesMaxIter + 1; num += iter) {
				int probesToSelect = num;
				if (num == 0) {
					probesToSelect = 1;
				}
				if (num > numProbesMax) {
					probesToSelect = numProbesMax;
				}
				double[] allCorrelations = new double[permutations];
				double[] allCorrelationSpearman = new double[permutations];
				for (int permutation = 0; permutation < permutations; permutation++) {
					String inputdirPerm = outputdir + probesToSelect + "-Probes/Permutation-" + permutation + "/CellTypeProxyFile.txt";
					DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(inputdirPerm);
					ds2.transposeDataset(); // samples on the column
					double[] y = new double[x.length];
					System.out.println("y: " + y.length);
					double[] ytmp = ds2.rawData[0];
//                    if (ytmp.length != x.length) {
//                        System.err.println("Error: " + y.length);
//                        System.exit(-1);
//                    } else {
					for (int col = 0; col < ds.nrCols; col++) {
						int otherCol = ds2.hashCols.get(ds.colObjects.get(col));
						y[col] = ytmp[otherCol];
					}
					double corr = JSci.maths.ArrayMath.correlation(x, y);
					System.out.println(num + "\t" + permutation + "\t" + corr);
					double spearman = sp.correlation(x, y);
					allCorrelations[permutation] = corr;
					allCorrelationSpearman[permutation] = spearman;
//                    }
				}
				// doe
				double meanP = JSci.maths.ArrayMath.mean(allCorrelations);
				double sdP = JSci.maths.ArrayMath.standardDeviation(allCorrelations);
				double meanSP = JSci.maths.ArrayMath.mean(allCorrelationSpearman);
				double sdSP = JSci.maths.ArrayMath.standardDeviation(allCorrelationSpearman);
				statsout.writeln(num + "\t" + meanP + "\t" + sdP + "\t" + meanSP + "\t" + sdSP);
			}

			statsout.close();

		}

	}

	public void correlate(String file1, boolean transpose1, String file2, boolean transpose2, boolean addN, String probeFilter, String outfile) throws IOException {

		// exp has individuals on columns, check whether this is also the case for the cel
		// otherwise, transpose, make index
		int ct = 0;

		Set<String> probesToFilter = null;
		if (probeFilter != null) {
			TextFile tf1 = new TextFile(probeFilter, TextFile.R);
			probesToFilter = tf1.readAsSet(0, TextFile.tab);
			tf1.close();
		}

		DoubleMatrixDataset<String, String> exp = new DoubleMatrixDataset<String, String>(file1, probesToFilter);
		if (transpose1) {
			exp.transposeDataset();
		}
		DoubleMatrixDataset<String, String> cel = new DoubleMatrixDataset<String, String>(file2);
		if (transpose2) {
			cel.transposeDataset();
		}

		List<String> colsCel = cel.colObjects;

		for (String s : colsCel) {
			if (exp.hashCols.containsKey(s)) {

				ct++;
			}
		}
		if (ct == 0) {
			cel.transposeDataset();
			colsCel = cel.colObjects;
			for (String s : colsCel) {
				if (exp.hashCols.containsKey(s)) {
					ct++;
				}
			}
		}

		int[] index = new int[exp.nrCols];
		for (int i = 0; i < exp.nrCols; i++) {
			Integer id = cel.hashCols.get(exp.colObjects.get(i));
			if (id == null) {
				index[i] = -1;
			} else {
				index[i] = id;
			}
		}
		System.out.println("Found " + ct + " matching columns");
		SpearmansCorrelation corr = new SpearmansCorrelation();
		TextFile outputfile1 = new TextFile(outfile + "-Spearman.txt", TextFile.W);
		TextFile outputfile2 = new TextFile(outfile + "-Pearson.txt", TextFile.W);

		String header = "-";
		for (int row2 = 0; row2 < cel.nrRows; row2++) {
			if (addN) {
				header += "\t" + cel.rowObjects.get(row2) + "\t" + cel.rowObjects.get(row2) + "-N";
			} else {
				header += "\t" + cel.rowObjects.get(row2);
			}
		}
		outputfile1.writeln(header);
		outputfile2.writeln(header);

		for (int row = 0; row < exp.nrRows; row++) {

			String spearmanout = exp.rowObjects.get(row);
			String pearsonout = exp.rowObjects.get(row);

			for (int row2 = 0; row2 < cel.nrRows; row2++) {

				int ct2 = 0;
				for (int col = 0; col < exp.nrCols; col++) {
					int id = index[col];
					if (id != -1) {
						// check whether value != null;
						if (!Double.isNaN(exp.rawData[row][col])
								&& !Double.isNaN(cel.rawData[row2][id])) {
							ct2++;
						}
					}
				}

				double[] valsx = new double[ct2];
				double[] valsy = new double[ct2];
				int z = 0;
				for (int col = 0; col < exp.nrCols; col++) {
					int id = index[col];
					if (id != -1) {
						// check whether value != null;
						if (!Double.isNaN(exp.rawData[row][col])
								&& !Double.isNaN(cel.rawData[row2][id])) {
							valsx[z] = exp.rawData[row][col];
							valsy[z] = cel.rawData[row2][id];
							z++;
						}
					}
				}
				valsx = ztransform(valsx);
				valsy = ztransform(valsy);
				double spearman = corr.correlation(valsx, valsy);
				double pearson = JSci.maths.ArrayMath.correlation(valsx, valsy);
				if (addN) {
					pearsonout += "\t" + pearson + "\t" + ct2;
					spearmanout += "\t" + spearman + "\t" + ct2;
				} else {
					pearsonout += "\t" + pearson;
					spearmanout += "\t" + spearman;
				}
			}
			outputfile1.writeln(spearmanout);
			outputfile2.writeln(pearsonout);

		}
		outputfile1.close();
		outputfile2.close();
	}

	public void runIterativeNormalizer(String f,
									   String expression,
									   String mdscomponents,
									   String cellcountfile,
									   String gte,
									   String ingt,
									   String snpprobecombos,
									   String inexpPCCorrected,
									   Integer threads,
									   String outdirectory) throws Exception {
		TextFile tf = new TextFile(f, TextFile.R);
		String[] data = tf.readAsArray();
		tf.close();
		outdirectory = Gpio.formatAsDirectory(outdirectory);
		Gpio.createDir(outdirectory);
		int q = 1;

		while (q < data.length) {

			String probeFileName = outdirectory + q + "probes.txt";
			TextFile out = new TextFile(probeFileName, TextFile.W);
			System.out.println(outdirectory + q + "probes.txt");
			for (int z = 0; z < q; z++) {
				out.writeln(data[z]);
			}
			out.close();

			// run normalizer
			System.out.println("Running normalizer");
			String normalizerOutputDir = outdirectory + q + "-Normalizer/";
			Double correlationthreshold = 0.9;

			runNormalizationMethod(expression,
					normalizerOutputDir,
					correlationthreshold,
					cellcountfile,
					mdscomponents,
					cellcountfile,
					gte,
					threads);

			// run the interaction model
			String outputCellCountFile = normalizerOutputDir + "CellTypeProxyFile.txt";
			String covariatelist = null;
			runInteractionModel(inexpPCCorrected, outputCellCountFile, ingt, gte, snpprobecombos, threads, normalizerOutputDir, covariatelist);

			q++;
		}

		if (cellcountfile != null) {

			// tfout.writeln("Correlation between actual cell counts and PC1 scores: " + r + "\tr2: " + (r * r) + "\tn: " + xArr.length
			TextFile listOfCorrelations = new TextFile(outdirectory + "correlationList.txt", TextFile.W);
			listOfCorrelations.writeln("#probes\tspearman\tpearson");
			for (q = 1; q < data.length + 1; q++) {
				String outdir = outdirectory + q + "-Normalizer/";
				String correlationFile = outdir + "ComparisonToCellCount.txt";
				tf = new TextFile(correlationFile, TextFile.R);
				tf.readLine();
				String ln = tf.readLine();
				String[] elems = ln.split("\t");
				String pearson = elems[0];
				String spearman = elems[1];
				tf.close();
				listOfCorrelations.writeln(q + "\t" + spearman + "\t" + pearson);
			}
			listOfCorrelations.close();
		}

	}

	public void runNormalizationMethod(String inexpraw,
									   String outdirectory,
									   Double correlationThreshold,
									   String celltypeSpecificProbeFile,
									   String mdsComponentFile,
									   String cellCountFile,
									   String gte,
									   Integer threads) throws IOException {
		InteractionAnalysisMultiThreaded iamt = new InteractionAnalysisMultiThreaded();
		iamt.prepareDataForCelltypeSpecificEQTLMapping(inexpraw, outdirectory, correlationThreshold, celltypeSpecificProbeFile, mdsComponentFile, cellCountFile, gte, threads);
	}

	private void runInteractionModel(String inExpPCCorrected, String covariateFile, String ingt,
									 String gte, String snpprobecombinationfile, Integer nrThreads, String out,
									 String covariateList) throws Exception {
		InteractionAnalysisMultiThreaded iamt = new InteractionAnalysisMultiThreaded();

        /*
		 String inExpPCCorrected, String covariateFile, String ingt,
         String gte, String snpprobecombinationfile, Integer nrThreads, String out,
         String covariateList
         */
		iamt.runInteractionAnalysis(inExpPCCorrected, covariateFile, ingt, gte, snpprobecombinationfile, nrThreads, out, covariateList);
	}

	private double[] ztransform(double[] valsx) {
		double mean = JSci.maths.ArrayMath.mean(valsx);
		double sd = JSci.maths.ArrayMath.standardDeviation(valsx);
		for (int i = 0; i < valsx.length; i++) {
			valsx[i] = (valsx[i] - mean) / sd;
		}
		return valsx;
	}

	public void prepareDataForCelltypeSpecificEQTLMapping(DoubleMatrixDataset<String, String> rawExpressionDataset, String inexpraw, String outdirectory, Double correlationThreshold, String celltypeSpecificProbeFile, String mdsComponentFile, String cellCountFile, String gte, Integer threads) throws IOException {
		String rawExpressionDataFile = inexpraw;
		// 7. select Cell type specific probes
		System.out.println("Loading list of cell type specific probes from: " + celltypeSpecificProbeFile);
		HashSet<String> cellTypeSpecificProbeSet = new HashSet<String>();
		TextFile cellSpecificProbeTF = new TextFile(celltypeSpecificProbeFile, TextFile.R);
		cellTypeSpecificProbeSet.addAll(cellSpecificProbeTF.readAsArrayList());
		cellSpecificProbeTF.close();

		if (cellTypeSpecificProbeSet.isEmpty()) {
			System.err.println("Error: " + celltypeSpecificProbeFile + " is empty!");
			System.exit(-1);
		} else {
			System.out.println(cellTypeSpecificProbeSet.size() + " cell type specific probes loaded.");
		}

		// 1. load gene expression data
		System.out.println("Loading gene expression data.");

		double[][] rawExpressionData = rawExpressionDataset.getRawData();

		// determine the number of cell type specific probes in this dataset
		int probeCounter = 0;
		List<String> probes = rawExpressionDataset.rowObjects;
		for (int i = 0; i < probes.size(); i++) {
			if (cellTypeSpecificProbeSet.contains(probes.get(i))) {
				probeCounter++;
			}
		}

		if (probeCounter == 0) {
			System.err.println("Error: none of the cell type specific probes defined in " + celltypeSpecificProbeFile + " are present in expression dataset: " + rawExpressionDataset.fileName);
			System.exit(-1);
		} else {
			System.out.println(probeCounter + " of the cell type specific probes are in your dataset.");
		}

		System.out.println("Now reloading the gene expression data for the samples that passed the QC.");
		// 6. Remove samples with r < 0.9 for PC1
		// reload expression file, include only samples that pass QC...
//        rawExpressionDataset = new DoubleMatrixDataset<String, String>(rawExpressionDataFile);
//        rawExpressionData = rawExpressionDataset.getRawData();

//        // quantile normalize, log2 transform again, because the number of samples might have been changed..
//        QuantileNormalization.quantilenormalize(rawExpressionData);
//        Log2Transform.log2transform(rawExpressionData);
		rawExpressionData = rawExpressionDataset.rawData;

		// collect data for cell type specific probes
		double[][] probeData = new double[probeCounter][rawExpressionDataset.colObjects.size()];
		probeCounter = 0;
		ArrayList<String> cellTypeSpecificProbeDatasetRowNames = new ArrayList<String>();
		for (int i = 0; i < probes.size(); i++) {
			if (cellTypeSpecificProbeSet.contains(probes.get(i))) {
				probeData[probeCounter] = rawExpressionData[i];
				cellTypeSpecificProbeDatasetRowNames.add(probes.get(i));
				probeCounter++;
			}
		}

		// initiate cell type specific probe correlation matrix
		double[][] celltypeSpecificCorrelationMatrix = new double[probeCounter][probeCounter];
		for (int i = 0; i < probeCounter; i++) {
			for (int j = i + 1; j < probeCounter; j++) {
				double r = Correlation.correlate(probeData[i], probeData[j]);
				celltypeSpecificCorrelationMatrix[i][j] = r;
				celltypeSpecificCorrelationMatrix[j][i] = r;
			}
			celltypeSpecificCorrelationMatrix[i][i] = 1;
		}

		// save the correlation matrix
		DoubleMatrixDataset<String, String> probeCorrelationMatrixOut = new DoubleMatrixDataset<String, String>();
		probeCorrelationMatrixOut.colObjects = cellTypeSpecificProbeDatasetRowNames;
		probeCorrelationMatrixOut.rowObjects = cellTypeSpecificProbeDatasetRowNames;
		probeCorrelationMatrixOut.rawData = celltypeSpecificCorrelationMatrix;
		probeCorrelationMatrixOut.recalculateHashMaps();
//        probeCorrelationMatrixOut.save(outdirectory + "CelltypeSpecificProbeCorrelationMatrix.txt.gz");

		// 9. PCA over cell specific probe correlation matrix
		DoubleMatrixDataset<String, String> cellTypeSpecificDataset = new DoubleMatrixDataset<String, String>(probeData);
		cellTypeSpecificDataset.colObjects = rawExpressionDataset.colObjects;
		cellTypeSpecificDataset.rowObjects = cellTypeSpecificProbeDatasetRowNames;
//        cellTypeSpecificDataset.save(expressionOutputDirectory + "CellTypeSpecificProbeExpression.txt.gz");
		cellTypeSpecificDataset.transposeDataset();
		Normalizer n = new Normalizer();
		// calculate first Principal Component over the cell type specific probe matrix...
		Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = n.calculatePCA(cellTypeSpecificDataset, celltypeSpecificCorrelationMatrix, outdirectory + "CellTypeSpecificProbePCA", 1);

		// 10. PC1 scores: cell specific proxy -- write to file for future use...
		DoubleMatrixDataset<String, String> cellSpecificPCScores = PCAResults.getLeft();

		//Ensure that the cellTypeSpecificPCScores correlate positively with the set of probes that we have used to determine this component:
		double[] pcScoresSamples = new double[cellSpecificPCScores.nrRows];
		for (int i = 0; i < cellSpecificPCScores.nrRows; i++) {
			pcScoresSamples[i] = cellSpecificPCScores.rawData[i][0];
		}
		cellTypeSpecificDataset.transposeDataset();
		int nrProbesCorrelatingPositively = 0;
		for (int i = 0; i < cellTypeSpecificDataset.rawData.length; i++) {
			double corr = JSci.maths.ArrayMath.correlation(pcScoresSamples, cellTypeSpecificDataset.rawData[i]);
			if (corr >= 0) {
				nrProbesCorrelatingPositively++;
			} else {
				nrProbesCorrelatingPositively--;
			}
		}
		if (nrProbesCorrelatingPositively < 0) {
			for (int i = 0; i < cellSpecificPCScores.nrRows; i++) {
				cellSpecificPCScores.rawData[i][0] = -cellSpecificPCScores.rawData[i][0];
			}
		}

		TextFile tfOutCellSpecific = new TextFile(outdirectory + "CellTypeProxyFile.txt", TextFile.W);
		tfOutCellSpecific.writeln("Sample\tCellCountProxyValue");
		for (int i = 0; i < cellSpecificPCScores.nrRows; i++) {
			tfOutCellSpecific.writeln(cellSpecificPCScores.rowObjects.get(i) + "\t" + cellSpecificPCScores.rawData[i][0]);
		}
		tfOutCellSpecific.close();

		File f = new File(outdirectory + "CellTypeSpecificProbePCA.PCAOverSamplesEigenvalues.txt.gz");
		f.delete();
		f = new File(outdirectory + "CellTypeSpecificProbePCA.PCAOverSamplesEigenvectors.txt.gz");
		f.delete();
		f = new File(outdirectory + "CellTypeSpecificProbePCA.PCAOverSamplesEigenvectorsTransposed.txt.gz");
		f.delete();
		f = new File(outdirectory + "CellTypeSpecificProbePCA.PCAOverSamplesPrincipalComponents.txt.gz");
		f.delete();

	}

	private void combine(String file1, String file2, String pheno, String outfile) throws IOException {
		DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>(file1);
		DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(file2);

		ds2.transposeDataset();

		HashMap<Integer, Integer> sharedSamples = new HashMap<Integer, Integer>();
		for (int i = 0; i < ds1.rowObjects.size(); i++) {
			Integer indexI = ds2.hashRows.get(ds1.rowObjects.get(i));
			if (indexI != null) {
				sharedSamples.put(i, indexI);
			}
		}
		System.out.println(sharedSamples.size() + " shared samples");

		double[][] data2 = new double[sharedSamples.size()][ds1.colObjects.size() + 1];
		ArrayList<String> newRows = new ArrayList<String>();
		int ctr = 0;
		Integer phenoIndex = ds2.hashCols.get(pheno);
		System.out.println("Merging col: " + phenoIndex);
		for (int i = 0; i < ds1.rowObjects.size(); i++) {
			Integer otherSample = sharedSamples.get(i);
			if (otherSample != null) {
				System.arraycopy(ds1.rawData[i], 0, data2[ctr], 0, ds1.rawData[i].length);

				data2[ctr][data2[ctr].length - 1] = ds2.rawData[otherSample][phenoIndex];
				newRows.add(ds1.rowObjects.get(i));
				ctr++;
			}
		}

		DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<>();
		dsout.rawData = data2;
		dsout.rowObjects = newRows;
		ds1.colObjects.add(pheno);
		dsout.colObjects = ds1.colObjects;
		dsout.recalculateHashMaps();
		dsout.save(outfile);

	}

	private void plotBoxPlots(String file1, String pheno1, String out) throws IOException {
		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(file1);

		Integer phenoId = ds.hashCols.get(pheno1);
		int nrOther = 0;
		for (int i = 0; i < ds.nrRows; i++) {
			if (ds.rawData[i][phenoId] > 1) {
				nrOther++;
			}
		}

		for (int col = 0; col < ds.nrCols; col++) {
			double[] x = new double[ds.nrRows - nrOther];
			double[] y = new double[nrOther];
			double[][][] container = new double[1][2][];
			String[] datasetNames = new String[]{ds.colObjects.get(col)};
			String[][] xlabels = new String[1][];
			xlabels[0] = new String[]{"1", "2"};
			String ylabel = ds.colObjects.get(col);

			container[0][0] = x;
			container[0][1] = y;

			int q1 = 0;
			int q2 = 0;
			for (int row = 0; row < ds.nrRows; row++) {
				if (ds.rawData[row][phenoId] > 1) {
					y[q2] = ds.rawData[row][col];
					q2++;
				} else {
					x[q1] = ds.rawData[row][col];
					q1++;
				}
			}
//            System.out.println((ds.nrRows - nrOther) + "\t" + q1 + "\t" + nrOther + "\t" + q2);
			WilcoxonMannWhitney wmw = new WilcoxonMannWhitney();
			double pw = wmw.returnWilcoxonMannWhitneyPValue(x, y);
			double pt = TTest.test(x, y);
			System.out.println(ds.colObjects.get(col) + "\tw: " + pw + "\tt: " + pt);

		}
	}

	private void annotate(String in, String annotation, String out) throws IOException {
		TextFile tf = new TextFile(annotation, TextFile.R);
		HashMap<String, String> probeToGene = new HashMap<String, String>();
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String probe = elems[1];
			String gene = elems[2];
			probeToGene.put(probe, gene);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(in, TextFile.R);
		TextFile outf = new TextFile(out, TextFile.W);
		String header = tf2.readLine() + "\tGene Symbol";
		outf.writeln(header);
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String probe = elems[0];

			String gene = probeToGene.get(probe);

			outf.writeln(Strings.concat(elems, Strings.tab) + "\t" + gene);
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		outf.close();

	}

	private void determineCorrelationFromZScore(String fileIn, String fileOut) throws IOException {
		TextFile tf = new TextFile(fileIn, TextFile.R);
		String[] header = tf.readLineElems(TextFile.tab);
		int sampleCol = 0;
		int zscorecol = 0;
		int fdrcol = 0;
		int probeCol = 4;
// EGCUT,SHIP_TREND,Groningen-HT12,-,Rotterdam,DILGOM,INCHIANTI,HVH-HT12v3,-
		HashMap<String, Integer> sampleSizes = new HashMap<String, Integer>();
        /*
         The order of datasets is as follows: EGCUT, SHIP-TREND, Fehrmann-HT12v3, Fehrmann-H8v2, Rotterdam Study, DILGOM, InChianti, HVH-HT12v3, HVH-HT12v4
         Their respective sample sizes are: 891, 963, 1240, 229, 762, 509, 611, 43, 63
         */
		sampleSizes.put("EGCUT", 891);
		sampleSizes.put("SHIP_TREND", 963);
		sampleSizes.put("Groningen-HT12", 1240);
		sampleSizes.put("Groningen-H8v2", 229);
		sampleSizes.put("Rotterdam", 762);
		sampleSizes.put("DILGOM", 509);
		sampleSizes.put("INCHIANTI", 611);
		sampleSizes.put("HVH-HT12v3", 43);
		sampleSizes.put("HVH-HT12v4", 63);

		for (int col = 0; col < header.length; col++) {
			String s = header[col];
			if (s.equals("DatasetsWhereSNPProbePairIsAvailableAndPassesQC")) {
				sampleCol = col;
			}
			if (s.equals("OverallZScore")) {
				zscorecol = col;
			}

			if (s.equals("FDR")) {
				fdrcol = col;
			}

			if (s.equals("ProbeName")) {
				probeCol = col;
			}
		}

		System.out.println("Sample " + sampleCol);
		System.out.println("FDR " + fdrcol);
		System.out.println("probe " + probeCol);
		System.out.println("Z " + zscorecol);

		String[] elems = tf.readLineElems(TextFile.tab);
		TextFile out = new TextFile(fileOut, TextFile.W);
		out.writeln("Pval\tSNP\tProbe\tZ\tSample\tFDR\tR");
		while (elems != null) {

			String snp = elems[1];
			String pval = elems[0];
			String samples = elems[sampleCol];
			double z = Double.parseDouble(elems[zscorecol]);
			double fdr = Double.parseDouble(elems[fdrcol]);
			String probe = elems[probeCol];

			samples = samples.replaceAll("\"", "");
			String[] sampleStr = samples.split(",");
			int sum = 0;
			for (String s : sampleStr) {
				if (!s.equals("-")) {
					Integer sampleSize = sampleSizes.get(s);
					if (sampleSize == null) {
						System.err.println("Could not find: " + s);

						sampleSize = 0;
					}
					sum += sampleSize;
				}
			}

			Double R = ZScores.zScoreToCorrelation(z, sum);
			out.writeln(pval + "\t" + snp + "\t" + probe + "\t" + z + "\t" + sum + "\t" + fdr + "\t" + R);

			elems = tf.readLineElems(TextFile.tab);
		}
		out.close();

		tf.close();

	}

	private void intersectEQTLFiles(String file1, String correlationFile, String fileout) throws IOException {
		HashMap<String, Double> eQTLToR = new HashMap<String, Double>();
		TextFile tf = new TextFile(correlationFile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[1];
			String probe = elems[2];
			Double r = Double.parseDouble(elems[6]);
			eQTLToR.put(snp + "-" + probe, r);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(file1, TextFile.R);

		TextFile out = new TextFile(fileout, TextFile.W);
		out.writeln(tf2.readLine() + "\tR");
		String[] elems2 = tf2.readLineElems(TextFile.tab);

		while (elems2 != null) {
			String snp = elems2[1];
			String probe = elems2[4];

			double r = eQTLToR.get(snp + "-" + probe);
			String output = Strings.concat(elems2, Strings.tab) + "\t" + r;
			out.writeln(output);
			elems2 = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		out.close();
	}

	public void merge(String normal1, String robust1, String normal2, String robust2, String out) throws IOException {

		HashMap<String, Double> eQTLToSE = loadSE(normal1);
		HashMap<String, Double> eQTLToSERobust = loadSE(robust1);

		HashMap<String, Double> eQTLToSE2 = loadSE(normal2);
		HashMap<String, Double> eQTLToSERobust2 = loadSE(robust2);

		Set<String> eqtls = eQTLToSE.keySet();
		TextFile tf = new TextFile(out, TextFile.W);
		tf.writeln("SNP-Probe\tSE1\tSE1Robust\tSE2\tSE2Robust");
		for (String s : eqtls) {
			if (eQTLToSE2.containsKey(s)) {

				String ln = s + "\t" + eQTLToSE.get(s) + "\t" + eQTLToSERobust.get(s)
						+ "\t" + eQTLToSE2.get(s) + "\t" + eQTLToSERobust2.get(s);
				tf.writeln(ln);
			}
		}
		tf.close();

	}

	public HashMap<String, Double> loadSE(String file) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, Double> output = new HashMap<String, Double>();
		while (elems != null) {
			String snp = elems[0];
			String probe = elems[1];
			String eQTL = snp + "-" + probe;
			Double d = Double.parseDouble(elems[elems.length - 2]);
			output.put(eQTL, d);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}



	public void spearmanSE(String normal1, String robust1) throws IOException {
		HashMap<String, Double> eQTLToSE = loadSE(normal1);
		HashMap<String, Double> eQTLToSERobust = loadSE(robust1);

		Set<String> eqtls = eQTLToSE.keySet();
		ArrayList<Double> x = new ArrayList<Double>();
		ArrayList<Double> y = new ArrayList<Double>();
		for (String s : eqtls) {

			x.add(eQTLToSE.get(s));
			y.add(eQTLToSERobust.get(s));

		}

		SpearmansCorrelation corr = new SpearmansCorrelation();
		double c = corr.correlation(Primitives.toPrimitiveArr(x.toArray(new Double[0])), Primitives.toPrimitiveArr(y.toArray(new Double[0])));
		System.out.println(c);
	}
}
