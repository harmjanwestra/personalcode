package hms.hwestra.interactionrebuttal2;

import hms.hwestra.interactionrebuttal.InteractionRebuttal;
import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.*;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author hwestra
 */
public class InteractionRebuttal2 {

	public static void main(String[] args) {
		InteractionRebuttal2 b = new InteractionRebuttal2();
		try {
//			String axisdir = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/0.1.PreinigerEtAlAxesOfVariation/";
//			String axisdir = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/0.1.PreinigerEtAlAxesOfVariationFORCENORMAL/";
//			String annot = "/Data/ProbeAnnotation/2015-03-23-HT12v3ILMNIdToArrayAddress.txt";
////			b.rewritePleinigerAxisToArrayAddress(axisdir, annot);
//			b.createProxysFromPleinigerAxis(axisdir);

//			String pcFile = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/3.3.AdditionalPCsInModelTop58Probes/EGCUTData/CellTypeSpecificProbePCA.PCAOverSamplesPrincipalComponents.txt";
//			String cellcountFile = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/EGCUTData/EGCUTValidSamplesNeutrosOnly.txt";
//			String pheno = "Neutrophils";
			//b.iterativelyIncreaseNumberOfPCsInCellCountPredictionModel(pcFile, cellcountFile, pheno);
//			InteractionRebuttal a = new InteractionRebuttal();
////
//			String normal1 = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/3.4.2.RobustSEPValues/EGCUT/Normal/InteractionResults.txt";
//			String robust1 = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/3.4.2.RobustSEPValues/EGCUT/Robust/InteractionResults.txt";
//			String out = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/3.4.2.RobustSEPValues/EGCUT.txt";
//			b.determineInteractionPvalueAndMerge(normal1, robust1, out);
//
//
//			String normal2 = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/3.4.2.RobustSEPValues/GRNG/Normal/InteractionResults.txt";
//			String robust2 = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/3.4.2.RobustSEPValues/GRNG/Robust/InteractionResults.txt";
//			out = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/3.4.2.RobustSEPValues/GRNG.txt";
//			b.determineInteractionPvalueAndMerge(normal2, robust2, out);


//			String fileIn = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/GeneExpressionData/InChianti/GSE48152_RAW.txt";
//			String fileOut = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/GeneExpressionData/InChianti/InChiantiData.txt";
//
//			b.rewriteExpressionMatrix(fileIn, fileOut, annot);
//			String[] filesIn = new String[]{"/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/GeneExpressionData/Kora/E-MTAB-1708.raw.1/",
//					"/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/GeneExpressionData/Kora/E-MTAB-1708.raw.2/",
//					"/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/GeneExpressionData/Kora/E-MTAB-1708.raw.3/"};
//			fileOut = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/GeneExpressionData/Kora/KoraData.txt";
//
//			b.rewriteMTabToMatrix(filesIn, fileOut, annot);

			String[] files = new String[10];
			files[0] = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/3.4.2.RobustSEPValues/Meta-forceNormal/output.txt";
			String[] fileNames = new String[10];
			fileNames[0] = "NeutrophilProxy";
			for (int i = 1; i < 10; i++) {
				files[i] = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/0.1.PreinigerEtAlAxesOfVariationFORCENORMAL/meta/Axis" + i + ".txt";
				fileNames[i] = "Axis-" + i;
			}

			String out = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/0.1.PreinigerEtAlAxesOfVariationFORCENORMAL/meta/merged-noforcednormalneutro.txt";

			String annot = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/ilmnToGene.txt";

			double threshold = 0.05 / 13037;
			System.out.println("P < " + threshold);
			b.mergeMetaFiles(files, fileNames, out, annot, threshold);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void mergeMetaFiles(String[] files, String[] fileNames, String out, String annot, double threshold) throws IOException {

		TextFile tf1 = new TextFile(annot, TextFile.R);
		Map<String, String> ilmnToArr = tf1.readAsHashMap(0, 1);
		tf1.close();

		ArrayList<HashMap<String, Double>> eqtls = new ArrayList<HashMap<String, Double>>();
		HashSet<String> uniqueEQTLs = new HashSet<String>();
		for (String file : files) {
			HashMap<String, Double> eqtl = loadInteractionMetaTableZScoreBlaat(file);
			Set<String> keyset = eqtl.keySet();
			for (String key : keyset) {
				uniqueEQTLs.add(key);
			}
			System.out.println(eqtl.size() + " loaded from " + file);
			eqtls.add(eqtl);
		}

		System.out.println(uniqueEQTLs.size() + " unique eQTLs");
		TextFile outfile = new TextFile(out, TextFile.W);
		String header = "SNP\tProbe\tGene";
		for (int i = 0; i < files.length; i++) {
			header += "\tZ-" + fileNames[i] + "\tP-" + fileNames[i];
		}
		outfile.writeln(header);

		int[] nrEQTLsBelowThreshold = new int[files.length];
		int[] nrTotalEQTLs = new int[files.length];
		for (String eqtl : uniqueEQTLs) {
			String[] eqtlelems = eqtl.split("-");
			String outln = eqtlelems[0] + "\t" + eqtlelems[1] + "\t" + ilmnToArr.get(eqtlelems[1]);
			for (int i = 0; i < files.length; i++) {
				HashMap<String, Double> map = eqtls.get(i);
				Double z = map.get(eqtl);

				if (z != null && !Double.isNaN(z)) {

					double p = ZScores.zToP(z);
					if (p < threshold) {
						nrEQTLsBelowThreshold[i]++;
					}
					nrTotalEQTLs[i]++;
					outln += "\t" + z + "\t" + p;
				} else {
					outln += "\tNaN\tNaN";
				}
			}
			outfile.writeln(outln);
		}

		outfile.close();

		System.out.println("pvals below threshold:");
		for (int i = 0; i < files.length; i++) {
			System.out.println(fileNames[i] + "\t" + nrEQTLsBelowThreshold[i] + "\t" + nrTotalEQTLs[i]);
		}

	}

	private HashMap<String, Double> loadInteractionMetaTableZScoreBlaat(String file) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, Double> output = new HashMap<String, Double>();
		boolean headerparsed = false;
		int col = -1;
		while (elems != null) {
			if (!elems[0].startsWith("#")) {
				if (!headerparsed) {
					for (int i = 0; i < elems.length; i++) {
						if (elems[i].equals("Meta_flipped_interaction_Z-score")) {
							col = i;
						}
					}
					headerparsed = true;
				} else {
					if (col == -1) {
						System.err.println("Error: Meta_flipped_interaction_Z-score not found");
						System.exit(-1);
					}
					String snp = elems[0];
					String probe = elems[1];
					String eQTL = snp + "-" + probe;
					Double d = Double.parseDouble(elems[col]);
					output.put(eQTL, d);
				}
			}

			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}

	private void rewriteMTabToMatrix(String[] filesIn, String fileOut, String annot) throws IOException {
		ArrayList<String> allFiles = new ArrayList<>();
		for (int i = 0; i < filesIn.length; i++) {
			String[] listOfFiles = Gpio.getListOfFiles(filesIn[i]);
			for (String file : listOfFiles) {
				allFiles.add(filesIn[i] + file);
			}

		}

		TextFile tf1 = new TextFile(annot, TextFile.R);
		Map<String, String> ilmnToArr = tf1.readAsHashMap(0, 1);
		tf1.close();

		String file = allFiles.get(0);
		TextFile tf = new TextFile(file, TextFile.R);
		String header = tf.readLine();
		String[] lnElems = tf.readLineElems(TextFile.tab);
		ArrayList<String> probes = new ArrayList<String>();
		HashMap<String, Integer> probeToInt = new HashMap<String, Integer>();
		int probectr = 0;
		while (lnElems != null) {
			String probe = lnElems[0];
			probes.add(ilmnToArr.get(probe));
			probeToInt.put(probe, probectr);
			probectr++;
			lnElems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		int nrSamples = allFiles.size();
		double[][] data = new double[probes.size()][nrSamples];
		ArrayList<String> samples = new ArrayList<String>();
		for (int f = 0; f < allFiles.size(); f++) {
			file = allFiles.get(f);
			tf = new TextFile(file, TextFile.R);
			String[] headerElems = tf.readLineElems(TextFile.tab);

			samples.add(headerElems[1].replaceAll("_AVG_Signal", ""));
			String[] elems = tf.readLineElems(TextFile.tab);

			while (elems != null) {
				double d = Double.parseDouble(elems[1]);
				String probe = elems[0];
				Integer probeId = probeToInt.get(probe);
				data[probeId][f] = d;
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

		}

		DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>();
		dsout.rawData = data;
		dsout.rowObjects = probes;
		dsout.colObjects = samples;
		dsout.recalculateHashMaps();

		dsout.save(fileOut);

	}

	private void rewriteExpressionMatrix(String fileIn, String fileOut, String annot) throws IOException {
		TextFile tf1 = new TextFile(annot, TextFile.R);
		Map<String, String> ilmnToArr = tf1.readAsHashMap(0, 1);
		tf1.close();

		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(fileIn);
		boolean[] keepColumn = new boolean[ds.nrCols];
		int nrToKeep = 0;
		for (int col = 0; col < ds.nrCols; col++) {
			if (!ds.colObjects.get(col).equals("Detection Pval")) {
				keepColumn[col] = true;
				nrToKeep++;
			}
		}

		ArrayList<String> newRows = new ArrayList<String>();
		for (String s : ds.rowObjects) {
			newRows.add(ilmnToArr.get(s));
		}

		double[][] newData = new double[ds.nrRows][nrToKeep];
		int colctr = 0;
		ArrayList<String> newCols = new ArrayList<>();
		for (int col = 0; col < ds.nrCols; col++) {
			if (keepColumn[col]) {
				for (int row = 0; row < ds.nrRows; row++) {
					newData[row][colctr] = ds.rawData[row][col];
				}
				colctr++;
				newCols.add(ds.colObjects.get(col));
			}
		}

		DoubleMatrixDataset<String, String> out = new DoubleMatrixDataset<String, String>();
		out.rawData = newData;
		out.colObjects = newCols;
		out.rowObjects = newRows;
		out.recalculateHashMaps();
		out.save(fileOut);


	}

	private void determineInteractionPvalueAndMerge(String normal1, String robust1, String out) throws IOException {

		InteractionRebuttal a = new InteractionRebuttal();
		HashMap<String, Double> normalSE = a.loadSE(normal1);
		HashMap<String, Double> robustSE = a.loadSE(robust1);

		HashMap<String, Double> normalZ = loadZ(normal1);
		HashMap<String, Double> robustZ = loadZ(robust1);

		Set<String> eqtls = normalZ.keySet();
		TextFile tf = new TextFile(out, TextFile.W);
		SpearmansCorrelation corr = new SpearmansCorrelation();
		tf.writeln("SNP-Probe\tZNorm\tZRobust\tPNorm\tPRobust\tSENorm\tSERobust");
		ArrayList<Double> x = new ArrayList<Double>();
		ArrayList<Double> y = new ArrayList<Double>();
		ArrayList<Double> x2 = new ArrayList<Double>();
		ArrayList<Double> y2 = new ArrayList<Double>();

		ArrayList<Double> x3 = new ArrayList<Double>();
		ArrayList<Double> y3 = new ArrayList<Double>();
		for (String s : eqtls) {
			if (robustZ.containsKey(s)) {
				if (!Double.isNaN(normalZ.get(s))) {
					x.add(normalZ.get(s));
					y.add(robustZ.get(s));
					x2.add(-Math.log10(ZScores.zToP(normalZ.get(s))));
					y2.add(-Math.log10(ZScores.zToP(robustZ.get(s))));
					x3.add(normalSE.get(s));
					y3.add(robustSE.get(s));
				}
				String ln = s + "\t" + normalZ.get(s) + "\t" + robustZ.get(s)
						+ "\t" + -Math.log10(ZScores.zToP(normalZ.get(s))) + "\t" + -Math.log10(ZScores.zToP(robustZ.get(s)))
						+ "\t" + normalSE.get(s) + "\t" + robustSE.get(s);
				tf.writeln(ln);
			}
		}
		tf.close();
		System.out.println(x.size());
		double c = corr.correlation(Primitives.toPrimitiveArr(x.toArray(new Double[0])), Primitives.toPrimitiveArr(y.toArray(new Double[0])));
		double c2 = corr.correlation(Primitives.toPrimitiveArr(x2.toArray(new Double[0])), Primitives.toPrimitiveArr(y2.toArray(new Double[0])));

		double c3 = Correlation.correlate(Primitives.toPrimitiveArr(x.toArray(new Double[0])), Primitives.toPrimitiveArr(y.toArray(new Double[0])));
		double c4 = Correlation.correlate(Primitives.toPrimitiveArr(x2.toArray(new Double[0])), Primitives.toPrimitiveArr(y2.toArray(new Double[0])));


		double c5 = corr.correlation(Primitives.toPrimitiveArr(x3.toArray(new Double[0])), Primitives.toPrimitiveArr(y3.toArray(new Double[0])));
		double c6 = Correlation.correlate(Primitives.toPrimitiveArr(x3.toArray(new Double[0])), Primitives.toPrimitiveArr(y3.toArray(new Double[0])));

		System.out.println("Z Spearman: " + c + " Pearson: " + c3);
		System.out.println("P Spearman: " + c2 + " Pearson: " + c4);
		System.out.println("SE Spearman: " + c5 + " Pearson: " + c6);
	}

	public HashMap<String, Double> loadZ(String file) throws IOException {
		TextFile tf = new TextFile(file, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, Double> output = new HashMap<String, Double>();
		while (elems != null) {
			String snp = elems[0];
			String probe = elems[1];
			String eQTL = snp + "-" + probe;
			Double d = Double.parseDouble(elems[7]);
			output.put(eQTL, d);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		return output;
	}

	private void iterativelyIncreaseNumberOfPCsInCellCountPredictionModel(String pcFile, String cellcountFile, String pheno) throws IOException {

		DoubleMatrixDataset<String, String> pcs = new DoubleMatrixDataset<String, String>(pcFile); // samples on rows, pcs on cols?
		DoubleMatrixDataset<String, String> cellcounts = new DoubleMatrixDataset<String, String>(cellcountFile); // samples on rows, celltype on cols


		Integer phenoId = cellcounts.hashCols.get(pheno);

		boolean[] includeRow = new boolean[pcs.nrRows];
		int shared = 0;
		for (int i = 0; i < pcs.nrRows; i++) {
			String sample = pcs.rowObjects.get(i);
			if (cellcounts.hashRows.containsKey(sample)) {
				shared++;
				includeRow[i] = true;
			}
		}


		// order the samples of the cell count in the order of the pcs
		double[] olsY = new double[shared]; //Ordinary least squares: cell count
		int ctr = 0;
		for (int i = 0; i < pcs.nrRows; i++) {
			String sample = pcs.rowObjects.get(i);
			Integer sampleId = cellcounts.hashRows.get(sample);
			if (sampleId != null) {
				olsY[ctr] = cellcounts.rawData[sampleId][phenoId];
				ctr++;
			}
		}

		org.apache.commons.math3.distribution.FDistribution fDist = null;
		cern.jet.random.tdouble.engine.DoubleRandomEngine randomEngine = null;
		cern.jet.random.tdouble.StudentT tDistColt = null;

		OLSMultipleLinearRegression previousFullModel = null;

		for (int col = 0; col < pcs.nrCols; col++) {
			OLSMultipleLinearRegression regressionFullModel = new OLSMultipleLinearRegression();
			OLSMultipleLinearRegression regressionOrigModel = new OLSMultipleLinearRegression();

			int nrPcs = col + 1;
			double[][] olsX = new double[shared][nrPcs];
			double[][] olsXN = new double[shared][1];
			for (int inc = 0; inc < col + 1; inc++) {
				ctr = 0;
				for (int i = 0; i < pcs.nrRows; i++) {
					if (includeRow[i]) {
						olsX[ctr][inc] = pcs.rawData[i][inc];
						ctr++;
					}
				}
			}

			double[] pc = new double[shared];
			ctr = 0;
			for (int i = 0; i < pcs.nrRows; i++) {
				if (includeRow[i]) {
					pc[ctr] = pcs.rawData[i][col];
					olsXN[ctr][0] = pcs.rawData[i][0];
					ctr++;
				}
			}


			double corr = JSci.maths.ArrayMath.correlation(pc, olsY);
			Correlation.correlationToZScore(olsY.length);
			double z = Correlation.convertCorrelationToZScore(olsY.length, corr);
			double p = ZScores.zToP(z);

			regressionFullModel.newSampleData(olsY, olsX);
			regressionOrigModel.newSampleData(olsY, olsXN);


			double rsquaredadj = regressionFullModel.calculateAdjustedRSquared();
			double rsquared = regressionFullModel.calculateRSquared();

			double rse = regressionOrigModel.estimateRegressionStandardError();
			double rsefull = regressionFullModel.estimateRegressionStandardError();


			double rss1 = regressionOrigModel.calculateResidualSumOfSquares();
			double rss2 = regressionFullModel.calculateResidualSumOfSquares();
			double F = ((rss1 - rss2) / (3 - 2)) / (rss2 / (olsY.length - 3));
			int numParams1 = 1; // regressor + intercept
			int numParams2 = nrPcs; // regressors + intercept
			if (nrPcs > 1) {


				double F2 = ((rss1 - rss2) / (numParams2 - numParams1)) / (rss2 / (olsY.length - numParams2));

				double rss3 = previousFullModel.calculateResidualSumOfSquares();
				int numParams3 = nrPcs - 1;
				double FPrevious = ((rss3 - rss2) / (numParams2 - numParams3)) / (rss2 / (olsY.length - numParams2));


				// pf(f, m1$df.residual-m2$df.residual, m2$df.residual, lower.tail = FALSE)
				// (double numeratorDegreesOfFreedom, double denominatorDegreesOfFreedom)
				fDist = new org.apache.commons.math3.distribution.FDistribution((numParams2 - numParams1), olsY.length - numParams2);
				FDistribution fDistPrev = new FDistribution((numParams2 - numParams3), olsY.length - numParams2);


				double anovaFTestP = -1;
				double anovaFTestP2 = -1;
				try {
					anovaFTestP = 1 - fDist.cumulativeProbability(F2);
					anovaFTestP2 = 1 - fDist.cumulativeProbability(FPrevious);
					if (anovaFTestP < 1E-160) {
						anovaFTestP = 1E-16;
					}

					if (anovaFTestP2 < 1E-160) {
						anovaFTestP2 = 1E-16;
					}
				} catch (Exception err) {
				}

				System.out.println(nrPcs + "\t" + corr + "\t" + z + "\t" + p + "\t" + rsquared + "\t" + numParams2 + "\t" + F2 + "\t" + FPrevious + "\t" + anovaFTestP + "\t" + anovaFTestP2);
			} else {
				System.out.println(nrPcs + "\t" + corr + "\t" + z + "\t" + p + "\t" + rsquared + "\t" + numParams1);
			}

			previousFullModel = regressionFullModel;

		}


		ArrayList<String> colNames = new ArrayList<String>();
		colNames.add("CellCount");
		double[][] data = new double[shared][pcs.nrCols + 1];
		for (int i = 0; i < olsY.length; i++) {
			data[i][0] = olsY[i];
		}


		ArrayList<String> rowNames = new ArrayList<String>();
		for (int col = 0; col < pcs.nrCols; col++) {
			ctr = 0;
			colNames.add(pcs.colObjects.get(col));
			for (int row = 0; row < pcs.nrRows; row++) {
				if (includeRow[row]) {
					data[ctr][col + 1] = pcs.rawData[row][col];
					ctr++;
				}

			}
		}

		for (int row = 0; row < pcs.nrRows; row++) {
			if (includeRow[row]) {
				rowNames.add("Sample_" + pcs.rowObjects.get(row));
			}
		}


		DoubleMatrixDataset<String, String> dsout = new DoubleMatrixDataset<String, String>();
		dsout.rawData = data;
		dsout.rowObjects = rowNames;
		dsout.colObjects = colNames;
		dsout.recalculateHashMaps();
		dsout.save(pcFile + "-mergedWCellCount.txt");


	}

	public DoubleMatrixDataset<String, String> loadDataset(String d, String gte) throws IOException {

		if (gte == null) {
			return new DoubleMatrixDataset<>(d);
		} else {
			TextFile tf = new TextFile(gte, TextFile.R);
			Set<String> set = tf.readAsSet(0, TextFile.tab);
			tf.close();

			return new DoubleMatrixDataset<String, String>(d, null, set);
		}
	}

	private void createProxysFromPleinigerAxis(String axisdir) throws IOException {
		// create proxy phenotypes
//		String egcut = "/Sync/AeroFS/cellTypeeQTL/2014-10-28-Rebuttal/EGCUTData/EGCUT-RawDataQNLog2.txt";
		String egcut = "/Data/tmp/EGCUTForceNormal/EGCUT-RawDataQNLog2.ForcedNormal.txt.gz";

		InteractionRebuttal b1 = new InteractionRebuttal();
		String gte = "/Sync/AeroFS/cellTypeeQTL/2015-01-31-Rebuttal2/EGCUTGTE.txt";
		DoubleMatrixDataset<String, String> ds = loadDataset(egcut, gte);

		String proxyOut = axisdir + "/EGCUT/";
		Gpio.createDir(proxyOut);
		for (int i = 1; i < 10; i++) {
			String inexpraw = egcut;
			String outdir = proxyOut + "Axis" + i + "/";
			Gpio.createDir(outdir);
			Double correlationthreshold = 0.8;
			String probefile = axisdir + "Axis" + i + "-arr.txt";
			b1.prepareDataForCelltypeSpecificEQTLMapping(ds, inexpraw, outdir, correlationthreshold, probefile, null, null, gte, 1);
		}

//		String grng = "/Volumes/Data/Datasets/GeneticalGenomicsDatasets/BloodHT12/2015-03-26-InteractionModel/ExpressionData/ExpressionDataRaw-QNormLog2Transformed.txt.gz";
		String grng = "/Data/tmp/GRNGForceNormal/ExpressionDataRaw-QNormLog2Transformed.ForcedNormal.txt.gz";
		gte = null;
		ds = loadDataset(grng, gte);
		proxyOut = axisdir + "/GRNG/";
		Gpio.createDir(proxyOut);
		for (int i = 1; i < 10; i++) {
			String inexpraw = grng;
			String outdir = proxyOut + "Axis" + i + "/";
			Gpio.createDir(outdir);
			Double correlationthreshold = 0.8;
			String probefile = axisdir + "Axis" + i + "-arr.txt";
			b1.prepareDataForCelltypeSpecificEQTLMapping(ds, inexpraw, outdir, correlationthreshold, probefile, null, null, gte, 1);
		}
	}

	public String[] readAsArray(String f, int col) throws IOException {
		TextFile tf = new TextFile(f, TextFile.R);
		ArrayList<String> s = new ArrayList<String>();
		String[] lnElems = tf.readLineElems(TextFile.tab);
		while (lnElems != null) {
			if (lnElems.length > col) {
				if (lnElems[col].trim().length() > 0) {
					s.add(lnElems[col]);
				}
			}
			lnElems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		return s.toArray(new String[0]);
	}

	public void rewritePleinigerAxisToArrayAddress(String dir, String annot) throws IOException {

		TextFile tf1 = new TextFile(annot, TextFile.R);
		Map<String, String> ilmnToArr = tf1.readAsHashMap(0, 1);
		tf1.close();
		tf1.open();
		Map<String, String> ilmnToG = tf1.readAsHashMap(0, 2);
		tf1.close();

		for (int i = 1; i < 10; i++) {

			String[] list = readAsArray(dir + "Axis" + i + ".txt", 0);
			String[] genes = readAsArray(dir + "Axis" + i + ".txt", 1);

			int notfound = 0;
			int equalG = 0;
			TextFile out = new TextFile(dir + "Axis" + i + "-arr.txt", TextFile.W);
			for (int q = 0; q < list.length; q++) {
				String s = list[q];
				String g = genes[q];
				String arr = ilmnToArr.get(s);
				String gen = ilmnToG.get(s);
				if (g.toLowerCase().equals(gen.toLowerCase())) {
					equalG++;
				} else {
					System.out.println(g + "\t" + gen + "\t" + s);

				}
				if (arr == null) {
					notfound++;
				}
				out.writeln(arr);
			}
			out.close();
			System.out.println(i + "\t" + notfound + "\t" + list.length + "\t" + equalG);

		}

	}
}
