/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.zscoretables;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class zScoreTableCorrelator {

    int distributionsize = 1000000;
    GWASCatalog catalog = new GWASCatalog();
    HashMap<String, String> snpToChr = new HashMap<String, String>();
    HashMap<String, String> snpToPos = new HashMap<String, String>();
    long[] r2PermutedFrequencyDistribution = new long[distributionsize];
    long[] r2PermutedCumulativeFrequencyDistribution = new long[distributionsize];
    long grandTotal = 0;
    double fdrcutoff = 0.05;
    String probeListFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt";
    String ensemblAnnotationFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-OnlyEnsemblAnnotation.txt";
    HashSet<String> allowedTraits = null;
    boolean nonparametric = true;

    public void run(String in, String pw, String out, String gwasCatalog, String dbSNP, int perm) throws IOException {

	catalog.read(gwasCatalog);
//
//	// "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-27-SNPMappings-dbSNP130.txt.gz"
	loadDBSNP(dbSNP);
//	// run correlation on permuted data
//
//
	ZScoreTableManipulation manipulator = new ZScoreTableManipulation();
//
	if (perm > 0) {


	    for (int p = 1; p < perm + 1; p++) {
		System.out.println("Running permutation " + p);
		String infile = in + "-Permutation" + p + ".txt.gz";
		String filename = infile;
		String finalFileName = filename + "-filtered.txt" + "-ens.txt" + "-collapsed.txt";
		if (!Gpio.exists(finalFileName)) {
		    if (!Gpio.exists(filename + "-filtered.txt")) {
			manipulator.filterZScoreTableForSetOfProbes(probeListFile, infile);
		    }

		    filename += "-filtered.txt";
		    if (!Gpio.exists(filename + "-ens.txt")) {
			manipulator.convertProbeNumbersToEnsemblId(ensemblAnnotationFile, filename);
		    }
		    filename += "-ens.txt";

		    if (!Gpio.exists(filename + "-collapsed.txt")) {
			manipulator.collapseDuplicateProbes(filename);
		    }
		}



		correlate(finalFileName, pw, out, true, null);

	    }

	    for (int i = 0; i < r2PermutedCumulativeFrequencyDistribution.length; i++) {
		if (i == 0) {
		    r2PermutedCumulativeFrequencyDistribution[i] = r2PermutedFrequencyDistribution[i];
		} else {
		    r2PermutedCumulativeFrequencyDistribution[i] = r2PermutedCumulativeFrequencyDistribution[i - 1] + r2PermutedFrequencyDistribution[i];

		}
		grandTotal += r2PermutedFrequencyDistribution[i];
	    }

	    // determine null distribution
//
	    TextFile permutedNullDist = new TextFile(out + "PermutedNullDistribution.txt.gz", TextFile.W);

	    for (int i = 0; i < r2PermutedCumulativeFrequencyDistribution.length; i++) {
		permutedNullDist.writeln(i + "\t" + (i * 1d / r2PermutedCumulativeFrequencyDistribution.length) + "\t" + r2PermutedFrequencyDistribution[i] + "\t" + r2PermutedCumulativeFrequencyDistribution[i] + "\t" + (r2PermutedFrequencyDistribution[i] / perm) + "\t" + (r2PermutedCumulativeFrequencyDistribution[i] / perm));
	    }
	    permutedNullDist.close();

	    System.out.println("Grand total: " + grandTotal);
	}

	// calculate effects for real data given the FDR threshold.


	String infile = in + ".txt.gz";


	String filename = infile;

	String finalFileName = filename + "-filtered.txt" + "-ens.txt" + "-collapsed.txt";
	if (!Gpio.exists(finalFileName)) {
	    if (!Gpio.exists(filename + "-filtered.txt")) {
		manipulator.filterZScoreTableForSetOfProbes(probeListFile, infile);
	    }

	    filename += "-filtered.txt";
	    if (!Gpio.exists(filename + "-ens.txt")) {
		manipulator.convertProbeNumbersToEnsemblId(ensemblAnnotationFile, filename);
	    }
	    filename += "-ens.txt";

	    if (!Gpio.exists(filename + "-collapsed.txt")) {
		manipulator.collapseDuplicateProbes(filename);
	    }
	    filename += "-collapsed.txt";
	}


//	correlate(finalFileName, pw, out, false, null);

//	double cutoff = determineFDRThreshold(out + "PermutedNullDistribution.txt.gz", out + "-corr-r2-dist.txt");



//	String filename = in + ".txt.gz-filtered.txt-ens.txt-collapsed.txt";
//	String traitselectionfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog-uniquetraits.txt";
//	loadTraitSelection(traitselectionfile);
//	correlate(finalFileName, pw, out, false, cutoff); // 0.0029969999999999997
//	correlate(finalFileName, pw, out, false, 0.0029969999999999997); // 0.0029969999999999997


    }

    private void loadDBSNP(String dbSNP) throws IOException {

	ArrayList<GWASSNP> gwassnps = new ArrayList<GWASSNP>();
	gwassnps.addAll(catalog.getSnps());

	HashSet<String> gwasSNPstr = new HashSet<String>();
	for (int i = 0; i < gwassnps.size(); i++) {
	    gwasSNPstr.add(gwassnps.get(i).getName());
	}

	TextFile dbsnptf = new TextFile(dbSNP, TextFile.R);
	String[] dbsnpelems = dbsnptf.readLineElems(TextFile.tab);

	while (dbsnpelems != null) {

	    String snpstr = dbsnpelems[2];
	    if (gwasSNPstr.contains(snpstr)) {
		snpToChr.put(snpstr, dbsnpelems[0]);
		snpToPos.put(snpstr, dbsnpelems[1]);
	    }
	    dbsnpelems = dbsnptf.readLineElems(TextFile.tab);
	}

	dbsnptf.close();
    }

    private void correlate(String inputData, String pathwayData, String outfile, boolean permutedData, Double r2Cutoff) throws IOException {

	DoubleMatrixDataset<String, String> eQTLZScoreData = new DoubleMatrixDataset<String, String>(inputData);
	DoubleMatrixDataset<String, String> pathwayZScoreData = new DoubleMatrixDataset<String, String>(pathwayData);

	ArrayList<String> eQTLZScoreColumns = (ArrayList<String>) eQTLZScoreData.colObjects;
	ArrayList<String> eQTLZScoreRows = (ArrayList<String>) eQTLZScoreData.rowObjects;

	ArrayList<String> pathwayZScoreColumns = (ArrayList<String>) pathwayZScoreData.colObjects;
	ArrayList<String> pathwayZScoreRows = (ArrayList<String>) pathwayZScoreData.rowObjects;

	int[] colToCol = new int[eQTLZScoreColumns.size()];
	int sharedgenes = 0;
	for (int i = 0; i < colToCol.length; i++) {
	    String transGene = eQTLZScoreColumns.get(i);
	    colToCol[i] = -1;
	    for (int j = 0; j < pathwayZScoreColumns.size(); j++) {
		if (pathwayZScoreColumns.get(j).equals(transGene)) {
		    colToCol[i] = j;
		    sharedgenes++;
		    break;
		}
	    }
	}
	System.out.println("Nr of shared genes:\t" + sharedgenes);

	int[] correlationdistribution = new int[10000];
	int[] correlationdistributionr2 = new int[r2PermutedFrequencyDistribution.length];

	TextFile out = null;
	if (!permutedData && r2Cutoff != null) {
	    String textoutfile = outfile + "PathWayAnnotation-FDR0.05-Cutoff" + r2Cutoff + ".txt";
	    out = new TextFile(textoutfile, TextFile.W);
	}
	ProgressBar pb = new ProgressBar(eQTLZScoreData.nrRows);
//	for (int i = 0; i < eQTLZScoreData.nrRows; i++) {


	int nrprocs = Runtime.getRuntime().availableProcessors();
	ExecutorService threadPool = Executors.newFixedThreadPool(nrprocs);

	CompletionService<Pair<Integer, Double>> pool = new ExecutorCompletionService<Pair<Integer, Double>>(threadPool);



//for(int i = 0; i &lt; 10; i++){
//
//   pool.submit(new StringTask());
//
//}
//
//
//
//for(int i = 0; i &lt; 10; i++){
//
//   String result = pool.take().get();
//
//
//
//   //Compute the result
//
//}




	for (int i = 0; i < eQTLZScoreData.nrRows; i++) {
	    String snp = eQTLZScoreRows.get(i);
	    String traits = null;

	    GWASSNP gwasnp = catalog.getSnpToObj().get(snp);

	    boolean testSNP = true;


	    if (gwasnp == null) {
		System.out.println(snp + " not present in catalog!");
	    } else {
		HashSet<GWASTrait> traithash = gwasnp.getAssociatedTraits();
		if (traithash == null) {
		    System.out.println(snp + "\t has no associated traits??");
		} else {
		    traits = traithash.size() + "";

		    GWASTrait[] traitra = new GWASTrait[traithash.size()];
		    traithash.toArray(traitra);

		    int traitcounter = 0;
		    for (int t = 0; t < traitra.length; t++) {
			if (allowedTraits != null) {
			    if (allowedTraits.contains(traitra[t].getName())) {
				if (gwasnp.getPValueAssociatedWithTrait(traitra[t]) != null && gwasnp.getPValueAssociatedWithTrait(traitra[t]) < 5E-8) {
				    traitcounter++;
				}
			    }
			}
			traits += "," + traitra[t].getName();
		    }

		    if (allowedTraits != null && traitcounter == 0) {
			testSNP = false;
		    }
		}
	    }

	    int nrMissing = 0;
	    for (int p = 0; p < eQTLZScoreData.nrCols; p++) {
		if (Double.isNaN(eQTLZScoreData.rawData[i][p])) {
		    nrMissing++;
		}
	    }

	    if (nrMissing == 0 && testSNP) {

		// distribute to threads
		for (int j = 0; j < pathwayZScoreData.nrRows; j++) {
		    ZScoreTableCorrelationTask c = new ZScoreTableCorrelationTask(pathwayZScoreData, sharedgenes, eQTLZScoreData, colToCol, i, j, r2PermutedCumulativeFrequencyDistribution.length, nonparametric);
		    pool.submit(c);
		}

		// return results :)
		int returnedResults = 0;
		Double[] results = new Double[pathwayZScoreData.nrRows];
		while (returnedResults < pathwayZScoreData.nrRows) {
		    try {

			Pair<Integer, Double> result = pool.take().get();

			results[result.getLeft()] = result.getRight();
			returnedResults++;

		    } catch (InterruptedException ex) {
			Logger.getLogger(zScoreTableCorrelator.class.getName()).log(Level.SEVERE, null, ex);
		    } catch (ExecutionException ex) {
			Logger.getLogger(zScoreTableCorrelator.class.getName()).log(Level.SEVERE, null, ex);
		    }


		}

		for (int q = 0; q < pathwayZScoreData.nrRows; q++) {

		    Double r2 = results[q];
		    if (r2 != null) {
			int actualRow = q;
			int corr2bin = (int) Math.floor((r2 / (1d / distributionsize)));
			if (corr2bin < 0) {
			    corr2bin = 0;
			}

			if (corr2bin >= distributionsize) {
			    corr2bin = distributionsize;
			}

			if (permutedData) {
			    r2PermutedFrequencyDistribution[corr2bin]++;
			}
			correlationdistributionr2[corr2bin]++;

			if (!permutedData && r2Cutoff != null && r2 >= r2Cutoff) {

			    String chr = snpToChr.get(eQTLZScoreRows.get(i)); // + chr +"\t" + chrpos + "\t"
			    String chrpos = snpToPos.get(eQTLZScoreRows.get(i));
			    System.out.println(i + "\t" + traits + "\t" + eQTLZScoreRows.get(i) + "\t" + chr + "\t" + chrpos + "\t" + pathwayZScoreRows.get(actualRow) + "\t" + Math.sqrt(r2) + "\t" + r2);
			    out.writeln(i + "\t" + traits + "\t" + eQTLZScoreRows.get(i) + "\t" + chr + "\t" + chrpos + "\t" + pathwayZScoreRows.get(actualRow) + "\t" + Math.sqrt(r2) + "\t" + r2);

			}
		    }

		}


	    }
	    pb.set(i);
	}
	threadPool.shutdown();

	pb.close();
	if (!permutedData && r2Cutoff != null) {
	    out.close();
	}

//	TextFile outcorrr = new TextFile(outfile + "-corrdist.txt", TextFile.W);
//
//	int sum = 0;
//	for (int i = 0; i < correlationdistribution.length; i++) {
//	    double b = 2d / correlationdistribution.length * i;
//	    sum += correlationdistribution[i];
//	    outcorrr.writeln(i + "\t" + b + "\t" + correlationdistribution[i] + "\t" + sum);
//	}
//
//	outcorrr.close();

	if (!permutedData) {
	    TextFile outcorrr2 = new TextFile(outfile + "-corr-r2-dist.txt", TextFile.W);

	    int sum = 0;
	    for (int i = 0; i < correlationdistributionr2.length; i++) {
		double b = (double) i * 1d / correlationdistributionr2.length;
		sum += correlationdistributionr2[i];
		outcorrr2.writeln(i + "\t" + b + "\t" + correlationdistributionr2[i] + "\t" + sum);

	    }
	    outcorrr2.close();
	}

    }

    //		for (int j = 0; j < pathwayZScoreData.nrRows; j++) {
    //		    double[] valsX = new double[sharedgenes];
    //		    double[] valsY = new double[sharedgenes];
    //		    int itr = 0;
    //		    for (int p = 0; p < eQTLZScoreData.nrCols; p++) {
    //			if (colToCol[p] != -1) {
    //			    valsX[itr] = Math.abs(eQTLZScoreData.rawData[i][p]);
    //			    valsY[itr] = Math.abs(pathwayZScoreData.rawData[j][colToCol[p]]);
    //			    itr++;
    //			}
    //		    }
    //
    //		    double correlation = JSci.maths.ArrayMath.correlation(valsX, valsY);
    //
    //		    int correlationbin = (int) Math.floor((1 + correlation) / 2d * correlationdistribution.length);
    //		    if (correlationbin >= correlationdistribution.length) {
    //			correlationbin = correlationdistribution.length - 1;
    //		    }
    //
    //		    if (correlationbin < 0) {
    //			correlation = 0;
    //		    }
    //
    //		    correlationdistribution[correlationbin]++;
    //
    //		    double r2 = correlation * correlation;
    //
    //		    int corr2bin = (int) Math.floor((r2 / (1d / correlationdistributionr2.length)));
    //		    if (corr2bin < 0) {
    //			corr2bin = 0;
    //		    }
    //
    //		    if (corr2bin >= correlationdistributionr2.length) {
    //			corr2bin = correlationdistributionr2.length;
    //		    }
    //
    ////		    if (r2 > 0.0001) {
    ////			System.out.println(r2 + "\t" + corr2bin + "\t" + (r2 / (1d / correlationdistributionr2.length)));
    ////		    }
    //		    if (permutedData) {
    //			r2PermutedFrequencyDistribution[corr2bin]++;
    //		    }
    //
    //		    correlationdistributionr2[corr2bin]++;
    //		    if (!permutedData) {
    //
    //
    //			double fdr = r2;
    //
    //
    //
    //			if (fdr > cutoff) {
    //			    String chr = snpToChr.get(eQTLZScoreRows.get(i)); // + chr +"\t" + chrpos + "\t"
    //			    String chrpos = snpToPos.get(eQTLZScoreRows.get(i));
    //			    System.out.println(i + "\t" + traits + "\t" + eQTLZScoreRows.get(i) + "\t" + chr + "\t" + chrpos + "\t" + pathwayZScoreRows.get(j) + "\t" + correlation + "\t" + r2);
    //			    out.writeln(i + "\t" + traits + "\t" + eQTLZScoreRows.get(i) + "\t" + chr + "\t" + chrpos + "\t" + pathwayZScoreRows.get(j) + "\t" + correlation + "\t" + r2);
    //			}
    //		    }
    //
    //		}
//		for (int j = 0; j < pathwayZScoreData.nrRows; j++) {
//		    double[] valsX = new double[sharedgenes];
//		    double[] valsY = new double[sharedgenes];
//		    int itr = 0;
//		    for (int p = 0; p < eQTLZScoreData.nrCols; p++) {
//			if (colToCol[p] != -1) {
//			    valsX[itr] = Math.abs(eQTLZScoreData.rawData[i][p]);
//			    valsY[itr] = Math.abs(pathwayZScoreData.rawData[j][colToCol[p]]);
//			    itr++;
//			}
//		    }
//
//		    double correlation = JSci.maths.ArrayMath.correlation(valsX, valsY);
//
//		    int correlationbin = (int) Math.floor((1 + correlation) / 2d * correlationdistribution.length);
//		    if (correlationbin >= correlationdistribution.length) {
//			correlationbin = correlationdistribution.length - 1;
//		    }
//
//		    if (correlationbin < 0) {
//			correlation = 0;
//		    }
//
//		    correlationdistribution[correlationbin]++;
//
//		    double r2 = correlation * correlation;
//
//		    int corr2bin = (int) Math.floor((r2 / (1d / correlationdistributionr2.length)));
//		    if (corr2bin < 0) {
//			corr2bin = 0;
//		    }
//
//		    if (corr2bin >= correlationdistributionr2.length) {
//			corr2bin = correlationdistributionr2.length;
//		    }
//
////		    if (r2 > 0.0001) {
////			System.out.println(r2 + "\t" + corr2bin + "\t" + (r2 / (1d / correlationdistributionr2.length)));
////		    }
//		    if (permutedData) {
//			r2PermutedFrequencyDistribution[corr2bin]++;
//		    }
//
//		    correlationdistributionr2[corr2bin]++;
//		    if (!permutedData) {
//
//
//			double fdr = r2;
//
//
//
//			if (fdr > cutoff) {
//			    String chr = snpToChr.get(eQTLZScoreRows.get(i)); // + chr +"\t" + chrpos + "\t"
//			    String chrpos = snpToPos.get(eQTLZScoreRows.get(i));
//			    System.out.println(i + "\t" + traits + "\t" + eQTLZScoreRows.get(i) + "\t" + chr + "\t" + chrpos + "\t" + pathwayZScoreRows.get(j) + "\t" + correlation + "\t" + r2);
//			    out.writeln(i + "\t" + traits + "\t" + eQTLZScoreRows.get(i) + "\t" + chr + "\t" + chrpos + "\t" + pathwayZScoreRows.get(j) + "\t" + correlation + "\t" + r2);
//			}
//		    }
//
//		}
    private void createtraitlist(String gwascatalog, String snpids) throws IOException {

	GWASCatalog c = new GWASCatalog();
	c.read(gwascatalog);

	TextFile dbsnptf = new TextFile(snpids, TextFile.R);

	String[] dbsnpelems = dbsnptf.readLineElems(TextFile.tab);

	while (dbsnpelems != null) {

	    String snp = dbsnpelems[2];
	    GWASSNP gwasnp = c.getSnpToObj().get(snp);
	    if (gwasnp != null) {
		String traits = "";

		HashSet<GWASTrait> traithash = gwasnp.getAssociatedTraits();
		if (traithash == null) {
		    System.out.println(snp + "\t has no associated traits??");
		} else {
		    traits = traithash.size() + "";

		    GWASTrait[] traitra = new GWASTrait[traithash.size()];
		    traithash.toArray(traitra);
		    for (int t = 0; t < traitra.length; t++) {
			traits += "," + traitra[t].getName();
		    }
		}

		System.out.println(snp + "\t" + dbsnpelems[0] + "\t" + dbsnpelems[1] + "\t" + traits);
	    }
	    dbsnpelems = dbsnptf.readLineElems(TextFile.tab);
	}

	dbsnptf.close();
    }

    private double determineFDRThreshold(String nullDistFile, String realDataFile) throws IOException {


	int[] permutedDist = new int[distributionsize];
	int[] permutedDistCumulative = new int[distributionsize];

	int[] realDist = new int[distributionsize];
	int[] realDistCumulative = new int[distributionsize];


	TextFile nulldist = new TextFile(nullDistFile, TextFile.R);

	String[] elems = nulldist.readLineElems(TextFile.tab);
	int i = 0;
	while (elems != null) {

	    permutedDist[i] = Integer.parseInt(elems[4]);
//	    permutedDistCumulative[i] = Integer.parseInt(elems[5]);

	    elems = nulldist.readLineElems(TextFile.tab);
	    i++;
	}

	nulldist.close();
	TextFile realdist = new TextFile(realDataFile, TextFile.R);

	elems = realdist.readLineElems(TextFile.tab);
	i = 0;
	while (elems != null) {

	    realDist[i] = Integer.parseInt(elems[2]);
//	    realDistCumulative[i] = Integer.parseInt(elems[3]);


	    elems = realdist.readLineElems(TextFile.tab);
	    i++;
	}

	realdist.close();


	for (i = distributionsize - 1; i > -1; i--) {
	    if (i == distributionsize - 1) {
		permutedDistCumulative[i] = permutedDist[i];
		realDistCumulative[i] = realDist[i];
	    } else {
		permutedDistCumulative[i] = permutedDistCumulative[i + 1] + permutedDist[i];
		realDistCumulative[i] = realDistCumulative[i + 1] + realDist[i];
	    }
	}


	// now for each bin, determine the FDR
	double r2 = 0;
	for (i = distributionsize - 1; i > -1; i--) {
	    int fp = permutedDistCumulative[i];
	    int tp = realDistCumulative[i];
	    double fdr = (double) fp / (tp + fp);
	    if (fdr > fdrcutoff) {
		break;
	    }
	    r2 = i * (1d / distributionsize);
//	    System.out.println(i+"\t"+r2+"\t"+tp+"\t"+fp+"\t"+fdr);
	}


	System.out.println("Detected " + r2 + " FDR threshold for cutoff " + fdrcutoff);
	return r2;


    }

    private void loadTraitSelection(String traitselectionfile) throws IOException {
	System.out.println("Loading trait selection from: " + traitselectionfile);
	TextFile tf = new TextFile(traitselectionfile, TextFile.R);
	allowedTraits = new HashSet<String>();
	String[] elems = tf.readLineElems(TextFile.tab);

	while (elems != null) {
	    if (elems.length > 1) {
		if (elems[1].equals("FALSE")) {
		    allowedTraits.add(elems[0]);
		}
	    }
	    elems = tf.readLineElems(TextFile.tab);
	}
	tf.close();

	System.out.println(allowedTraits.size() + " traits selected");
    }
}
