package nl.harmjanwestra.playground.methylation;

import com.google.common.util.concurrent.AtomicDouble;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class QuantileNormalizeDiskBased450K {

	private static final RankingAlgorithm COV_RANKER_TIE = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);

	public void run(String in, String out) throws IOException {


		// note this expects probes on columns

		int maxNonNAvalues = Integer.MIN_VALUE;


		DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(in);
		// expect genes/probes on columns
		AtomicDouble[] rankedMeanAtomic = new AtomicDouble[it.getNrCols()];
		AtomicInteger[] rankedNrOfInds = new AtomicInteger[it.getNrCols()];
		for (int d = 0; d < rankedMeanAtomic.length; d++) {
			rankedMeanAtomic[d] = new AtomicDouble(0d);
		}

		int probeCount = it.getNrCols();
		int sampleCount = it.getNrRows();

		// samples on rows
		ProgressBar pb = new ProgressBar(it.getNrRows(),"Determining means.");
		for (double[] sample : it) {

			double[] x = new double[probeCount];
			for (int probeID = 0; probeID < probeCount; probeID++) {
				x[probeID] = sample[probeID]; // data.getElementQuick(probeID, sampleID); //rawData[probeID][sampleID];
			}
			java.util.Arrays.parallelSort(x);
			for (int probeID = 0; probeID < probeCount; probeID++) {

				if (!Double.isNaN(x[probeID])) {
					rankedNrOfInds[probeID].getAndIncrement();
					rankedMeanAtomic[probeID].getAndAdd(x[probeID]);
				}
			}
			pb.iterate();
		}
		pb.close();

		it.close();
		double[] rankedMean = new double[probeCount];
		for (int probeID = 0; probeID < probeCount; probeID++) {
			if (rankedNrOfInds[probeID].get() > 0) {
				rankedMean[probeID] = rankedMeanAtomic[probeID].get() / rankedNrOfInds[probeID].get();
			}
		}

		double[] rankedMeanClasses = new double[probeCount - 1];
		for (int probeID = 0; probeID < (probeCount - 1); probeID++) {
			rankedMeanClasses[probeID] = ((rankedMean[probeID] + rankedMean[probeID + 1]) / 2);
		}

		//Iterate through each sample:
		it = new DoubleMatrixDatasetRowIterable(in);
		DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(new ArrayList<String>(it.getCols()), out);
		ArrayList<String> samplenames = new ArrayList<String>(it.getNrRows());
		org.apache.commons.math3.stat.correlation.SpearmansCorrelation spearman = new org.apache.commons.math3.stat.correlation.SpearmansCorrelation();
		int rowctr = 0;
		for (double[] sample : it) {
			double[] probes = new double[probeCount];
			for (int p = 0; p < probeCount; p++) {
				probes[p] = sample[p];
			}
			double[] probesRanked = COV_RANKER_TIE.rank(probes);
			double[] probesQuantileNormalized = new double[probeCount];

			for (int p = 0; p < probeCount; p++) {

				if ((probesRanked[p] % 1) != 0) {
					probesQuantileNormalized[p] = rankedMeanClasses[(int) Math.floor((probesRanked[p] - 1))];
				} else {
					probesQuantileNormalized[p] = rankedMean[(int) (probesRanked[p] - 1)];
				}
			}
			String samplename = samplenames.get(rowctr);
			writer.append(probesQuantileNormalized, samplename);
			System.out.println("Normalized sample:\t" + samplename + "\tPearson correlation original data and ranked data:\t" + JSci.maths.ArrayMath.correlation(probes, probesRanked) + "\tSpearman correlation original data and quantile normalized data:\t" + spearman.correlation(probes, probesQuantileNormalized));
			rowctr++;
		}

		it.close();

	}
}
