package nl.harmjanwestra.playground.methylation;

import com.google.common.util.concurrent.AtomicDouble;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessReader;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class QuantileNormalizeDiskBased450K {

//	private static final RankingAlgorithm COV_RANKER_TIE = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);

	public void run(String in, String out) throws IOException {


		// note this expects probes on columns

		int maxNonNAvalues = Integer.MIN_VALUE;


		DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(in);
		System.out.println("Input: " + in);
		System.out.println("Output: " + out);
		System.out.println("Input is size: " + it.getNrRows() + "\t" + it.getNrCols());
		int probeCount = it.getNrCols();
		int sampleCount = it.getNrRows();

		// expect genes/probes on columns
		double[] rankedMeanAtomic = new double[probeCount];
		int[] rankedNrOfInds = new int[probeCount];


		// samples on rows
		ProgressBar pb = new ProgressBar(it.getNrRows(), "Determining means.");

		int ctr = 0;
		for (double[] sample : it) {
			java.util.Arrays.parallelSort(sample);
			IntStream.range(0, probeCount).parallel().forEach(probeID -> {
				if (!Double.isNaN(sample[probeID])) {
					rankedNrOfInds[probeID]++;
					rankedMeanAtomic[probeID] += sample[probeID];
				}
			});
			ctr++;
			pb.set(ctr);

		}
		pb.close();

		it.close();
		double[] rankedMean = new double[probeCount];
		for (int probeID = 0; probeID < probeCount; probeID++) {
			if (rankedNrOfInds[probeID] > 0) {
				rankedMean[probeID] = rankedMeanAtomic[probeID] / rankedNrOfInds[probeID];
			}
		}

		double[] rankedMeanClasses = new double[probeCount - 1];
		for (int probeID = 0; probeID < (probeCount - 1); probeID++) {
			rankedMeanClasses[probeID] = ((rankedMean[probeID] + rankedMean[probeID + 1]) / 2);
		}

		//Iterate through each sample:
		it = new DoubleMatrixDatasetRowIterable(in);
		DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(new ArrayList<String>(it.getCols()), out);
		ArrayList<String> samplenames = new ArrayList<String>(it.getRows());
		org.apache.commons.math3.stat.correlation.SpearmansCorrelation spearman = new org.apache.commons.math3.stat.correlation.SpearmansCorrelation();
		int rowctr = 0;

//		DoubleMatrixDatasetRandomAccessReader reader = new DoubleMatrixDatasetRandomAccessReader(in);
		ProgressBar pb2 = new ProgressBar(it.getNrRows(), "Quantile normalizing...");
		for (double[] sample : it) {
			RankingAlgorithm COV_RANKER_TIE = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);
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
//			try {
//				synchronized (writer) {
			writer.append(probesQuantileNormalized, samplename);
//				}
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
			pb.set(rowctr);
			rowctr++;
			//System.out.println("Normalized sample:\t" + rowctr + "/" + sampleCount + " " + samplename + "\tPearson correlation original data and ranked data:\t" + JSci.maths.ArrayMath.correlation(probes, probesRanked) + "\tSpearman correlation original data and quantile normalized data:\t" + spearman.correlation(probes, probesQuantileNormalized));

		}
		pb2.close();
//		for (double[] sample : it) {
//
//			rowctr++;
//		}

		it.close();
		writer.close();

	}
}
