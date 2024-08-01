package nl.harmjanwestra.methylation;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.Marker;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessReader;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.util.RunTimer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.stream.IntStream;

public class CorrelateDiskBased450K {

	int clBlockSize = 128;


	boolean ocl = false;

	public static void main(String[] args) {


		CorrelateDiskBased450K c = new CorrelateDiskBased450K();

		try {
			c.runblocks("E:\\Methylation\\T1MAndT2M-mval-qqnorm",
					"E:\\Methylation\\T1MAndT2M-mval-qqnorm-correl",
					1000);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String in, String out) throws IOException {
		// correlate rows



		DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(in);
		ArrayList<String> rowids = new ArrayList<String>(it.getRows());
		DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(rowids, out);
		int row = 0;
		for (double[] sample : it) {
			try {
				DoubleMatrixDatasetRandomAccessReader reader2 = new DoubleMatrixDatasetRandomAccessReader(in);
				double[] correlations = new double[it.getNrRows()];
				IntStream.range(row, rowids.size()).parallel().forEachOrdered(sample2id -> {
					try {
						RankingAlgorithm COV_RANKER_TIE = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);
						SpearmansCorrelation correlation = new SpearmansCorrelation(COV_RANKER_TIE);
						double[] sample2 = null;
						synchronized (reader2) {
							sample2 = reader2.getRow(sample2id);
						}
						double cor = correlation.correlation(sample, sample2);
						correlations[sample2id] = cor;
					} catch (IOException e) {
						e.printStackTrace();
					}
				});

				writer.append(correlations, rowids.get(row));

			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		it.close();
		writer.close();

	}


	public void runblocks(String in, String out, int rowsPerBlock) throws IOException {
		// correlate rows

		DoubleMatrixDatasetRandomAccessReader reader = new DoubleMatrixDatasetRandomAccessReader(in);
		ArrayList<String> rowids = new ArrayList<String>(reader.getRowObjects());
		DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(rowids, out);

		int nrRows = reader.getRowObjects().size();
		int remainingblocks = nrRows % rowsPerBlock;
		System.out.println(remainingblocks + " blocks will remain.");
		int rowsPerBlockA = rowsPerBlock;
		int rowsPerBlockB = rowsPerBlock;
		int row = 0;
		RunTimer timer = new RunTimer();
		timer.start();
		double[][] blockA = new double[rowsPerBlock][];
		double[] meanA = new double[rowsPerBlock];
		double[] varA = new double[rowsPerBlock];

		double[][] blockB = new double[rowsPerBlock][];
		double[] meanB = new double[rowsPerBlock];
		double[] varB = new double[rowsPerBlock];

		double[][] correlations = new double[rowsPerBlockA][nrRows];


		while (row < nrRows) {

			if (row + rowsPerBlockA > nrRows) {
				rowsPerBlockA = remainingblocks;
				blockA = new double[rowsPerBlockA][];
				meanA = new double[rowsPerBlockA];
				varA = new double[rowsPerBlockA];
				correlations = new double[rowsPerBlockA][nrRows];
			}

			// read first block
			readBlock(row, rowsPerBlockA, reader, blockA);
			rank(blockA);
			calcmean(blockA, meanA, varA);
			System.gc();


			// correlate within block
			correlate(blockA, meanA, varA, 0, blockA, meanA, varA, 0, correlations, true);

			// read future blocks until end
			int rowB = row + rowsPerBlockB;

			RunTimer timer2 = new RunTimer();

			while (rowB < nrRows) {
				if (rowB + rowsPerBlockB >= nrRows) {
					rowsPerBlockB = remainingblocks;
					blockB = new double[rowsPerBlockB][];
					meanB = new double[rowsPerBlockB];
					varB = new double[rowsPerBlockB];
				}
				blockB = readBlock(rowB, rowsPerBlockB, reader, blockB);
				rank(blockB);
				calcmean(blockB, meanB, varB);
				System.gc();
//				correlate(blockA, meanA, varA, 0, blockB, meanB, varB, rowB, correlations, false);
				correlate(blockA, meanA, varA, 0, blockB, meanB, varB, rowB, correlations, false);
				rowB += rowsPerBlockB;

				int nrRowsB = nrRows - row;
				long diff = timer2.getTimeDiff() / 1000000000L;
				double timePerIter = (double) diff / (double) rowB;
				double timeLeft = timePerIter * (double) (nrRowsB - rowB);
				String strTimeLeft = timer2.getTimeDesc((long) timeLeft * 1000000000L);

				System.out.println(row + " / " + nrRows + " | " + rowB + " / " + nrRows + " (blockB) rows done. " + timer2.getTimeDesc() + " time spent. ETA: " + strTimeLeft);
			}

			// save block to the disk
			ProgressBar pb = new ProgressBar(correlations.length, "Writing correlations to disk");
			for (int i = 0; i < correlations.length; i++) {
				String rowId = rowids.get(row + i);
				writer.append(correlations[i], rowId);
				pb.iterate();
			}
			pb.close();

			row += rowsPerBlockA;
			rowsPerBlockB = rowsPerBlock;

			// reset block B data
			varB = new double[rowsPerBlockB];
			meanB = new double[rowsPerBlockB];
			blockB = new double[rowsPerBlockB][];

			long diff = timer.getTimeDiff() / 1000000000L;
			double timePerIter = (double) diff / (double) row;
			double timeLeft = timePerIter * (double) (nrRows - row);
			String strTimeLeft = timer.getTimeDesc((long) timeLeft * 1000000000L);

			System.out.println(row + " / " + nrRows + " rows done. " + timer.getTimeDesc() + " time spent. ETA: " + strTimeLeft);
		}

		System.out.println("Total time: " + timer.getTimeDesc());
		writer.close();

	}

	private void calcmean(double[][] blockB, double[] meanB, double[] varB) {
		ProgressBar pb = new ProgressBar(blockB.length, "Calculating means");
		IntStream.range(0, blockB.length).parallel().forEach(i -> {
			meanB[i] = Descriptives.mean(blockB[i]);
			varB[i] = Descriptives.variance(blockB[i]);
			pb.iterateSynched();
		});
		pb.close();
	}

	public double[][] readBlock(int row, int rowsPerBlock, DoubleMatrixDatasetRandomAccessReader reader, double[][] blockA) throws IOException {

		ProgressBar pb1 = new ProgressBar(rowsPerBlock, "Reading block: " + row + " - " + (row + rowsPerBlock));
		int endrow = row + rowsPerBlock;
		for (int rowA = row; rowA < endrow; rowA++) {
			blockA[rowA - row] = reader.getRow(rowA);
			pb1.iterate();
		}
		pb1.close();
		return blockA;
	}

	private void correlate(final double[][] blockA, final double[] meanA, final double[] varA, int rowA, final double[][] blockB, final double[] meanB, final double[] varB, int rowB, double[][] correlations, boolean triangular) {


		ProgressBar pb = new ProgressBar(blockA.length, "Performing correlations for row " + rowA + " - " + (rowA + blockA.length) + " blockA and " + rowB + " - " + (rowB + blockB.length) + " blockB");
		IntStream.range(0, blockA.length).parallel().forEach(i -> {
			if (triangular) {
				for (int j = i; j < blockB.length; j++) {
					correlations[i][rowB + j] = Correlation.correlate(blockA[i], blockB[j], meanA[i], meanB[j], varA[i], varB[j]);

				}
			} else {
				for (int j = 0; j < blockB.length; j++) {
					correlations[i][rowB + j] = Correlation.correlate(blockA[i], blockB[j], meanA[i], meanB[j], varA[i], varB[j]);
				}
			}
			pb.iterateSynched();
		});
		pb.close();
	}

	private double[][] rank(double[][] blockA) {
		ProgressBar pb = new ProgressBar(blockA.length, "Ranking..");
		IntStream.range(0, blockA.length).parallel().forEach(i -> {
			RankingAlgorithm COV_RANKER_TIE = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);
			blockA[i] = COV_RANKER_TIE.rank(blockA[i]);
			pb.iterateSynched();
		});
		pb.close();
		return blockA;
	}


	private static final Logger LOG = LogManager.getLogger(CorrelateDiskBased450K.class);

	private static String humanReadableByteCount(long bytes, boolean si) {
		final int unit = si ? 1000 : 1024;
		if (bytes < unit) {
			return bytes + " B";
		}
		final int exp = (int) (Math.log(bytes) / Math.log(unit));
		final String pre = (si ? "kMGTPE" : "KMGTPE").charAt(exp - 1) + (si ? "" : "i");

		return String.format("%.1f %sB", bytes / Math.pow(unit, exp), pre);
	}




	private void destripe(float[] results, double[][] correlations, int blockSizeA, int startA, int blockSizeB, int startB) {

//		System.out.println("Destriping: " + results.length + "\t" + blockSizeA + "\t" + blockSizeB + "\t" + startA + "\t" + startB + "\t" + correlations.length + "\t" + correlations[0].length);
		IntStream.range(0, blockSizeA).parallel().forEach(i -> {
			for (int j = 0; j < blockSizeB; j++) {
//				System.out.println(i + "\t" + j + "\t" + correlations.length + "\t" + correlations[0].length + "\t" + results.length + "\t" + (startA + i) + "\t" + (startB + j) + "\t" + (startA + i) + "\t" + ((i * blockSizeA) + j));
				correlations[startA + i][startB + j] = results[(i * blockSizeB) + j];
			}
		});
	}




	private float[] initStripe(int nrRows, int nrcols) {
		int nrvals = (nrRows) * nrcols;
		return new float[nrvals];
	}

	private float[] initStripe(int nrrows) {
		return new float[nrrows];
	}

	private void stripe(double[][] blockA, int startrow, int endrow, float[] output) {
		int nrcols = blockA[0].length;
		int nrvals = (endrow - startrow) * nrcols;

		int ctr = 0;
		for (int i = startrow; i < endrow; i++) {
			for (int j = 0; j < nrcols; j++) {
				output[ctr] = (float) blockA[i][j];
				ctr++;
			}
		}
//		IntStream.range(startrow, endrow).parallel().forEach(i -> {
//			int ctr = 0;
//			int start = i * nrcols;
//
//		});

	}

	private void stripe(double[] arrA, int startrow, int endrow, float[] output) {
		int ctr = 0;
		for (int i = startrow; i < endrow; i++) {
			output[ctr] = (float) arrA[i];
			ctr++;
		}
	}
}
