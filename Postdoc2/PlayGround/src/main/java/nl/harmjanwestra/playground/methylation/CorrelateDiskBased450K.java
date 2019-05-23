package nl.harmjanwestra.playground.methylation;

import com.aparapi.Kernel;
import com.aparapi.Range;
import com.aparapi.device.Device;
import com.aparapi.device.OpenCLDevice;
import com.aparapi.internal.kernel.KernelManager;
import org.apache.commons.math3.exception.MathInternalError;
import org.apache.commons.math3.exception.NotANumberException;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessReader;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.util.RunTimer;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class CorrelateDiskBased450K {

	int clBlockSize = 128;
	private LinkedHashSet<Device> devices;
	private OpenCLDevice device;
	private KernelManager manager;

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
		if (manager == null && ocl) {
			initOCL();
			System.exit(9);
		}
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


	private static final Logger LOG = Logger.getLogger(CorrelateDiskBased450K.class);

	private static String humanReadableByteCount(long bytes, boolean si) {
		final int unit = si ? 1000 : 1024;
		if (bytes < unit) {
			return bytes + " B";
		}
		final int exp = (int) (Math.log(bytes) / Math.log(unit));
		final String pre = (si ? "kMGTPE" : "KMGTPE").charAt(exp - 1) + (si ? "" : "i");

		return String.format("%.1f %sB", bytes / Math.pow(unit, exp), pre);
	}

	public void initOCL() {
		manager = KernelManager.instance();

		List<OpenCLDevice> alldevices = OpenCLDevice.listDevices(Device.TYPE.GPU);

		device = alldevices.get(0);
		devices = new LinkedHashSet<>();
		devices.add(device);
		final long globalMemSize = device.getGlobalMemSize();
		// final long maxMemAllocSize = Math.max((globalMemSize/4), 128*1024*1024);
		final long maxMemAllocSize = device.getMaxMemAllocSize();
		final long localmemsize = device.getLocalMemSize();

		int maxcompute = device.getMaxComputeUnits();
		LOG.debug("Available OpenCL globalMemSize: " + humanReadableByteCount(globalMemSize, true));
		LOG.debug("Available OpenCL maxMemAllocSize: " + humanReadableByteCount(maxMemAllocSize, true));
		LOG.debug("Available OpenCL localMemSize: " + humanReadableByteCount(localmemsize, true));
		LOG.debug("Available OpenCL compute units: " + maxcompute);


	}

	private void correlateOCL(final double[][] blockA, final double[] meanA, final double[] var1, int rowA, final double[][] blockB, final double[] meanB, final double[] var2, int rowB, double[][] correlations, boolean triangular) {

		if (manager == null) {
			initOCL();
		}

		int nrcols = blockA[0].length;
		int nrRowsB = blockB.length;
		int bytesToAlloc = ((blockA[0].length * 4) * 2) + 12; // 4xrowlen*2 (blockA+B) + mean + var + result

		LOG.debug("Bytes to alloc per comparison: " + bytesToAlloc);
		LOG.debug("Max GPU block size (per matrix): " + clBlockSize);

		int blockSizeA = clBlockSize;
		int blockSizeB = clBlockSize;

		if (blockSizeA > blockA.length) {
			blockSizeA = blockA.length;
		}
		if (blockSizeB > blockB.length) {
			blockSizeB = blockB.length;
		}

		LOG.debug("Max GPU block size (per matrix, adjusted for blockA size): " + humanReadableByteCount(blockSizeA, true) + " and " + humanReadableByteCount(blockSizeB, true));
		LOG.debug("Eventual size of array: " + humanReadableByteCount((blockSizeA * nrcols * 4), true) + " and " + humanReadableByteCount((blockSizeB * nrcols * 4), true));

//		if ((blockSizeA * nrcols * 4) < 0 || (blockSizeB * nrcols * 4) < 0 ||
//				((blockSizeA * nrcols * 4) + (blockSizeB * nrcols * 4)) > 2 * 1048576 * 1024) {
//
//			System.out.println("Memory usage exceeding 2Gb");
//
//			System.exit(-1);
//		}


		int row = 0;


		float[] blockAs = initStripe(blockSizeA, nrcols);
		float[] meanAs = initStripe(blockSizeA);
		float[] varAs = initStripe(blockSizeB);

		float[] blockBs = initStripe(blockSizeB, nrcols);
		float[] meanBs = initStripe(blockSizeB);
		float[] varBs = initStripe(blockSizeB);

		float[] results = initStripe(blockSizeA, blockSizeB);

		OCLCorrelKernel kernel = new OCLCorrelKernel(blockAs, blockSizeA, meanAs, varAs, blockBs, blockSizeB, meanBs, varBs, nrcols, results);
		Range range = Range.create2D(device, blockSizeA, blockSizeB);
		kernel.setExplicit(false);
		manager.setPreferredDevices(kernel, devices);
		ProgressBar pb = new ProgressBar(blockA.length, "Calculating correlations...");
		while (row < blockA.length) {

			if (row + clBlockSize > blockA.length) {
				blockSizeA = blockA.length % clBlockSize;


				blockAs = initStripe(blockSizeA, nrcols);
				meanAs = initStripe(blockSizeA);
				varAs = initStripe(blockSizeB);
				results = initStripe(blockSizeA, blockSizeB);

				kernel.dispose();
				kernel = new OCLCorrelKernel(blockAs, blockSizeA, meanAs, varAs, blockBs, blockSizeB, meanBs, varBs, nrcols, results);
				range = Range.create2D(device, blockSizeA, blockSizeB);
				kernel.setExplicit(false);
				manager.setPreferredDevices(kernel, devices);

			}
			int endrow = row + blockSizeA;

			stripe(blockA, row, endrow, blockAs);
			stripe(meanA, row, endrow, meanAs);
			stripe(var1, row, endrow, varAs);

			kernel.put(blockAs);
			kernel.put(meanAs);
			kernel.put(varAs);

			int row2 = 0;

			while (row2 < blockB.length) {
				if (row2 + clBlockSize > blockB.length) {
					blockSizeB = blockB.length % clBlockSize;

					blockBs = initStripe(blockSizeB, nrcols);
					meanBs = initStripe(blockSizeB);
					varBs = initStripe(blockSizeB);
					results = initStripe(blockSizeA, blockSizeB);

					kernel.dispose();
					kernel = new OCLCorrelKernel(blockAs, blockSizeA, meanAs, varAs, blockBs, blockSizeB, meanBs, varBs, nrcols, results);
					range = Range.create2D(device, blockSizeA, blockSizeB);
					kernel.setExplicit(false);
					manager.setPreferredDevices(kernel, devices);
					kernel.put(results);
					kernel.put(blockAs);
					kernel.put(meanAs);
					kernel.put(varAs);
				}

				int endrow2 = row2 + blockSizeB;
				stripe(blockB, row2, endrow2, blockBs);
				stripe(meanB, row2, endrow2, meanBs);
				stripe(var2, row2, endrow2, varBs);


				// correlate
				kernel.put(blockBs);

				kernel.put(meanBs);
				kernel.put(varBs);

				kernel.execute(range);

				destripe(results, correlations, blockSizeA, rowA + row, blockSizeB, rowB + row2);
				row2 += blockSizeB;


			}


			// reinit block sizes
			blockSizeB = clBlockSize;
			blockBs = initStripe(blockSizeB, nrcols);
			meanBs = initStripe(blockSizeB);
			varBs = initStripe(blockSizeB);
			results = initStripe(blockSizeA, blockSizeB);

			row += blockSizeA;
			pb.set(row);
		}
		pb.close();

		StringBuilder builder = new StringBuilder();
		manager.reportDeviceUsage(builder, true);
		System.out.println(builder.toString());
		kernel.dispose();


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


	public class OCLCorrelKernel extends Kernel {
		final float[] matrixA;
		final float[] matrixB;
		final int nrCols;

		final int nrRowsA;
		final int nrRowsB;
		float[] results;
		float[] varA;
		float[] varB;
		float[] meanA;
		float[] meanB;

		public OCLCorrelKernel(float[] matrixA, int nrRowsA, float[] meanA, float[] varA, float[] matrixB, int nrRowsB, float[] meanB, float[] varB, int nrCols, float[] results) {
			this.matrixA = matrixA;
			this.matrixB = matrixB;
			this.nrCols = nrCols;
			this.nrRowsA = nrRowsA;
			this.nrRowsB = nrRowsB;
			this.meanA = meanA;
			this.meanB = meanB;
			this.varA = varA;
			this.varB = varB;
			this.results = results;
		}

		@Override
		public void run() {
			final int i = this.getGlobalId(0);

			if (i < nrRowsA) {
				final int j = this.getGlobalId(1);
				if (j < nrRowsB) {
					results[(i * nrRowsB) + j] = correlate(i, j);
				}
			}
		}

		private float correlate(int i, int j) {
			float ans = 0f;
			float denom = (float) Math.sqrt(varA[i] * varB[j]);
			if (denom == 0) {
				if (varA[i] == 0 && varB[j] == 0) {
					return 0;
				} else {
					return 1;
				}
			} else {
				for (int q = 0; q < nrCols; q++) {
					ans += (matrixA[(i * nrCols) + q] - meanA[i]) * (matrixB[(j * nrCols) + q] - meanB[j]);
				}
				ans /= (nrCols - 1);
				ans /= denom;
				return ans;
			}
		}

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
