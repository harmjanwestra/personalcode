package nl.harmjanwestra.playground.methylation;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessReader;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;

public class CorrelateDiskBased450K {


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
}
