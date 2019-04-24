package nl.harmjanwestra.playground.methylation;

import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.stream.IntStream;

public class Log2TransformDiskBased450K {

	public void run(String in, String out) throws IOException {



		double minValue = Double.MAX_VALUE;

		DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(in);

		int probeCount = it.getNrCols();
		int sampleCount = it.getNrRows();

		for (double[] row : it) {
			for (double d : row) {
				if (!Double.isNaN(d) && d < minValue) {
					minValue = d;
				}
			}
		}
		it.close();


		System.out.println("\nLog2 transforming data: Absolute minimum Expression Value:\t" + minValue);
		double multiplier = 1.0d / Math.log10(2.0d);
		double finalMinValue = minValue;

		it = new DoubleMatrixDatasetRowIterable(in);
		ArrayList<String> rows = new ArrayList<String>(it.getRows());
		DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(new ArrayList<>(it.getCols()), out);
		int rctr = 0;
		double[] rowoutput = new double[it.getNrCols()];
		for (double[] row : it) {

			for (int i = 0; i < it.getNrCols(); i++) {
				if (Double.isNaN(row[i])) {
					rowoutput[i] = Double.NaN;
				} else {
					if (minValue <= 0) {
						rowoutput[i] = Math.log10(row[i] - finalMinValue + 1) * multiplier;
					} else {
						rowoutput[i] = Math.log10(row[i] + 1) * multiplier;
					}
				}
			}
			writer.append(rowoutput, rows.get(rctr));
			rctr++;
		}
		writer.close();
		it.close();

	}
}
