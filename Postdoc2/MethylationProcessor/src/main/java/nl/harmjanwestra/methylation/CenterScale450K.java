package nl.harmjanwestra.methylation;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;
import umcg.genetica.math.stats.Descriptives;

import java.io.IOException;
import java.util.ArrayList;

public class CenterScale450K {


	public void run(String in, String out) throws IOException {
		DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(in);

		DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(new ArrayList<String>(it.getCols()), out);

		ArrayList<String> rowids = new ArrayList<>(it.getRows());
		ProgressBar pb = new ProgressBar(it.getNrRows(), "Centering and scaling...");
		int rowid = 0;
		for (double[] row : it) {

			double mean = Descriptives.mean(row);
			double sd = Math.sqrt(Descriptives.variance(row));

			for (int d = 0; d < row.length; d++) {
				row[d] = (row[d] - mean) / sd;
			}

			writer.append(row, rowids.get(rowid));
			pb.iterate();
			rowid++;
		}
		pb.close();

		it.close();
		writer.close();

	}
}
