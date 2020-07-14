package nl.harmjanwestra.methylation;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;

import java.io.IOException;
import java.util.ArrayList;

public class DetermineNAs450K {

	public void run(String in, String out) throws IOException {

		TextFile outrow = new TextFile(out + "-NaNInRows.txt", TextFile.W);
		TextFile outcol = new TextFile(out + "-NaNInCols.txt", TextFile.W);


		DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(in);
		ArrayList<String> rowIds = new ArrayList<String>(it.getRows());
		ArrayList<String> colIds = new ArrayList<String>(it.getCols());

		int[] nrNaNCols = new int[it.getNrCols()];

		int rowId = 0;
		ProgressBar pb = new ProgressBar(rowIds.size(), "Counting NaNs in " + in);
		for (double[] row : it) {
			int nrNan = 0;
			for (int col = 0; col < row.length; col++) {
				if (Double.isNaN(row[col])) {
					nrNan++;
					nrNaNCols[col]++;
				}
			}
			outrow.writeln(rowIds.get(rowId) + "\t" + nrNan + "\t" + ((double) nrNan / row.length));
			pb.iterate();
			rowId++;
		}
		pb.close();

		outrow.close();
		for (int col = 0; col < nrNaNCols.length; col++) {
			outcol.writeln(colIds.get(col) + "\t" + nrNaNCols[col] + "\t" + ((double) nrNaNCols[col] / it.getNrRows()));
		}
		outcol.close();

	}
}
