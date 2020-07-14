package nl.harmjanwestra.methylation;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class DatasetFilter {


	public enum FILTERMODE {
		INCLUDE,
		EXCLUDE;

		public static FILTERMODE parse(String colFilterModeStr) {
			if (colFilterModeStr.toLowerCase().equals("exclude")) {
				return EXCLUDE;
			} else {
				return INCLUDE;
			}

		}
	}

	;

	public void run(String in, String out, String rowfilter, String rowFilterModeStr, String colfilter, String colFilterModeStr) throws IOException {

		FILTERMODE rowFilterMode = FILTERMODE.parse(rowFilterModeStr);
		FILTERMODE colFilterMode = FILTERMODE.parse(colFilterModeStr);


		DoubleMatrixDatasetRowIterable it = new DoubleMatrixDatasetRowIterable(in);

		HashSet<String> rowset = null;
		ArrayList<String> rowIds = new ArrayList<>(it.getRows());
		if (rowfilter != null && !rowfilter.equals("null")) {
			System.out.println("Will filter rows on: " + rowfilter + " with mode " + rowFilterMode);
			rowset = readset(rowfilter);
		}

		boolean[] includecols = null;

		ArrayList<String> newColIds = null;
		if (colfilter != null && !colfilter.equals("null")) {
			System.out.println("Will filter cols on: " + colfilter + " with mode " + colFilterMode);
			HashSet<String> colset = readset(colfilter);
			includecols = new boolean[it.getNrRows()];
			int cctr = 0;
			for (String s : it.getCols()) {
				if (colFilterMode.equals(FILTERMODE.INCLUDE) && colset.contains(s)
						|| colFilterMode.equals(FILTERMODE.EXCLUDE) && !colset.contains(s)) {
					includecols[cctr] = true;
					newColIds.add(s);
				}
			}
		} else {
			newColIds = new ArrayList<>(it.getCols());
		}

		int nrcols = newColIds.size();
		DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(newColIds, out);

		int rctr = 0;
		for (double[] row : it) {
			String rowId = rowIds.get(rctr);

			if (rowset == null ||
					rowFilterMode.equals(FILTERMODE.INCLUDE) && rowset.contains(rowId) ||
					rowFilterMode.equals(FILTERMODE.EXCLUDE) && !rowset.contains(rowId)) {
				double[] output = filterRow(row, includecols, nrcols);
				writer.append(output, rowId);
			}
			rctr++;
		}

	}

	private double[] filterRow(double[] row, boolean[] includecols, int nrrows) {
		double[] output = new double[nrrows];
		int c = 0;
		for (int d = 0; d < row.length; d++) {
			if (includecols[d]) {
				output[c] = row[d];
			}
			c++;
		}
		return output;
	}

	private HashSet<String> readset(String rowfilter) throws IOException {
		HashSet<String> set = new HashSet<String>();
		TextFile tf = new TextFile(rowfilter, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			set.add(elems[0]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(rowfilter + " has " + set.size() + " unique elements to filter.");
		return set;
	}
}
