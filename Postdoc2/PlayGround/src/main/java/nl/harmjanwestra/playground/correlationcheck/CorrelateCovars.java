package nl.harmjanwestra.playground.correlationcheck;

import org.checkerframework.checker.units.qual.C;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.util.Primitives;

import java.util.ArrayList;
import java.util.HashSet;

public class CorrelateCovars {

	public static void main(String[] args) {
		CorrelateCovars c = new CorrelateCovars();
		if(args.length<4){
			System.out.println("Usage: f1 f2 samplelist outfile");
		} else {
			try {
				c.run(args[0], args[1], args[2], args[3]);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	public void run(String f1, String f2, String sampleLimitFile, String outfile) throws Exception {
		TextFile tf = new TextFile(sampleLimitFile, TextFile.R);
		HashSet<String> requestedRows = new HashSet<>();
		String ln = tf.readLine();
		while (ln != null) {
			requestedRows.add(ln);
			ln = tf.readLine();
		}
		tf.close();

		// samples on the rows
		DoubleMatrixDataset<String, String> d1 = DoubleMatrixDataset.loadSubsetOfTextDoubleData(f1, '\t', requestedRows, null);
		DoubleMatrixDataset<String, String> d2 = DoubleMatrixDataset.loadSubsetOfTextDoubleData(f2, '\t', requestedRows, null);

		TextFile output = new TextFile(outfile, TextFile.W);
		for (int c1 = 0; c1 < d1.columns(); c1++) {
			for (int c2 = 0; c2 < d2.columns(); c2++) {
				ArrayList<Double> x = new ArrayList<>();
				ArrayList<Double> y = new ArrayList<>();
				for (int r1 = 0; r1 < d1.rows(); r1++) {
					String id1 = d1.getRowObjects().get(r1);
					Integer r2 = d2.getRowIndex(id1);
					if (r2 != null && r2 > -1) {
						x.add(d1.getElement(r1, c1));
						y.add(d2.getElement(r2, c2));
					}
				}
				double corr = Correlation.correlate(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
				System.out.println();
				output.writeln(d1.getColObjects().get(c1) + "\t" + d2.getColObjects().get(c2) + "\t" + corr);
			}
		}
		output.close();

	}
}
