package nl.harmjanwestra.methylation;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;

public class DeTriangularize {

	public void run(String in, String out) throws IOException {
		System.out.println("Detriangle:");
		System.out.println(in);

		DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleBinaryData(in);

		ProgressBar pb = new ProgressBar(ds.rows());
		for (int i = 0; i < ds.rows(); i++) {
			for (int j = 0; j < ds.columns(); j++) {
				ds.setElementQuick(j, i, ds.getElementQuick(i, j));
			}
			ds.setElementQuick(i, i, 1d);
			pb.iterateSynched();
		}
		pb.close();

		System.out.println("Saving to: " + out);
		ds.save(out);


	}
}
