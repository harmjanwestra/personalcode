package nl.harmjanwestra.playground.methylation;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixConverter;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessTranspose;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;

public class TestTools {


	public static void main(String[] args) {

		TestTools t = new TestTools();
		try {
//			t.createMatrix("D:\\tmp\\MatrixA.txt", 2250, 0, 1500, 0);
//			t.createMatrix("D:\\tmp\\MatrixB.txt", 1000, 2250, 1500, 0);

			DoubleMatrixConverter.TextToBinary("D:\\tmp\\MatrixA.txt", "D:\\tmp\\MatrixA");
			DoubleMatrixConverter.TextToBinary("D:\\tmp\\MatrixB.txt", "D:\\tmp\\MatrixB");

			Merge450KDiskBased d = new Merge450KDiskBased();
			d.run(null, "D:\\tmp\\MatrixA", "D:\\tmp\\MatrixB", "d:\\tmp\\MatrixC");

//			compare("D:\\tmp\\MatrixA", "d:\\tmp\\MatrixC");
//			compare("D:\\tmp\\MatrixB", "d:\\tmp\\MatrixC");

			testTranspose("d:\\tmp\\MatrixC");

//			DoubleMatrixConverter.BinaryToText("d:\\tmp\\MatrixC", "d:\\tmp\\MatrixC.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void testTranspose(String matrixAloc) throws IOException {

		DoubleMatrixDataset<String, String> matrixA = DoubleMatrixDataset.loadDoubleBinaryData(matrixAloc);
		matrixA = matrixA.viewDice();
		DoubleMatrixDatasetRandomAccessTranspose t = new DoubleMatrixDatasetRandomAccessTranspose();
		t.transposeLargeMatrix(matrixAloc, matrixAloc + "T", 55);

		DoubleMatrixDataset<String, String> matrixB = DoubleMatrixDataset.loadDoubleBinaryData(matrixAloc+"T");

		ProgressBar pb = new ProgressBar(matrixA.rows(), "Correlating");
		int[] bins = new int[10];
		for (int r = 0; r < matrixA.rows(); r++) {
			double[] a = matrixA.getRow(r).toArray();
			String name = matrixA.getRowObjects().get(r);

			Integer id = matrixB.getHashRows().get(name);
			if (id != null) {
				double[] btmp = matrixB.getRow(id).toArray();
				double[] b = new double[matrixA.columns()];
				for (int c = 0; c < matrixA.columns(); c++) {
					String namecol = matrixA.getColObjects().get(c);
					Integer idcol = matrixB.getHashCols().get(namecol);
					if (idcol != null) {
						b[c] = btmp[idcol];
					}
				}

				double cor = Correlation.correlate(a, b);
				int bin = (int) Math.floor(10 * cor);
				if (bin == 10) {
					bin = 9;
				}
				bins[bin]++;
				if (cor < 0.95) {
					System.out.println(name + " poor correlation " + cor + " between " + matrixAloc + " and disk based transpose.");
				}

			} else {
				System.out.println("row " + name + " not present in " + matrixB);
			}

			pb.iterate();
			;
		}
		pb.close();

		System.out.println("rows tested: " + matrixA.rows());
		for (int i = 0; i < bins.length; i++) {
			System.out.println(i + "\t" + bins[i] + "\t" + ((double) bins[i] / matrixA.rows()));
		}
	}

	private static void compare(String matrixAloc, String matrixBloc) throws IOException {
		System.out.println("testing " + matrixAloc + " vs " + matrixBloc);
		DoubleMatrixDataset<String, String> matrixA = DoubleMatrixDataset.loadDoubleBinaryData(matrixAloc);
		DoubleMatrixDataset<String, String> matrixB = DoubleMatrixDataset.loadDoubleBinaryData(matrixBloc);

		ProgressBar pb = new ProgressBar(matrixA.rows(), "Correlating");
		int[] bins = new int[10];
		for (int r = 0; r < matrixA.rows(); r++) {
			double[] a = matrixA.getRow(r).toArray();
			String name = matrixA.getRowObjects().get(r);

			Integer id = matrixB.getHashRows().get(name);
			if (id != null) {
				double[] btmp = matrixB.getRow(id).toArray();
				double[] b = new double[matrixA.columns()];
				for (int c = 0; c < matrixA.columns(); c++) {
					String namecol = matrixA.getColObjects().get(c);
					Integer idcol = matrixB.getHashCols().get(namecol);
					if (idcol != null) {
						b[c] = btmp[idcol];
					}
				}

				double cor = Correlation.correlate(a, b);
				int bin = (int) Math.floor(10 * cor);
				if (bin == 10) {
					bin = 9;
				}
				bins[bin]++;
				if (cor < 0.95) {
					System.out.println(name + " poor correlation " + cor + " between " + matrixAloc + " and " + matrixBloc);
				}

			} else {
				System.out.println("row " + name + " not present in " + matrixB);
			}

			pb.iterate();
			;
		}
		pb.close();

		System.out.println("rows tested: " + matrixA.rows());
		for (int i = 0; i < bins.length; i++) {
			System.out.println(i + "\t" + bins[i] + "\t" + ((double) bins[i] / matrixA.rows()));
		}


	}

	public void createMatrix(String output, int rows, int rowStartId, int cols, int colStartId) throws IOException {

		TextFile tf = new TextFile(output, TextFile.W);
		String[] colheader = new String[cols];
		for (int i = 0; i < cols; i++) {
			colheader[i] = "col" + (i + colStartId);
		}

		tf.writeln("-\t" + Strings.concat(colheader, Strings.tab));
		double[] outputln = new double[cols];
		for (int r = 0; r < rows; r++) {
			String rowheader = "row" + (r + rowStartId);
			for (int c = 0; c < cols; c++) {
				outputln[c] = Math.random();
			}
			tf.writeln(rowheader + "\t" + Strings.concat(outputln, Strings.tab));
		}
		tf.close();

	}

}
