package nl.harmjanwestra.playground.mat;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;

import java.io.IOException;
import java.util.HashSet;

public class MeanAndVariance {
	
	public static void main(String[] args) {
		if (args.length < 4) {
			System.out.println("Usage: input output rowselect colselect");
		} else {
			MeanAndVariance v = new MeanAndVariance();
			try {
				v.run(args[0], args[1], args[2], args[3]);
			} catch (IOException e) {
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	
	public void run(String input, String output, String rowsToSelect, String colsToSelect) throws Exception {
		
		HashSet<String> rowselect = null;
		if (rowsToSelect != null) {
			rowselect = new HashSet<String>();
			TextFile tf = new TextFile(rowsToSelect, TextFile.R);
			rowselect.addAll(tf.readAsArrayList());
			tf.close();
			
		}
		
		HashSet<String> colselect = null;
		if (colsToSelect != null) {
			colselect = new HashSet<>();
			TextFile tf = new TextFile(colsToSelect, TextFile.R);
			colselect.addAll(tf.readAsArrayList());
			tf.close();
		}
		
		DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadSubsetOfTextDoubleData(input, '\t', rowselect, colselect);
		
		// over rows
		TextFile out1 = new TextFile(output + "-MeanAndVarianceRows.txt", TextFile.W);
		out1.writeln("ID\tN\tMean\tVar");
		for (int r = 0; r < ds.rows(); r++) {
			DoubleMatrix1D row = ds.getMatrix().viewRow(r);
			double[] rowarr = row.toArray();
			String lnout = ds.getRowObjects().get(r) + "\t" + rowarr.length + "\t" + Descriptives.mean(rowarr) + "\t" + Descriptives.variance(rowarr);
			out1.writeln(lnout);
		}
		out1.close();
		
		// over rows
		TextFile out2 = new TextFile(output + "-MeanAndVarianceCols.txt", TextFile.W);
		out2.writeln("ID\tN\tMean\tVar");
		for (int r = 0; r < ds.columns(); r++) {
			DoubleMatrix1D col = ds.getMatrix().viewColumn(r);
			double[] colarr = col.toArray();
			String lnout = ds.getColObjects().get(r) + "\t" + colarr.length + "\t" + Descriptives.mean(colarr) + "\t" + Descriptives.variance(colarr);
			out2.writeln(lnout);
		}
		out2.close();
		
		
	}
}
