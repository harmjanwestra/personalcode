package nl.harmjanwestra.playground.conv;

import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class CleanExpTable {
	
	public static void main(String[] args) {
		CleanExpTable c = new CleanExpTable();
		try {
			c.runFeatureCounts(args[0], args[1]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void run(String in, String out) throws IOException {
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tf2 = new TextFile(out, TextFile.W);
		String[] elems = tf.readLineElems(Strings.tab);
		while (elems != null) {
			elems[1] = elems[1].split("\\.")[0];
			tf2.writeln(Strings.concat(elems, Strings.tab, 1, elems.length));
			elems = tf.readLineElems(Strings.tab);
		}
		tf2.close();
		tf.close();
	}
	
	public static void runFeatureCounts(String in, String out) throws IOException {
		
		String[] datasets = in.split(",");
		
		ArrayList<DoubleMatrixDataset<String, String>> list = new ArrayList<DoubleMatrixDataset<String, String>>();
		for (int i = 0; i < datasets.length; i++) {
			System.out.println("parsing: " + datasets[i]);
			
			TextFile tf = new TextFile(datasets[i], TextFile.R);
			TextFile tf2 = new TextFile(datasets[i] + "-filtered.txt.gz", TextFile.W);
			
			tf2.writeln(tf.readLine()); // header
			tf.readLine(); //  unmapped
			tf.readLine(); //  multimapped
			tf.readLine(); // nofeature
			tf.readLine(); // ambiguous
			String[] elems = tf.readLineElems(Strings.tab); // ensg.1
			int written = 0;
			int read = 0;
			while (elems != null) {
				elems[0] = elems[0].split("\\.")[0];
				int nrzeros = 0;
				for (int q = 1; q < elems.length; q++) {
					if (elems[q].equals("0")) {
						nrzeros++;
					}
				}
				if (nrzeros < 0.5 * (elems.length - 1)) {
					tf2.writeln(Strings.concat(elems, Strings.tab));
				}
				
				
				elems = tf.readLineElems(Strings.tab);
			}
			tf2.close();
			tf.close();
			
			
			DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(datasets[i] + "-filtered.txt.gz");
			list.add(ds);
		}
		
		// now merge datasets
		HashMap<String, Integer> inds = new HashMap<>();
		ArrayList<String> indlist = new ArrayList<>();
		HashMap<String, Integer> genes = new HashMap<String, Integer>();
		ArrayList<String> geneList = new ArrayList<>();
		HashMap<String, Integer> genectr = new HashMap<>();
		int indctr = 0;
		int gctr = 0;
		for (DoubleMatrixDataset<String, String> ds : list) {
			ArrayList<String> cols = ds.getColObjects();
			for (String col : cols) {
				if (inds.containsKey(col)) {
					System.out.println("ERROR: duplicate ind id: " + col);
				} else {
					inds.put(col, indctr);
					indlist.add(col);
					indctr++;
				}
			}
			ArrayList<String> rows = ds.getRowObjects();
			for (String row : rows) {
				if (!genes.containsKey(row)) {
					genes.put(row, gctr);
					geneList.add(row);
					
					gctr++;
				}
				Integer ctr = genectr.get(row);
				if (ctr == null) {
					ctr = 1;
				} else {
					ctr++;
				}
				genectr.put(row, ctr);
				
			}
		}
		System.out.println(geneList.size() + " genes total ");
		HashMap<String, Integer> finalgenes = new HashMap<>();
		ArrayList<String> finalGenelis = new ArrayList<>();
		int ctr2 = 0;
		for (String gene : geneList) {
			if (genectr.get(gene).equals(datasets.length)) {
				finalGenelis.add(gene);
				finalgenes.put(gene, ctr2);
				ctr2++;
			}
		}
		geneList = finalGenelis;
		genes = finalgenes;
		System.out.println(geneList.size() + " present in all three datasets.");
		// very inefficient :/
		double[][] data = new double[genes.size()][inds.size()];
		int dctr = 0;
		int pos = 0;
		for (DoubleMatrixDataset<String, String> ds : list) {
			System.out.println("merging " + dctr);
			
			for (int row = 0; row < ds.rows(); row++) {
				String gene = ds.getRowObjects().get(row);
				double[] datasetRow = ds.getRow(row).toArray();
				Integer geneId = genes.get(gene);
				if (geneId != null) {
					System.arraycopy(datasetRow, 0, data[geneId], pos, datasetRow.length);
				}
			}
			dctr++;
			pos += ds.columns();
			System.out.println(pos);
		}
		
		// remove nulls
		TextFile outf = new TextFile(out, TextFile.W);
		String header = "-\t" + Strings.concat(indlist, Strings.tab);
		outf.writeln(header);
		double threshold = 0.5 * indlist.size();
		int written = 0;
		for (int row = 0; row < data.length; row++) {
			int nrzeros = 0;
			String ln = geneList.get(row);
			for (int col = 0; col < data[row].length; col++) {
				if (data[row][col] == 0) {
					nrzeros++;
				}
		
			}
			
			if (nrzeros < threshold) {
				ln += "\t" + Strings.concat(data[row], Strings.tab);
				outf.writeln(ln);
				written++;
			}
			System.out.print(row + "\t" + written + "\r");
		}
		outf.close();
		System.out.println();
	}
}
