package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class MergeRNASeq {
	
	public static void main(String[] args) {
		MergeRNASeq m = new MergeRNASeq();
		if (args.length < 4) {
			System.out.println("Usage: exp samplestokeep mapfile out");
		} else {
			try {
				m.run(args[0], args[1], args[2], args[3]);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	public void run(String exp, String samplestokeep, String mapfile, String out) throws Exception {
		
		// load the data
		DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>().loadDoubleData(exp);
		
		HashSet<String> keep = loadHash(samplestokeep);
		
		HashMap<String, HashSet<String>> barcodesperind = new HashMap<>();
		TextFile tf = new TextFile(mapfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.comma);
		while (elems != null) {
			String ind = elems[0];
			String bar = "X" + elems[1];
			System.out.println(ind + " --> " + bar);
			if (keep == null || keep.contains(bar) && ds.getHashCols().containsKey(bar)) {
				HashSet<String> select = barcodesperind.get(ind);
				if (select == null) {
					select = new HashSet<>();
				}
				select.add(bar);
				barcodesperind.put(ind, select);
			}
			elems = tf.readLineElems(TextFile.comma);
		}
		tf.close();
		System.out.println(barcodesperind.size() + " individuals total");
		
		double[][] data = new double[ds.rows()][barcodesperind.size()];
		ArrayList<String> samples = new ArrayList<>();
		samples.addAll(barcodesperind.keySet());
		Collections.sort(samples);
		for (int i = 0; i < samples.size(); i++) {
			HashSet<String> select = barcodesperind.get(samples.get(i));
			double[] rows = new double[ds.rows()];
			for (String rnasample : select) {
				Integer rnaindid = ds.getHashCols().get(rnasample);
				for (int row = 0; row < ds.rows(); row++) {
					data[row][i] += ds.getMatrix().getQuick(row, rnaindid);
				}
				for (int row = 0; row < ds.rows(); row++) {
					data[row][i] /= select.size();
				}
			}
		}
		
		DoubleMatrixDataset<String, String> outds = new DoubleMatrixDataset<>();
		
		outds.setMatrix(data);
		outds.setRowObjects(ds.getRowObjects());
		outds.setColObjects(samples);
		outds.save(out);
	}
	
	private HashSet<String> loadHash(String rnasamplestokeep) throws IOException {
		HashSet<String> rnakeephash = null;
		if (rnasamplestokeep != null) {
			TextFile tfk1 = new TextFile(rnasamplestokeep, TextFile.R);
			ArrayList<String> rnasamplestokeeparr = tfk1.readAsArrayList();
			tfk1.close();
			
			rnakeephash = new HashSet<>();
			rnakeephash.addAll(rnasamplestokeeparr);
		}
		return rnakeephash;
	}
}
