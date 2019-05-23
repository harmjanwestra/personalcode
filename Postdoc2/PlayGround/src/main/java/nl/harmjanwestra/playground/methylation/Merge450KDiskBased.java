package nl.harmjanwestra.playground.methylation;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessReader;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

public class Merge450KDiskBased {

	public void run(String probelist, String inT1, String inT2, String out) throws IOException {

		System.out.println("Merging:");
		System.out.println(inT1);
		System.out.println(inT2);

		DoubleMatrixDatasetRowIterable t2iterator = new DoubleMatrixDatasetRowIterable(inT2);
		DoubleMatrixDatasetRandomAccessReader t1reader = new DoubleMatrixDatasetRandomAccessReader(inT1);
		HashMap<String, Integer> probeMap = new HashMap<String, Integer>();
		ArrayList<String> probenames = new ArrayList<>();

		if (probelist != null) {
			TextFile tf = new TextFile(probelist, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			int ctr = 0;
			while (elems != null) {
				probenames.add(elems[0]);
				probeMap.put(elems[0], ctr);
				ctr++;
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
		} else {
			Set<String> probesT2 = t2iterator.getCols();
			int ctr = 0;
			for (String s : probesT2) {
				probenames.add(s);
				probeMap.put(s, ctr);
				ctr++;
			}
			ArrayList<String> probesT1 = t1reader.getColObjects();
			for (String s : probesT1) {
				probenames.add(s);
				probeMap.put(s, ctr);
				ctr++;
			}

		}

		System.out.println("Col/probe map has " + probeMap.size() + " elements.");


		ArrayList<String> probesT2 = new ArrayList<String>(t2iterator.getCols());
		int[] t2colmap = new int[t2iterator.getNrCols()];
		int nrcolsfound = 0;
		for (int d = 0; d < t2colmap.length; d++) {
			String probe = probesT2.get(d);
			Integer id = probeMap.get(probe);
			if (id != null) {
				t2colmap[d] = id;
				nrcolsfound++;
			} else {
				t2colmap[d] = -1;
			}
		}

		ArrayList<String> probesT1 = new ArrayList<String>(t1reader.getColObjects());
		int[] t1colmap = new int[probesT1.size()];
		for (int d = 0; d < t1colmap.length; d++) {
			String probe = probesT1.get(d);
			Integer id = probeMap.get(probe);
			if (id != null) {
				t1colmap[d] = id;
				nrcolsfound++;
			} else {
				t1colmap[d] = -1;
			}
		}
		System.out.println(nrcolsfound + " cols can be found in both datasets.");

		int[] samplemap = new int[t2iterator.getNrRows()];
		ArrayList<String> samplesT2 = new ArrayList<String>(t2iterator.getRows());
		ArrayList<String> samplesT1 = new ArrayList<String>(t2iterator.getRows());
		HashMap<String, Integer> s1samplemap = new HashMap<>();
		for (int i = 0; i < samplesT1.size(); i++) {
			s1samplemap.put(samplesT1.get(i), i);
		}

		for (int i = 0; i < samplesT2.size(); i++) {
			String sample = samplesT2.get(i);
			Integer id = s1samplemap.get(sample);
			samplemap[i] = id;
		}


		// probes are on columns?
		DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(probenames, out);
		double[] output = new double[probenames.size()];
		int row = 0;
		ProgressBar pb = new ProgressBar(samplesT2.size(), "Merging rows..");
		for (double[] d : t2iterator) {

			Arrays.fill(output, Double.NaN);
			String sampleT2 = samplesT2.get(row);

			// merge T2 probes
			for (int i = 0; i < d.length; i++) {
				int id = t2colmap[i];
				if (id > -1) {
					output[id] = d[i];
				}
			}


			// get T1 sample row
			Integer t1rowid = samplemap[row];
			double[] d2 = t1reader.getRow(t1rowid);

			// merge T1 probes
			for (int i = 0; i < d2.length; i++) {
				int id = t1colmap[i];
				if (id > -1) {
					output[id] = d[i];
				}
			}

			writer.append(output, sampleT2);
			row++;
			pb.set(row);

		}
		pb.close();
		writer.close();
		t2iterator.close();
		t1reader.close();

	}
}
