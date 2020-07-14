package nl.harmjanwestra.methylation;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRandomAccessReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class Merge450KDiskBased {


	public void run(String probelist, String inT1, String inT2, String out) throws IOException {

		System.out.println("Merging:");
		System.out.println(inT1);
		System.out.println(inT2);

		DoubleMatrixDatasetRandomAccessReader t2reader = new DoubleMatrixDatasetRandomAccessReader(inT2);
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
			int ctr = 0;
			ArrayList<String> probesT1 = t1reader.getColObjects();
			for (String s : probesT1) {
				probenames.add(s);
				probeMap.put(s, ctr);
				ctr++;
			}
			System.out.println(probeMap.size() + " unique probes in d1");
			ArrayList<String> probesT2 = t2reader.getColObjects();
			for (String s : probesT2) {
				if (!probeMap.containsKey(s)) {
					probenames.add(s);
					probeMap.put(s, ctr);
					ctr++;
				}
			}
			System.out.println(probeMap.size() + " unique probes after merging in d2");
		}

		System.out.println("Col/probe map has " + probeMap.size() + " elements.");


		// create probe index
		ArrayList<String> probesT2 = new ArrayList<String>(t2reader.getColObjects());
		int[] t2colmap = new int[t2reader.cols()];
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

		// index samples
		ArrayList<String> samples = new ArrayList<String>();
		HashMap<String, Integer> sampleMap = new HashMap<String, Integer>();

		samples.addAll(t1reader.getRowObjects());
		int sctr = 0;
		for (int s = 0; s < samples.size(); s++) {
			sampleMap.put(samples.get(s), sctr);
			sctr++;
		}

		ArrayList<String> samples2 = new ArrayList<String>(t2reader.getRowObjects());
		for (int s = 0; s < samples2.size(); s++) {
			if (!sampleMap.containsKey(samples2.get(s))) {
				samples.add(samples2.get(s));
				sampleMap.put(samples2.get(s), sctr);
				sctr++;
			}
		}

		// create index
		int[] d1sampleindex = new int[samples.size()];
		int[] d2sampleindex = new int[samples.size()];

		for (int s = 0; s < samples.size(); s++) {
			String samplename = samples.get(s);
			Integer id1 = (Integer) t1reader.getHashRows().get(samplename);
			Integer id2 = (Integer) t2reader.getHashRows().get(samplename);
			if (id1 == null) {
				d1sampleindex[s] = -1;
			} else {
				d1sampleindex[s] = id1;
			}
			if (id2 == null) {
				d2sampleindex[s] = -1;
			} else {
				d2sampleindex[s] = id2;
			}


		}

		System.out.println("Merged matrix will be " + samples.size() + " x " + probenames.size());
		DoubleMatrixDatasetAppendableWriter writer = new DoubleMatrixDatasetAppendableWriter(probenames, out);
		double[] output = new double[probenames.size()];
		ProgressBar pb = new ProgressBar(samples.size(), "Merging rows..");
		for (int s = 0; s < samples.size(); s++) {
			Arrays.fill(output, Double.NaN);


			int rowid1 = d1sampleindex[s];
			if (rowid1 > -1) {

				// read the row
				double[] row = t1reader.getRow(rowid1);

				// merge T1 probes
				for (int i = 0; i < row.length; i++) {
					int id = t1colmap[i];
					if (id > -1) {
						output[id] = row[i];
					}
				}
			}

			int rowid2 = d2sampleindex[s];
			if (rowid2 > -1) {

				// read the row
				double[] row = t2reader.getRow(rowid2);

				// merge T1 probes
				for (int i = 0; i < row.length; i++) {
					int id = t2colmap[i];
					if (id > -1) {
						output[id] = row[i];
					}
				}
			}
			writer.append(output, samples.get(s));
			pb.set(s);
		}

		pb.close();
		writer.close();
		t2reader.close();
		t1reader.close();

	}
}
