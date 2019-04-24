package nl.harmjanwestra.playground.biogen.freeze2;

import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class GeneticSimFilter {


	public static void main(String[] args) {

		String geneticsimfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\geneticsimfiles.txt";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\duplicates\\";
		try {
			TextFile tf = new TextFile(geneticsimfile, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String ds = elems[0];
				String dsname = elems[1];
				GeneticSimFilter f = new GeneticSimFilter();
				f.runplinkgenome(ds, out + "Duplicates-" + dsname + ".txt", 0.125);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
//
//		GeneticSimFilter s = new GeneticSimFilter();
//		try {
//			s.runplinkgenome("D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-genome.genome.gz",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\MSBB-dups.txt", 0.125);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
	}


	public void runplinkgenome(String file, String fileout, double threshold) throws IOException {
		System.out.println(file);
		TextFile tf = new TextFile(file, TextFile.R);
		String ln = tf.readLine();
		int lnctr = 0;
		int picol = -1;
		int sa1col = -1;
		int sa2col = -1;
		ArrayList<Pair<String, String>> dups = new ArrayList<>();
		while (ln != null) {
			while (ln.contains("  ")) {
				ln = ln.replaceAll("  ", " ");
			}
			while (ln.startsWith(" ")) {
				ln = ln.substring(1);
			}
			String[] elems = ln.split(" ");

			if (lnctr == 0) {
				for (int i = 0; i < elems.length; i++) {
					String e = elems[i];
					if (e.equals("IID1")) {
						sa1col = i;
					} else if (e.equals("IID2")) {
						sa2col = i;
					} else if (e.equals("PI_HAT")) {
						picol = i;
					}
				}
				System.out.println("SA1: " + sa1col);
				System.out.println("SA2: " + sa2col);
				System.out.println("PIH: " + picol);

			} else {
				String s1 = elems[sa1col];
				String s2 = elems[sa2col];
				String pi = elems[picol];

				double pid = Double.parseDouble(pi);
				if (pid > threshold && !s1.equals(s2)) {
					dups.add(new Pair<String, String>(s1, s2));
				}
			}
			lnctr++;
			ln = tf.readLine();
		}
		tf.close();


		writedups(fileout, dups);
	}

	public void runGeneticSim(String file, String fileout, double threshold) throws IOException {


		TextFile tf = new TextFile(file, TextFile.R);

		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);

		ArrayList<Pair<String, String>> dups = new ArrayList<>();

		while (elems != null) {

			String sample1 = elems[0];
			String sample2 = elems[1];

			if (!sample1.equals(sample2)) {
				double d = Double.parseDouble(elems[3]);
				if (d >= threshold) {
					dups.add(new Pair<String, String>(sample1, sample2));
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		writedups(fileout, dups);
	}

	public void writedups(String fileout, ArrayList<Pair<String, String>> dups) throws IOException {
		HashMap<String, HashSet<String>> dupset = new HashMap<String, HashSet<String>>();
		for (Pair<String, String> sample : dups) {
			String sample1 = sample.getLeft();
			String sample2 = sample.getRight();
			HashSet<String> dupsforsample1 = dupset.get(sample1);
			HashSet<String> dupsforsample2 = dupset.get(sample2);
			if (dupsforsample1 == null && dupsforsample2 == null) {
				dupsforsample1 = new HashSet<>();
				dupsforsample1.add(sample2);
				dupset.put(sample1, dupsforsample1);
			} else if (dupsforsample1 == null && dupsforsample2 != null) {

				dupsforsample2.add(sample1);
				dupset.put(sample2, dupsforsample2);
			} else if (dupsforsample1 != null && dupsforsample2 == null) {
				dupsforsample1.add(sample2);
				dupset.put(sample1, dupsforsample1);
			}
		}

		if (!dupset.isEmpty()) {
			HashSet<String> alreadyWritten = new HashSet<>();
			TextFile tfo = new TextFile(fileout, TextFile.W);
			for (String key : dupset.keySet()) {
				HashSet<String> duplicates = dupset.get(key);
				for (String s : duplicates) {
					if (!alreadyWritten.contains(s)) {
						tfo.writeln(s);
						alreadyWritten.add(s);
					}
				}
			}
			tfo.close();
		}
	}
}
