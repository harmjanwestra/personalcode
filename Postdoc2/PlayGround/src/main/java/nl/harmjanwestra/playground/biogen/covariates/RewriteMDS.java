package nl.harmjanwestra.playground.biogen.covariates;

import nl.harmjanwestra.playground.methylation.PCA;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class RewriteMDS {


	public static void main(String[] args) {

		String[] linksfiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MAYO-V2.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-MSBB-V2.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-AMPAD-ROSMAP-V2.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-Braineac.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-BrainGVEX.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-CMC.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-CMC_HBCC_set2.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-ENA.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-GTEx.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-GVEX.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-LIBD_1M.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-NABEC-H610.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-NABEC-WES.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\links\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt-TargetALS.txt"
		};
		String[] mdsfiles = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\AMPAD-MAYO-V2\\AMPAD-MAYO-V2-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\AMPAD-MSBB-V2\\AMPAD-MSBB-V2-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\AMPAD-ROSMAP-V2\\AMPAD-ROSMAP-V2-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\Braineac\\Braineac-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\BrainGVEX\\BrainGVEX-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\CMC\\CMC-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\CMC_HBCC_set2\\CMC_HBCC_set2-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\ENA\\ENA-postqc-mds-ibd.mds.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\GTEx\\GTEx-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\GVEX\\GVEX-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\LIBD_1M\\LIBD_1M-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\NABEC-H610\\NABEC-H610-mds-ibd.mds.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\NABEC-WES\\NABEC-WES-mds.txt",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\TargetALS\\TargetALS-mds-ibd.mds.gz"
		};

		String[] names = new String[]{
				"AMPAD-MAYO",
				"AMPAD-MSBB",
				"AMPAD-ROSMAP",
				"Braineac",
				"BrainGVEX",
				"CMC",
				"CMC_HBCC_set2",
				"ENA",
				"GTEx",
				"GVEX",
				"LIBD_1M",
				"NABEC-H610",
				"NABEC-WES",
				"TargetALS"
		};

		String outfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\genotypeqc\\qcfiles\\2019-06-01-CovarsAndMDS-Cortex.txt";

		String covars = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-11-Freeze2RNAQc\\2019-04-11-Freeze2.TMM.Covariates-Numeric-Top10Covariates.txt";

		RewriteMDS m = new RewriteMDS();

		try {
			m.run(linksfiles, names, mdsfiles, covars, outfile);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


	public void run(String[] linkfiles, String[] datasetnames, String[] mdsfiles, String covarfile, String output) throws Exception {

		DoubleMatrixDataset<String, String> covars = null;
		if (covarfile != null) {
			covars = DoubleMatrixDataset.loadDoubleData(covarfile);
		}

		int[] datasets = new int[linkfiles.length - 1];
		TextFile tfout = new TextFile(output, TextFile.W);

		String header = "Sample\tMDS1\tMDS2\tMDS3\tMDS4";
		if (covars != null) {
			String covarString = Strings.concat(covars.getColObjects(), Strings.tab);
			header = "Sample\t" + covarString + "\tMDS1\tMDS2\tMDS3\tMDS4";
		}
		for (int d = 0; d < datasets.length; d++) {
			header += "\t" + datasetnames[d];
		}
		tfout.writeln(header);
		int total = 0;
		int totalcovars = 0;
		for (int i = 0; i < linkfiles.length; i++) {
			datasets = new int[linkfiles.length - 1];
			HashMap<String, HashSet<String>> dnaToRNA = loadLinkFile(linkfiles[i]);
			if (i != linkfiles.length - 1) {
				datasets[i] = 1;
			}

			int matched = 0;
			int hascovars = 0;
			// parse MDS file
			TextFile tf = new TextFile(mdsfiles[i], TextFile.R);
			tf.readLine(); // skip header

			String line = tf.readLine();
			HashSet<String> foundActualDNAs = new HashSet<>();
			while (line != null) {

				while (line.startsWith(" ")) {
					line = line.substring(1);
				}
				while (line.contains("  ")) {
					line = line.replaceAll("  ", " ");
				}

				String[] elems = line.split(" ");
				String sample = elems[1];
				HashSet<String> idlist = dnaToRNA.get(sample);
				if (idlist == null) {
					idlist = dnaToRNA.get("0_" + sample + "_" + sample);
				}


				boolean idhascovars = false;
				if (idlist == null && datasetnames[i].equals("ENA")) {
//					System.out.println("Could not find sample " + sample + " for dataset: " + datasetnames[i]);
				} else if (idlist != null) {

					foundActualDNAs.add(sample);
					for (String id : idlist) {

						if (covars != null) {
							Integer covarId = covars.getHashRows().get(id);
							double[] covarvals = new double[covars.getColObjects().size()];
							if (covarId != null) {
								covarvals = covars.getRow(covarId).toArray();
								idhascovars = true;
							}
							String outln = id + "\t" + Strings.concat(covarvals, Strings.tab) + "\t" + Strings.concat(elems, Strings.tab, 3, 7) + "\t" + Strings.concat(datasets, Strings.tab);
							tfout.writeln(outln);
						} else {
							String outln = id + "\t" + Strings.concat(elems, Strings.tab, 3, 7) + "\t" + Strings.concat(datasets, Strings.tab);
							tfout.writeln(outln);
						}

					}
					matched++;
					if (idhascovars) {
						hascovars++;
					}
				}
				line = tf.readLine();
			}
			totalcovars += hascovars;
			total += matched;

			tf.close();

			for (String s : dnaToRNA.keySet()) {
				if (!foundActualDNAs.contains(s)) {
					System.out.println(datasetnames[i] + "\t" + s + "\tnot found.");
				}
			}

			System.out.println(datasetnames[i] + "\t" + matched + " matched out of " + dnaToRNA.size() + ", " + hascovars + " with covars");
		}

		tfout.close();

		System.out.println(total + "/" + totalcovars);


	}

	private HashMap<String, HashSet<String>> loadLinkFile(String linkfile) throws IOException {

		HashMap<String, HashSet<String>> dnaToRNA = new HashMap<String, HashSet<String>>();
		TextFile tf = new TextFile(linkfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);

		while (elems != null) {
			HashSet<String> list = dnaToRNA.get(elems[0]);
			if (list == null) {
				list = new HashSet<>();
			}
			list.add(elems[1]);
			dnaToRNA.put(elems[0], list);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		return dnaToRNA;

	}


}
