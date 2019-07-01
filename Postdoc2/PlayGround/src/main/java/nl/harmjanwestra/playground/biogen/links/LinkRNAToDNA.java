package nl.harmjanwestra.playground.biogen.links;


import org.apache.poi.ss.formula.functions.T;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class LinkRNAToDNA {

	public static void main(String[] args) {
		LinkRNAToDNA l = new LinkRNAToDNA();

		String origlink = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\2019-04-13-RNASeqIdToGenotypeID.txt";
		String individualfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\genotypeIndividuals\\2019-06-18-individualfiles.txt";

		String freeze1links = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\Freeze1Links.txt";

		String psychencodelinks = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\PsychEncodeIndividualToGT.txt";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksNoIntegrative\\links-";

		try {
			l.run(individualfile, origlink, freeze1links, psychencodelinks, output);
		} catch (IOException e) {
			e.printStackTrace();
		}

		String file = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\links-FoundInds.txt";
		String defupfileout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\links-FoundInds-dedup.txt";
		try {
			l.dedupRNA(file, defupfileout);


			String tissuemapfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-rnaqc\\2019-04-11\\2019-04-14-FullSampleAssignment.txt";
			String splittissuefolder = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\splitpertissue\\";
			l.splitBasedOnTissueType(defupfileout, tissuemapfile, splittissuefolder);

			String[] files = new String[]{
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\splitpertissue\\Cerebellum.txt",
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\splitpertissue\\Cortex.txt",
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\splitpertissue\\Spinal cord.txt",
			};
			for (String linkfile : files) {
				l.dedupDNA(linkfile);
			}
			files = new String[]{
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\splitpertissue\\Cerebellum.txt-dedup-gte.txt",
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\splitpertissue\\Cortex.txt-dedup-gte.txt",
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\splitpertissue\\Spinal cord.txt-dedup-gte.txt",
			};
			String popfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\populationassignment\\allGTSamplesPopulationAssignment.txt";
			for (String linkfile : files) {
				l.splitPerPopulation(linkfile, popfile);
			}


			String pathstrippeddedupfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\links-FoundInds-dedup-pathstripped.txt";
			files = new String[]{
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt",
					"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt"
			};

			for (String linkfile : files) {
				l.splitPerDataset(linkfile, pathstrippeddedupfile);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
//		String[] inds = new String[]{
//
//		}

//		try {
//			l.run(inds, origlink, output);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
	}

	private void splitPerDataset(String linkfile, String defupfileout) throws IOException {
		HashMap<String, String> samplePerDataset = new HashMap<>();
		HashMap<String, TextFile> dsout = new HashMap<>();

		TextFile tf = new TextFile(defupfileout, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);


		while (elems != null) {

			String sample = elems[0];
			String dataset = elems[4];
			TextFile tfo = dsout.get(dataset);
			if (tfo == null) {
				dsout.put(dataset, new TextFile(linkfile + "-" + dataset + ".txt", TextFile.W));

			}
			samplePerDataset.put(sample, dataset);

			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(linkfile, TextFile.R);

		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {

			String sample = elems[1];
			String ds = samplePerDataset.get(sample);
			if (ds == null) {
				System.out.println("Could not find ds for " + sample);

			} else {
				dsout.get(ds).writeln(elems[0] + "\t" + elems[1]);
			}

			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		for (String key : dsout.keySet()) {
			dsout.get(key).close();
		}

	}

	private void splitPerPopulation(String file, String popfile) throws IOException {
		HashMap<String, String> popPerSample = new HashMap<>();
		HashMap<String, TextFile> popout = new HashMap<>();
		TextFile tf = new TextFile(popfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			popPerSample.put(elems[0], elems[1]);
			TextFile tfout = popout.get(elems[1]);
			if (tfout == null) {
				popout.put(elems[1], new TextFile(file + "-" + elems[1] + ".txt", TextFile.W));
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(file, TextFile.R);

		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {

			String gt = elems[0];
			String pop = popPerSample.get(gt);
			if (gt.startsWith("SAMN") || gt.startsWith("SAMEA")) {
				pop = "EUR";
			}

			if (pop == null) {
				System.out.println("Could not find sample: " + gt);
			} else {
				popout.get(pop).writeln(elems[0] + "\t" + elems[1]);
			}

			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		for (String key : popout.keySet()) {
			popout.get(key).close();
		}
	}

	private void dedupDNA(String file) throws IOException {

		HashSet<String> visitedInd = new HashSet<>();
		TextFile tf = new TextFile(file, TextFile.R);
		TextFile tfo = new TextFile(file + "-dedup-gte.txt", TextFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (!visitedInd.contains(elems[1])) {
				tfo.writeln(elems[1] + "\t" + elems[0]);
				visitedInd.add(elems[1]);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfo.close();


	}


	private void splitBasedOnTissueType(String defupfileout, String tissuemapfile, String splittissuefolder) throws IOException {

		HashMap<String, String> tissuemap = new HashMap<>();
		HashMap<String, TextFile> tissueout = new HashMap<>();
		TextFile tf = new TextFile(tissuemapfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			String sample = elems[0];
			String tissue = elems[1];

			tissuemap.put(sample, tissue);
			TextFile tissuetf = tissueout.get(tissue);
			if (tissuetf == null) {
				tissueout.put(tissue, new TextFile(splittissuefolder + tissue + ".txt", TextFile.W));
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(defupfileout, TextFile.R);
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {

			String sample = elems[0];
			String tissue = tissuemap.get(sample);
			tissueout.get(tissue).writeln(sample + "\t" + elems[1]);
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		for (String key : tissueout.keySet()) {
			tissueout.get(key).close();
		}

	}

	private void dedupRNA(String file, String fileout) throws IOException {
		HashSet<String> visitedRNA = new HashSet<>();
		TextFile tf = new TextFile(file, TextFile.R);
		TextFile tfo = new TextFile(fileout, TextFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {

			if (!visitedRNA.contains(elems[0])) {
				tfo.writeln(Strings.concat(elems, Strings.tab));
				visitedRNA.add(elems[0]);
			}

			elems = tf.readLineElems(TextFile.tab);
		}
		tfo.close();
		tf.close();
	}

	public class Quad {
		String rna;
		String dna;
		String ds;
		String meta;
	}

	public HashMap<String, String> readAsMap(String file) throws IOException {
		HashMap<String, String> map = new HashMap<>();
		TextFile tfq = new TextFile(file, TextFile.R);
		String[] elems = tfq.readLineElems(TextFile.tab);
		while (elems != null) {
			if (!elems[1].equals("NotMapped")) {
				map.put(elems[0], elems[1]);
			}
			elems = tfq.readLineElems(TextFile.tab);

		}
		tfq.close();
		return map;
	}

	public void run(String indfile, String origlink, String freeze1links, String psychencodelinks, String output) throws IOException {

		ArrayList<Quad> data = new ArrayList<>();

		TextFile tf1 = new TextFile(origlink, TextFile.R);
		String[] elems1 = tf1.readLineElems(TextFile.tab);

		while (elems1 != null) {
			Quad q = new Quad();
			q.dna = elems1[1];
			q.rna = elems1[0];
			q.meta = elems1[2];
			q.ds = elems1[3];
			data.add(q);

			elems1 = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();

		HashMap<String, String> freeze1map = readAsMap(freeze1links);
		HashMap<String, String> psychmap = readAsMap(psychencodelinks);

		TextFile tf = new TextFile(indfile, TextFile.R);

		String ln = tf.readLine();
		HashSet<String> foundinds = new HashSet<>();
		TextFile allindOut = new TextFile(output + "AllInds.txt", TextFile.W);
		TextFile foundindOut = new TextFile(output + "FoundInds.txt", TextFile.W);
		while (ln != null) {

			System.out.println(ln);
			//
			HashSet<String> inds = new HashSet<>();
			TextFile tf3 = new TextFile(ln, TextFile.R);
			String ln2 = tf3.readLine();
			while (ln2 != null) {
				inds.add(ln2);
				allindOut.writeln(ln2 + "\t" + ln);
				ln2 = tf3.readLine();
			}
			tf3.close();


			boolean found = false;
			for (Quad q : data) {
				found = false;
				String rna = q.rna;
				String f1 = freeze1map.get(rna);

				String dna = q.dna;
				if (inds.contains(dna)) {
					foundinds.add(f1);
					foundinds.add(dna);
					foundindOut.writeln(rna + "\t" + dna + "\t" + q.ds + "\t" + q.meta + "\t" + ln);
					found = true;
				}
				if (!found && f1 != null) {

					// try freeze 1 mapping
					if (inds.contains(f1)) {
						foundinds.add(f1);
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + ln);
						found = true;
					}
				}

				// targetals
				if (!found) {
					f1 = dna + "-b38";
					if (inds.contains(f1)) {
						foundinds.add(f1);
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + ln);
						found = true;
					}
				}

				// NABEC
				if (!found) {
					f1 = "0_" + dna;
					if (inds.contains(f1)) {
						foundinds.add(f1);
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + ln);
						found = true;
					}
				}

				if (!found) {
					if (dna.contains("UMARY")) {
//						System.out.println("Try this!");
					}
					f1 = "0_" + dna.replaceAll("-", "_");
					if (inds.contains(f1)) {
						foundinds.add(f1);
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + ln);
						found = true;
					}
				}

				// HBCC
				if (!found) {

					for (String ind : inds) {
						String[] indelems = ind.split("_");
						if (indelems.length == 3) {
							String actual = indelems[1] + "_" + indelems[2];
							if (actual.equals(dna)) {

								foundinds.add(f1);
								foundinds.add(dna);
								foundindOut.writeln(rna + "\t" + ind + "\t" + q.ds + "\t" + q.meta + "\t" + ln);
								found = true;
								break;
							}
						}
					}
				}
				// braineac
				if (!found) {
					f1 = "0_" + dna + "_1";
					if (inds.contains(f1)) {
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + ln);
						found = true;
					}
				}

				// NABEC-H610
				if (!found) {
					// SH-07-73 --> 0_SH-07-73_SH-07-73
					f1 = "0_" + dna + "_" + dna;
					if (inds.contains(f1)) {
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + ln);
						found = true;
					}
				}

				if (!found) {
					if (rna.equals("Br1157_R2836") && ln.contains("CMC_HBCC_set1")) {
						System.out.println("Found it");
					}
					String psychencodeind = psychmap.get(dna);
					if (psychencodeind != null) {
						for (String ind : inds) {
							String[] indelems = ind.split("_");
							if (indelems.length == 3) {
								String actual = indelems[1] + "_" + indelems[2];
								if (actual.equals(psychencodeind)) {
									foundinds.add(psychencodeind);
									foundinds.add(dna);
									foundindOut.writeln(rna + "\t" + ind + "\t" + q.ds + "\t" + q.meta + "\t" + ln);
									found = true;
									break;
								}
							}
						}

					}
				}

			}

			System.out.println(foundinds.size() + " ids found sofar..");
			ln = tf.readLine();
		}
		tf.close();
		allindOut.close();
		foundindOut.close();

		TextFile notfoundOut = new TextFile(output + "NotFoundInds.txt", TextFile.W);
		int nrnoutfound = 0;
		for (Quad d : data) {
			if (!foundinds.contains(d.dna)) {
				notfoundOut.writeln(d.rna + "\t" + d.dna + "\t" + d.ds + "\t" + d.meta);
				nrnoutfound++;
			}
		}
		notfoundOut.close();
		System.out.println(nrnoutfound + " not found out of " + data.size());


	}

	public void run(String[] individualFiles, String origlink, String output) throws IOException {
		HashMap<String, String> rnaToGenotype = new HashMap<>();
		TextFile tf = new TextFile(origlink, TextFile.R);
		tf.readLine();

		HashSet<String> visitedRNA = new HashSet<>();

		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 2) {
				String id = elems[0] + "_" + elems[1];
				if (visitedRNA.contains(elems[1])) {
					System.out.println("Dup RNA: " + elems[1]);
				}
				visitedRNA.add(elems[1]);
				String gt = elems[2];
				if (rnaToGenotype.containsKey(id)) {
					System.out.println("Dup: " + id + "\t" + elems[0] + "\t" + elems[1]);
				} else {
					rnaToGenotype.put(id, gt);
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


	}
}
