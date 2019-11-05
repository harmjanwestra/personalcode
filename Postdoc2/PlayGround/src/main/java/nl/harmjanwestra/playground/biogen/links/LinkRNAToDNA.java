package nl.harmjanwestra.playground.biogen.links;


import it.unimi.dsi.fastutil.Hash;
import org.apache.poi.ss.formula.functions.T;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
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
		String individualfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\2019-07-01-individualfiles.txt";

		String freeze1links = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\Freeze1Links.txt";

		String psychencodelinks = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-14-linksIntegrativeFirst\\PsychEncodeIndividualToGT.txt";
		String outdir = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-01-linksNoIntegrative\\";
		String output = outdir + "links-";

		String rnaseqoutliers = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-rnaqc\\2019-04-11\\2019-04-11-AllRNASeqOutliers.txt";

		try {
			l.run(individualfile, origlink, freeze1links, psychencodelinks, rnaseqoutliers, output);
		} catch (IOException e) {
			e.printStackTrace();
		}
//		System.exit(0);

		boolean preferresequencedsamples = true;
		String file = outdir + "links-FoundInds.txt";
		String defupfileout = outdir + "links-FoundInds-dedup.txt";
		try {
			// remove duplicate genotype samples
			l.dedupRNA(file, defupfileout);
			String pathstrippeddedupfile = outdir + "links-FoundInds-dedup-pathstripped.txt";
			l.stripPaths(defupfileout, pathstrippeddedupfile);

//			System.exit(-1);

			String tissuemapfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-12-rnaqc\\2019-04-11\\2019-04-14-FullSampleAssignment.txt";
			String splittissuefolder = outdir + "splitpertissue\\";
			l.splitBasedOnTissueType(pathstrippeddedupfile, tissuemapfile, splittissuefolder);

			String[] files = new String[]{
					outdir + "splitpertissue\\Cerebellum.txt",
					outdir + "splitpertissue\\Cortex.txt",
					outdir + "splitpertissue\\Spinal cord.txt",
			};

			for (String linkfile : files) {
				l.dedupDNA(linkfile, preferresequencedsamples);
			}

			files = new String[]{
					outdir + "splitpertissue\\Cerebellum.txt-dedup-gte.txt",
					outdir + "splitpertissue\\Cortex.txt-dedup-gte.txt",
					outdir + "splitpertissue\\Spinal cord.txt-dedup-gte.txt",
			};

			String popfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-04-13-genotypeqc\\populationassignment\\allGTSamplesPopulationAssignment.txt";
			for (String linkfile : files) {
				l.splitPerPopulation(linkfile, popfile);
			}

			files = new String[]{
//					outdir + "splitpertissue\\Cortex.txt-dedup-gte.txt-EUR.txt",
//					outdir + "splitpertissue\\Cortex.txt-dedup-gte.txt-AFR.txt"
					outdir + "splitpertissue\\Cerebellum.txt-dedup-gte.txt-EUR.txt",
					outdir + "splitpertissue\\Cerebellum.txt-dedup-gte.txt-AFR.txt"
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

	private void removeDuplicateAndRelatedDNA(String defupfileout, String duplicateDNAIds, String pathstrippeddedupfilewithoutrelated) throws IOException {
		HashSet<String> idsToExclude = new HashSet<>();
		TextFile tf = new TextFile(duplicateDNAIds, TextFile.R);
		idsToExclude.addAll(tf.readAsArrayList());
		tf.close();
		TextFile tf2 = new TextFile(defupfileout, TextFile.R);
		TextFile outf = new TextFile(pathstrippeddedupfilewithoutrelated, TextFile.W);

		String[] elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {

			String dna = elems[1];
			if (!idsToExclude.contains(dna)) {
				outf.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = tf2.readLineElems(TextFile.tab);
		}

		outf.close();
		tf2.close();

	}

	private void stripPaths(String defupfileout, String pathstrippeddedupfile) throws IOException {

		TextFile in = new TextFile(defupfileout, TextFile.R);
		TextFile out = new TextFile(pathstrippeddedupfile, TextFile.W);

		String[] elems = in.readLineElems(TextFile.tab);
		while (elems != null) {

			String path = elems[4];
			String[] pathelems = path.split("\\\\");
			String name = pathelems[pathelems.length - 1];
			name = name.replace("-Individuals.txt", "");
			elems[4] = name;

			out.writeln(Strings.concat(elems, Strings.tab));
			elems = in.readLineElems(TextFile.tab);
		}
		in.close();
		out.close();
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
		TextFile sampleToDs = new TextFile(linkfile + "-SampleToDataset.txt", TextFile.W);
		while (elems != null) {

			String sample = elems[1];
			String ds = samplePerDataset.get(sample);
			sampleToDs.writeln(sample + "\t" + ds);
			if (ds == null) {
				System.out.println("Could not find ds for " + sample);

			} else {
				dsout.get(ds).writeln(elems[0] + "\t" + elems[1]);
			}

			elems = tf2.readLineElems(TextFile.tab);
		}
		sampleToDs.close();
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

	private void dedupDNA(String file, boolean preferreseqenced) throws IOException {

		if (preferreseqenced) {


			HashMap<String, HashSet<String>> links = new HashMap<>();
			TextFile tf = new TextFile(file, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				HashSet<String> rnas = links.get(elems[1]);
				if (rnas == null) {
					rnas = new HashSet<>();
				}
				rnas.add(elems[0]);
				links.put(elems[1], rnas);

				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();

			TextFile tfo = new TextFile(file + "-dedup-gte.txt", TextFile.W);
			for (String key : links.keySet()) {
				HashSet<String> rnas = links.get(key);
				String reseqsample = null;
				String othersample = null;
				for (String s : rnas) {
					if (s.contains("reseq")) {
						reseqsample = s;
					} else {
						othersample = s;
					}
				}

				if (reseqsample == null) {
					// pick a random sample
					tfo.writeln(key + "\t" + othersample);
				} else {
					tfo.writeln(key + "\t" + reseqsample);
				}

			}
			tfo.close();

		} else {
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


	}


	private void splitBasedOnTissueType(String defupfileout, String tissuemapfile, String splittissuefolder) throws IOException {
		Gpio.createDir(splittissuefolder);
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

	public void run(String indfile, String origlink, String freeze1links, String psychencodelinks, String rnseqoutlierfile, String output) throws IOException {

		HashSet<String> rnaseqoutliers = null;
		if (rnseqoutlierfile != null) {
			rnaseqoutliers = new HashSet<>();
			TextFile tf = new TextFile(rnseqoutlierfile, TextFile.R);
			rnaseqoutliers.addAll(tf.readAsArrayList());
			tf.close();
		}

		ArrayList<Quad> data = new ArrayList<>();

		TextFile tf1 = new TextFile(origlink, TextFile.R);
		String[] elems1 = tf1.readLineElems(TextFile.tab);

		while (elems1 != null) {
			Quad q = new Quad();
			q.dna = elems1[1];
			q.rna = elems1[0];
			q.meta = elems1[2];
			q.ds = elems1[3];
			if (rnaseqoutliers == null || !rnaseqoutliers.contains(q.rna)) {
				data.add(q);
			}

			elems1 = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();

		HashMap<String, String> freeze1map = readAsMap(freeze1links);
		HashMap<String, String> psychmap = readAsMap(psychencodelinks);

		TextFile tf = new TextFile(indfile, TextFile.R);

		String[] lnelems = tf.readLineElems(TextFile.tab);
		HashSet<String> foundinds = new HashSet<>();
		TextFile allindOut = new TextFile(output + "AllInds.txt", TextFile.W);
		TextFile foundindOut = new TextFile(output + "FoundInds.txt", TextFile.W);

		while (lnelems != null) {

			System.out.println(lnelems[0]);
			HashSet<String> excludeInds = new HashSet<String>();
			if (lnelems.length > 1) {
				System.out.println(lnelems[1]);
				TextFile tf2 = new TextFile(lnelems[1], TextFile.R);
				String[] dupelems = tf2.readLineElems(TextFile.tab);
				while (dupelems != null) {
					excludeInds.add(dupelems[3]);
					dupelems = tf2.readLineElems(TextFile.tab);
				}
				tf2.close();
			}


			//
			HashSet<String> inds = new HashSet<>();
			String individualFile = lnelems[0];
			TextFile tf3 = new TextFile(individualFile, TextFile.R);
			String ln2 = tf3.readLine();
			int included = 0;
			int excluded = 0;
			while (ln2 != null) {
				if (excludeInds.contains(ln2)) {
					excluded++;
				} else {
					inds.add(ln2);
					included++;
					allindOut.writeln(ln2 + "\t" + lnelems[0]);
				}

				ln2 = tf3.readLine();
			}
			tf3.close();

			System.out.println("File has: " + included + " included, " + excluded + " excluded.");

			boolean found = false;
			for (Quad q : data) {
				found = false;
				String rna = q.rna;
				String f1 = freeze1map.get(rna);

				String dna = q.dna;
				if (inds.contains(dna)) {
					foundinds.add(f1);
					foundinds.add(dna);
					foundindOut.writeln(rna + "\t" + dna + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
					found = true;
				}
				if (!found && f1 != null) {

					// try freeze 1 mapping
					if (inds.contains(f1)) {
						foundinds.add(f1);
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
						found = true;
					}
				}

				// targetals
				if (!found) {
					f1 = dna + "-b38";
					if (inds.contains(f1)) {
						foundinds.add(f1);
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
						found = true;
					}
				}

				// NABEC
				if (!found) {
					f1 = "0_" + dna;
					if (inds.contains(f1)) {
						foundinds.add(f1);
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
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
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
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
								foundindOut.writeln(rna + "\t" + ind + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
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
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
						found = true;
					}
				}

				// NABEC-H610
				if (!found) {
					// SH-07-73 --> 0_SH-07-73_SH-07-73
					f1 = "0_" + dna + "_" + dna;
					if (inds.contains(f1)) {
						foundinds.add(dna);
						foundindOut.writeln(rna + "\t" + f1 + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
						found = true;
					}
				}

				if (!found) {
					if (rna.equals("Br1157_R2836") && individualFile.contains("CMC_HBCC_set1")) {
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
									foundindOut.writeln(rna + "\t" + ind + "\t" + q.ds + "\t" + q.meta + "\t" + individualFile);
									found = true;
									break;
								}
							}
						}

					}
				}

			}

			System.out.println(foundinds.size() + " ids found sofar..");
			lnelems = tf.readLineElems(TextFile.tab);
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
