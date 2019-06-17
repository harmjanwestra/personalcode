package nl.harmjanwestra.playground.methylation;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixConverter;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetAppendableWriter;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class Process450KData {
	public static void main(String[] args) {
		Process450KData m = new Process450KData();
//
//		String[] exclusion = new String[]{
//				"Z:\\projects\\2018-methylation\\code\\450K_DataProcessing\\ADDITIONAL_INFO\\ProbeFiltering\\freq5percent\\probeToFilter_450K_1000G_omni2.5.hg19.EUR_alleleFreq5percent_50bp_wInterroSite.txt",
//				"Z:\\projects\\2018-methylation\\code\\450K_DataProcessing\\ADDITIONAL_INFO\\ProbeFiltering\\ProbesBindingNonOptimal\\Source&BSProbesMappingMultipleTimesOrNotBothToBSandNormalGenome.txt"
//		};
//
//		String infolder = "Z:\\projects\\2018-methylation\\GPL13534_450k_OUT\\";
//		String probelistout = "Z:\\projects\\2018-methylation\\probelists\\GPL13534_450k";
//		try {
//			m.makeprobelistPerType(infolder, probelistout, exclusion);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		try {
			if (args.length == 0) {
				printUsage();
			} else if (args[0].equals("makeprobelist")) {
				if (args.length < 4) {
					System.out.println("Usage: makeprobelistfile infolder probeoutdir exclusionlist");
					System.out.println("Exclusion list comma separated list of files");
				} else {
					m.makeprobelistPerType(args[1], args[2], args[3].split(";"));
				}
			} else if (args[0].equals("mergetables")) {
				if (args.length < 5) {
					System.out.println("Usage: mergetables probelist1 probelist2 infolder outfolder");
				} else {
					m.runappendIndependentMatrices(args[1], args[2], args[3], args[4]);
				}
			} else if (args[0].equals("qqnorm")) {
				if (args.length < 3) {
					System.out.println("Usage: qqnorm input output");
				} else {
					QuantileNormalizeDiskBased450K q = new QuantileNormalizeDiskBased450K();
					q.run(args[1], args[2]);
				}
			} else if (args[0].equals("logtransform")) {
				if (args.length < 3) {
					System.out.println("Usage: logtransform input output");
				} else {
					Log2TransformDiskBased450K q = new Log2TransformDiskBased450K();
					q.run(args[1], args[2]);

				}
			} else if (args[0].equals("countnans")) {
				if (args.length < 3) {
					System.out.println("Usage: countnans input output");
				} else {
					DetermineNAs450K d = new DetermineNAs450K();
					d.run(args[1], args[2]);
				}
			} else if (args[0].equals("calcmvals")) {

				if (args.length < 4) {
					System.out.println("Usage: calcmvals inU inM out");
				} else {
					CalculateMValuesDiskBased450K d = new CalculateMValuesDiskBased450K();
					d.run(args[1], args[2], args[3]);
				}
			} else if (args[0].equals("mergebinarymatrices")) {
				Merge450KDiskBased d = new Merge450KDiskBased();
				if (args.length < 4) {
					System.out.println("Usage: inT1 inT2 out [probelist]");
				} else {
					String pblist = null;
					if (args.length >= 5) {
						pblist = args[4];
					}
					d.run(pblist, args[1], args[2], args[3]);
				}
			} else if (args[0].equals("correl")) {
				CorrelateDiskBased450K d = new CorrelateDiskBased450K();
				if (args.length < 4) {
					System.out.println("Usage: correl in out nrrowsperblock");
				} else {
					d.runblocks(args[1], args[2], Integer.parseInt(args[3]));
				}
			} else if (args[0].equals("detriangle")) {
				if (args.length < 3) {
					System.out.println("Usage: input output");
				} else {
					DeTriangularize d = new DeTriangularize();
					d.run(args[1], args[2]);
				}
			} else if (args[0].equals("pca")) {
				if (args.length < 3) {
					System.out.println("Usage: in nrpcs");
				} else {
					PCA p = new PCA();
					p.run(args[1], Integer.parseInt(args[2]));
				}
			} else if (args[0].equals("filter")) {
				if (args.length < 7) {
					System.out.println("Usage: input output (rowset|null) (include|exclude) (colset|null) (include|exclude)");
				} else {
					DatasetFilter filter = new DatasetFilter();
					filter.run(args[1], args[2], args[3], args[4], args[5], args[6]);
				}
			} else if (args[0].equals("centerscale")) {
				CenterScale450K c = new CenterScale450K();
				if (args.length < 3) {
					System.out.println("Usage; input output");
				} else {
					c.run(args[1], args[2]);
				}
			} else if (args[0].equals("txttobin")) {
				if (args.length < 3) {
					System.out.println("Usage: input output");
				} else {
					DoubleMatrixConverter.TextToBinary(args[1], args[2]);
				}
			} else if (args[0].equals("bintotxt")) {
				if (args.length < 3) {
					System.out.println("Usage: input output");
				} else {
					DoubleMatrixConverter.BinaryToText(args[1], args[2]);
				}
			} else {
				printUsage();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private static void printUsage() {
		System.out.println("Usage:\n" +
				"txttobin\tConvert text to binary matrix\n" +
				"bintotxt\tConvert binary to text matrix\n\n" +
				"makeprobelist\tlist probes in output (probes expected on rows, samples on cols)\n" +
				"mergetables\tmerge 450K pipeline folder output\n" +
				"calcmvals\tcalculate m-values\n" +
				"qqnorm\tquantile normalize (expect samples on rows)\n" +
				"logtransform\tlog2transform (expect samples on rows)\n" +
				"correl\tcorrelations (expect samples on rows); outputs upper triangle correlation matrix over rows\n" +
				"countnans\tcount the number of nans in a matrix\n" +
				"mergebinarymatrices\tmerge binary matrices\n" +
				"detriangle\tfill in the lower triangle of the matrix\n" +
				"pca\tcalculate eigenvectors\n" +
				"centerscale\tcenter and scale data (expect samples on rows)\n" +
				"filter\tfilter rows or columns");

	}

	private void makeprobelist(String infolder, String probelistfile, String[] exclusion) throws IOException {
		HashSet<String> exclusionset = new HashSet<>();
		for (String s : exclusion) {
			TextFile tf = new TextFile(s, TextFile.R);
			ArrayList<String> list = tf.readAsArrayList();
			exclusionset.addAll(list);
			tf.close();
		}
		System.out.println(exclusionset.size() + " probes to exclude");

		HashSet<String> probes = new HashSet<>();
		Files.walk(Paths.get(infolder))
				.filter(Files::isDirectory)
				.forEach(v -> {
							if (probes.size() < 422801) {
								File[] list = v.toFile().listFiles();
								for (File f : list) {
									try {
										if (f.getName().endsWith("T1_Grn_M_Signal.txt.gz")) {
											addprobes(f, probes, exclusionset);
										} else if (f.getName().endsWith("T1_Red_M_Signal.txt.gz")) {
											addprobes(f, probes, exclusionset);
										} else if (f.getName().endsWith("T2_M_Signal.txt.gz")) {
											addprobes(f, probes, exclusionset);
										} else if (f.getName().endsWith("T2_U_Signal.txt.gz")) {
											addprobes(f, probes, exclusionset);
										}
									} catch (IOException e) {
										e.printStackTrace();
									}
								}
							}
						}
				);
		System.out.println(probes.size() + " probes total.");
		TextFile out = new TextFile(probelistfile, TextFile.W);
		for (String s : probes) {
			out.writeln(s);
		}
		out.close();
	}

	private void addprobes(File v, HashSet<String> probes, HashSet<String> exclusionset) throws IOException {
		System.out.println(v.getName() + "\t" + probes.size() + " probes sofar");
		TextFile tf = new TextFile(v, TextFile.R);
		tf.readLine();
		String ln = tf.readLine();

		while (ln != null) {
			String[] sub = Strings.subsplit(ln, Strings.tab, 0, 1);
			if (!exclusionset.contains(sub[0])) {
				probes.add(sub[0]);
			}
			ln = tf.readLine();
		}
		tf.close();
		System.out.println(v.getName() + " done. \t" + probes.size() + " probes sofar");
	}

	public void runappendIndependentMatrices(String probelistfileT1, String probelistfileT2, String
			inFolder, String outfolder) throws IOException {

		Pair<HashMap<String, Integer>, ArrayList<String>> probeT1 = loadAndIndexProbelist(probelistfileT1);
		Pair<HashMap<String, Integer>, ArrayList<String>> probeT2 = loadAndIndexProbelist(probelistfileT2);


		DoubleMatrixDatasetAppendableWriter outT1U = new DoubleMatrixDatasetAppendableWriter(probeT1.getRight(), outfolder + "T1U");
		DoubleMatrixDatasetAppendableWriter outT1M = new DoubleMatrixDatasetAppendableWriter(probeT1.getRight(), outfolder + "T1M");
		DoubleMatrixDatasetAppendableWriter outT2U = new DoubleMatrixDatasetAppendableWriter(probeT2.getRight(), outfolder + "T2U");
		DoubleMatrixDatasetAppendableWriter outT2M = new DoubleMatrixDatasetAppendableWriter(probeT2.getRight(), outfolder + "T2M");

		File[] list = new File(inFolder).listFiles();
		int ctr = 0;
		for (File f : list) {
			if (f.isDirectory()) {
				try {
					System.out.println("Folder item: " + ctr + " / " + list.length);
					appendGSEProjectToIndependentMatrices(f, outT1M, outT1U, outT2M, outT2U, probeT1.getLeft(), probeT2.getLeft());
					System.gc();

				} catch (IOException e) {
					e.printStackTrace();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			ctr++;
		}

		outT1U.close();
		outT1M.close();
		outT2M.close();
		outT2U.close();
	}

	private Pair<HashMap<String, Integer>, ArrayList<String>> loadAndIndexProbelist(String probelistfileT1) throws
			IOException {
		TextFile lf = new TextFile(probelistfileT1, TextFile.R);
		ArrayList<String> probelist = lf.readAsArrayList();
		Collections.sort(probelist);
		lf.close();

		HashMap<String, Integer> probeHashs = new HashMap<String, Integer>();
		int c = 0;
		for (String s : probelist) {
			probeHashs.put(s, c);
			c++;
		}
		return new Pair<>(probeHashs, probelist);
	}

	public void runappendGSEProject(String probelistfile, String inFolder, String outfile) throws IOException {

		Pair<HashMap<String, Integer>, ArrayList<String>> probes = loadAndIndexProbelist(probelistfile);

		TextFile out = new TextFile(outfile, TextFile.W);
		String header = "-\t" + Strings.concat(probes.getRight(), Strings.tab);
		out.writeln(header);
		File[] list = new File(inFolder).listFiles();
		for (File f : list) {
			if (f.isDirectory()) {
				try {
					appendGSEProject(f, out, probes.getLeft());
					System.gc();
				} catch (IOException e) {
					e.printStackTrace();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
//		Files.walk(Paths.get(inFolder))
//				.filter(Files::isDirectory)
//				.forEach(v -> {
//
//				});
//
		out.close();
	}


	public static String T1PGRNM = "T1_Grn_M_Signal.txt.gz";
	public static String T1PGRNU = "T1_Grn_U_Signal.txt.gz";
	public static String T1PREDM = "T1_Red_M_Signal.txt.gz";
	public static String T1PREDU = "T1_Red_U_Signal.txt.gz";
	public static String T2PM = "T2_M_Signal.txt.gz";
	public static String T2PU = "T2_U_Signal.txt.gz";

	public static String T1M = "T1M";
	public static String T1U = "T1U";
	public static String T2M = "T2M";
	public static String T2U = "T2U";


	public void appendGSEProject(File v, TextFile out, HashMap<String, Integer> probeHashs) throws Exception {

		File[] list = v.listFiles();
		if (list.length > 5) {
			DoubleMatrixDataset<String, String> m1red = null;
			DoubleMatrixDataset<String, String> m1grn = null;
			DoubleMatrixDataset<String, String> u1red = null;
			DoubleMatrixDataset<String, String> u1grn = null;
			DoubleMatrixDataset<String, String> m2 = null;
			DoubleMatrixDataset<String, String> u2 = null;

			// GSE100134_T1_Grn_M_Signal.txt.gz  GSE100134_T1_Grn_U_Signal.txt.gz  GSE100134_T1_Red_M_Signal.txt.gz  GSE100134_T1_Red_U_Signal.txt.gz  GSE100134_T2_M_Signal.txt.gz  GSE100134_T2_U_Signal.txt.gz

			String gse = "";

			HashSet<String> uniqueSamples = new HashSet<String>();
			for (File f : list) {
				if (f.getName().endsWith(T1PGRNM)) {
					m1grn = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(m1grn.getColObjects());
					gse = f.getName().replaceAll("_" + T1PGRNM, "");
				} else if (f.getName().endsWith(T1PGRNU)) {
					u1grn = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(u1grn.getColObjects());
				} else if (f.getName().endsWith(T1PREDM)) {
					m1red = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(m1red.getColObjects());
				} else if (f.getName().endsWith(T1PREDU)) {
					u1red = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(u1red.getColObjects());
				} else if (f.getName().endsWith(T2PM)) {
					m2 = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(m2.getColObjects());
				} else if (f.getName().endsWith(T2PU)) {
					u2 = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(u2.getColObjects());
				}

			}

			// get a list of samples
			ArrayList<String> allsamples = new ArrayList<>();
			allsamples.addAll(uniqueSamples);
			Collections.sort(allsamples);
			System.out.println(uniqueSamples.size() + " samples.");

			System.out.println("Indexing GRN");

			int[] mtougrn = new int[u1grn.getRowObjects().size()];
			int[] grn1ind = new int[probeHashs.size()];
			{
				System.out.println("starting index");
				ArrayList<String> rowobs = m1grn.getRowObjects();
				LinkedHashMap<String, Integer> rowobs2 = u1grn.getHashRows();
				IntStream.range(0, rowobs.size()).parallel().forEach(i -> {
							String probe = rowobs.get(i);
							Integer index = probeHashs.get(probe);
							if (index != null) {
								Integer uind = rowobs2.get(probe);
								mtougrn[i] = uind;
								grn1ind[i] = index;
							} else {
								grn1ind[i] = -1;
							}
						}
				);
			}

			System.out.println("Indexing RED");
			int[] mtoured = new int[m1red.getRowObjects().size()];
			int[] red1ind = new int[probeHashs.size()];
			{
				ArrayList<String> rowobs = m1red.getRowObjects();
				LinkedHashMap<String, Integer> rowobs2 = u1red.getHashRows();
				IntStream.range(0, rowobs.size()).parallel().forEach(i -> {
					String probe = rowobs.get(i);
					Integer index = probeHashs.get(probe);
					if (index != null) {
						Integer uind = rowobs2.get(probe);
						mtoured[i] = uind;
						red1ind[i] = index;
					} else {
						red1ind[i] = -1;
					}
				});
			}

			System.out.println("Indexing T2");
			int[] m2tou = new int[m2.getRowObjects().size()];
			int[] m2ind = new int[probeHashs.size()];
			{
				ArrayList<String> rowobs = m2.getRowObjects();
				LinkedHashMap<String, Integer> rowobs2 = u2.getHashRows();
				IntStream.range(0, rowobs.size()).parallel().forEach(i -> {
					String probe = rowobs.get(i);
					Integer index = probeHashs.get(probe);
					if (index != null) {
						Integer uind = rowobs2.get(probe);

						m2tou[i] = uind;
						m2ind[i] = index;
					} else {
						m2ind[i] = -1;
					}
				});

			}


			DoubleMatrixDataset<String, String> finalM1grn = m1grn;
			DoubleMatrixDataset<String, String> finalU1grn = u1grn;
			DoubleMatrixDataset<String, String> finalM1red = m1red;
			DoubleMatrixDataset<String, String> finalU1red = u1red;
			DoubleMatrixDataset<String, String> finalM = m2;
			DoubleMatrixDataset<String, String> finalU = u2;
			AtomicInteger ctr = new AtomicInteger();
			String finalGse = gse;
			ProgressBar pb = new ProgressBar(allsamples.size());
			allsamples.stream().parallel().forEach(sample ->
					{
						try {
							Integer vm1grn = finalM1grn.getHashCols().get(sample);
							Integer vu1grn = finalU1grn.getHashCols().get(sample);
							Integer vm1red = finalM1red.getHashCols().get(sample);
							Integer vu1red = finalU1red.getHashCols().get(sample);
							Integer vm2 = finalM.getHashCols().get(sample);
							Integer vu2 = finalU.getHashCols().get(sample);

//							System.out.println("Sample: " + sample);
							// now iterate probes
//							System.out.println(sample + "\tinit");
							double[] output = new double[probeHashs.size()];
							Arrays.fill(output, Double.NaN);
//							System.out.println(sample + "\tinit done");
							int nrt1grn = finalM1grn.getRowObjects().size();
							for (int i = 0; i < nrt1grn; i++) {
								int index = grn1ind[i];
								if (index != -1) {
									int uind = mtougrn[i];
									double mval = Math.log(finalM1grn.getElementQuick(i, vm1grn) / finalU1grn.getElementQuick(uind, vu1grn));
									output[index] = mval;
								}
							}

//							System.out.println(sample + "\tdone init m1 grn");
							int nrt1red = finalM1red.getRowObjects().size();
							for (int i = 0; i < nrt1red; i++) {
								int index = red1ind[i];
								if (index != -1) {
									int uind = mtoured[i];
									double mval = Math.log(finalM1red.getElementQuick(i, vm1red) / finalU1red.getElementQuick(uind, vu1red));
									output[index] = mval;
								}
							}
//							System.out.println(sample + "\tdone init m1 red");
							int nrt2 = finalM.getRowObjects().size();
							for (int i = 0; i < nrt2; i++) {
								int index = m2ind[i];
								if (index != -1) {
									int uind = m2tou[i];
									double mval = Math.log(finalM.getElementQuick(i, vm2) / finalU.getElementQuick(uind, vu2));
									output[index] = mval;
								}
							}

//							System.out.println(sample + "\tdone init m2");
							String mergedsample = finalGse + "_" + sample;
							out.writelnsynced(mergedsample + "\t" + Strings.concat(output, Strings.tab));


							pb.set(ctr.getAndIncrement());
						} catch (NoSuchElementException e) {
							System.out.println("Sample " + sample + " not found for all probe types");
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
			);
			pb.close();
		}
	}


	private void makeprobelistPerType(String infolder, String probelistfile, String[] exclusion) throws IOException {
		HashSet<String> exclusionset = new HashSet<>();
		for (String s : exclusion) {
			TextFile tf = new TextFile(s, TextFile.R);
			ArrayList<String> list = tf.readAsArrayList();
			exclusionset.addAll(list);
			tf.close();
		}
		System.out.println(exclusionset.size() + " probes to exclude");

		HashMap<String, HashSet<String>> probesPerDatatype = new HashMap<>();

		Files.walk(Paths.get(infolder))
				.filter(Files::isDirectory)
				.parallel()
				.forEach(v -> {

							File[] list = v.toFile().listFiles();
							for (File f : list) {
								try {
									String probeset = null;
									if (f.getName().endsWith(T1PGRNM) || f.getName().endsWith(T1PREDM)) {
										probeset = T1M;
									} else if (f.getName().endsWith(T1PREDU) || f.getName().endsWith(T1PGRNU)) {
										probeset = T1U;
									} else if (f.getName().endsWith(T2PM)) {
										probeset = T2M;
									} else if (f.getName().endsWith(T2PU)) {
										probeset = T2U;
									}

									if (probeset != null) {
										HashSet<String> probes = new HashSet<>();

										System.out.println(f.getName() + "\t" + probes.size() + " probes sofar");
										TextFile tf = new TextFile(f, TextFile.R);
										tf.readLine();
										String ln = tf.readLine();

										while (ln != null) {
											String[] sub = Strings.subsplit(ln, Strings.tab, 0, 1);
											if (!exclusionset.contains(sub[0])) {
												probes.add(sub[0]);
											}
											ln = tf.readLine();
										}
										tf.close();
										System.out.println(f.getName() + " done. \t" + probes.size() + " probes sofar");
										synchronized (probesPerDatatype) {
											HashSet<String> tmpprobes = probesPerDatatype.get(probeset);
											if (tmpprobes == null) {
												probesPerDatatype.put(probeset, probes);
											} else {
												tmpprobes.addAll(probes);
												probesPerDatatype.put(probeset, tmpprobes);
											}
										}
									}
								} catch (IOException e) {
									e.printStackTrace();
								}
							}

						}
				);

		for (String key : probesPerDatatype.keySet()) {
			HashSet<String> probes = probesPerDatatype.get(key);
			System.out.println(probes.size() + " probes total.");
			TextFile out = new TextFile(probelistfile + "-" + key + ".txt.gz", TextFile.W);
			for (String s : probes) {
				out.writeln(s);
			}
			out.close();
		}

	}


	private class FolderObj {
		File m1redfile = null;
		File m1grnfile = null;

		File u1redfile = null;
		File u1grnfile = null;

		File u2file = null;
		File m2file = null;

		String gse = null;

		public FolderObj(File v) {
			System.out.println("Scanning folder: " + v.getName());
			File[] list = v.listFiles();
			IntStream.range(0, list.length).parallel().forEach(fi -> {
				File f = list[fi];
				if (f.getName().endsWith(T1PGRNM)) {
					m1grnfile = f;
					gse = f.getName().replaceAll("_" + T1PGRNM, "");
				} else if (f.getName().endsWith(T1PGRNU)) {
					u1grnfile = f;
				} else if (f.getName().endsWith(T1PREDM)) {
					m1redfile = f;
				} else if (f.getName().endsWith(T1PREDU)) {
					u1redfile = f;
				} else if (f.getName().endsWith(T2PM)) {
					m2file = f;
				} else if (f.getName().endsWith(T2PU)) {
					u2file = f;
				}
			});
		}

		public boolean allfound() {
			return (m1redfile != null &&
					m1grnfile != null &&
					u1redfile != null &&
					u1grnfile != null &&
					m2file != null &&
					u2file != null);
		}

		public ArrayList<String> getUniqueSamples() throws IOException {

			HashSet<String> samples = new HashSet<>();
			if (allfound()) {
				TextFile tf = new TextFile(m1grnfile, TextFile.R);
				String[] elems = tf.readLineElems(TextFile.tab);
				for (int i = 1; i < elems.length; i++) {
					samples.add(elems[i]);
				}
				tf.close();
			}
			ArrayList<String> sampleList = new ArrayList<>();
			sampleList.addAll(samples);
			Collections.sort(sampleList);
			System.out.println(sampleList.size() + " samples.");
			return sampleList;
		}

	}

	public void appendGSEProjectToIndependentMatrices(File v,
													  DoubleMatrixDatasetAppendableWriter t1Mout,
													  DoubleMatrixDatasetAppendableWriter t1Uout,
													  DoubleMatrixDatasetAppendableWriter t2Mout,
													  DoubleMatrixDatasetAppendableWriter t2Uout,
													  HashMap<String, Integer> probeHashsT1,
													  HashMap<String, Integer> probeHashsT2) throws Exception {
		FolderObj folder = new FolderObj(v);
		if (folder.allfound()) {
			// get a list of samples
			ArrayList<String> allsamples = folder.getUniqueSamples();
			{
				DoubleMatrixDataset<String, String> m1red = DoubleMatrixDataset.loadDoubleData(folder.m1redfile.getAbsolutePath());
				DoubleMatrixDataset<String, String> m1grn = DoubleMatrixDataset.loadDoubleData(folder.m1grnfile.getAbsolutePath());

				int[] t1redmIndex = indexprobes(m1red, probeHashsT1);
				int[] t1grnmIndex = indexprobes(m1grn, probeHashsT1);

				ProgressBar pb = new ProgressBar(allsamples.size(), "Merging T1M");
				allsamples.stream().parallel().forEach(sample ->
						{
							try {
								Integer vm1grn = m1grn.getHashCols().get(sample);
								Integer vm1red = m1red.getHashCols().get(sample);

								String mergedsample = folder.gse + "_" + sample;
								mergeandWriteT1(m1grn, m1red, t1grnmIndex, t1redmIndex, mergedsample, vm1grn, vm1red, probeHashsT1.size(), t1Mout);


								pb.iterateSynched();
							} catch (NoSuchElementException e) {
								System.out.println("Sample " + sample + " not found for all probe types");
							} catch (IOException e) {
								e.printStackTrace();
							}
						}
				);
				pb.close();
			}
			{
				DoubleMatrixDataset<String, String> u1red = DoubleMatrixDataset.loadDoubleData(folder.u1redfile.getAbsolutePath());
				DoubleMatrixDataset<String, String> u1grn = DoubleMatrixDataset.loadDoubleData(folder.u1grnfile.getAbsolutePath());
				int[] t1redmIndex = indexprobes(u1red, probeHashsT1);
				int[] t1grnmIndex = indexprobes(u1grn, probeHashsT1);
				ProgressBar pb = new ProgressBar(allsamples.size(), "Merging T1U");
				allsamples.stream().parallel().forEach(sample ->
						{
							try {

								Integer vu1grn = u1grn.getHashCols().get(sample);
								Integer vu1red = u1red.getHashCols().get(sample);

//							System.out.println("Sample: " + sample);
								// now iterate probes
//							System.out.println(sample + "\tinit");
								String mergedsample = folder.gse + "_" + sample;
								mergeandWriteT1(u1grn, u1red, t1grnmIndex, t1redmIndex, mergedsample, vu1grn, vu1red, probeHashsT1.size(), t1Uout);
								pb.iterateSynched();
							} catch (NoSuchElementException e) {
								System.out.println("Sample " + sample + " not found for all probe types");
							} catch (IOException e) {
								e.printStackTrace();
							}
						}
				);
				pb.close();
			}

			{
				DoubleMatrixDataset<String, String> m2 = DoubleMatrixDataset.loadDoubleData(folder.m2file.getAbsolutePath());
				DoubleMatrixDataset<String, String> u2 = DoubleMatrixDataset.loadDoubleData(folder.u2file.getAbsolutePath());

				// index probes
				int[] t2uIndex = indexprobes(u2, probeHashsT2);
				int[] t2mIndex = indexprobes(m2, probeHashsT2);

				// merge U matrices for type 1
				ProgressBar pb = new ProgressBar(allsamples.size());
				allsamples.stream().parallel().forEach(sample ->
						{
							try {
								Integer vm2 = m2.getHashCols().get(sample);
								Integer vu2 = u2.getHashCols().get(sample);

//							System.out.println("Sample: " + sample);
								// now iterate probes
//							System.out.println(sample + "\tinit");
								String mergedsample = folder.gse + "_" + sample;
								mergeandWriteT2(m2, t2mIndex, mergedsample, vm2, probeHashsT2.size(), t2Mout);
								mergeandWriteT2(u2, t2uIndex, mergedsample, vu2, probeHashsT2.size(), t2Uout);


								pb.iterateSynched();
							} catch (NoSuchElementException e) {
								System.out.println("Sample " + sample + " not found for all probe types");
							} catch (IOException e) {
								e.printStackTrace();
							}
						}
				);
				pb.close();
			}

		}
	}

	private void mergeandWriteT2(DoubleMatrixDataset<String, String> m2, int[] t2mIndex, String
			mergedsample, Integer vm2, int size, DoubleMatrixDatasetAppendableWriter t2Mout) throws IOException {

		double[] output = new double[size];
		Arrays.fill(output, Double.NaN);
		for (int i = 0; i < t2mIndex.length; i++) {
			int index = t2mIndex[i];
			if (index != -1) {
//				int uind = mtoured[i];
				output[index] = m2.getElementQuick(i, vm2);
			}
		}

		t2Mout.append(output, mergedsample);

	}

	private void mergeandWriteT1(DoubleMatrixDataset<String, String> green,
								 DoubleMatrixDataset<String, String> red,
								 int[] t1grnmIndex,
								 int[] t1redmIndex,
								 String sample,
								 Integer vm1grn,
								 Integer vm1red,
								 int size,
								 DoubleMatrixDatasetAppendableWriter t1Mout) throws IOException {

		double[] output = new double[size];

		Arrays.fill(output, Double.NaN);

		// copy green
		for (int i = 0; i < t1grnmIndex.length; i++) {
			int index = t1grnmIndex[i];
			if (index != -1) {
				output[index] = green.getElementQuick(i, vm1grn);
			}
		}

		// copy red
		for (int i = 0; i < t1redmIndex.length; i++) {
			int index = t1redmIndex[i];
			if (index != -1) {
				output[index] = red.getElementQuick(i, vm1red);
			}
		}


		t1Mout.append(output, sample);

	}

	private int[] indexprobes(DoubleMatrixDataset<String, String> matrix, HashMap<String, Integer> probeHashs) {

		ArrayList<String> rowobs = matrix.getRowObjects();

		int[] indexarr = new int[rowobs.size()];

		IntStream.range(0, rowobs.size()).parallel().forEach(i -> {
					String probe = rowobs.get(i);
					Integer index = probeHashs.get(probe);
					if (index != null) {
						indexarr[i] = index;
					} else {
						indexarr[i] = -1;
					}
				}
		);

		return indexarr;
	}
}
