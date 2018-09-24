package nl.harmjanwestra.playground.methylation;

import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class Merge450KData {
	public static void main(String[] args) {
		Merge450KData m = new Merge450KData();
		
		String[] exclusion = new String[]{
				"S:\\projects\\2018-methylation\\code\\450K_DataProcessing\\ADDITIONAL_INFO\\ProbeFiltering\\freq5percent\\probeToFilter_450K_1000G_omni2.5.hg19.EUR_alleleFreq5percent_50bp_wInterroSite.txt",
				"S:\\projects\\2018-methylation\\code\\450K_DataProcessing\\ADDITIONAL_INFO\\ProbeFiltering\\ProbesBindingNonOptimal\\Source&BSProbesMappingMultipleTimesOrNotBothToBSandNormalGenome.txt"
		};

//		String infolder = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT\\";
		String probelistfile = args[0]; //"S:\\projects\\2018-methylation\\GPL13534_450k_probelist.txt";
//		try {
//			m.makeprobelist(infolder, probelistfile, exclusion);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

//		String probelistfile = "";
		String infolder = args[1]; //"S:\\projects\\2018-methylation\\GPL13534_450k_OUT\\";
		String outfile = args[2]; // "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged\\FromIdat.txt.gz";
		
		
		try {
			m.run(probelistfile, infolder, outfile);
		} catch (IOException e) {
			e.printStackTrace();
		}
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
	
	public void run(String probelistfile, String inFolder, String outfile) throws IOException {
		
		TextFile lf = new TextFile(probelistfile, TextFile.R);
		ArrayList<String> probelist = lf.readAsArrayList();
		Collections.sort(probelist);
		lf.close();
		
		HashMap<String, Integer> probeHashs = new HashMap<String, Integer>();
		int c = 0;
		for (String s : probelist) {
			probeHashs.put(s, c);
			c++;
		}
		
		TextFile out = new TextFile(outfile, TextFile.W);
		String header = "-\t" + Strings.concat(probelist, Strings.tab);
		out.writeln(header);
		File[] list = new File(inFolder).listFiles();
		for (File f : list) {
			if (f.isDirectory()) {
				try {
					append(f, out, probeHashs);
					System.gc();
				} catch (IOException e) {
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
	
	public void append(File v, TextFile out, HashMap<String, Integer> probeHashs) throws IOException {
		
		File[] list = v.listFiles();
		if (list.length > 5) {
			DoubleMatrixDataset<String, String> m1red = null;
			DoubleMatrixDataset<String, String> m1grn = null;
			DoubleMatrixDataset<String, String> u1red = null;
			DoubleMatrixDataset<String, String> u1grn = null;
			DoubleMatrixDataset<String, String> m2 = null;
			DoubleMatrixDataset<String, String> u2 = null;
			
			// GSE100134_T1_Grn_M_Signal.txt.gz  GSE100134_T1_Grn_U_Signal.txt.gz  GSE100134_T1_Red_M_Signal.txt.gz  GSE100134_T1_Red_U_Signal.txt.gz  GSE100134_T2_M_Signal.txt.gz  GSE100134_T2_U_Signal.txt.gz
			
			HashSet<String> uniqueSamples = new HashSet<String>();
			for (File f : list) {
				if (f.getName().endsWith("T1_Grn_M_Signal.txt.gz")) {
					m1grn = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(m1grn.getColObjects());
				} else if (f.getName().endsWith("T1_Grn_U_Signal.txt.gz")) {
					u1grn = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(u1grn.getColObjects());
				} else if (f.getName().endsWith("T1_Red_M_Signal.txt.gz")) {
					m1red = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(m1red.getColObjects());
				} else if (f.getName().endsWith("T1_Red_U_Signal.txt.gz")) {
					u1red = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(u1red.getColObjects());
				} else if (f.getName().endsWith("T2_M_Signal.txt.gz")) {
					m2 = DoubleMatrixDataset.loadDoubleData(f.getPath());
					uniqueSamples.addAll(m2.getColObjects());
				} else if (f.getName().endsWith("T2_U_Signal.txt.gz")) {
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
			
			
			ProgressBar pb = new ProgressBar(allsamples.size());
			
			DoubleMatrixDataset<String, String> finalM1grn = m1grn;
			DoubleMatrixDataset<String, String> finalU1grn = u1grn;
			DoubleMatrixDataset<String, String> finalM1red = m1red;
			DoubleMatrixDataset<String, String> finalU1red = u1red;
			DoubleMatrixDataset<String, String> finalM = m2;
			DoubleMatrixDataset<String, String> finalU = u2;
			AtomicInteger ctr = new AtomicInteger();
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
							out.writelnsynced(sample + "\t" + Strings.concat(output, Strings.tab));
							
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
}
