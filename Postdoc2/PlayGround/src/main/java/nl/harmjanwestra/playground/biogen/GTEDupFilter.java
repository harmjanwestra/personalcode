package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class GTEDupFilter {
	public static void main(String[] args) {
		GTEDupFilter f = new GTEDupFilter();
		
		String targetalsgenotypes = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-SampleAssignment.txt";
		String targetalsoriggte = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-GTE-orig.txt";
		String additionaldnaids = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-WGSToExternal-orig.txt";
		String additionalrnaids = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-RNAToExternal-orig.txt";
		String rnaseqsamples = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS.geneCounts.samples.txt";
		String linksout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-GTE-unique.txt";
		
		GTEDupFilter t = new GTEDupFilter();
		try {
//			t.targetals(rnaseqsamples, targetalsgenotypes, targetalsoriggte, additionaldnaids, additionalrnaids, linksout);
			
			String rnaToTissue = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\RNAToTissueType.txt";
			String linkfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-GTE-unique.txt-allcombos.txt";
			String outdir = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\gtepertissue\\";
			t.targetalsSplitPerTissue(rnaToTissue, linkfile, outdir);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void targetalsSplitPerTissue(String rnaToTissue, String linkfile, String outdir) throws IOException {
		
		HashSet<String> uniqueGroups = new HashSet<String>();
		HashMap<String, TextFile> groupToTextFile = new HashMap<>();
		
		HashMap<String, String> sampleToGroup = new HashMap<String, String>();
		TextFile tf = new TextFile(rnaToTissue, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String sample = elems[0];
			String group = elems[2];
			
			if (elems[0].contains("01218")) {
				System.out.println("Ping!");
			}
			sampleToGroup.put(sample, group);
			
			if (!groupToTextFile.containsKey(group)) {
				TextFile o = new TextFile(outdir + "GTE-" + elems[2] + ".txt", TextFile.W);
				groupToTextFile.put(elems[2], o);
			}
			uniqueGroups.add(elems[2]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		
		HashMap<String, HashSet<String>> dnaswrittenpergroup = new HashMap<>();
		
		TextFile tf2 = new TextFile(linkfile, TextFile.R);
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String dna = elems[0];
			String rna = elems[1];
			String shrna = elems[3];
			String group = sampleToGroup.get(shrna);
			HashSet<String> list = dnaswrittenpergroup.get(group);
			if (list == null) {
				list = new HashSet<>();
			}
			if (!list.contains(dna)) {
				// write
				
				TextFile tfo = groupToTextFile.get(group);
				if (tfo == null) {
					System.out.println("Could not find tfo for group: " + group + " for sample " + shrna);
				} else {
					tfo.writeln(dna + "\t" + rna);
				}
				
				
				list.add(dna);
			}
			dnaswrittenpergroup.put(group, list);
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		for (String group : uniqueGroups) {
			groupToTextFile.get(group).close();
		}
	}
	
	public void targetals(String rnaseqsamples, String genotypesamples, String links, String additionalWGSIds, String additionalRNAIds, String linkout) throws IOException {
		
		HashMap<String, String> additionalWGSIdSet = new HashMap<String, String>();
		{
			TextFile tf5 = new TextFile(additionalWGSIds, TextFile.R);
			String[] elems5 = tf5.readLineElems(TextFile.tab);
			while (elems5 != null) {
				if (elems5.length > 1) {
					String[] meh = elems5[0].split("/");
					for (String s : meh) {
						additionalWGSIdSet.put(s, elems5[1]);
					}
					
				}
				elems5 = tf5.readLineElems(TextFile.tab);
			}
			tf5.close();
		}
		HashMap<String, String> additionalRNAIdSet = new HashMap<String, String>();
		TextFile tf5 = new TextFile(additionalRNAIds, TextFile.R);
		String[] elems5 = tf5.readLineElems(TextFile.tab);
		while (elems5 != null) {
			if (elems5.length > 1) {
				String[] meh = elems5[0].split("/");
				for (String s : meh) {
					additionalRNAIdSet.put(s, elems5[1]);
				}
				
			}
			elems5 = tf5.readLineElems(TextFile.tab);
		}
		tf5.close();
		
		TextFile tf = new TextFile(rnaseqsamples, TextFile.R);
		String lnr = tf.readLine();
		HashSet<String> uniqueRNAseq = new HashSet<String>();
		while (lnr != null) {
			if (lnr.length() > 0) {
				uniqueRNAseq.add(lnr);
			}
			lnr = tf.readLine();
		}
		tf.close();
		
		System.out.println(uniqueRNAseq.size() + " unique RNAs");
		
		TextFile tf2 = new TextFile(genotypesamples, TextFile.R);
		HashSet<String> uniqueDNAseq = new HashSet<>();
		String[] elems2 = tf2.readLineElems(TextFile.tab);
		while (elems2 != null) {
			uniqueDNAseq.add(elems2[0]);
			elems2 = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		System.out.println(uniqueDNAseq.size() + " unique DNAs");
		
		TextFile tf3 = new TextFile(links, TextFile.R);
		HashSet<String> dnaidsinlinkfile = new HashSet<>();
		HashSet<String> rnaidsinlinkfile = new HashSet<>();
		String[] elems = tf3.readLineElems(TextFile.tab);
		
		HashSet<String> dnaswritten = new HashSet<>();
		HashSet<String> dnasfound = new HashSet<>();
		HashSet<String> rnasfound = new HashSet<>();
		TextFile tfout = new TextFile(linkout, TextFile.W);
		TextFile tfout3 = new TextFile(linkout + "-allcombos.txt", TextFile.W);
		TextFile tfout2 = new TextFile(linkout + "-dnasnotfound.txt", TextFile.W);
		int nrdnasfound = 0;
		while (elems != null) {
			if (elems.length > 1) {
				String[] dnaelems = elems[0].split("/");
				boolean written = false;
				for (String dna : dnaelems) {
					
					String origS = new String(dna);
					if (!dna.equals("Not Applicable")) {
						
						if (dna.contains("01696")) {
							System.out.println("ping!");
						}
						
						dna = dna.replaceAll("_", "-");
						if (!uniqueDNAseq.contains(dna)) {
							dna = dna + "-b38";
						}
						
						
						if (!uniqueDNAseq.contains(dna)) {
							String otherId = additionalWGSIdSet.get(origS);
							if (!uniqueDNAseq.contains(dna)) {
								otherId += "-b38";
							}
							if (otherId != null && uniqueDNAseq.contains(otherId)) {
								dna = otherId;
//								System.out.println("replacing " + s + " with " + otherId);
							}
						}
						
						dnaidsinlinkfile.add(dna);
						if (!uniqueDNAseq.contains(dna)) {
//							System.out.println("Could not find DNA: " + s);
							
							tfout2.writeln(dna + "\t" + elems[1]);
						} else {
							nrdnasfound++;
							
							// now match the rnas
							String rna = elems[1];
							dnasfound.add(dna);
//							if (!uniqueRNAseq.contains(rna)) {
							// try combo with external id
							String origrna = new String(rna);
							String match = null;
							for (String rnaseqs : uniqueRNAseq) {
								if (rnaseqs.contains(rna)) {
									match = rnaseqs;
								}
							}
							if (match == null) {
								System.out.println("Could not find rna: " + origrna);
							} else {
								rna = match;
							}


//							}
							if (uniqueRNAseq.contains(rna)) {
								rnasfound.add(rna);
								String out = dna + "\t" + rna + "\t" + elems[0] + "\t" + elems[1];
								tfout3.writeln(out);
								if (!dnaswritten.contains(dna)) {
									tfout.writeln(out);
									dnaswritten.add(dna);
									written = true;
								}
							} else {
//								System.out.println(rna + " RNA not found");
							}
						}
					}
				}
				rnaidsinlinkfile.add(elems[1]);
			}
			elems = tf3.readLineElems(TextFile.tab);
		}
		
		tf3.close();
		tfout.close();
		tfout2.close();
		tfout3.close();
		
		System.out.println(dnaidsinlinkfile.size() + " DNAs in link file");
		System.out.println(rnaidsinlinkfile.size() + " RNAs in link file");
		
		System.out.println(nrdnasfound + " DNAs from link file found in WGS data");
		
		
		System.out.println();
		System.out.println();
		System.out.println("DNAs not found");
		int ctr = 0;
		for (String s : uniqueDNAseq) {
			if (!dnasfound.contains(s)) {
				System.out.println(s);
				ctr++;
			}
		}
		System.out.println("Total DNAs missing: " + ctr);
		
		System.out.println();
		System.out.println("RNAs not found: ");
		ctr = 0;
		for (String s : uniqueRNAseq) {
			if (!rnasfound.contains(s)) {
				System.out.println(s);
				ctr++;
			}
		}
		System.out.println("Total RNAs missing: " + ctr);
	}
}
