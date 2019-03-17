package nl.harmjanwestra.playground.biogen.datasets.gtex;


import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class SplitTissues {
	
	public static void main(String[] args) {
		
		String annotfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\SraRunTable.txt";
		String samplefile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\rnasamples\\rnaseqsamples.txt";
		String rnaseqwithtissues = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\rnasamples\\rnasamplesWithTissue.txt";
		
		SplitTissues t = new SplitTissues();
		try {
			t.determineTissuePerSample(annotfile, samplefile, rnaseqwithtissues);
			
			String abbreviationFile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\rnasamples\\tissueAbbreviations.txt";
			String outdir = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\initialGTEs\\";
			String europeanGTs = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-09-27-GTEx\\genotypepca\\GTEx-Europeans.txt";
			
			t.writeGTEs(annotfile, samplefile, abbreviationFile, outdir, europeanGTs);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void writeGTEs(String annotfile, String sampleFile, String abbreviationfile, String outdir, String europeangts) throws IOException {
		
		TextFile tf3 = new TextFile(europeangts, TextFile.R);
		ArrayList<String> eurlist = tf3.readAsArrayList();
		HashSet<String> eurset = new HashSet<String>();
		eurset.addAll(eurlist);
		tf3.close();
		
		HashMap<String, String> tissueToAbbrev = new HashMap<>();
		TextFile tfq = new TextFile(abbreviationfile, TextFile.R);
		String[] aelems = tfq.readLineElems(TextFile.tab);
		while (aelems != null) {
			tissueToAbbrev.put(aelems[1], aelems[0]);
			aelems = tfq.readLineElems(TextFile.tab);
		}
		tfq.close();
		
		TextFile annot = new TextFile(annotfile, TextFile.R);
		String[] elems = annot.readLineElems(TextFile.tab);
		HashMap<String, SRASample> sraToSample = new HashMap<String, SRASample>();
		while (elems != null) {
			String sample = elems[17];
			String gt = elems[33];
			String tissue = elems[22];
			String hist = elems[24];
			SRASample sra = new SRASample();
			sra.id = sample;
			sra.gt = gt;
			sra.tissue = tissue;
			sra.hist = hist;
			sraToSample.put(sample + "_" + sample, sra);
			elems = annot.readLineElems(TextFile.tab);
		}
		annot.close();
		
		TextFile tf = new TextFile(sampleFile, TextFile.R);
		ArrayList<String> list = tf.readAsArrayList();
		tf.close();
		
		HashMap<String, TextFile> abbrevToTXT = new HashMap<>();
		HashMap<String, HashSet<String>> abbrevToSampleSet = new HashMap<>();
		
		for (String l : list) {
			SRASample s = sraToSample.get(l);
			if (s != null) {
				
				String abbrev = tissueToAbbrev.get(s.tissue);
				if (abbrev != null) {
					HashSet<String> seenGTs = abbrevToSampleSet.get(abbrev);
					if (seenGTs == null || !seenGTs.contains(s.gt) && eurset.contains(s.gt)) {
						if (seenGTs == null) {
							seenGTs = new HashSet<>();
						}
						TextFile tfo = abbrevToTXT.get(abbrev);
						if (tfo == null) {
							tfo = new TextFile(outdir + "GTEx-" + abbrev + ".txt", TextFile.W);
						}
						tfo.writeln(s.gt + "\t" + l);
						seenGTs.add(s.gt);
						abbrevToSampleSet.put(abbrev, seenGTs);
						abbrevToTXT.put(abbrev, tfo);
					}
				}
			}
		}
		
		// close text files
		for (String key : abbrevToTXT.keySet()) {
			
			abbrevToTXT.get(key).close();
		}
		
	}
	
	
	class SRASample {
		
		String id;
		String gt;
		String tissue;
		String hist;
		
	}
	
	public void determineTissuePerSample(String annotfile, String sampleFile, String out) throws IOException {
		
		TextFile annot = new TextFile(annotfile, TextFile.R);
		String[] elems = annot.readLineElems(TextFile.tab);
		HashMap<String, SRASample> sraToSample = new HashMap<String, SRASample>();
		while (elems != null) {
			String sample = elems[17];
			String gt = elems[33];
			String tissue = elems[22];
			String hist = elems[24];
			SRASample sra = new SRASample();
			sra.id = sample;
			sra.gt = gt;
			sra.tissue = tissue;
			sra.hist = hist;
			
			sraToSample.put(sample + "_" + sample, sra);
			
			elems = annot.readLineElems(TextFile.tab);
		}
		annot.close();
		
		TextFile tf = new TextFile(sampleFile, TextFile.R);
		ArrayList<String> list = tf.readAsArrayList();
		tf.close();
		
		TextFile outf = new TextFile(out, TextFile.W);
		TextFile outf2 = new TextFile(out + "-withoutGenotypeDups.txt", TextFile.W);
		outf.writeln("Gt\tSample\tHistologicalSite\tTissue");
		HashSet<String> seenGTs = new HashSet<>();
		for (String l : list) {
			SRASample s = sraToSample.get(l);
			if (s != null) {
				outf.writeln(s.gt + "\t" + l + "\t" + s.hist + "\t" + s.tissue);
				if (!seenGTs.contains(s.gt)) {
					outf2.writeln(s.gt + "\t" + l + "\t" + s.hist + "\t" + s.tissue);
					seenGTs.add(s.gt);
				}
			} else {
				System.out.println("Could not find file: " + l);
			}
		}
		outf.close();
		outf2.close();
	}
}
