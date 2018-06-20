package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class SampleLinker {
	
	public static void main(String[] args) {
		
		
		SampleLinker s = new SampleLinker();
		try {
//			s.mayo("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\CBE-samples.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\samplestokeep-cbe.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\Individuals.txt",
//					null,
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\MAYO-DNAToRNA-CBE-samples.txt");
//
//			s.mayo("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\TCX-samples.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\samplestokeep-tcx.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\Individuals.txt",
//					null,
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMayo\\MAYO-DNAToRNA-TCX-samples.txt");
//
//			s.msbb("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMSBB\\avgsamples.txt",
//					null,
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMSBB\\Individuals.txt",
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMSBB\\AMP-AD_MSBB_WGS__sample_barcode_brainBankID.csv",
//					null,
//					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADMSBB\\MSBB-DNAToRNA.txt"
//			);
			
			s.rosmap("D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADRosmap\\samples.txt",
					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADRosmap\\samplestokeep.txt",
					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADRosmap\\Individuals.txt",
					null,
					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADRosmap\\ROSMAP_samplecoupling.txt",
					"D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-06-17-initialAnalysis\\0-Couplings\\AMPADRosmap\\ROSMAP-DNAToRNA.txt");
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void rosmap(String rna, String rnasamplestokeep, String dna, String dnasamplestokeep, String map, String outputfile) throws IOException {
		HashSet<String> rnakeephash = loadHash(rnasamplestokeep);
		HashSet<String> dnakeephash = loadHash(dnasamplestokeep);
		HashSet<String> rnasamples = loadHash(rna);
		HashMap<String, String> rnasamplesproc = new HashMap<String, String>();
		for (String sample : rnasamples) {
			if (rnakeephash == null || rnakeephash.contains(sample)) {
				String clean = sample.replaceAll("X", "");
				int lastunderscore = clean.lastIndexOf("_");
				clean = clean.substring(0, lastunderscore);
				rnasamplesproc.put(clean, sample);
			}
		}
		
		System.out.println(rnasamplesproc.size() + " rna samples");
		HashSet<String> dnasamples = loadHash(dna);
		
		int written = 0;
		TextFile out = new TextFile(outputfile, TextFile.W);
		TextFile tf = new TextFile(map, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> visiteddna = new HashSet<>();
		HashSet<String> visitedrna = new HashSet<>();
		while (elems != null) {
			String dnasample = elems[0];
			String rnasample = elems[1];
			
			if (dnasamples.contains(dnasample) && (dnakeephash == null || dnakeephash.contains(dnasample))) {
				String exprna = rnasamplesproc.get(rnasample);
				if (exprna != null && rnasamples.contains(exprna)) {
					if (!visiteddna.contains(dnasample) && !visitedrna.contains(rnasample)) {
						out.writeln(dnasample + "\t" + exprna);
						written++;
						visiteddna.add(dnasample);
						visitedrna.add(rnasample);
					}
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		out.close();
		tf.close();
		System.out.println(written + " written");
	}
	
	private void msbb(String rna, String rnasamplestokeep, String dna, String dnabarcodes, String dnasamplestokeep, String outputfile) throws IOException {
		HashSet<String> rnakeephash = loadHash(rnasamplestokeep);
		HashSet<String> dnakeephash = loadHash(dnasamplestokeep);
		HashSet<String> rnasamples = loadHash(rna);
		HashSet<String> dnasamples = loadHash(dna);

//		// match on the basis of ampad sample ids
//		HashMap<String, String> rnaToMSBB = new HashMap<String, String>();
//		TextFile tf = new TextFile(rnabarcode, TextFile.R);
//		tf.readLine();
//		String[] elems = tf.readLineElems(Strings.comma);
//		while (elems != null) {
//			rnaToMSBB.put("X" + elems[1], elems[0]);
//
//			elems = tf.readLineElems(Strings.comma);
//		}
//		tf.close();
		
		// match on the basis of ampad sample ids
		HashMap<String, String> MSBBToDNA = new HashMap<String, String>();
		TextFile tf = new TextFile(dnabarcodes, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(Strings.comma);
		while (elems != null) {
			MSBBToDNA.put(elems[0], elems[1]);
			elems = tf.readLineElems(Strings.comma);
		}
		tf.close();
		
		System.out.println(rnasamples.size() + " rna ids");
		System.out.println(dnasamples.size() + " dna ids");
		int written = 0;
		TextFile out = new TextFile(outputfile, TextFile.W);
		for (String rnasample : rnasamples) {
			if (rnakeephash == null || rnakeephash.contains(rnasample)) {
				
				String msbb = rnasample;//rnaToMSBB.get(rnasample);
				if (msbb != null) {
					String dnasample = MSBBToDNA.get(msbb);
					if (dnasample != null && dnasamples.contains(dnasample)) {
						if (dnakeephash == null || dnakeephash.contains(dnasample)) {
							out.writeln(dnasample + "\t" + rnasample);
							written++;
						}
					}
				}
			}
		}
		out.close();
		System.out.println(written + " written ");
	}
	
	public void mayo(String rna, String rnasamplestokeep, String dna, String dnasamplestokeep, String out) throws IOException {
		
		
		HashSet<String> rnakeephash = loadHash(rnasamplestokeep);
		HashSet<String> dnakeephash = loadHash(dnasamplestokeep);
		HashSet<String> rnasamples = loadHash(rna);
		HashSet<String> dnasamples = loadHash(dna);
		
		TextFile outf = new TextFile(out, TextFile.W);
		int nrsamples = rnasamples.size();
		int nrmatched = 0;
		for (String s : rnasamples) {
			if (rnakeephash == null || rnakeephash.contains(s)) {
				String clean = s.replaceAll("X", "");
				clean = clean.replaceAll("_CER", "");
				clean = clean.replaceAll("_TC", "");
				if ((dnakeephash != null && dnakeephash.contains(clean)) || dnasamples.contains(clean)) {
					outf.writeln(clean + "\t" + s);
					nrmatched++;
				}
			}
		}
		outf.close();
		System.out.println(nrmatched + " matched out of " + nrsamples);
		
	}
	
	private HashSet<String> loadHash(String rnasamplestokeep) throws IOException {
		HashSet<String> rnakeephash = null;
		if (rnasamplestokeep != null) {
			TextFile tfk1 = new TextFile(rnasamplestokeep, TextFile.R);
			ArrayList<String> rnasamplestokeeparr = tfk1.readAsArrayList();
			tfk1.close();
			
			rnakeephash = new HashSet<>();
			rnakeephash.addAll(rnasamplestokeeparr);
		}
		return rnakeephash;
	}
}
