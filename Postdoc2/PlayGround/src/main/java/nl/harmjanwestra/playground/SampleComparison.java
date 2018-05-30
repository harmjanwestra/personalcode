package nl.harmjanwestra.playground;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

public class SampleComparison {
	public static void main(String[] args) {
		SampleComparison c = new SampleComparison();

//		try {
////			c.ntr("C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\SampleSheets\\LudeMatch.txt",
////					"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GTE\\NTR_Methylation_GTE.txt",
////					"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\NTR_Methylation_GTE.txt"
////			);
////			c.sampleSheet("C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\SampleSheets\\sample_sheet_NTR.txt",
////					"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\SampleSheets\\sample_sheet_NTR_genotypeToMethylation.txt");
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//
		
		
		String[] gtefiles = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GTE\\CODAM_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GTE\\LLDeep_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GTE\\LLS_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GTE\\LLS_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GTE\\RS_Methylation_GTE.txt"
		};
		
		String[] genotypeInds = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GT\\CODAMIndividuals.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GT\\LLIndividuals.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GT\\LL_660QIndividuals.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GT\\LLS_OmniExprIndividuals.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\GT\\RS_OmniExprIndividuals.txt"
		};
		
		String[] outfiles = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\CODAM_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\LLDeep_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\LLS_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\LLS2_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\RS_Methylation_GTE.txt"
		};


//		SampleComparison c = new SampleComparison();
//		try {
//			c.combineEPRSWithCisAndTransFX(gtefiles[0], genotypeInds[0], outfiles[0]);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//
//		for (int i = 0; i < gtefiles.length; i++) {
//			try {
////				c.combineEPRSWithCisAndTransFX(gtefiles[i], genotypeInds[i], outfiles[i]);
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}
		
		String[] infilesmeth = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\CODAM_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\LLDeep_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\LLS_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\LLS2_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\RS_Methylation_GTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newGTE\\NTR_Methylation_GTE.txt"
		};
		
		String[] infilesge = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\gteExp\\GTE_LLDEEP_and_BIOS_last_related_removed_110417.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\gteExp\\GTE_LLDEEP_and_BIOS_last_related_removed_110417.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\gteExp\\GTE_LLDEEP_and_BIOS_last_related_removed_110417.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\gteExp\\GTE_LLDEEP_and_BIOS_last_related_removed_110417.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\gteExp\\GTE_LLDEEP_and_BIOS_last_related_removed_110417.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\gteExp\\GTE_LLDEEP_and_BIOS_last_related_removed_110417.txt",
		};
		
		String[] mtefiles = new String[]{
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newMTE\\CODAM_MTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newMTE\\LLDeep_MTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newMTE\\LLS_660Q_MTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newMTE\\LLS_Omni_MTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newMTE\\RS_MTE.txt",
				"C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\newMTE\\NTR.txt",
		};
		
		String sampleList = "C:\\Sync\\OneDrive\\Postdoc2\\2017-11-eQTLMeta\\meQTL\\expressionSamples.txt";
		sampleList = null;
		
		try {
			c.makeMethylationToExpressionFiles(infilesmeth, infilesge, sampleList, mtefiles);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void makeMethylationToExpressionFiles(String[] infilesmeth, String[] infilesge, String expressionSampleFile, String[] mtefiles) throws IOException {
		
		// go from G->M
		// G->E .. E->G->M
		
		HashSet<String> allowedExp = null;
		if (expressionSampleFile != null) {
			
			TextFile tf = new TextFile(expressionSampleFile, TextFile.R);
			ArrayList<String> samples = tf.readAsArrayList();
			HashSet<String> s = new HashSet<String>();
			s.addAll(samples);
			tf.close();
			allowedExp = s;
			System.out.println(allowedExp.size() + " samples in expression sample list");
			if (allowedExp.contains("BD1NYRACXX-5-3")) {
				System.out.println("Found it!");
			}
			
		}
		
		int sum = 0;
		for (int i = 0; i < infilesmeth.length; i++) {
			TextFile tf = new TextFile(infilesmeth[i], TextFile.R);
			HashMap<String, String> gtm = new HashMap<String, String>();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length > 1) {
					gtm.put(elems[0], elems[1]);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			TextFile out = new TextFile(mtefiles[i], TextFile.W);
			TextFile gte = new TextFile(infilesge[i], TextFile.R);
			elems = gte.readLineElems(TextFile.tab);
			int samples = 0;
			while (elems != null) {
				if (elems.length > 1) {
					String gt = elems[0];
					String exp = elems[1];
					if (allowedExp == null || allowedExp.contains(exp)) {
						String meth = gtm.get(gt);
						if (meth != null) {
							out.writeln(meth + "\t" + exp);
							samples++;
						}
						
					}
				}
				elems = gte.readLineElems(TextFile.tab);
			}
			gte.close();
			out.close();
			
			System.out.println(samples + "\t" + mtefiles[i]);
			
			sum += samples;
		}
		System.out.println(sum + " total samples");
		
	}
	
	public void ntr(String lmatch, String gte, String out) throws IOException {
		TextFile tf = new TextFile(lmatch, TextFile.R);
		HashMap<String, String> sampleMap = new HashMap<String, String>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String samplenew = elems[0];
			String sampleOld = elems[1];
			if (sampleOld.contains("NTR")) {
				sampleOld = sampleOld.replaceAll("NTR-", "");
				sampleMap.put(sampleOld, samplenew);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(sampleMap.size() + " loaded from: " + lmatch);
		
		TextFile tf3 = new TextFile(out, TextFile.W);
		TextFile tf2 = new TextFile(gte, TextFile.R);
		elems = tf2.readLineElems(TextFile.tab);
		int found = 0;
		int searched = 0;
		while (elems != null) {
			if (elems.length > 1) {
				String gt = elems[0];
				gt = gt.substring(1);
				String meth = elems[1];
				String newGt = sampleMap.get(gt);
				if (newGt != null) {
					tf3.writeln(newGt + "\t" + meth);
					found++;
				} else {
					System.out.println("Could not find: " + gt + "\t" + meth);
				}
				searched++;
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf3.close();
		tf2.close();
		System.out.println(searched + " searched and " + found + " found");
	}
	
	private void sampleSheet(String s, String outf) throws IOException {
		TextFile tf = new TextFile(s, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.comma);
		TextFile out = new TextFile(outf, TextFile.W);
		while (elems != null) {
			
			String samples = elems[0];
			String meth = elems[3];
			
			String[] sampleElems = samples.split("-");
			out.writeln(sampleElems[sampleElems.length - 1] + "\t" + meth + "\t" + samples);
			
			elems = tf.readLineElems(TextFile.comma);
		}
		out.close();
		tf.close();
		
		
	}
	
	public void run(String sampleset1, String sampleSet2, String out) throws IOException {
		System.out.flush();
		System.err.flush();
		
		ArrayList<String> list1 = getSamples(sampleset1);
		System.out.println(list1.size() + " from " + sampleset1);
		ArrayList<String> list2 = getSamples(sampleSet2);
		System.out.println(list2.size() + " from " + sampleSet2);
		
		HashSet<String> list1hash = new HashSet<String>();
		list1hash.addAll(list1);
		HashMap<String, String> sampleToSample = new HashMap<String, String>();
		for (String sample : list2) {
			if (list1hash.contains(sample)) {
				sampleToSample.put(sample, sample);
			} else {
				// try splitting
				String[] sampleElems = sample.split("_");
				if (sampleElems.length > 1) {
//					System.out.println("Trying ");
					String sampleSplit = sampleElems[1];
					if (list1hash.contains(sampleSplit)) {
						sampleToSample.put(sampleSplit, sample);
					} else {
						// try a combo of
						sampleSplit = Strings.concat(sampleElems, Pattern.compile("_"), 1, sampleElems.length);
						if (list1hash.contains(sampleSplit)) {
							sampleToSample.put(sampleSplit, sample);
						} else {
							sampleSplit = "F" + sampleElems[1];
							if (list1hash.contains(sampleSplit)) {
								sampleToSample.put(sampleSplit, sample);
							}
						}
					}
				}
			}
		}
		System.out.println("Found: " + sampleToSample.size() + "\tout of\t" + list1.size() + "\tor " + list1hash.size() + " samples from: " + sampleset1);
		for (String sample : list1) {
			if (!sampleToSample.containsKey(sample)) {
				System.out.println("Could not find:\t" + sample);
			}
		}
		System.out.flush();
		System.err.flush();
		
		TextFile outf = new TextFile(out, TextFile.W);
		TextFile in = new TextFile(sampleset1, TextFile.R);
		String[] elems = in.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 1) {
				String sample = sampleToSample.get(elems[0]);
				if (sample != null) {
					outf.writeln(sample + "\t" + elems[1]);
				}
			}
			elems = in.readLineElems(TextFile.tab);
		}
		outf.close();
		System.out.println();
	}
	
	public ArrayList<String> getSamples(String list) throws IOException {
		ArrayList<String> samples = new ArrayList<String>();
		TextFile tf = new TextFile(list, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			ln = ln.trim();
			String[] elms = ln.split("\t");
			ln = elms[0];
			
			if (ln.length() > 0) {
				samples.add(ln);
			}
			ln = tf.readLine();
		}
		tf.close();
		return samples;
	}
}
