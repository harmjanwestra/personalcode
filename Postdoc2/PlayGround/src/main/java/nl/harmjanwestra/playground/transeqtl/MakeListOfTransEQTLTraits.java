package nl.harmjanwestra.playground.transeqtl;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class MakeListOfTransEQTLTraits {
	
	
	public static void main(String[] args) {
		String transfx = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\trans\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
		String traitfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\genetic_risk_factors_cleaned_traits_standardized_tested_in_trans_20170908.txt";
		String rowtraitfile = null;//"";
		String coltraitfile = null; //"";
		boolean usejaccard = true;
		String outfilename = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-30-comparisons\\transcomparisons.txt";
		
		MakeListOfTransEQTLTraits t = new MakeListOfTransEQTLTraits();
		try {
//			t.runTrans(transfx, traitfile, rowtraitfile, coltraitfile, usejaccard, outfilename);
			
			
			String prsfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\eprs\\eQTLsFDR-Significant-0.05.txt.gz";
			rowtraitfile = null;
			coltraitfile = null;
			usejaccard = true;
			outfilename = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-30-comparisons\\prscomparisons.txt";
//			t.runPRS(prsfile, rowtraitfile, coltraitfile, usejaccard, outfilename);
			
			prsfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\eprs\\eQTLsFDR-Significant-0.05.txt.gz";
			rowtraitfile = null;
			coltraitfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\prstraits\\2018-05-23-immuneprs.txt";
			usejaccard = true;
			outfilename = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-30-comparisons\\prscomparisons-selectedimmunetraits.txt";
//			t.runPRS(prsfile, rowtraitfile, coltraitfile, usejaccard, outfilename);
			
			transfx = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\trans\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
			traitfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\genetic_risk_factors_cleaned_traits_standardized_tested_in_trans_20170908.txt";
			prsfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\eprs\\eQTLsFDR-Significant-0.05.txt.gz";
			String transtraitfile = null;
			String prstraitfile = null;
			usejaccard = true;
			outfilename = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-30-comparisons\\prscomparisonswithtrans.txt";
			t.prsVersusTrans(prsfile, prstraitfile, traitfile, transfx, transtraitfile, usejaccard, outfilename);
			
			transfx = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\trans\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
			traitfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\genetic_risk_factors_cleaned_traits_standardized_tested_in_trans_20170908.txt";
			prsfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\eprs\\eQTLsFDR-Significant-0.05.txt.gz";
			transtraitfile = null;
			prstraitfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\prstraits\\2018-05-23-immuneprs.txt";
			usejaccard = true;
			outfilename = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-30-comparisons\\prscomparisonswithtrans-selectedimmunetraits.txt";
			t.prsVersusTrans(prsfile, prstraitfile, traitfile, transfx, transtraitfile, usejaccard, outfilename);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void prsVersusTrans(String prsfile, String prstraitfile, String gwastraitfile, String transfx, String transtraitfile, boolean usejaccard, String outfilename) throws IOException {
		QTLTextFile tf = new QTLTextFile(prsfile, TextFile.R);
		EQTL[] eprs = tf.read();
		tf.close();
		
		HashMap<String, HashSet<String>> prsgenespertrait = new HashMap<>();
		for (EQTL e : eprs) {
			String prs = e.getRsName();
			int lastunderscore = prs.lastIndexOf("_");
			String prsstripped = prs.substring(0, lastunderscore);
			HashSet<String> genes = prsgenespertrait.get(prsstripped);
			if (genes == null) {
				genes = new HashSet<>();
			}
			genes.add(e.getProbe());
			prsgenespertrait.put(prsstripped, genes);
		}
		
		TextFile tf2 = new TextFile(gwastraitfile, TextFile.R);
		tf2.readLine();
		HashMap<String, HashSet<String>> snptotrait = new HashMap<>();
		String[] elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			String traits = elems[5];
			String[] traitelems = traits.split("; ");
			
			HashSet<String> snpset = snptotrait.get(snp);
			if (snpset == null) {
				snpset = new HashSet<>();
			}
			for (String trait : traitelems) {
				if (trait.trim().length() > 0) {
					snpset.add(trait);
				}
			}
			snptotrait.put(snp, snpset);
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		System.out.println(snptotrait.size() + " trait snps ");
		
		// load transfx
		QTLTextFile t = new QTLTextFile(transfx, TextFile.R);
		EQTL[] eqtls = t.read();
		t.close();
		
		HashMap<String, ArrayList<EQTL>> eqtlspertrait = new HashMap<String, ArrayList<EQTL>>();
		HashSet<String> selectedTraits = new HashSet<>();
		HashSet<String> genes = new HashSet<String>();
		for (EQTL e : eqtls) {
			HashSet<String> traits = snptotrait.get(e.getRsName());
			genes.add(e.getProbe());
			if (traits != null) {
				for (String trait : traits) {
					ArrayList<EQTL> traiteqtls = eqtlspertrait.get(trait);
					if (traiteqtls == null) {
						traiteqtls = new ArrayList<>();
					}
					traiteqtls.add(e);
					selectedTraits.add(trait);
					eqtlspertrait.put(trait, traiteqtls);
				}
			}
		}
		System.out.println(eqtlspertrait.size() + " traits with transfx");
		System.out.println(genes.size() + " total genes");
		
		
		ArrayList<String> prstraitAr = new ArrayList<>();
		prstraitAr.addAll(prsgenespertrait.keySet());
		Collections.sort(prstraitAr);
		
		ArrayList<String> transtraitarr = new ArrayList<>();
		transtraitarr.addAll(selectedTraits);
		Collections.sort(transtraitarr);
		
		ArrayList<String> rowtraits;
		ArrayList<String> coltraits;
		
		if (prstraitfile != null) {
			coltraits = determinetraitlist(prstraitfile, prstraitAr);
		} else {
			coltraits = prstraitAr;
		}
		
		if (transtraitfile != null) {
			rowtraits = determinetraitlist(transtraitfile, transtraitarr);
		} else {
			rowtraits = transtraitarr;
		}
		
		// make header
		String header0 = "-\t-";
		for (String s : coltraits) {
			header0 += "\t" + s;
		}
		String header1 = "TraitA\tNumberOfGenes";
		for (String s : coltraits) {
			HashSet<String> genesA = prsgenespertrait.get(s);
			header1 += "\t" + genesA.size();
		}
		
		TextFile outf = new TextFile(outfilename, TextFile.W);
		
		outf.writeln(header0);
		outf.writeln(header1);
		
		for (int i = 0; i < rowtraits.size(); i++) {
			String traitA = rowtraits.get(i);
			ArrayList<EQTL> eqtlsfortraitA = eqtlspertrait.get(traitA);
			HashSet<String> genesA = getGenes(eqtlsfortraitA);
			String ln = traitA + "\t" + genesA.size();
			for (int j = 0; j < coltraits.size(); j++) {
				String traitB = coltraits.get(j);
				HashSet<String> genesB = prsgenespertrait.get(traitB);
				int intersect = 0;
				for (String e : genesA) {
					if (genesB.contains(e)) {
						intersect++;
					}
				}
				if (usejaccard) {
					HashSet<String> union = new HashSet<>();
					union.addAll(genesA);
					union.addAll(genesB);
					double jaccard = (double) intersect / union.size();
					ln += "\t" + jaccard;
				} else {
					ln += "\t" + intersect;
				}
			}
			outf.writeln(ln);
		}
		outf.close();
		
		
	}
	
	public void runPRS(String prsfile, String rowtraitfile, String coltraitfile, boolean usejaccard, String outfilename) throws IOException {
		QTLTextFile tf = new QTLTextFile(prsfile, TextFile.R);
		EQTL[] eprs = tf.read();
		tf.close();
		
		HashMap<String, HashSet<String>> genespertrait = new HashMap<>();
		for (EQTL e : eprs) {
			String prs = e.getRsName();
			int lastunderscore = prs.lastIndexOf("_");
			String prsstripped = prs.substring(0, lastunderscore);
			HashSet<String> genes = genespertrait.get(prsstripped);
			if (genes == null) {
				genes = new HashSet<>();
			}
			genes.add(e.getProbe());
			genespertrait.put(prsstripped, genes);
		}
		
		ArrayList<String> traitAr = new ArrayList<>();
		traitAr.addAll(genespertrait.keySet());
		Collections.sort(traitAr);
		
		ArrayList<String> rowtraits;
		ArrayList<String> coltraits;
		if (rowtraitfile != null) {
			rowtraits = determinetraitlist(rowtraitfile, traitAr);
		} else {
			rowtraits = traitAr;
		}
		if (coltraitfile != null) {
			coltraits = determinetraitlist(coltraitfile, traitAr);
		} else {
			coltraits = traitAr;
		}
		
		// make header
		String header0 = "-\t-";
		for (String s : coltraits) {
			header0 += "\t" + s;
		}
		String header1 = "TraitA\tNumberOfGenes";
		for (String s : coltraits) {
			HashSet<String> genesA = genespertrait.get(s);
			header1 += "\t" + genesA.size();
		}
		
		TextFile outf = new TextFile(outfilename, TextFile.W);
		
		outf.writeln(header0);
		outf.writeln(header1);
		
		for (int i = 0; i < rowtraits.size(); i++) {
			String traitA = rowtraits.get(i);
			HashSet<String> genesA = genespertrait.get(traitA);
			String ln = traitA + "\t" + genesA.size();
			for (int j = 0; j < coltraits.size(); j++) {
				String traitB = coltraits.get(j);
				HashSet<String> genesB = genespertrait.get(traitB);
				int intersect = 0;
				for (String e : genesA) {
					if (genesB.contains(e)) {
						intersect++;
					}
				}
				if (usejaccard) {
					HashSet<String> union = new HashSet<>();
					union.addAll(genesA);
					union.addAll(genesB);
					double jaccard = (double) intersect / union.size();
					ln += "\t" + jaccard;
				} else {
					ln += "\t" + intersect;
				}
			}
			outf.writeln(ln);
		}
		outf.close();
		
		
	}
	
	public void runTrans(String transfx, String traitfile, String rowtraitfile, String coltraitfile, boolean usejaccard, String outfilename) throws IOException {
		
		TextFile tf2 = new TextFile(traitfile, TextFile.R);
		tf2.readLine();
		HashMap<String, HashSet<String>> snptotrait = new HashMap<>();
		String[] elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			String traits = elems[5];
			String[] traitelems = traits.split("; ");
			
			HashSet<String> snpset = snptotrait.get(snp);
			if (snpset == null) {
				snpset = new HashSet<>();
			}
			for (String trait : traitelems) {
				if (trait.trim().length() > 0) {
					snpset.add(trait);
				}
			}
			snptotrait.put(snp, snpset);
			
			
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		System.out.println(snptotrait.size() + " trait snps ");
		
		// parse transfx
		QTLTextFile t = new QTLTextFile(transfx, TextFile.R);
		EQTL[] eqtls = t.read();
		t.close();
		
		HashMap<String, ArrayList<EQTL>> eqtlspertrait = new HashMap<String, ArrayList<EQTL>>();
		HashSet<String> selectedTraits = new HashSet<>();
		HashSet<String> genes = new HashSet<String>();
		for (EQTL e : eqtls) {
			HashSet<String> traits = snptotrait.get(e.getRsName());
			genes.add(e.getProbe());
			if (traits != null) {
				for (String trait : traits) {
					ArrayList<EQTL> traiteqtls = eqtlspertrait.get(trait);
					if (traiteqtls == null) {
						traiteqtls = new ArrayList<>();
					}
					traiteqtls.add(e);
					selectedTraits.add(trait);
					eqtlspertrait.put(trait, traiteqtls);
				}
			}
		}
		System.out.println(eqtlspertrait.size() + " traits with transfx");
		System.out.println(genes.size() + " total genes");
		ArrayList<String> traitAr = new ArrayList<>();
		traitAr.addAll(selectedTraits);
		Collections.sort(traitAr);
		
		ArrayList<String> rowtraits;
		ArrayList<String> coltraits;
		if (rowtraitfile != null) {
			rowtraits = determinetraitlist(rowtraitfile, traitAr);
		} else {
			rowtraits = traitAr;
		}
		if (coltraitfile != null) {
			coltraits = determinetraitlist(coltraitfile, traitAr);
		} else {
			coltraits = traitAr;
		}
		
		
		// make header
		String header0 = "-\t-";
		for (String s : coltraits) {
			header0 += "\t" + s;
		}
		String header1 = "TraitA\tNumberOfGenes";
		for (String s : coltraits) {
			HashSet<String> genesA = getGenes(eqtlspertrait.get(s));
			header1 += "\t" + genesA.size();
		}
		
		TextFile outf = new TextFile(outfilename, TextFile.W);
		
		outf.writeln(header0);
		outf.writeln(header1);
		
		for (int i = 0; i < rowtraits.size(); i++) {
			String traitA = rowtraits.get(i);
			ArrayList<EQTL> eqtlsfortraitA = eqtlspertrait.get(traitA);
			HashSet<String> genesA = getGenes(eqtlsfortraitA);
//			HashMap<String, EQTL> eqtlhash = hasheqtl(eqtlsfortraitA);
			String ln = traitA + "\t" + genesA.size();
			for (int j = 0; j < coltraits.size(); j++) {
				String traitB = coltraits.get(j);
				ArrayList<EQTL> eqtlsfortraitB = eqtlspertrait.get(traitB);
				HashSet<String> genesB = getGenes(eqtlsfortraitB);


//				HashMap<String, EQTL> eqtlhash = hasheqtl(eqtlsfortraitB);
				int intersect = 0;
				for (String e : genesA) {
					if (genesB.contains(e)) {
						intersect++;
					}
				}
				if (usejaccard) {
					HashSet<String> union = new HashSet<>();
					union.addAll(genesA);
					union.addAll(genesB);
					double jaccard = (double) intersect / union.size();
					ln += "\t" + jaccard;
				} else {
					ln += "\t" + intersect;
				}
			}
			outf.writeln(ln);
		}
		outf.close();
		
	}
	
	private HashSet<String> getGenes(ArrayList<EQTL> eqtlsfortraitB) {
		HashSet<String> genes = new HashSet<>();
		
		for (EQTL e : eqtlsfortraitB) {
			genes.add(e.getProbe());
		}
		return genes;
		
	}
	
	private ArrayList<String> determinetraitlist(String rowtraitfile, ArrayList<String> traitAr) throws IOException {
		TextFile tf = new TextFile(rowtraitfile, TextFile.R);
		HashSet<String> set = new HashSet<String>();
		set.addAll(tf.readAsArrayList());
		tf.close();
		
		ArrayList<String> output = new ArrayList<>();
		for (String s : traitAr) {
			if (set.contains(s)) {
				output.add(s);
			}
		}
		return output;
	}
	
	private HashMap<String, EQTL> hasheqtl(ArrayList<EQTL> eqtlsfortraitA) {
		HashMap<String, EQTL> hash = new HashMap<>();
		for (EQTL e : eqtlsfortraitA) {
			hash.put(e.getRsName() + "_" + e.getProbe(), e);
		}
		return hash;
	}
	
}
