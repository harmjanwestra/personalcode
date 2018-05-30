package nl.harmjanwestra.playground.transeqtl;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class TraitsWithTrans {
	
	public static void main(String[] args) {
		
		
		String immuneprsfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-23-immuneprs.txt";
		String allprsfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-23-prsOverlapWithCisAndTransFX.txt";
		String transfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\trans\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
		String immunesnps = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-23-immunesnps.txt";
		String gwassnps = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\genetic_risk_factors_cleaned_traits_standardized_tested_in_trans_20170908.txt";
		String nonimmuneoverlapout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-23-nonimmunetraits.txt";
		String immunetraitsout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-23-immunetraits.txt";
		String nonimmuneprstraitout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-23-nonimmuneprstraits.txt";
		
		TraitsWithTrans t = new TraitsWithTrans();
		try {
			t.run(immuneprsfile, allprsfile, transfile, immunesnps, gwassnps, nonimmuneoverlapout, immunetraitsout, nonimmuneprstraitout);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void run(String immuneprsfile, String nonimmuneprsfile, String transfile, String immunesnps, String gwassnps, String nonimmuneoverlapout, String immuneoverlapout, String nonimmuneprstraitout) throws IOException {
		
		// merge prs per gene
		TextFile prsf = new TextFile(immuneprsfile, TextFile.R);
		HashMap<String, HashSet<String>> immuneprsset = new HashMap<>();
		HashSet<String> immuneprstraits = new HashSet<>();
		prsf.readLine();
		String[] elems = prsf.readLineElems(TextFile.tab);
		while (elems != null) {
			String p = elems[0];
			String gene = elems[1];
			HashSet<String> prsss = immuneprsset.get(gene);
			if (prsss == null) {
				prsss = new HashSet<>();
			}
			int lastunderscore = p.lastIndexOf("_");
			String prsstripped = p.substring(0, lastunderscore);
			prsss.add(prsstripped);
			immuneprstraits.add(prsstripped);
			immuneprsset.put(gene, prsss);
			elems = prsf.readLineElems(TextFile.tab);
		}
		prsf.close();
		
		// merge prs per gene
		TextFile prsf2 = new TextFile(nonimmuneprsfile, TextFile.R);
		HashMap<String, HashSet<String>> nonimmuneprsset = new HashMap<>();
		prsf2.readLine();
		String[] elems2 = prsf2.readLineElems(TextFile.tab);
		while (elems2 != null) {
			String p = elems2[0];
			String gene = elems2[1];
			HashSet<String> prsss = nonimmuneprsset.get(gene);
			if (prsss == null) {
				prsss = new HashSet<>();
			}
			int lastunderscore = p.lastIndexOf("_");
			String prsstripped = p.substring(0, lastunderscore);
			prsss.add(prsstripped);
			nonimmuneprsset.put(gene, prsss);
			elems2 = prsf2.readLineElems(TextFile.tab);
		}
		prsf2.close();
		
		// merge genes with immune trans fx
		HashMap<String, HashSet<String>> immunesnpset = new HashMap<>();
		TextFile tf1 = new TextFile(immunesnps, TextFile.R);
		tf1.readLine();
		elems = tf1.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length >= 5) {
				String snp = elems[0];
				String traits = elems[5];
				String[] traitelems = traits.split("; ");
				HashSet<String> snpset = immunesnpset.get(snp);
				if (snpset == null) {
					snpset = new HashSet<>();
				}
				for (String trait : traitelems) {
					snpset.add(trait);
				}
				immunesnpset.put(snp, snpset);
			}
			elems = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();
		System.out.println(immunesnpset.size() + " immune snps ");
		// load other traits
		HashMap<String, HashSet<String>> nonimmunesnpset = new HashMap<>();
		TextFile tf2 = new TextFile(gwassnps, TextFile.R);
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[0];
			String traits = elems[5];
			String[] traitelems = traits.split("; ");
			if (!immunesnpset.containsKey(snp)) {
				HashSet<String> snpset = nonimmunesnpset.get(snp);
				if (snpset == null) {
					snpset = new HashSet<>();
				}
				for (String trait : traitelems) {
					snpset.add(trait);
				}
				nonimmunesnpset.put(snp, snpset);
			}
			
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		// load trans snps per gene
		TextFile tf3 = new TextFile(transfile, TextFile.R);
		tf3.readLine();
		String[] telems = tf3.readLineElems(TextFile.tab);
		HashMap<String, HashSet<String>> prsgenetransfx = new HashMap<>();
		
		while (telems != null) {
			String snp = telems[1];
			String gene = telems[4];
			if (nonimmuneprsset.containsKey(gene)) {
				HashSet<String> transset = prsgenetransfx.get(gene);
				if (transset == null) {
					transset = new HashSet<>();
				}
				transset.add(snp);
				prsgenetransfx.put(gene, transset);
			}
			telems = tf3.readLineElems(TextFile.tab);
		}
		tf3.close();
		
		
		HashSet<String> outimmunetraits = new HashSet<>();
		HashSet<String> outnonimmunetraits = new HashSet<>();
		for (String gene : nonimmuneprsset.keySet()) {
			HashSet<String> transsnps = prsgenetransfx.get(gene);
			if (transsnps != null) {
				for (String snp : transsnps) {
					HashSet<String> nonummunetraits = nonimmunesnpset.get(snp);
					HashSet<String> immunetraits = immunesnpset.get(snp);
					if (nonummunetraits != null) {
						outnonimmunetraits.addAll(nonummunetraits);
					}
					if (immunetraits != null) {
						outimmunetraits.addAll(immunetraits);
					}
				}
			}
		}
		
		
		TextFile traitout = new TextFile(nonimmuneoverlapout, TextFile.W);
		for (String t : outnonimmunetraits) {
			traitout.writeln(t);
		}
		traitout.close();
		traitout = new TextFile(immuneoverlapout, TextFile.W);
		for (String t : outimmunetraits) {
			traitout.writeln(t);
		}
		traitout.close();
		
		HashSet<String> nonimmuneprstraits = new HashSet<>();
		for (String gene : immuneprsset.keySet()) {
			if (nonimmuneprsset.containsKey(gene)) {
				nonimmuneprstraits.addAll(nonimmuneprsset.get(gene));
			}
		}
		traitout = new TextFile(nonimmuneprstraitout, TextFile.W);
		for (String t : nonimmuneprstraits) {
			if (!immuneprstraits.contains(t)) {
				traitout.writeln(t);
			}
		}
		traitout.close();
	}
	
	public void mergePRS(String immuneprsfile, String allprsfile, String bloodcellprsfile) {
		
		// gene immunetrait bloodcelltrait othertrait
	}
}
