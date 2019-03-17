package nl.harmjanwestra.playground.eprs;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.*;

public class EPRSConvergence {
	
	
	public static void main(String[] args) {
		
		
		String prsfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\eprs\\eQTLsFDR-Significant-0.05.txt.gz";
		String prsqcfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-22-eprs-regionsoverlappinggenes.txt";
		String prssnpfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-22-overlapoutput.txtsnpsPerPRS.txt";
		String transeqtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\trans\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
		String ciseqtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz";
		String gwasCatalog = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\genetic_risk_factors_cleaned_traits_standardized_tested_in_trans_20170908.txt";
		String outputfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-23-prsOverlapWithCisAndTransFX.txt";
		EPRSConvergence c = new EPRSConvergence();
		
		try {
			c.combineEPRSWithCisAndTransFX(prsfile, prsqcfile, prssnpfile, transeqtlfile, ciseqtlfile, gwasCatalog, outputfile);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void combineEPRSWithCisAndTransFX(String prsfile, String prsqcfile, String prssnpfile, String transeqtlfile, String ciseqtlfile, String gwasCatalog, String outputfile) throws IOException {
		HashSet<String> prsGenes = new HashSet<>();
		TextFile tf0 = new TextFile(prsqcfile, TextFile.R);
		tf0.readLine();
		String[] prselems = tf0.readLineElems(TextFile.tab);
		
		while (prselems != null) {
			String gene = prselems[1];
			prsGenes.add(gene);
			
			prselems = tf0.readLineElems(TextFile.tab);
		}
		tf0.close();
		
		System.out.println(prsGenes.size() + " prs genes ");
		// PRS --> snps
		HashMap<String, HashSet<String>> prssnps = new HashMap<>();
		TextFile tfp = new TextFile(prssnpfile, TextFile.R);
		tfp.readLine();
		String[] tfpelems = tfp.readLineElems(TextFile.tab);
		while (tfpelems != null) {
			
			String prs = tfpelems[0];
			HashSet<String> set = new HashSet<String>();
			set.addAll(Arrays.asList(tfpelems[2].split(";")));
			prssnps.put(prs, set);
			tfpelems = tfp.readLineElems(TextFile.tab);
		}
		tfp.close();
		
		System.out.println(prssnps.size() + " sets of prssnps");
		
		// diseases --> snps
		HashMap<String, HashSet<String>> snpToDisease = new HashMap<>();
		TextFile tf = new TextFile(gwasCatalog, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> diseaseSNPs = new HashSet<>();
		while (elems != null) {
			String snp = elems[0];
			String disease = elems[5];
			diseaseSNPs.add(snp);
			HashSet<String> diseases = snpToDisease.get(snp);
			if (diseases == null) {
				diseases = new HashSet<>();
			}
			diseases.add(disease);
			snpToDisease.put(snp, diseases);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(snpToDisease.size() + " disease snps");
		
		HashMap<String, ArrayList<EQTL>> transeqtls = loadGeneToEQTL(transeqtlfile, null, null);
		System.out.println(transeqtls.size() + " trans eqtl genes");
		HashMap<String, ArrayList<EQTL>> prsqtls = loadGeneToEQTL(prsfile, null, null);
		System.out.println(prsqtls.size() + " prs genes");
		HashMap<String, ArrayList<EQTL>> ciseqtls = loadGeneToEQTL(ciseqtlfile, null, prsGenes);
		System.out.println(ciseqtls.size() + " cis eqtl genes");
		
		
		// TODO: if there's overlap with cis or trans-fx, what is the LD with the index snps?
		
		// parse prs once more
		TextFile tfout = new TextFile(outputfile, TextFile.W);
		String header = "PRS\tGene\tPRSZ\tPRSCisRegionSNP\tCisSNPs\tCisZ\tCisDisease\tPRSTransSNPs\tTransSNPs\tTransZ\tTransDisease";
		tfout.writeln(header);
		tf0 = new TextFile(prsqcfile, TextFile.R);
		tf0.readLineElems(TextFile.tab);
		prselems = tf0.readLineElems(TextFile.tab);
		while (prselems != null) {
			// PRS	Gene	GeneCoordinates	#PRSSNPs	#SNPRegionsOverlap	TssDistances	OverlappingSNPPos	OverlappingSNPIds
			
			String prs = prselems[0];
			String gene = prselems[1];
			
			String overlapsnp = "";
			if (prselems.length >= 8) {
				overlapsnp = prselems[7];
			} else {
				System.out.println(prselems.length);
			}
			String[] overlapsnps = overlapsnp.split(";");
			HashSet<String> overlapsnphash = new HashSet<>();
			
			overlapsnphash.addAll(Arrays.asList(overlapsnps));
			
			ArrayList<EQTL> trans = transeqtls.get(gene);
			if (trans == null) {
				trans = new ArrayList<>();
			}
			ArrayList<EQTL> cis = ciseqtls.get(gene);
			if (cis == null) {
				cis = new ArrayList<>();
			}
			ArrayList<EQTL> eprs = prsqtls.get(gene);
			
			// link prs gene to prs
			EQTL prsobj = null;
			for (EQTL e : eprs) {
				if (e.getRsName().equals(prs)) {
					prsobj = e;
				}
			}
			
			// link prs gene to trans-eqtls
			HashSet<String> snpsinprs = prssnps.get(prs);
			String[] transsnps = new String[trans.size()];
			String[] transZ = new String[trans.size()];
			String[] transDisease = new String[trans.size()];
			ArrayList<String> transsnpsinprs = new ArrayList<>();
			for (int i = 0; i < trans.size(); i++) {
				transsnps[i] = trans.get(i).getRsName();
				if (snpsinprs.contains(transsnps[i])) {
					transsnpsinprs.add(transsnps[i]);
				}
				transZ[i] = "" + trans.get(i).getZscore();
				HashSet<String> disease = snpToDisease.get(transsnps[i]);
				transDisease[i] = Strings.concat(disease.toArray(new String[0]), Strings.semicolon);
			}
			
			// is the overlapping snp also a cis-fx?
			ArrayList<EQTL> cisfx = new ArrayList<>();
			for (int i = 0; i < cis.size(); i++) {
				EQTL c = cis.get(i);
				if (overlapsnphash.contains(c.getRsName())) {
					cisfx.add(c);
				}
			}
			
			// link prs gene to trans-eqtls
			String[] cissnps = new String[cisfx.size()];
			String[] cisZ = new String[cisfx.size()];
			String[] cisDisease = new String[cisfx.size()];
			for (int i = 0; i < cisfx.size(); i++) {
				cissnps[i] = cisfx.get(i).getRsName();
				cisZ[i] = "" + cisfx.get(i).getZscore();
				HashSet<String> disease = snpToDisease.get(cissnps[i]);
				if (disease != null) {
					cisDisease[i] = Strings.concat(disease.toArray(new String[0]), Strings.semicolon);
				} else {
					cisDisease[i] = "-";
				}
			}
			
			// PRS	Gene	GeneCoordinates	#PRSSNPs	#SNPRegionsOverlap	TssDistances	OverlappingSNPPos	OverlappingSNPIds
			
			String output = prs
					+ "\t" + gene
					+ "\t" + prsobj.getZscore()
					+ "\t" + overlapsnp
					+ "\t" + Strings.concat(cissnps, Strings.semicolon)
					+ "\t" + Strings.concat(cisZ, Strings.semicolon)
					+ "\t" + Strings.concat(cisDisease, Strings.semicolon)
					+ "\t" + Strings.concat(transsnpsinprs, Strings.semicolon)
					+ "\t" + Strings.concat(transsnps, Strings.semicolon)
					+ "\t" + Strings.concat(transZ, Strings.semicolon)
					+ "\t" + Strings.concat(transDisease, Strings.semicolon);
			
			tfout.writeln(output);
			prselems = tf0.readLineElems(TextFile.tab);
		}
		tfout.close();
		tf0.close();
		
		
	}
	
	private HashMap<String, ArrayList<EQTL>> loadGeneToEQTL(String transeqtlfile, HashSet<String> snplimit, HashSet<String> genelimit) throws IOException {
		HashMap<String, ArrayList<EQTL>> output = new HashMap<>();
		QTLTextFile tf2 = new QTLTextFile(transeqtlfile, TextFile.R);
		Iterator<EQTL> transiterator = tf2.getEQtlIterator();
		
		while (transiterator.hasNext()) {
			EQTL e = transiterator.next();
			
			String snp = e.getRsName();
			String gene = e.getProbe();
			if (snplimit == null || snplimit.contains(snp)) {
				if (genelimit == null || genelimit.contains(gene)) {
					ArrayList<EQTL> eqtls = output.get(gene);
					if (eqtls == null) {
						eqtls = new ArrayList<>();
					}
					eqtls.add(e);
					output.put(gene, eqtls);
				}
			}
		}
		tf2.close();
		return output;
	}
	
	
}
