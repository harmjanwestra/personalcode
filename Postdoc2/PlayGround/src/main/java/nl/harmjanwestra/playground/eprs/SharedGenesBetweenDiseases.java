package nl.harmjanwestra.playground.eprs;

import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class SharedGenesBetweenDiseases {
	
	public static void main(String[] arg) {
		
		String prsfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\eprs\\eQTLsFDR-Significant-0.05.txt.gz";
		String prstraitfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\prstraits\\2018-05-23-immuneprs.txt";
		String ensembl = "D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
		String outputfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-eprsqc\\2018-05-23-immunetrait-genecomp.txt";
		SharedGenesBetweenDiseases g = new SharedGenesBetweenDiseases();
		try {
			g.run(prsfile, prstraitfile, ensembl, outputfile);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String prsfile, String prstraitfile, String ensembl, String outputfile) throws IOException {
		
		GTFAnnotation annot = new GTFAnnotation(ensembl);
		
		HashSet<String> selectedtraits = new HashSet<>();
		TextFile tf1 = new TextFile(prstraitfile, TextFile.R);
		selectedtraits.addAll(tf1.readAsArrayList());
		tf1.close();
		
		QTLTextFile tf = new QTLTextFile(prsfile, TextFile.R);
		EQTL[] eprs = tf.read();
		tf.close();
		
		HashMap<String, HashSet<String>> prsgenespertrait = new HashMap<>();
		HashMap<String, ArrayList<EQTL>> eprspertrait = new HashMap<>();
		
		for (EQTL e : eprs) {
			String prs = e.getRsName();
			int lastunderscore = prs.lastIndexOf("_");
			String prsstripped = prs.substring(0, lastunderscore);
			
			if (selectedtraits.contains(prsstripped)) {
				int index = prs.indexOf("_");
				String traitname = prsstripped.substring(0, index);
				HashSet<String> genes = prsgenespertrait.get(traitname);
				if (genes == null) {
					genes = new HashSet<>();
				}
				genes.add(e.getProbe());
				prsgenespertrait.put(traitname, genes);
				
				ArrayList<EQTL> eqtls = eprspertrait.get(prsstripped);
				if (eqtls == null) {
					eqtls = new ArrayList<>();
				}
				eqtls.add(e);
				eprspertrait.put(prsstripped, eqtls);
			}
		}
		
		ArrayList<String> alldisease = new ArrayList<String>();
		alldisease.addAll(prsgenespertrait.keySet());
		Collections.sort(alldisease);
		
		System.out.println(alldisease.size() + " traits");
		HashSet<String> allgenes = new HashSet<>();
		for (String key : alldisease) {
			System.out.println(key + "\t" + prsgenespertrait.get(key).size());
			allgenes.addAll(prsgenespertrait.get(key));
		}
		
		System.out.println();
		System.out.println("unique genes: " + allgenes.size());
		System.out.println();
		
		for (int i = 0; i < alldisease.size(); i++) {
			String diseasA = alldisease.get(i);
			HashSet<String> setA = prsgenespertrait.get(diseasA);
			for (int j = i + 1; j < alldisease.size(); j++) {
				String diseasB = alldisease.get(j);
				HashSet<String> setB = prsgenespertrait.get(diseasB);
				int intersect = 0;
				ArrayList<String> shared = new ArrayList<>();
				for (String s : setA) {
					if (setB.contains(s)) {
						intersect++;
						Gene g = annot.getStrToGene().get(s);
						shared.add(g.getGeneSymbol());
					}
				}
				if (intersect > 0) {
					System.out.println(diseasA + "\t" + diseasB + "\t" + intersect + "\t" + Strings.concat(shared, Strings.semicolon));
				}
			}
		}
		
		// output directions per gene
		EQTL[][] zscoremat = new EQTL[allgenes.size()][eprspertrait.size()];
		ArrayList<String> genelist = new ArrayList<>();
		genelist.addAll(allgenes);
		Collections.sort(genelist);
		HashMap<String, Integer> geneindex = new HashMap<>();
		for (int i = 0; i < genelist.size(); i++) {
			geneindex.put(genelist.get(i), i);
		}
		
		ArrayList<String> traitlist = new ArrayList<>();
		traitlist.addAll(eprspertrait.keySet());
		Collections.sort(traitlist);
		HashMap<String, Integer> traitindex = new HashMap<>();
		for (int i = 0; i < traitlist.size(); i++) {
			traitindex.put(traitlist.get(i), i);
		}
		
		for (String key : eprspertrait.keySet()) {
			Integer id = traitindex.get(key);
			ArrayList<EQTL> eqtls = eprspertrait.get(key);
			for (EQTL e : eqtls) {
				String gene = e.getProbe();
				Integer geneid = geneindex.get(gene);
				zscoremat[geneid][id] = e;
			}
		}
		
		TextFile output = new TextFile(outputfile, TextFile.W);
		String header = "-";
		for (String d : traitlist) {
			header += "\t" + d;
		}
		output.writeln(header);
		for (int i = 0; i < genelist.size(); i++) {
			String gene = genelist.get(i);
			Gene g = annot.getStrToGene().get(gene);
			String ln = g.getGeneSymbol();
			EQTL ref = null;
			for (int j = 0; j < traitlist.size(); j++) {
				EQTL e = zscoremat[i][j];
				if (e != null) {
					if (ref == null) {
						ref = e;
						ln += "\t" + e.getZscore();
					} else {
						Boolean flip = BaseAnnot.flipalleles(ref.getAlleles(), ref.getAlleleAssessed(), e.getAlleles(), e.getAlleleAssessed());
						if (flip) {
							ln += "\t" + (-e.getZscore());
						} else {
							ln += "\t" + e.getZscore();
						}
					}
					
				} else {
					ln += "\t-";
				}
			}
			output.writeln(ln);
		}
		output.close();
	}
	
}
