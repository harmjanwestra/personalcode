package nl.harmjanwestra.playground.transeqtl;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.FeatureComparator;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class SNPProbePairSorter {
	
	public static void main(String[] args) {
		
		
		String snpannot = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz";
		String gtf = "D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
		String snpprobe = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\allsnpprobecombos.txt.gz";
		String probesnp = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\sortedGeneSNPCombos.txt.gz";
		
		SNPProbePairSorter p = new SNPProbePairSorter();
		try {
			p.run(snpprobe, gtf, snpannot, probesnp);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	public void run(String snpprobe, String gtf, String snpannot, String out) throws IOException {
		
		HashSet<String> snps = new HashSet<>();
		
		HashMap<String, ArrayList<String>> snpspergene = new HashMap<>();
		TextFile tf = new TextFile(snpprobe, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		
		int lnctr = 0;
		while (elems != null) {
			if (elems.length == 2) {
				String gene = Strings.cache(elems[1]);
				String snp = Strings.cache(elems[0]);
				
				snps.add(snp);
				ArrayList<String> set = snpspergene.get(gene);
				if (set == null) {
					set = new ArrayList<>();
				}
				set.add(snp);
				snpspergene.put(gene, set);
			}
			lnctr++;
			if (lnctr % 100 == 0) {
				System.out.print("\r"+lnctr + "\tsnpprobe lines read. "+snpspergene.size()+" genes sofar.");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(snpspergene.size() + " genes");
		
		
		ArrayList<String> sortedGenes = new ArrayList<>(snpspergene.keySet().size());
		{
			GTFAnnotation g = new GTFAnnotation(gtf);
			ArrayList<Gene> geneArrayList = new ArrayList<>();
			for (String gene : snpspergene.keySet()) {
				Gene gobj = g.getStrToGene().get(gene);
				if (gobj != null) {
					geneArrayList.add(gobj);
				}
			}
			
			Collections.sort(geneArrayList, new FeatureComparator(false));
			for (Gene gobj : geneArrayList) {
				sortedGenes.add(gobj.getName());
			}
		}
		
		System.out.println(sortedGenes.size() + " sorted genes");

//		ArrayList<String> sortedSNPs = new ArrayList<>(snps.size());
		HashMap<String, Integer> snpToInt = new HashMap<>();
		{
			TextFile tf2 = new TextFile(snpannot, TextFile.R);
			elems = tf2.readLineElems(TextFile.tab);
			int ctr = 0;
			lnctr = 0;
			while (elems != null) {
				if (snps.contains(elems[2])) {
					snpToInt.put(Strings.cache(elems[2]), ctr);
					ctr++;
				}
			
				lnctr++;
				if (lnctr % 10000 == 0) {
					System.out.print(lnctr + "\tsnp lines read\r");
				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
		}
		System.out.println(snpToInt.size() + " indexed snps");
		
		
		class SNPIndexSortObj implements Comparable<SNPIndexSortObj> {
			
			int id;
			String snp;
			
			@Override
			public boolean equals(Object o) {
				if (this == o) return true;
				if (o == null || getClass() != o.getClass()) return false;
				
				SNPIndexSortObj that = (SNPIndexSortObj) o;
				
				return id == that.id;
			}
			
			@Override
			public int hashCode() {
				return id;
			}
			
			SNPIndexSortObj(Integer id, String snp) {
				this.id = id;
				this.snp = snp;
			}
			
			
			@Override
			public int compareTo(SNPIndexSortObj o) {
				if (this.equals(o)) {
					return 0;
				} else if (this.id > o.id) {
					return 1;
				} else {
					return -1;
				}
			}
		}
		
		// now we kindof know the sorting for snp probe combos?
		TextFile outf = new TextFile(out, TextFile.W);
		TextFile outf2 = new TextFile(out+"-genes.txt", TextFile.W);
		
		for (String gene : sortedGenes) {
			ArrayList<String> set = snpspergene.get(gene);
			ArrayList<SNPIndexSortObj> sort = new ArrayList<>();
			for (String snp : set) {
				Integer id = snpToInt.get(snp);
				if (id != null) {
					sort.add(new SNPIndexSortObj(id, snp));
				}
			}
			Collections.sort(sort);
			for (SNPIndexSortObj s : sort) {
				outf.writeln(gene + "\t" + s.snp);
			}
			outf2.writeln(gene);
		}
		outf.close();
		outf2.close();
		
		
	}
}
