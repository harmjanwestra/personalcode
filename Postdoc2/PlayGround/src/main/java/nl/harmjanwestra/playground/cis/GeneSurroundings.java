package nl.harmjanwestra.playground.cis;


import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.FeatureComparator;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.*;

public class GeneSurroundings {
	
	public static void main(String[] args) {
		
		String genelistfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\sortedGeneSNPCombos.txt.gz-genes.txt";
		String geneAnnotation = "D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
		String querygene = "ENSG00000169047";
		
		GeneSurroundings s = new GeneSurroundings();
		try {
			s.run(null, geneAnnotation, genelistfile, querygene, null);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	
	public void run(String eqtlfile, String geneAnnotation, String testedGeneList, String queryGene, String out) throws IOException {
		
		// which genes are located left of the querygene
		
		TextFile tf = new TextFile(testedGeneList, TextFile.R);
		ArrayList<String> listOfGenes = tf.readAsArrayList();
		HashSet<String> geneHash = new HashSet<String>();
		geneHash.addAll(listOfGenes);
		tf.close();
		
		GTFAnnotation gtf = new GTFAnnotation(geneAnnotation);
		Collection<Gene> allgenes = gtf.getGenes();
		ArrayList<Gene> allgenelist = new ArrayList<>();
		allgenelist.addAll(allgenes);
		Collections.sort(allgenelist, new FeatureComparator(true));
		
		HashMap<String, Integer> geneToInt = new HashMap<>();
		
		int ctr = 0;
		for (Gene g : allgenelist) {
			geneToInt.put(g.getName(), ctr);
			ctr++;
		}
		
		Integer geneIndex = geneToInt.get(queryGene);
		if (geneIndex != null) {
			Gene qg = allgenelist.get(geneIndex);
			int midpoint = (qg.getStop() + qg.getStop()) / 2;
			
			System.out.println(qg.getName() + "\t" + qg.getChromosome() + "\t" + midpoint);
			System.out.println(geneIndex + "\t" + allgenelist.size());
			System.out.println("------");
			
			int indexdown = geneIndex;
			int distance = 0;
			while (indexdown > -1 && distance < 2000000) {
				Gene g = allgenelist.get(indexdown);
				int midpoint2 = (g.getStop() + g.getStop()) / 2;
				distance = Math.abs(midpoint - midpoint2);
				System.out.println(g.getName() + "\t" + g.getChromosome() + "\t" + midpoint2 + "\t" + distance + "\t" + geneHash.contains(g.getName()));
				indexdown--;
			}
			
			
			
			int indexup = geneIndex;
			distance = 0;
			while (indexup < allgenelist.size() && distance < 2000000) {
				Gene g = allgenelist.get(indexup);
				int midpoint2 = (g.getStop() + g.getStop()) / 2;
				distance = Math.abs(midpoint - midpoint2);
				System.out.println(g.getName() + "\t" + midpoint2 + "\t" + distance + "\t" + geneHash.contains(g.getName()));
				indexup++;
			}
			
			
		}
		
		
	}
}
