package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class EnrichmentForLargeLDBlocks {
	public static void main(String[] args) {
	
	}
	
	public void run(String in) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		HashMap<String, Integer> snpspergene = new HashMap<String, Integer>();
		HashSet<String> significantGenes = new HashSet<String>();
		
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			String gene = elems[4];
			Integer ctr = snpspergene.get(gene);
			if (ctr == null) {
				ctr = 0;
			}
			ctr++;
			
			snpspergene.put(gene, ctr);
			
			double fdr = Double.parseDouble(elems[elems.length-1]);
			if(fdr<0.05){
				significantGenes.add(gene);
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		
		// 
	}
}
