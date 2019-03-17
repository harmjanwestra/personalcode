package nl.harmjanwestra.playground.transeqtl;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

public class MakeGeneTranslation {
	
	public static void main(String[] args) {
		String ref = "D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
		String out = "D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71-HGNCIdsPerGene.txt.gz";
		
		MakeGeneTranslation t = new MakeGeneTranslation();
		try {
			t.run(ref, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String in, String out) throws IOException {
		GTFAnnotation a = new GTFAnnotation(in);
		Collection<Gene> genes = a.getGenes();
		
		HashMap<String, HashSet<String>> ensgToHGNC = new HashMap<String, HashSet<String>>();
		for (Gene g : genes) {
			HashSet<String> name = ensgToHGNC.get(g.getName());
			if (g.getGeneSymbol() != null) {
				if (name == null) {
					name = new HashSet<>();
				}
				name.add(g.getGeneSymbol());
				ensgToHGNC.put(g.getName(), name);
			}
		}
		
		
		TextFile outf = new TextFile(out, TextFile.W);
		for (String k : ensgToHGNC.keySet()) {
			HashSet<String> name = ensgToHGNC.get(k);
			String[] names = name.toArray(new String[0]);
			String gout = "-";
			if (names.length >= 0) {
				gout = Strings.concat(names, Strings.semicolon);
			}
			outf.writeln(k + "\t" + gout);
		}
		
		outf.close();
		
		
	}
}
