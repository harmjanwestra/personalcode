package nl.harmjanwestra.playground;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.enums.Strand;
import umcg.genetica.features.Gene;
import umcg.genetica.features.Transcript;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

public class EnsemblReWrite {
	
	
	public static void main(String[] args) {
		
		try {
			
			TextFile tf2 = new TextFile("D:\\Work\\TMP\\listGenes.txt", TextFile.R);
			ArrayList<String> ensg = tf2.readAsArrayList();
			tf2.close();
			
			GTFAnnotation f = new GTFAnnotation("D:\\Work\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz");
			
			Collection<Gene> genes = f.getGenes();
			TextFile tf = new TextFile("D:\\Work\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71-rewrite.txt.gz", TextFile.W);
			
			
			HashSet<String> s = new HashSet<>();
			for (Gene g : genes) {
				s.add(g.getName());
			}
			
			int ctr = 0;
			HashSet<String> meh = new HashSet<>();
			for (String g : ensg) {
				if (!s.contains(g)) {
					ctr++;
					meh.add(g);
				}
			}
			
			System.out.println(ctr + " eqtms unaccounted for");
			System.out.println(meh.size() + " missing ENSG");
			
			for (Gene g : genes) {
				ArrayList<Transcript> transcripts = g.getTranscripts();
				int minpos = Integer.MAX_VALUE;
				int maxpos = 0;
				for (Transcript t : transcripts) {
					if (t.getStart() < minpos) {
						minpos = t.getStart();
					}
					if (t.getStop() > maxpos) {
						maxpos = t.getStop();
					}
				}
				
				if (g.getStrand().equals(Strand.POS)) {
					tf.writeln(g.getName() + "\t" + minpos);
				} else {
					tf.writeln(g.getName() + "\t" + maxpos);
				}
			}
			
			tf.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
}
