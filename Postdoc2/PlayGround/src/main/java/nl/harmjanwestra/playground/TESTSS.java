package nl.harmjanwestra.playground;

import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Transcript;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;

public class TESTSS {
	public static void main(String[] args) {
		
		try {
			GTFAnnotation g = new GTFAnnotation("D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz");
			TextFile out = new TextFile("d:\\Homo_sapiens.GRCh37.71-tes.txt", TextFile.W);
			for (Gene gene : g.getGenes()) {
				out.writeln(gene.getName() + "\t" + gene.getStart() + "\t" + gene.getStop() + "\t" + gene.getStrand().toString());
			}
			out.close();
			
			out = new TextFile("d:\\Homo_sapiens.GRCh37.71-shortestTranscript-tes.txt", TextFile.W);
			for (Gene gene : g.getGenes()) {
				ArrayList<Transcript> transcripts = gene.getTranscripts();
				Transcript shortest = null;
				int len = Integer.MAX_VALUE;
				for (Transcript t : transcripts) {
					int size = t.getSize();
					if (size < len) {
						shortest = t;
						len = size;
					}
				}
				
				out.writeln(gene.getName() + "\t" + shortest.getName() + "\t" + shortest.getStart() + "\t" + shortest.getStop() + "\t" + shortest.getStrand().toString());
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
