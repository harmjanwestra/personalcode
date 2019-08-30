import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class SimpliFyGTF {

	public static void main(String[] args) {
		try {
			GTFAnnotation annot = new GTFAnnotation("D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz");
			TextFile tf = new TextFile("D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71-genes.txt.gz", TextFile.W);
			for (Gene g : annot.getGenes()) {
				int midpoint = g.getStart() + ((g.getStop() - g.getStart()) / 2);
				tf.writeln(g.getName() + "\t" + g.getChromosome().getNumber() + "\t" + g.getStart() + "\t" + g.getStop() + "\t" + midpoint);
			}
			tf.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
