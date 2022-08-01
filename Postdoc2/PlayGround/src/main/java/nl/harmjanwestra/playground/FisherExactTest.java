package nl.harmjanwestra.playground;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.enums.Strand;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.Collection;

public class FisherExactTest {

	public static void main(String[] args) {
		int ldp = 1980;
		int nonldp = 13915;
		int ldnp = 881;
		int nonldnp = 16219;
//		umcg.genetica.math.stats.FisherExactTest t = new umcg.genetica.math.stats.FisherExactTest();
//		System.out.println(t.getFisherPValue(ldp, nonldp, ldnp, nonldnp));

//		String gtf = "D:\\Sync\\SyncThing\\Data\\Ref\\gencode\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";
		String gtf = "D:\\Homo_sapiens.GRCh37.75.gtf";
		try {
			GTFAnnotation a = new GTFAnnotation(gtf);
			Collection<Gene> genes = a.getGenes();
//			TextFile out = new TextFile("D:\\Sync\\SyncThing\\Data\\Ref\\gencode\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf-genePos.txt.gz", TextFile.W);
			TextFile out = new TextFile("D:\\Homo_sapiens.GRCh37.75.gtf-genepos.txt", TextFile.W);
			out.writeln("Gene\tSymbol\tStrand\tStart\tStop\tmidpoint\tTSS");
			for (Gene g : genes) {
				int midpoint = g.getStart() + ((g.getStop() - g.getStart()) / 2);
				int tss = g.getStart();
				if (g.getStrand().equals(Strand.NEG)) {
					tss = g.getStop();
				}
				out.writeln(g.getName() + "\t" + g.getGeneSymbol() + "\t" + g.getStrand() + "\t" + g.getStart() + "\t" + g.getStop() + "\t" + midpoint + "\t" + tss);
			}
			out.close();

		} catch (IOException e) {
			throw new RuntimeException(e);
		}

//		int cisindextrans = 150;
//		int cistrans = 311;
//		int gwasvariantsindextrans = 15745;
//		int gwasvariantscis = 109172;
//		System.out.println(t.getFisherPValue(cisindextrans, cistrans, gwasvariantsindextrans, gwasvariantscis));


	}
}
