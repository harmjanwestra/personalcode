package nl.harmjanwestra.playground.biogen.locusdebug;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;

public class FindGenesOnMultipleChr {

	public static void main(String[] args) {
		String query = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\2019-09-26-GenesOfInterest.txt";
		String table = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\2019-09-26-patch_genes.txt";

		FindGenesOnMultipleChr c = new FindGenesOnMultipleChr();

		String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\annotation\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";
		try {
			c.run(query, table, gtf);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String query, String table, String gtf) throws IOException {
		HashSet<String> querygenes = new HashSet<String>();
		TextFile tf = new TextFile(query, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			querygenes.add(elems[0]);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println("Overlapping:");
		TextFile tf2 = new TextFile(table, TextFile.R);
		tf2.readLine();
		String[] elems2 = tf2.readLineElems(TextFile.tab);
		HashSet<String> affectedENSGs = new HashSet<String>();
		while (elems2 != null) {
			String gene = elems2[0];
			for (int q = 1; q < elems2.length; q++) {
				if (!elems2[q].equals("-")) {
					String ensg = elems2[q].split("\\.")[0];
					affectedENSGs.add(ensg);
				}

			}
			if (querygenes.contains(gene)) {
				System.out.println(gene);
			}
			elems2 = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		System.out.println();
		GTFAnnotation annot = new GTFAnnotation(gtf);
		Collection<Gene> genes = annot.getGenes();
		HashSet<String> foundGenes = new HashSet<>();
		System.out.println("Found genes");

		for (Gene g : genes) {
			if (querygenes.contains(g.getGeneSymbol())) {
				foundGenes.add(g.getGeneSymbol());
				String ensg = g.getName().split("\\.")[0];
				if (affectedENSGs.contains(ensg)) {
					System.out.println("AFFECTED: " + g.getGeneSymbol() + "\t" + g.getName());
				} else {
					System.out.println("NOT AFFECTED: " + g.getGeneSymbol() + "\t" + g.getName());
				}
			}
		}

		System.out.println();
		System.out.println("Not found:");
		for (String g : querygenes) {
			if (!foundGenes.contains(g)) {
				System.out.println(g);
			}
		}
	}
}
