package nl.harmjanwestra.playground.biogen.locusdebug;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class DetermineNonWantedSequences {


	public static void main(String[] args) {
		String patchgenes = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\2019-09-26-patch_genes.txt";
		String patchgenesout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-GeneQueries\\2019-09-26-patch_genes-problematicsetcomparedtogtex.txt";
		String gtexdict = "D:\\Sync\\SyncThing\\Data\\HumanGenome\\Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict";

		DetermineNonWantedSequences d = new DetermineNonWantedSequences();
		try {
			d.run(gtexdict, patchgenes, patchgenesout);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public void run(String trimmeddict, String table, String output) throws IOException {
		TextFile tf = new TextFile(trimmeddict, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);

		HashSet<String> allowedNames = new HashSet<String>();
		while (elems != null) {
			if (elems[0].startsWith("@SQ")) {
				String name = elems[1];
				name = name.split(":")[1];
				name = name.split("v")[0];
				String[] nameelems = name.split("_");
				if (nameelems.length > 1) {
					name = nameelems[1];
				}

				allowedNames.add(name);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(table, TextFile.R);
		String[] header = tf2.readLineElems(TextFile.tab);
		boolean[] seqok = new boolean[header.length];
		for (int i = 0; i < header.length; i++) {
			String seqname = header[i].split("\\.")[0];
			if (i > 330) {
				System.out.println("Barbapapalekukeleku KI270706");
			}
			if (allowedNames.contains(seqname)) {
				seqok[i] = true;
			} else {
				System.out.println("Could not find: " + seqname);
			}
		}
		elems = tf2.readLineElems(TextFile.tab);
		HashSet<String> problemeticAutosomal = new HashSet<String>();
		int ctr = 0;
		int ctr2 = 0;
		while (elems != null) {
			String gene = elems[0];
			ArrayList<String> geneIdsOnAllowedChr = new ArrayList<String>();
			boolean problematicgene = false;
			for (int i = 1; i < elems.length; i++) {
				if (!elems[i].equals("-")) {
					if (seqok[i]) {
						String[] geneelems = elems[i].split(";");
						for (String g : geneelems) {
							geneIdsOnAllowedChr.add(g);
						}

					} else {
						problematicgene = true;
					}

				}
			}
			if (problematicgene) {
				ctr++;
				problemeticAutosomal.addAll(geneIdsOnAllowedChr);
			}
			ctr2++;
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();

		TextFile tfout = new TextFile(output, TextFile.W);
		for (String s : problemeticAutosomal) {
			tfout.writeln(s);
		}
		tfout.close();

		System.out.println(ctr + " / " + ctr2 + " problematic genes.");

	}
}
