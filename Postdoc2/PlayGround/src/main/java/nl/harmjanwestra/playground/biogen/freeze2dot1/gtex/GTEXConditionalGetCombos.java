package nl.harmjanwestra.playground.biogen.freeze2dot1.gtex;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.Gene;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public class GTEXConditionalGetCombos {

	public static void main(String[] args) {
		String indir = "U:\\2020-GTExV8\\GTEx_Analysis_v8_QTLs\\GTEx_Analysis_v8_eQTL_independent\\";
		GTEXConditionalGetCombos gc = new GTEXConditionalGetCombos();
		String output = "U:\\2020-GTExV8\\GTEx_Analysis_v8_QTLs\\GTEx_Analysis_v8_eQTL_independent-combosWoBrain.txt";
		String filelistfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-06-02-SNPIds\\lsof.txt";
		String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz";
		try {
//			gc.run(indir, output);
			String outputrewrite = "U:\\2020-GTExV8\\GTEx_Analysis_v8_QTLs\\GTEx_Analysis_v8_eQTL_independent-combosWoBrain-MetaBrainDot2IDs.txt";
			gc.rewriteIds(filelistfile, gtf, 2, 3, output, outputrewrite);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String indir, String outfile) throws IOException {
		String[] files = Gpio.getListOfFiles(indir);
		TextFile out = new TextFile(outfile, TextFile.W);
		out.writeln("Tissue\tRank\tSNP\tGene");
		for (String file : files) {
			if (file.endsWith("v8.independent_eqtls.txt.gz")) {
				String tissue = file.replaceAll(".v8.independent_eqtls.txt.gz", "");
				if (!tissue.startsWith("Brain")) {
					TextFile tf = new TextFile(indir + file, TextFile.R);
					tf.readLine();
					String[] elems = tf.readLineElems(TextFile.tab);
					while (elems != null) {
						String gene = elems[0];
						String[] snp = elems[6].split("_");
						String rank = elems[elems.length - 1];
						String outsnp = snp[0].replaceAll("chr", "") + ":" + snp[1];
						String outln = tissue + "\t" + rank + "\t" + outsnp + "\t" + gene;
						out.writeln(outln);
						elems = tf.readLineElems(TextFile.tab);
					}
					tf.close();
				}
			}
		}
		out.close();
	}

	public void rewriteIds(String filelistfile, String gtf, int snpcol, int genecol, String input, String output) throws IOException {

		GTFAnnotation g = new GTFAnnotation(gtf);
		Collection<Gene> genes = g.getGenes();
		HashMap<String, String> genemap = new HashMap<String, String>();
		for (Gene gene : genes) {
			String id = gene.getName().split("\\.")[0];
			genemap.put(id, gene.getName());
		}

		HashMap<String, String> rsToId = new HashMap<String, String>();
		if (snpcol > -1) {
			TextFile tf = new TextFile(input, TextFile.R);
			String[] elems = tf.readLineElems(TextFile.tab);
			int ctr = 0;
			while (elems != null) {
				String snpid = elems[snpcol];
				rsToId.put(snpid, "NotInDs-" + snpid);
				elems = tf.readLineElems(TextFile.tab);
				ctr++;
				if (ctr % 100000 == 0) {
					System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
				}
			}
			tf.close();
			System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
			System.out.println();
			TextFile tf1 = new TextFile(filelistfile, TextFile.R);
			ArrayList<String> list = tf1.readAsArrayList();
			tf1.close();


			for (String s : list) {
				TextFile tf2 = new TextFile(s, TextFile.R);

				String ln = tf2.readLine();

				while (ln != null) {
					elems = ln.split(":");
					String query = elems[0] + ":" + elems[1];
					if (rsToId.containsKey(query)) {
						rsToId.put(query, ln);
					}
					ln = tf2.readLine();
				}
				tf2.close();
				System.out.println(rsToId.size() + " ids after reading: " + s);
			}

			tf = new TextFile(input, TextFile.R);
			TextFile out = new TextFile(output, TextFile.W);
			out.writeln(tf.readLine());
			elems = tf.readLineElems(TextFile.tab);

			while (elems != null) {
				String snpid = elems[snpcol];
				String gene = elems[genecol];
				String id = gene.split("\\.")[0];
				elems[genecol] = genemap.get(id);
				String replacement = rsToId.get(snpid);
				elems[snpcol] = replacement;

				out.writeln(Strings.concat(elems, Strings.tab));
				elems = tf.readLineElems(TextFile.tab);
				ctr++;
				if (ctr % 100000 == 0) {
					System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
				}
			}
			out.close();
			tf.close();

		}
	}
}
