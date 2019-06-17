package nl.harmjanwestra.playground.biogen;

import nl.harmjanwestra.playground.transeqtl.MakeTranscriptome;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class MatchFreeze2SNPIds {


	public static void main(String[] args) {

		String filelistfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-06-02-SNPIds\\lsof.txt";

		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\annotation\\gwas_catalog_v1.0.2-associations_e96_r2019-05-03.SNPs.tsv";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\annotation\\gwas_catalog_v1.0.2-associations_e96_r2019-05-03.SNPs-MetaBrainFreeze2Ids.tsv";

		MatchFreeze2SNPIds s = new MatchFreeze2SNPIds();

		s.snplist(filelistfile, in, out);

	}

	public void snplist(String filelistfile, String input, String output) {
		String gwascatalog = input; // "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2018-01-31-eqtlgen-cis-snps.txt.gz";
		String gwascatalogout = output; //"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2018-01-31-eqtlgen-cis-snps-MetaBrainFreeze2Ids.txt.gz";
		try {


			TextFile tf = new TextFile(filelistfile, TextFile.R);
			ArrayList<String> list = tf.readAsArrayList();
			tf.close();

			HashMap<String, String> rsToId = new HashMap<String, String>();

			TextFile gwas = new TextFile(gwascatalog, TextFile.R);
			String[] elems = gwas.readLineElems(TextFile.tab);
			while (elems != null) {
				String rs = elems[0];
				rsToId.put(rs, null);
				elems = gwas.readLineElems(TextFile.tab);
			}
			gwas.close();

			for (String s : list) {
				TextFile tf2 = new TextFile(s, TextFile.R);

				String ln = tf2.readLine();

				while (ln != null) {
					elems = ln.split(":");
					if (rsToId.containsKey(elems[2])) {
						rsToId.put(elems[2], ln);
					}
					ln = tf2.readLine();
				}
				tf2.close();
				System.out.println(rsToId.size() + " ids after reading: " + s);
			}


			gwas.open();
			TextFile gwasout = new TextFile(gwascatalogout, TextFile.W);
			elems = gwas.readLineElems(TextFile.tab);
			int found = 0;
			int total = 0;


			while (elems != null) {
				String rs = elems[0];
				if (rsToId.containsKey(rs)) {
					elems[0] = rsToId.get(rs);

					gwasout.writeln(Strings.concat(elems, Strings.tab));
					found++;
				}
				total++;
				elems = gwas.readLineElems(TextFile.tab);
			}
			gwas.close();
			gwasout.close();
			System.out.println(found + "\t" + total + " GWAS snps found.");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
