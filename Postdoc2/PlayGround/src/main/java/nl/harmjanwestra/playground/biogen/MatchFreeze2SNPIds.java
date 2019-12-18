package nl.harmjanwestra.playground.biogen;

import nl.harmjanwestra.playground.transeqtl.MakeTranscriptome;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class MatchFreeze2SNPIds {


	public static void main(String[] args) {

		String filelistfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-06-02-SNPIds\\lsof.txt";

		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\annotation\\gwas_catalog_v1.0.2-associations_e96_r2019-05-03.SNPs.tsv";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\annotation\\gwas_catalog_v1.0.2-associations_e96_r2019-05-03.SNPs-MetaBrainFreeze2Ids.tsv";

		in = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e96_r2019-10-14.tsv.gz";
		String snpsout = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e96_r2019-10-14-snps.txt.gz";
		out = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e96_r2019-10-14-snps-MetaBrainFreeze2Ids.txt.gz";

//		try {
//			GWASCatalog c = new GWASCatalog(in);
//			HashSet<GWASSNP> snps = c.getSnps();
//			TextFile snpsoutf = new TextFile(snpsout, TextFile.W);
//			for (GWASSNP snp : snps) {
//				snpsoutf.writeln(snp.getName());
//			}
//			snpsoutf.close();
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		MatchFreeze2SNPIds s = new MatchFreeze2SNPIds();

		snpsout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-11-20-LPGAT\\PD_SNPS_2019_Nov.txt";
		out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-11-20-LPGAT\\PD_SNPS_2019_Nov-metabrainfreeze2ids.txt";
//		s.snplist(filelistfile, snpsout, out);

		snpsout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-12-18-BenteSNPs\\Alzheimer.txt";
		out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-12-18-BenteSNPs\\Alzheimer-metabrainfreeze2ids.txt";
		s.snplist(filelistfile, snpsout, out);
		snpsout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-12-18-BenteSNPs\\Depression.txt";
		out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-12-18-BenteSNPs\\Depression-metabrainfreeze2ids.txt";
		s.snplist(filelistfile, snpsout, out);


//		String snpTrait = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\annotation\\gwas_catalog_v1.0.2-associations_e93_r2018-10-29_SNP_TRAIT.txt.gz";
//		String freeze2efile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Trans\\eQTLsFDR0.05.txt.gz";
//		String freeze2Efileout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Trans\\eQTLsFDR0.05-traits.txt.gz";
//		try {
//			s.annotateEQTLFile(snpTrait, freeze2efile, freeze2Efileout);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
	}

	private void annotateEQTLFile(String snpTrait, String freeze2efile, String freeze2Efileout) throws IOException {
		HashMap<String, String> snptotrait = new HashMap<>();
		TextFile tf = new TextFile(snpTrait, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 1) {
				snptotrait.put(elems[0], elems[1]);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tfout = new TextFile(freeze2Efileout, TextFile.W);
		TextFile tfin = new TextFile(freeze2efile, TextFile.R);


		HashMap<String, Integer> snpToCt = new HashMap<>();

		tfout.writeln(tfin.readLine() + "\tTrait");
		String[] elems2 = tfin.readLineElems(TextFile.tab);
		while (elems2 != null) {

			String snp = elems2[1];
			String snprs = snp.split(":")[2];
			String traits = snptotrait.get(snprs);
			if (traits != null) {
				String[] telems = traits.split(",");
				for (String t : telems) {
					while (t.startsWith(" ")) {
						t = t.substring(1);
					}

					Integer ct = snpToCt.get(t);
					if (ct == null) {
						ct = 0;
					}
					ct++;
					snpToCt.put(t, ct);
				}
			}
			tfout.writeln(Strings.concat(elems2, Strings.tab) + "\t" + snptotrait.get(snprs));

			elems2 = tfin.readLineElems(TextFile.tab);
		}
		tfout.close();
		tfin.close();

		for (String key : snpToCt.keySet()) {
			System.out.println(key + "\t" + snpToCt.get(key));
		}

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
