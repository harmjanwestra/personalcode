package nl.harmjanwestra.playground.biogen.freeze2dot1;

import eqtlmappingpipeline.binarymeta.meta.ZScoreComparison;
import eqtlmappingpipeline.util.QTLFileCompare;
import eqtlmappingpipeline.util.QTLFileSorter;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class ReplicateMetabrainInEQTLgen {


	public static void main(String[] args) {

		// cis
		String ciseqtlfilemb = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\cis-cortex-eurandafr-iterative\\Iteration1-eQTLProbesFDR0.05-ProbeLevel.txt.gz";
		String eqtlgenloc = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-MetaBrain2IDs.txt.gz";
		String eqtlgenout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\cis-cortex-eurandafr-replicationineqtlgen\\2019-12-11-cis-eQTLsFDR-ProbeLevel-eqtlgen.txt.gz";

		BHFDR fdr = new BHFDR();
		QTLFileCompare c = new QTLFileCompare();
		ReplicateMetabrainInEQTLgen r = new ReplicateMetabrainInEQTLgen();
//		try {
////			r.runCis(ciseqtlfilemb, eqtlgenloc, eqtlgenout);
//
//			fdr.run(eqtlgenout, 0.05, eqtlgenout + "fdr.txt.gz");
//
//			String compout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\cis-cortex-eurandafr-replicationineqtlgen\\comp";
//			c.compareOverlapAndZScoreDirectionTwoEQTLFiles(ciseqtlfilemb, "metabrain",
//					eqtlgenout + "fdr.txt.gz", "eqtlgen",
//					compout, 1d, false, false, false);
//			compout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\cis-cortex-eurandafr-replicationineqtlgen\\compFDR";
//			c.compareOverlapAndZScoreDirectionTwoEQTLFiles(ciseqtlfilemb, "metabrain",
//					eqtlgenout + "fdr.txt.gz", "eqtlgen",
//					compout, 0.05, false, false, false);
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}

		String metabraintrans = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR\\eQTLsFDR0.05.txt.gz";
		String eqtlgenzscores = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\binary\\ZScoreMatrix-MetaBrain2dot1IDs.txt.gz";
		eqtlgenout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR-replicationineqtlgen\\eqtlgen-trans.txt.gz";
		QTLFileSorter s = new QTLFileSorter();
		try {
			r.runTrans(metabraintrans, eqtlgenzscores, eqtlgenout);
			s.run(eqtlgenout, eqtlgenout + "-sort.txt.gz", QTLFileSorter.SORTBY.P);
			fdr.run(eqtlgenout + "-sort.txt.gz", 0.05, eqtlgenout + "-sort-FDR.txt.gz");
			String compout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR-replicationineqtlgen\\comp";

			c.compareOverlapAndZScoreDirectionTwoEQTLFiles(metabraintrans, "metabrain",
					eqtlgenout + "-sort-FDR.txt.gz", "eqtlgen",
					compout, -1, true, true, false);
			compout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\trans-cortex-EURandAFR-replicationineqtlgen\\compFDR";
			c.compareOverlapAndZScoreDirectionTwoEQTLFiles(metabraintrans, "metabrain",
					eqtlgenout + "-sort-FDR.txt.gz", "eqtlgen",
					compout, 0.05, true, true, false);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}


	}

	private void runTrans(String metabraintrans, String eqtlgenzscores, String eqtlgenout) throws IOException {
		HashMap<String, minimaleqtl> set1 = new HashMap<>();
		TextFile tf = new TextFile(metabraintrans, TextFile.R);
		String[] mheader = tf.readLineElems(TextFile.tab);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> genes = new HashSet<String>();
		HashMap<String, String> genepos = new HashMap<>();
		HashMap<String, String> snppos = new HashMap<>();
		HashMap<String, String> hugo = new HashMap<>();
		HashMap<String, String> genemap = new HashMap<>();
		while (elems != null) {
			String snp = elems[1];
			snppos.put(snp, elems[2] + ";" + elems[3]);
			String gene = elems[4];
			String geneid = gene.split("\\.")[0];
			genepos.put(gene, elems[5] + ";" + elems[6]);
			genes.add(geneid);
			genemap.put(geneid, gene);
			hugo.put(gene, elems[QTLTextFile.HUGO]);
			Float z = Float.parseFloat(elems[10]);
			minimaleqtl m = new minimaleqtl();

			m.alleles = elems[8];
			m.assessed = elems[9];
			m.fdr = Float.parseFloat(elems[elems.length - 1]);

			m.snp = snp;
			m.z = z;
			set1.put(snp + "-" + gene, m);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(set1.size() + " eqtls loaded");
		System.out.println(genes.size() + " genes");

		tf = new TextFile(eqtlgenzscores, TextFile.R);
		String[] header = tf.readLineElems(TextFile.tab);
		int genesfound = 0;
		System.out.println((header.length - 3) + " genes present. ");

		HashSet<String> foundgenes = new HashSet<>();
		for (int i = 3; i < header.length; i++) {
			String gene = header[i].split("_")[0];

			if (genes.contains(gene)) {
				header[i] = genemap.get(gene);
				foundgenes.add(gene);
				genesfound++;
			}
		}
		System.out.println(genesfound + " genes found.");

		TextFile outfg = new TextFile(eqtlgenout + "genesnotfound.txt", TextFile.W);
		for (String gene : genes) {
			if (!foundgenes.contains(gene)) {
				outfg.writeln(gene);
			}
		}
		outfg.close();

		TextFile outf = new TextFile(eqtlgenout, TextFile.W);
		outf.writeln(Strings.concat(mheader, Strings.tab));
		elems = tf.readLineElems(TextFile.tab);
		int snpsfound = 0;
		while (elems != null) {
			String snp = elems[0];
			if (snppos.containsKey(snp)) {
				snpsfound++;
			}
			for (int i = 3; i < elems.length; i++) {

				if (set1.containsKey(snp + "-" + header[i])) {
					String[] outputln = new String[mheader.length];
					double z = Double.parseDouble(elems[i]);
					if (!Double.isNaN(z)) {
						outputln[0] = "" + ZScores.zToP(z);

						String[] genep = genepos.get(header[i]).split(";");
						String[] snpp = snppos.get(snp).split(";");


						outputln[1] = snp;
						outputln[2] = snpp[0];
						outputln[3] = snpp[1];
						outputln[4] = header[i];
						outputln[5] = genep[0];
						outputln[6] = genep[1];
						outputln[10] = elems[i];
						outputln[9] = elems[2];
						outputln[8] = elems[1];
						outputln[7] = "trans";
						outputln[QTLTextFile.HUGO] = hugo.get(header[i]);
						outf.writeln(Strings.concat(outputln, Strings.tab));
					}
				}

			}
			elems = tf.readLineElems(TextFile.tab);
		}
		outf.close();
		System.out.println(snpsfound + " snps found out of " + snppos.size());
	}


	public void runCis(String metabrain, String eqtlgen, String eqtlgenout) throws IOException {
		HashMap<String, minimaleqtl> set1 = new HashMap<>();
		TextFile tf = new TextFile(metabrain, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[1];
			String gene = elems[4];
			Float z = Float.parseFloat(elems[10]);
			minimaleqtl m = new minimaleqtl();

			m.alleles = elems[8];
			m.assessed = elems[9];
			m.fdr = Float.parseFloat(elems[elems.length - 1]);

			m.snp = snp;
			m.z = z;
			set1.put(snp + "-" + gene, m);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(set1.size() + " eqtls loaded");

		int ctr = 0;
		int written = 0;
		TextFile tf2 = new TextFile(eqtlgen, TextFile.R);
		TextFile outf = new TextFile(eqtlgenout, TextFile.W);
		outf.writeln(tf2.readLine());
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[1];
			String gene = elems[4];
			minimaleqtl m = set1.get(snp + "-" + gene);
			if (m != null) {
				outf.writeln(Strings.concat(elems, Strings.tab));
				written++;
			}
			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(ctr + " read, " + written + " written.");
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		outf.close();


	}
}
