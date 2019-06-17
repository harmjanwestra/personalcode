package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class eQTLGeneOverlap {

	public static void main(String[] args) {
		String eqtlgenprobes = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLProbesFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz";
		String eqtlgeneqtlsfdr = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz";
		String metabrainprobes = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-06-04-results\\cis-proteincoding\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
		String metabraineqtlsfdr = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-06-04-results\\cis-proteincoding\\eQTLsFDR0.05-ProbeLevel.txt.gz";

//		metabraineqtlsfdr = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-06-04-results\\cis\\2019-06-01-EUR-Cis-1mb-cortex-eQTLsFDR0.05-ProbeLevel.txt.gz";
//		metabrainprobes = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-06-04-results\\cis\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";

		String metabrainToEQTLGenout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-06-04-results\\cis-proteincoding\\eQTLsFDR0.05-ProbeLevel-metabraintoeqtlgen.txt";
		String eqtlgentometabrainout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-06-04-results\\cis-proteincoding\\eQTLsFDR0.05-ProbeLevel-eqtlgentometabrain.txt";

		eQTLGeneOverlap e = new eQTLGeneOverlap();
		e.eqtlgentometabrain(eqtlgenprobes, metabraineqtlsfdr, eqtlgentometabrainout);
		e.metabrainToEQTLGen(eqtlgeneqtlsfdr, metabrainprobes, metabrainToEQTLGenout);
	}


	public void eqtlgentometabrain(String eqtlgen, String metabrain, String f3) {
		try {
			HashSet<String> l1 = new HashSet<String>();

			QTLTextFile t = new QTLTextFile(eqtlgen, QTLTextFile.R);
			ArrayList<EQTL> eqtllist = t.readList();
			t.close();

			HashMap<String, EQTL> eqtlhash = new HashMap<>();
			for (EQTL e : eqtllist) {
				l1.add(e.getProbe());
				eqtlhash.put(e.getRsName() + "-" + e.getProbe(), e);
			}
			System.out.println(l1.size() + " genes in list 1");


			int opposite = 0;
			int concordant = 0;
			int incompat = 0;
			int shared = 0;

			TextFile tf = new TextFile(metabrain, TextFile.R);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			HashSet<String> ovlerap = new HashSet<>();
			TextFile tfo = new TextFile(f3, TextFile.W);
			tfo.writeln("eQTL\teqtlgenz\tmetabrainz\tflip\tmetabraintopsnp");
			HashSet<String> visitedGenes = new HashSet<String>();


			while (elems != null) {
				String gene = Strings.dot.split(elems[4])[0];
				String rs = elems[1].split(":")[2];
				boolean topsnp = false;
				if (!visitedGenes.contains(gene)) {
					topsnp = true;
					visitedGenes.add(gene);
				}

				EQTL e = eqtlhash.get(rs + "-" + gene);
				if (e != null) {
					shared++;
					String alleles = elems[8];
					String assessed = elems[9];
					Double z = Double.parseDouble(elems[10]);
					Boolean flip = BaseAnnot.flipalleles(e.getAlleles(), e.getAlleleAssessed(), alleles, assessed);

					if (flip != null) {


						if (flip) {
							z *= -1;
						}

//						System.out.println(rs + "-" + gene + "\t" + e.getZscore() + "\t" + z);
						tfo.writeln(rs + "-" + gene + "\t" + e.getZscore() + "\t" + z + "\t" + flip + "\t" + topsnp);
						if (z * e.getZscore() < 0) {
							opposite++;

						} else {
							concordant++;
						}

					} else {
						incompat++;
					}
				}
				if (l1.contains(gene)) {
					ovlerap.add(gene);
				}


				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			tfo.close();

			System.out.println();
			System.out.println("eqtlgen to metabrain");
			System.out.println(ovlerap.size() + " genes overlap");
			System.out.println(ovlerap.size() + " eqtls overlap");
			System.out.println(shared + " eqtls overlap");
			System.out.println(opposite + " eqtls overlap discordant, " + ((double) opposite / shared));
			System.out.println(concordant + " eqtls overlap concordant");
			System.out.println(incompat + " eqtls overlap incompatible");

		} catch (IOException e) {
			e.printStackTrace();
		}


	}

	public void metabrainToEQTLGen(String eqtlgen, String metabrain, String output) {
		try {

			QTLTextFile t = new QTLTextFile(metabrain, QTLTextFile.R);
			ArrayList<EQTL> eqtllist = t.readList();
			t.close();

			for (EQTL e : eqtllist) {
				e.setRsName(Strings.colon.split(e.getRsName())[2]);
				e.setProbe(Strings.dot.split(e.getProbe())[0]);
			}

			HashMap<String, EQTL> eqtlhash = new HashMap<>();
			for (EQTL e : eqtllist) {
				eqtlhash.put(e.getRsName() + "-" + e.getProbe(), e);
			}


			int opposite = 0;
			int concordant = 0;
			int incompat = 0;
			int shared = 0;

			TextFile tf = new TextFile(eqtlgen, TextFile.R);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			HashSet<String> ovlerap = new HashSet<>();
			TextFile tfo = new TextFile(output, TextFile.W);

			while (elems != null) {
				String gene = Strings.dot.split(elems[4])[0];
				String rs = elems[1]; // .split(":")[2];

				EQTL e = eqtlhash.get(rs + "-" + gene);
				if (e != null) {
					shared++;
					String alleles = elems[8];
					String assessed = elems[9];
					Double z = Double.parseDouble(elems[10]);
					Boolean flip = BaseAnnot.flipalleles(e.getAlleles(), e.getAlleleAssessed(), alleles, assessed);


					if (flip != null) {


						if (flip) {
							z *= -1;
						}

//						System.out.println(rs + "-" + gene + "\t" + e.getZscore() + "\t" + z);
						tfo.writeln(rs + "-" + gene + "\t" + e.getZscore() + "\t" + z + "\t" + flip);
						if (z * e.getZscore() < 0) {
							opposite++;
						} else {
							concordant++;
						}

					} else {
						incompat++;
					}
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			tfo.close();

			System.out.println();
			System.out.println("metabrain to eqtlgen");
			System.out.println(ovlerap.size() + " genes overlap");
			System.out.println(shared + " eqtls overlap");
			System.out.println(opposite + " eqtls overlap discordant, " + ((double) opposite / shared));
			System.out.println(concordant + " eqtls overlap concordant");
			System.out.println(incompat + " eqtls overlap incompatible");


		} catch (IOException e) {
			e.printStackTrace();
		}


	}

}
