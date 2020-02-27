package nl.harmjanwestra.playground.cis.GWAS;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class MergeWithGWAS {

	public static void main(String[] args) {
		String efile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\eQTLDump-sort-FDR-stripped.txt.gz";
		String efile2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\eQTLDump-sort-FDR-stripped-significantgenes.txt.gz";

		String snpstats = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\2019-07-08-Cortex-Cis-SNPQCLog-MAF-sort.txt.gz";

		String outfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\eqtlmerge\\";
		String[] gwas = new String[]{
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\ms2_se-hg38.fullgwas.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\als-hg38.fullgwas.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\al2-hg38.fullgwas.gz",
				"D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\pkd-hg38.fullgwas.gz"
		};

		String[] gwasnames = new String[]{
				"MS",
				"ALS",
				"ALZ",
				"PKD"
		};

		MergeWithGWAS g = new MergeWithGWAS();
		try {
			g.filterGenes(efile, efile2);
//			g.merge(efile2, gwas, gwasnames, snpstats, outfile);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void filterGenes(String efile, String efile2) throws IOException {
		TextFile tf = new TextFile(efile, TextFile.R);
		HashSet<String> significantGenes = new HashSet<String>();
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {

			double fdr = Double.parseDouble(elems[elems.length - 1]);
			if (fdr < 0.05) {
				significantGenes.add(elems[5]);
			}
			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(ctr + " read, " + significantGenes.size() + " significant");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();


		tf.open();
		TextFile tf2 = new TextFile(efile2, TextFile.W);
		tf2.writeln(tf.readLine());
		elems = tf.readLineElems(TextFile.tab);
		ctr = 0;
		int written = 0;
		while (elems != null) {

			if (significantGenes.contains(elems[5])) {
				tf2.writeln(Strings.concat(elems, Strings.tab));
				written++;
			}

			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(ctr + " read, " + written + " written");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf2.close();
		tf.close();
	}

	private class BetaPair {
		public String snp;
		public double fdr;
		double beta1;
		double beta1se;
		double p1;
		double beta2;
		double beta2se;
		double p2;
		boolean beta1sig;
		boolean beta2sig;
	}

	public void merge(String eqtl, String gwasfiles, String gwasnames, String snpstats, String output) throws IOException {
		merge(eqtl, new String[]{gwasfiles}, new String[]{gwasnames}, snpstats, output);
	}

	public void merge(String eqtl, String[] gwasfiles, String[] gwasnames, String snpstats, String output) throws IOException {

		HashMap<String, Double> mafdata = readSNPStats(snpstats);
		for (int i = 0; i < gwasfiles.length; i++) {
			String gwas = gwasfiles[i];
			String gwasname = gwasnames[i];
			HashMap<String, GWASHit> gwashits = loadGWAS(gwas);

			TextFile efile = new TextFile(eqtl, TextFile.R);

			efile.readLine(); // skip header.

			String[] elems = efile.readLineElems(TextFile.tab);

			int nwritten = 0;

			int read = 0;


			HashMap<String, ArrayList<BetaPair>> data = new HashMap<String, ArrayList<BetaPair>>();

			int nrskipped = 0;
			int nrSNPsNotFound = 0;
			while (elems != null) {
				String snp = elems[2];
				Double maf = mafdata.get(snp);
				if (maf == null || Double.isNaN(maf)) {
					// System.out.println("Could not find SNP in MAF data: " + snp);
					nrSNPsNotFound++;
				} else {

					String[] snpelems = snp.split(":");
					String rs = snpelems[2]; // match on rsid
					rs = snpelems[0] + ":" + snpelems[1]; // match on position
					GWASHit hit = gwashits.get(rs);
					String eqtlalleles = elems[3];
					String eqtlalleleassessed = elems[4];

					if (hit != null) {
						Boolean b = BaseAnnot.flipalleles(eqtlalleles, eqtlalleleassessed, hit.alleles, hit.assessed);
						if (b != null) {
							double beta2 = hit.beta;
							if (b) {
								beta2 *= -1;
							}

							double p1 = Double.parseDouble(elems[0]);

							double beta1 = Double.parseDouble(elems[6]);
							double beta1se = Double.parseDouble(elems[7]);
							int n1 = Integer.parseInt(elems[8]);
							double fdr = Double.parseDouble(elems[9]);


							BetaPair betaPair = new BetaPair();
							betaPair.snp = snp;

//						betaPair.n1 = n1;
							betaPair.beta1 = beta1;
							betaPair.beta1se = beta1se;
							betaPair.fdr = fdr;
							if (fdr < 0.05) {
								betaPair.beta1sig = true;
							}
							betaPair.p1 = p1;
							betaPair.p2 = hit.p;

							betaPair.beta2 = beta2;
							betaPair.beta2se = hit.se;
							if (hit.p < 5e-8) {
								betaPair.beta2sig = true;
							}


							// check for NA
							if (Double.isNaN(beta1) || Double.isNaN(beta2) || Double.isNaN(hit.se) || Double.isNaN(beta1se) ||
									beta1 == 0 || beta2 == 0 || hit.se == 0 || beta1se == 0) {
								nrskipped++;
							} else {
								String gene = elems[5];
								ArrayList<BetaPair> sets = data.get(gene);
								if (sets == null) {
									sets = new ArrayList<>();
								}

								sets.add(betaPair);
								data.put(gene, sets);
							}
							nwritten++;

						}
					}
				}
				read++;
				if (read % 1000000 == 0) {
					System.out.print(read + " read, " + nwritten + " added. " + nrSNPsNotFound + " have no MAF? " + nrskipped + " out of range.\r");
				}
				elems = efile.readLineElems(TextFile.tab);
			}
			efile.close();

			String dirout = output + gwasname + "/";
			if (gwasname != null) {
				dirout = output + gwasname + "/";
			}

			Gpio.createDir(dirout);
			for (String gene : data.keySet()) {
				ArrayList<BetaPair> sets = data.get(gene);

				// optional: determine if there are any significant hits in the locus

				boolean locussig = false;
				boolean eqtlsig = false;
				for (BetaPair p : sets) {
					if (p.beta1sig) {
						eqtlsig = true;
					}
					if (p.beta2sig) {
						locussig = true;
					}
				}

				if (sets.size() > 30 && locussig && eqtlsig) {
					gene = gene.split("\\.")[0];
					String fileout = dirout + gene + ".txt.gz";
					TextFile outf = new TextFile(fileout, TextFile.W);
					outf.writeln("SNP\tb1\tse1\tp1\tfdr\tb2\tse2\tp2\tmaf");
					int geneswritten = 0;
					for (BetaPair p : sets) {
						Double maf = mafdata.get(p.snp);
						if (maf != null) {
							outf.writeln(p.snp + "\t" + p.beta1 + "\t" + p.beta1se + "\t" + p.p1 + "\t" + p.fdr + "\t" + p.beta2 + "\t" + p.beta2se + "\t" + p.p2 + "\t" + maf);

						}

					}
					geneswritten++;
					if (geneswritten % 250 == 0) {
						System.out.print(geneswritten + "genes written sofar.\r");
					}
					outf.close();
				}
			}
		}

	}

	private HashMap<String, Double> readSNPStats(String snpstats) throws IOException {

		TextFile tf = new TextFile(snpstats, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, Double> snpToMaf = new HashMap<>();
		int ctr = 0;
		while (elems != null) {
			String snp = elems[0];
			double maf = Double.parseDouble(elems[elems.length - 1]);
			snpToMaf.put(snp, maf);
			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.print(ctr + " maf annotations read.\r");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(snpToMaf.size() + " MAF annotations loaded.");
		return snpToMaf;

	}

	private class GWASHit {
		String alleles;
		String assessed;
		double beta;
		double se;
		double p;
	}

	public HashMap<String, GWASHit> loadGWAS(String gwasFile) throws IOException {

		HashMap<String, GWASHit> map = new HashMap<>();
		TextFile tf = new TextFile(gwasFile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(Strings.whitespace);
		int nrNA = 0;
		int nrNA2 = 0;
		while (elems != null) {
			GWASHit hit = new GWASHit();

			hit.alleles = elems[1] + "/" + elems[2];
			hit.assessed = elems[2];
			try {
				hit.beta = Double.parseDouble(elems[3]);
				hit.se = Double.parseDouble(elems[4]);
				hit.p = Double.parseDouble(elems[5]);

				if (!Double.isNaN(hit.se) && !Double.isNaN(hit.beta)) {
					map.put(elems[0], hit);
				} else {
					nrNA2++;
				}


			} catch (NumberFormatException e) {
				nrNA++;
//				System.out.println("Error parsing " + elems[0] + ":\t" + elems[3] + "\t" + elems[4] + "\t" + elems[5]);
			}
			elems = tf.readLineElems(Strings.whitespace);
		}
		tf.close();

		System.out.println(map.size() + " gwas hits in " + gwasFile + ", " + nrNA + " with NA, " + nrNA2 + " with NA SE?");
		return map;
	}
}
