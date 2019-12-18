import nl.harmjanwestra.playground.legacy.vcf.VCFTabix;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.IntStream;

public class PowellRewrite {


	public static void main(String[] args) {


		String infile = "D:\\powell\\BMem.txt";
		String indir = "D:\\powell\\";
		String outdir = "D:\\powellemp\\";
		String outfile1 = "D:\\powellemp\\BMem.txt";
		String outfile2 = "D:\\powellemp\\BMem-comp.txt";
		String[] lsof = Gpio.getListOfFiles(indir);
		String efile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\trans\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
		IntStream.range(0, lsof.length).parallel().forEach(v -> {
					String file = lsof[v];
					PowellRewrite p = new PowellRewrite();
					try {
						p.processFile(indir + file, outdir + file);
						p.compare(efile, outdir + file, outdir + file + "-comp.txt");
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
		);


//
//

//


	}

	private class EQTL {
		String p;
		String z;
		String alleles;
		String assessed;
	}

	public void compare(String infile1, String infile2, String outfile) throws IOException {

		HashMap<String, EQTL> eqtls = new HashMap<String, EQTL>();

		TextFile tf1 = new TextFile(infile1, TextFile.R);
		tf1.readLine();
		String[] elems = tf1.readLineElems(TextFile.tab);
		while (elems != null) {

			String p = elems[0];
			String id = elems[1] + "-" + elems[4];
			String z = elems[10];
			String alleles = elems[8];
			String assessed = elems[9];

			EQTL e = new EQTL();
			e.alleles = alleles;
			e.assessed = assessed;
			e.p = p;
			e.z = z;

			eqtls.put(id, e);


			elems = tf1.readLineElems(TextFile.tab);
		}
		tf1.close();

		TextFile tf2 = new TextFile(infile2, TextFile.R);
		TextFile outf = new TextFile(outfile, TextFile.W);
		outf.writeln("SNP\tGene\tAllelesPowell\tAssessedPowell\tPPowell\tZPowell\tAllelesEQTLGen\tAssessedEQTLGen\tPEQTLGen\tZEQTLGen\tConcordant");
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);

		int zctr1 = 0;
		int zctr1concordant = 0;
		int zctr2concordant = 0;
		int zctr2 = 0;

		HashSet<String> snpcounted = new HashSet<String>();
		int nrPos = 0;
		int nrTotal = 0;
		while (elems != null) {


			String pstr = elems[0];
			String id = elems[1] + "-" + elems[4];
			String zstr = elems[10];
			String alleles = elems[8];
			String assessed = elems[9];

			EQTL e = eqtls.get(id);

			Boolean flipAllele = BaseAnnot.flipalleles(e.alleles, e.assessed, alleles, assessed);
			if (flipAllele != null) {
				Double z = Double.parseDouble(zstr);
				Double z2 = Double.parseDouble(e.z);
				if (flipAllele) {
					z *= -1;
				}
				int concordant = 1;
				if (z * z2 < 0) {
					concordant = 0;
				}

				if (z2 > 0) {
					nrPos++;
				}
				nrTotal++;

				double p = Double.parseDouble(pstr);
				String snp = elems[1];
//				if (!snpcounted.contains(snp)) {

					if (p < 0.001) {

						if (z2 >= 0) {
							zctr1++;
							if (concordant == 1) {
								zctr1concordant++;
							}
						} else {
							zctr2++;
							if (concordant == 1) {
								zctr2concordant++;
							}
						}
					}
					snpcounted.add(snp);
//				}


				String outln = elems[1] + "\t" + elems[4]
						+ "\t" + alleles
						+ "\t" + assessed
						+ "\t" + p
						+ "\t" + z
						+ "\t" + e.alleles
						+ "\t" + e.assessed
						+ "\t" + e.p
						+ "\t" + e.z
						+ "\t" + concordant;

				outf.writeln(outln);
			}

			eqtls.put(id, e);


			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		outf.close();
		System.out.println(infile2 + "\t" + zctr1 + "\t" + zctr1concordant + "\t" + ((double) zctr1concordant / zctr1)
				+ "\t" + zctr2 + "\t" + zctr2concordant + "\t" + ((double) zctr2concordant / zctr2) + "\t" + nrPos + "\t" + nrTotal + "\t" + ((double) nrPos / nrTotal));

	}

	public void processFile(String infile, String outfile) throws IOException {
		HashMap<String, String> positions = new HashMap<String, String>();
		TextFile tf = new TextFile(infile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snpid = elems[0].split("_")[0];
			String[] snpelems = snpid.split(":");

			if (!elems[2].equals("NA")) {
				positions.put(snpelems[0] + ":" + snpelems[1], null);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(positions.size() + " positions in file: " + infile);


		VCFTabix tabix = new VCFTabix("D:\\Work\\Data\\Ref\\1kg\\ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz");
		int nrread = 0;
		int updated = 0;
		for (String key : positions.keySet()) {
			String[] snpelems = key.split(":");
			SNPFeature feature = new SNPFeature();
			feature.setChromosome(Chromosome.parseChr(snpelems[0]));
			feature.setStart(Integer.parseInt(snpelems[1]));
			feature.setStop(Integer.parseInt(snpelems[1]) + 1);
			VCFVariant v = tabix.getVariant(feature, null);
			if (v != null && v.getAlleles().length > 1 && v.getAlleles().length < 3) {
				String pos = v.getChr() + ":" + v.getPos();
				if (positions.containsKey(pos)) {
					String[] alleles = v.getAlleles();
					positions.put(pos, alleles[0] + "/" + alleles[1] + "\t" + alleles[1]);
					updated++;
				}
			}
			nrread++;
//			if (nrread % 1000000 == 0) {
//			System.out.println(infile + "\t" + nrread + " read " + updated + " updated");
//			}
		}

//		TextFile tf2 = new TextFile("D:\\Work\\Data\\Ref\\1kg\\ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz", TextFile.R);
//
//		String ln = tf2.readLine();
//		int nrread = 0;
//		int updated = 0;
//		while (ln != null) {
//			if (!ln.startsWith("#")) {
//				VCFVariant v = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
//
//				if (v.isBiallelic()) {
//					String pos = v.getChr() + ":" + v.getPos();
//					if (positions.containsKey(pos)) {
//						String[] alleles = v.getAlleles();
//						positions.put(pos, alleles[0] + "/" + alleles[1] + "\t" + alleles[1]);
//						updated++;
//					}
//				}
//				nrread++;
//				if (nrread % 1000000 == 0) {
//					System.out.println(infile + "\t" + nrread + " read " + updated + " updated");
//				}
//			}
//			ln = tf2.readLine();
//		}
//		tf2.close();

		TextFile out = new TextFile(outfile, TextFile.W);
		out.writeln("PValue"
				+ "\tSNPName"
				+ "\tSNPChr"
				+ "\tSNPChrPos"
				+ "\tProbeName"
				+ "\tProbeChr"
				+ "\tProbeCenterChrPos"
				+ "\tCisTrans"
				+ "\tSNPType"
				+ "\tAlleleAssessed"
				+ "\tOverallZScore"
				+ "\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC"
				+ "\tDatasetsZScores"
				+ "\tDatasetsNrSamples"
				+ "\tIncludedDatasetsMeanProbeExpression"
				+ "\tIncludedDatasetsProbeExpressionVariance"
				+ "\tHGNCName"
				+ "\tIncludedDatasetsCorrelationCoefficient"
				+ "\tMeta-Beta (SE)"
				+ "\tBeta (SE)"
				+ "\tFoldChange"
				+ "\tFDR");

		tf = new TextFile(infile, TextFile.R);
		tf.readLine();
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			if (!elems[2].equals("NA")) {
				String snpid = elems[0].split("_")[0];
				String[] snpelems = snpid.split(":");

				String alleles = positions.get(snpelems[0] + ":" + snpelems[1]);
				if (alleles != null) {
					String outln = elems[4]
							+ "\t" + elems[6]
							+ "\t" + elems[7]
							+ "\t" + elems[8]
							+ "\t" + elems[12]
							+ "\t" + elems[13]
							+ "\t" + elems[14]
							+ "\ttrans"
							+ "\t" + alleles
//						+ "\t" + assessed
							+ "\t" + elems[3]
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-"
							+ "\t-";
					out.writeln(outln);
				}

			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		out.close();
	}
}
