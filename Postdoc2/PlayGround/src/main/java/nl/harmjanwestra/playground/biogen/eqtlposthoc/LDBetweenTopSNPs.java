package nl.harmjanwestra.playground.biogen.eqtlposthoc;

import nl.harmjanwestra.playground.legacy.vcf.DetermineLD;
import nl.harmjanwestra.playground.legacy.vcf.VCFTabix;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.IOException;
import java.util.HashMap;

public class LDBetweenTopSNPs {


	public static void main(String[] args) {
		String reftemplate = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-b38\\ALL.chrCHR.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz";
		String sampleset = null;

		String eqtlfile1 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-24-ReplicationInAFR\\ldcomparison\\EUR-iteration1\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
		String eqtlfile2 = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-24-ReplicationInAFR\\ldcomparison\\AFR-eQTLProbesFDR0.05-ProbeLevel.txt.gz";

		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-24-ReplicationInAFR\\ldcomparison\\EurVsAFR-TopSNPLD.txt";
		LDBetweenTopSNPs ld = new LDBetweenTopSNPs();

		try {
			ld.ldForTopSNPs(eqtlfile1, eqtlfile2, reftemplate, sampleset, output);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private class EQTL {
		String snp;
		String gene;
		String alleles;
		String assessed;
		int snppos;
		int chr;
		double p;
		double z;
		double beta;
		double se;

		public SNPFeature asFeature() {
			return new SNPFeature(Chromosome.parseChr("" + chr), snppos, snppos + 1);
		}
	}

	public void ldForTopSNPs(String eqtlfile1, String eqtlfile2, String tabixtemplate, String tabixsampleselect, String output) throws IOException {

		HashMap<String, EQTL> eqtlset1 = loadEQTLs(eqtlfile1);
		HashMap<String, EQTL> eqtlset2 = loadEQTLs(eqtlfile2);

		TextFile outf = new TextFile(output, TextFile.W);
		String header = "gene" +
				"\tchr" +
				"\tsnp1" +
				"\tsnp1pos" +
				"\tsnp1alleles" +
				"\tsnp1assessed" +
				"\tsnp1Z" +
				"\tsnp1P" +
				"\tsnp2" +
				"\tsnp2pos" +
				"\tsnp2alleles" +
				"\tsnp2assessed" +
				"\tsnp2Z" +
				"\tsnp2P" +
				"\tdprime" +
				"\trsquared" +
				"\tdistance" +
				"\tconcordant" +

				"\trefsnp1" +
				"\trefsnp1pos" +
				"\trefsnp1alleles" +
				"\trefsnp1assessed\tflip1" +
				"\trefsnp2" +
				"\trefsnp2pos" +
				"\trefsnp2alleles" +
				"\trefsnp2assessed\tflip2";
		outf.writeln(header);
		int overlap = 0;
		int checkable = 0;
		int failedcheck = 0;
		int ctr = 0;
		for (String gene : eqtlset1.keySet()) {
			EQTL e1 = eqtlset1.get(gene);
			EQTL e2 = eqtlset2.get(gene);

			if (e2 != null) {

				overlap++;
				System.out.println(ctr + "/" + overlap + " -- Comparing gene " + gene);
				if (e1.snppos == e2.snppos) {

					// perfect LD.. check for flippage of fx
					Boolean flip = BaseAnnot.flipalleles(e1.alleles, e1.assessed, e2.alleles, e2.assessed);
					if (flip == null) {
						failedcheck++;
					} else {
						checkable++;
						double z2 = e2.z;

						if (flip) {
							z2 *= -1;
						}
						boolean concordant = false;
						if (z2 * e1.z >= 0) {
							concordant = true;
						}

						double dprime = 1;
						double rsquared = 1;
						int distance = 0;
						String outln = gene
								+ "\t" + e1.chr
								+ "\t" + e1.snp + "\t" + e1.snppos + "\t" + e1.alleles + "\t" + e1.assessed
								+ "\t" + e1.z + "\t" + e1.p
								+ "\t" + e2.snp + "\t" + e2.snppos + "\t" + e2.alleles + "\t" + e2.assessed
								+ "\t" + z2 + "\t" + e2.p
								+ "\t" + dprime + "\t" + rsquared + "\t" + distance
								+ "\t" + concordant
								+ "\t" + e1.snp + "\t" + e1.snppos + "\t" + e1.alleles + "\t" + e1.assessed + "\t" + false
								+ "\t" + e2.snp + "\t" + e2.snppos + "\t" + e2.alleles + "\t" + e2.assessed + "\t" + flip;
						outf.writeln(outln);
					}
				} else {
					// check LD
					String tpl = tabixtemplate.replaceAll("CHR", "" + e1.chr);
					VCFTabix tabix = new VCFTabix(tpl);

					boolean[] includesample = null;
					if (tabixsampleselect != null) {
						includesample = tabix.getSampleFilter(tabixsampleselect);
					}
					VCFVariant var1 = tabix.getVariant(e1.asFeature(), includesample);
					VCFVariant var2 = tabix.getVariant(e2.asFeature(), includesample);
					if (var1 != null && var2 != null) {
						checkable++;
						DetermineLD ld = new DetermineLD();
						Pair<Double, Double> ldvals = ld.getLD(var1, var2);
						double dprime = ldvals.getLeft();
						double rsquared = ldvals.getRight();

						int distance = var1.getPos() - var2.getPos();


						// flip to minor allele
						String[] refalleles1 = getAlleleDesc(var1);
						Boolean flip1 = BaseAnnot.flipalleles(refalleles1[0], refalleles1[1], e1.alleles, e1.assessed);
						String[] refalleles2 = getAlleleDesc(var2);
						Boolean flip2 = BaseAnnot.flipalleles(refalleles2[0], refalleles2[1], e2.alleles, e2.assessed);

						if (flip1 == null || flip2 == null) {
							failedcheck++;
						} else {
							double z1 = e1.z;
							if (flip1) {
								z1 *= -1;
							}
							double z2 = e2.z;
							if (flip2) {
								z2 *= -1;
							}

							boolean concordant = false;
							if (z1 * z2 >= 0) {
								concordant = true;
							}

							String outln = gene
									+ "\t" + e1.chr
									+ "\t" + e1.snp + "\t" + e1.snppos + "\t" + e1.alleles + "\t" + e1.assessed
									+ "\t" + z1 + "\t" + e1.p
									+ "\t" + e2.snp + "\t" + e2.snppos + "\t" + e2.alleles + "\t" + e2.assessed
									+ "\t" + z2 + "\t" + e2.p
									+ "\t" + dprime + "\t" + rsquared + "\t" + distance
									+ "\t" + concordant
									+ "\t" + var1.toString() + "\t" + var1.getPos() + "\t" + refalleles1[0] + "\t" + refalleles1[1]
									+ "\t" + flip1
									+ "\t" + var2.toString() + "\t" + var2.getPos() + "\t" + refalleles2[0] + "\t" + refalleles2[1]
									+ "\t" + flip2;
							outf.writeln(outln);
						}


					}


				}
			}
			ctr++;
		}

		System.out.println(overlap + " overlap\t " + checkable + " checkable \t" + failedcheck + " failed check");

	}

	private String[] getAlleleDesc(VCFVariant var1) {
		String[] desc = new String[]{
				var1.getAlleles()[0] + "/" + var1.getAlleles()[1],
				var1.getMinorAllele()
		};
		return desc;

	}

	private HashMap<String, EQTL> loadEQTLs(String eqtlfile1) throws IOException {
		HashMap<String, EQTL> output = new HashMap<>();
		TextFile tf = new TextFile(eqtlfile1, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[1];
			int chr = Integer.parseInt(elems[2]);
			int pos = Integer.parseInt(elems[3]);
			String gene = elems[4];
			String alleles = elems[8];
			String assessed = elems[9];
			double z = Double.parseDouble(elems[10]);
			double p = Double.parseDouble(elems[0]);
			String betastr = elems[18];
			String[] betastrelems = betastr.split(" ");
			double beta = Double.parseDouble(betastrelems[0]);
			double se = Double.parseDouble(betastrelems[1].replaceAll("\\)", "").replaceAll("\\(", ""));

			EQTL e = new EQTL();
			e.gene = gene;
			e.alleles = alleles;
			e.assessed = assessed;
			e.snp = snp;
			e.snppos = pos;
			e.chr = chr;
			e.p = p;
			e.z = z;
			e.beta = beta;
			e.se = se;

			output.put(gene, e);


			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		System.out.println(output.size() + " eQTLs from " + eqtlfile1);
		return output;
	}
}
