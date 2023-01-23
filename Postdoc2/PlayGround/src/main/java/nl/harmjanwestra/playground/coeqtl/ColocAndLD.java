package nl.harmjanwestra.playground.coeqtl;

import nl.harmjanwestra.playground.legacy.vcf.DetermineLD;
import nl.harmjanwestra.playground.legacy.vcf.VCFTabix;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class ColocAndLD {


	public static void main(String[] args) {
		ColocAndLD r = new ColocAndLD();
		try {
			r.run();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public class ColocDataObj {
		String snp;
		int pos;
		Chromosome chr;
		double eqtlp;
		double gwasp;
	}

	public void run() throws IOException {
		String sumstats = "D:\\Sync\\TMP\\coeqtlcoloc\\COEQTL_example_input_Type_1_Diabetes.csv";
		String vcf = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-b38\\ALL.chr12.shapeit2_integrated_v1a.GRCh38.20181129.phased.nochr.rsid.vcf.gz";
		String vcfsamples = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p35va-europeans.txt.gz";
		String out = "D:\\Sync\\TMP\\coeqtlcoloc\\COEQTL_example_input_Type_1_Diabetes-wABP-wLD.txt";


		TextFile tf = new TextFile(sumstats, TextFile.R);
		int snpcol = -1;
		int snpposcol = -1;
		int snpchrcol = -1;
		int eqtlpcol = -1;
		int eqtlbetacol = -1;
		int eqtlsecol = -1;
		int gwaspcol = -1;
		int gwasbetacol = -1;
		int gwassecol = -1;

		String[] header = tf.readLineElems(TextFile.comma);
		for (int i = 0; i < header.length; i++) {
			header[i] = header[i].replaceAll("\"", "");
			if (header[i].equals("SNP")) {
				snpcol = i;
			}
			if (header[i].equals("SNPChr")) {
				snpchrcol = i;
			}
			if (header[i].equals("SNPPos")) {
				snpposcol = i;
			}
			if (header[i].equals("MetaP")) {
				eqtlpcol = i;
			}
			if (header[i].equals("MetaBeta")) {
				eqtlbetacol = i;
			}
			if (header[i].equals("MetaSE")) {
				eqtlsecol = i;
			}
			if (header[i].equals("gwas_pvalue")) {
				gwaspcol = i;
			}
			if (header[i].equals("gwas_effect_size")) {
				gwasbetacol = i;
			}
			if (header[i].equals("gwas_standard_error")) {
				gwassecol = i;
			}
		}

		ArrayList<AssociationResult> eQTLassociationResults = new ArrayList<>();
		ArrayList<AssociationResult> GWASassociationResults = new ArrayList<>();
		String[] elems = tf.readLineElems(TextFile.comma);

		AssociationResult indexAssocGWAS = null;
		AssociationResult indexAssocEQTL = null;

		HashSet<String> snps = new HashSet<>();
		int minpos = Integer.MAX_VALUE;
		int maxpos = 0;
		Chromosome chromosome = null;
		while (elems != null) {
			String snp = elems[snpcol].replaceAll("\"", "");
			int pos = Integer.parseInt(elems[snpposcol].replaceAll("\"", ""));

			if (pos > maxpos) {
				maxpos = pos;
			}
			if (pos < minpos) {
				minpos = pos;
			}

			Chromosome chr = Chromosome.parseChr(elems[snpchrcol].replaceAll("\"", ""));
			chromosome = chr;

			double eqtlp = Double.parseDouble(elems[eqtlpcol].replaceAll("\"", ""));
			double eqtlbeta = Double.parseDouble(elems[eqtlbetacol].replaceAll("\"", ""));
			double eqtlse = Double.parseDouble(elems[eqtlsecol].replaceAll("\"", ""));

			double gwasp = Double.parseDouble(elems[gwaspcol].replaceAll("\"", ""));
			double gwasbeta = Double.parseDouble(elems[gwasbetacol].replaceAll("\"", ""));
			double gwasse = Double.parseDouble(elems[gwassecol].replaceAll("\"", ""));

			snps.add(snp);

			AssociationResult eqtlobj = new AssociationResult();
			eqtlobj.snp = new SNPFeature(chr, pos, pos);
			eqtlobj.snp.setName(snp);
			eqtlobj.setBeta(eqtlbeta);
			eqtlobj.setSe(eqtlse);
			eqtlobj.setPval(eqtlp);
			eQTLassociationResults.add(eqtlobj);

			AssociationResult gwasobj = new AssociationResult();
			gwasobj.snp = new SNPFeature(chr, pos, pos);
			gwasobj.snp.setName(snp);
			gwasobj.setBeta(gwasbeta);
			gwasobj.setSe(gwasse);
			gwasobj.setPval(gwasp);
			GWASassociationResults.add(gwasobj);

			if (indexAssocGWAS == null || gwasp < indexAssocGWAS.getPval()) {
				indexAssocGWAS = gwasobj;
			}
			if (indexAssocEQTL == null || eqtlp < indexAssocEQTL.getPval()) {
				indexAssocEQTL = eqtlobj;
			}

			elems = tf.readLineElems(TextFile.comma);
		}

		tf.close();

		System.out.println("Top GWAS SNP: " + indexAssocGWAS.getSnp().getName() + "\t" + indexAssocGWAS.getPval());
		System.out.println("Top EQTL SNP: " + indexAssocEQTL.getSnp().getName() + "\t" + indexAssocEQTL.getPval());

		System.out.println(eQTLassociationResults.size() + " associations in " + sumstats);

		ApproximateBayesPosterior abp = new ApproximateBayesPosterior();
		abp.calculatePosterior(eQTLassociationResults);
		abp.calculatePosterior(GWASassociationResults);

		VCFTabix tabix = new VCFTabix(vcf);
		boolean[] filter = tabix.getSampleFilter(vcfsamples);
		maxpos += 1;
		minpos -= 1;
		Feature window = new Feature(chromosome, minpos, maxpos);

		ArrayList<VCFVariant> variants = tabix.getAllVariants(window, filter);
		System.out.println(variants.size() + " variants in region: " + window);
		VCFVariant indexVariantGWAS = null;
		VCFVariant indexVariantEQTL = null;
		HashMap<String, VCFVariant> variantMap = new HashMap<>();
		for (VCFVariant v : variants) {
			variantMap.put(v.getId(), v);
			if (v.getId().equals(indexAssocGWAS.getSnp().getName())) {
				indexVariantGWAS = v;
			}
			if (v.getId().equals(indexAssocEQTL.getSnp().getName())) {
				indexVariantEQTL = v;
			}
		}

		if (indexVariantEQTL != null && indexVariantGWAS != null) {

			DetermineLD ldcalc = new DetermineLD();

			TextFile tfo = new TextFile(out, TextFile.W);

			tfo.writeln("SNP\tChr\tPos\tEQTL-P\tEQTL-ABP\tGWAS-P\tGWAS-ABP\tRsq-TopEQTLSNP(" + indexAssocEQTL.getSnp().getName() + ")\tRsqTopGWASSNP(" + indexAssocGWAS.getSnp().getName() + ")");

			ArrayList<AssociationResult> eqtlspruned = new ArrayList<>();
			ArrayList<AssociationResult> gwasspruned = new ArrayList<>();
			AssociationResult newTopEQTLResult = null;
			AssociationResult newTopGWASResult = null;
			VCFVariant newTopVariantGWAS = null;
			VCFVariant newTopVariantEQTL = null;

			for (int i = 0; i < eQTLassociationResults.size(); i++) {
				AssociationResult eqtlresult = eQTLassociationResults.get(i);
				AssociationResult gwasresult = GWASassociationResults.get(i);
				if (eqtlresult.getSnp().getName().equals("rs1701704")) {
					System.out.println("ping!");
				}
				VCFVariant v = variantMap.get(eqtlresult.getSnp().getName());
				double rve = Double.NaN;
				double rvg = Double.NaN;
				if (v != null) {
					Pair<Double, Double> rg = ldcalc.getLD(v, indexVariantGWAS);
					if (rg != null) {
						rvg = rg.getRight();
					}
					Pair<Double, Double> re = ldcalc.getLD(v, indexVariantEQTL);
					if (re != null) {
						rve = re.getRight();
						if (rve < 0.7) {
							if (newTopEQTLResult == null || eqtlresult.getPval() < newTopEQTLResult.getPval()) {
								newTopEQTLResult = eqtlresult;
								newTopVariantEQTL = v;
							}
							if (newTopGWASResult == null || gwasresult.getPval() < newTopGWASResult.getPval()) {
								newTopGWASResult = gwasresult;
								newTopVariantGWAS = v;
							}
							eqtlspruned.add(eqtlresult);
							gwasspruned.add(gwasresult);
						}
					}
				}
				String outln = eqtlresult.getSnp().getName() + "\t" + eqtlresult.getSnp().getChromosome().getNumber() + "\t" + eqtlresult.getSnp().getStart()
						+ "\t" + eqtlresult.getPval() + "\t" + eqtlresult.getPosterior()
						+ "\t" + gwasresult.getPval() + "\t" + gwasresult.getPosterior()
						+ "\t" + rve + "\t" + rvg;
				tfo.writeln(outln);
			}
			tfo.close();
			System.out.println("Done with all snps");
			// recalculate posteriors after removing top results

			System.out.println(eqtlspruned.size() + " snps remain after pruning");
			System.out.println("Top GWAS SNP: " + newTopGWASResult.getSnp().getName() + "\t" + newTopGWASResult.getPval());
			System.out.println("Top EQTL SNP: " + newTopEQTLResult.getSnp().getName() + "\t" + newTopEQTLResult.getPval());

			abp.calculatePosterior(eqtlspruned);
			abp.calculatePosterior(gwasspruned);

			System.out.println("repeating on pruned snps");
			TextFile tfo2 = new TextFile(out + "-woTopEQTLLDSNPs.txt", TextFile.W);
			tfo2.writeln("SNP\tChr\tPos\tEQTL-P\tEQTL-ABP\tGWAS-P\tGWAS-ABP\tRsq-TopEQTLSNP(" + indexAssocEQTL.getSnp().getName() + ")\tRsqTopGWASSNP(" + indexAssocGWAS.getSnp().getName() + ")");
			for (int i = 0; i < eqtlspruned.size(); i++) {
				AssociationResult eqtlresult = eqtlspruned.get(i);
				AssociationResult gwasresult = gwasspruned.get(i);
				VCFVariant v = variantMap.get(eqtlresult.getSnp().getName());
				double rve = Double.NaN;
				double rvg = Double.NaN;
				if (v != null) {
					Pair<Double, Double> rg = ldcalc.getLD(v, newTopVariantGWAS);
					if (rg != null) {
						rvg = rg.getRight();
					}
					Pair<Double, Double> re = ldcalc.getLD(v, newTopVariantEQTL);
					if (re != null) {
						rve = re.getRight();
					}
				}
				String outln = eqtlresult.getSnp().getName() + "\t" + eqtlresult.getSnp().getChromosome().getNumber() + "\t" + eqtlresult.getSnp().getStart()
						+ "\t" + eqtlresult.getPval() + "\t" + eqtlresult.getPosterior()
						+ "\t" + gwasresult.getPval() + "\t" + gwasresult.getPosterior()
						+ "\t" + rve + "\t" + rvg;
				tfo2.writeln(outln);
			}
			tfo2.close();


		} else {
			System.out.println("Error: snps not found.");
		}

	}

}
