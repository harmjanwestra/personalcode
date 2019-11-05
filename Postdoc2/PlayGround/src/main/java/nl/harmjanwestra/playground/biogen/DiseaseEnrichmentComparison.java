package nl.harmjanwestra.playground.biogen;

import nl.harmjanwestra.playground.legacy.vcf.DetermineLD;
import nl.harmjanwestra.playground.legacy.vcf.VCFTabix;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.IntStream;

public class DiseaseEnrichmentComparison {


	public static void main(String[] args) {
		DiseaseEnrichmentComparison c = new DiseaseEnrichmentComparison();
		String gwascatalog = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e96_r2019-10-14.tsv.gz";
		String snpmap = "D:\\Sync\\SyncThing\\Data\\Ref\\dbsnp\\SNPMappings-dbsnp151-b38.txt.gz";
		String gwascatalogpositions = "D:\\Sync\\SyncThing\\Data\\Ref\\gwascatalog\\gwas_catalog_v1.0-associations_e96_r2019-10-14-b38positions.txt.gz";


		try {
			/* procedure:
		       0. read gwas catalog
		       1. map snps to positions (b38)
		    */
//			c.mapGWASCatalogToSNPPositions(gwascatalog, gwascatalogpositions, snpmap, 5e-8);

			String eqtlgengenelist = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-23-DiseaseEnrichmentComparison\\eqtlgen-annotation_meta_all_2017-08-17.txt";
			String metabraingenelist = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-23-DiseaseEnrichmentComparison\\2019-06-05-protein_coding_genes_gencode24.txt";
			String overlapout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-23-DiseaseEnrichmentComparison\\eqtlgen-metabrain-geneoverlap.txt";
//			c.overlapGeneLists(eqtlgengenelist, metabraingenelist, overlapout);


			String metabrainb38eqtls = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-23-DiseaseEnrichmentComparison\\2019-09-12-Cis-EurAndAFR-Iteration1-eQTLProbesFDR0.05-ProbeLevel.txt.gz";
			String eqtlgenb38eqtls = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2018-01-31-eqtlgen-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38-topfxpergene.txt.gz";
			String allowedgenesfile = overlapout;
			String gwassnpfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-23-DiseaseEnrichmentComparison\\gwas_catalog_v1.0-associations_e96_r2019-10-14-b38positions.txt.gz";
			int minnrsnps = 10;

			int ldDistanceThreshold = 1000000;
			double ldRSqThreshold = 0.2;

			String vcfTabixSampleLimitFile = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p35va-europeans.txt";
			String vcfTabixTemplate = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-b38\\ALL.chrCHR.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz";
			String outputfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-10-23-DiseaseEnrichmentComparison\\2019-10-23-enrichmentTest.txt";

			c.determineEnrichment(metabrainb38eqtls, eqtlgenb38eqtls, allowedgenesfile, gwassnpfile,
					minnrsnps, ldDistanceThreshold, ldRSqThreshold,
					vcfTabixSampleLimitFile, vcfTabixTemplate, outputfile);

		} catch (IOException e) {
			e.printStackTrace();
		}


	}


	private class EQTL {
		String gene;
		String snpchr;
		Integer snppos;
		String snprs;

	}

	HashMap<Chromosome, VCFCache> cache = new HashMap<>();

	public class VCFCache {
		HashMap<String, VCFVariant> variantcache = new HashMap();

		public void addVariant(String snp, VCFVariant variant) {
			variantcache.put(snp, variant);
		}

		public VCFVariant get(String querySNP) {
			return variantcache.get(querySNP);
		}
	}

	HashMap<String, VCFVariant> vcfCache = new HashMap<>();
	HashSet<String> nonExistingVariants = new HashSet<String>();

	public void determineEnrichment(String metabrainb38eqtls, String eqtlgenb38eqtls, String allowedgenesfile, String gwassnpfile, int minnrsnps,
									int ldDistanceThreshold, double ldRSqThreshold,
									String vcfTabixSampleLimitFile, String vcfTabixTemplate, String outputfile) throws IOException {

		// read list of allowed genes
		HashSet<String> allowedGenes = new HashSet<String>();
		TextFile tf = new TextFile(allowedgenesfile, TextFile.R);
		allowedGenes.addAll(tf.readAsArrayList());
		tf.close();

		System.out.println("Allowed genes: " + allowedGenes.size());

		// read eQTLs
		HashMap<Chromosome, ArrayList<EQTL>> metabraineqtls = readEQTLs(metabrainb38eqtls, allowedGenes);
		HashMap<Chromosome, ArrayList<EQTL>> eqtlgeneqtls = readEQTLs(eqtlgenb38eqtls, allowedGenes);


		HashSet<String> uniqueDiseaseSNPs = getUniqueDiseaseSNPs(gwassnpfile, minnrsnps);
		System.out.println("There are " + uniqueDiseaseSNPs.size() + " unique disease snps.");

		// iterate diseases
		TextFile tf2 = new TextFile(gwassnpfile, TextFile.R);
		tf2.readLine();
		String[] elems = tf2.readLineElems(TextFile.tab);
		TextFile output = new TextFile(outputfile, TextFile.W);
		String header = "Disease" +
				"\tDiseaseSNPs" +
				"\tNonDiseaseSNPs" +
				"\tdiseaseSNPsLinkedMetaBrain" +
				"\tnonDiseaseSNPsLinkedMetaBrain" +
				"\tMetaBrainFisherExactP" +
				"\tMetaBrainFisherExactZ" +
				"\tdiseaseSNPsLinkedEqtlgen" +
				"\tnonDiseaseSNPsLinkedEqtlgen" +
				"\tEQTLGenFisherExactP" +
				"\tEQTLGenFisherExactZ" +
				"\tZDifference" +
				"\tZDifferenceP";
		output.writeln(header);
		System.out.println(header);
		while (elems != null) {

			String disease = elems[0];
			String[] snps = elems[2].split(";");
			HashSet<String> diseaseSNPs = new HashSet<String>();
			for (String snp : snps) {
				String[] snpelems = snp.split(":");
				// remove non-autosomal disease snps...
				Chromosome chr = Chromosome.parseChr(snpelems[0]);
				if (chr.isAutosome()) {
					diseaseSNPs.add(snp);
				}
			}

			if (diseaseSNPs.size() > minnrsnps) {
				// get set of non-disease snps
				HashSet<String> nonDiseaseSNPs = subtract(diseaseSNPs, uniqueDiseaseSNPs);

				// group per chromosome
				HashMap<Chromosome, ArrayList<String>> diseasSNPsPerChr = groupSNPsPerChromosome(diseaseSNPs);
				HashMap<Chromosome, ArrayList<String>> nonDiseasSNPsPerChr = groupSNPsPerChromosome(nonDiseaseSNPs);

				// remove linked variants from nonDiseaseSet
				System.out.println(disease + " unlinking snps.");
				nonDiseasSNPsPerChr = unlink(diseasSNPsPerChr, nonDiseasSNPsPerChr, ldDistanceThreshold, ldRSqThreshold, vcfTabixSampleLimitFile, vcfTabixTemplate);

				// prune the remaining sets of variants
				System.out.println(disease + " pruning non disease snps");
				nonDiseasSNPsPerChr = prune(nonDiseasSNPsPerChr, ldDistanceThreshold, ldRSqThreshold, vcfTabixSampleLimitFile, vcfTabixTemplate);
				System.out.println(disease + " pruning non disease snps");
				diseasSNPsPerChr = prune(diseasSNPsPerChr, ldDistanceThreshold, ldRSqThreshold, vcfTabixSampleLimitFile, vcfTabixTemplate);

				// test snps for linkage with eQTLs
				diseaseSNPs = asSet(diseasSNPsPerChr);
				nonDiseaseSNPs = asSet(nonDiseasSNPsPerChr);
				System.out.println(disease + " linking disease snps to metabrain");
				HashSet<String> linkedDiseaseSNPsMetaBrain = link(diseasSNPsPerChr, metabraineqtls, ldDistanceThreshold, ldRSqThreshold, vcfTabixSampleLimitFile, vcfTabixTemplate);
				System.out.println(disease + " linking non disease snps to metabrain");
				HashSet<String> linkedNonDiseaseSNPsMetaBrain = link(nonDiseasSNPsPerChr, metabraineqtls, ldDistanceThreshold, ldRSqThreshold, vcfTabixSampleLimitFile, vcfTabixTemplate);

				System.out.println(disease + " linking disease snps to eqtlgen");
				HashSet<String> linkedDiseaseSNPsEQTLgen = link(diseasSNPsPerChr, eqtlgeneqtls, ldDistanceThreshold, ldRSqThreshold, vcfTabixSampleLimitFile, vcfTabixTemplate);

				System.out.println(disease + " linking non disease snps to eqtlgen");
				HashSet<String> linkedNonDiseaseSNPsEQTLgen = link(nonDiseasSNPsPerChr, eqtlgeneqtls, ldDistanceThreshold, ldRSqThreshold, vcfTabixSampleLimitFile, vcfTabixTemplate);

				FisherExactTest fet = new FisherExactTest();
				double metabrainFisherExactP = fet.getFisherPValue(linkedDiseaseSNPsMetaBrain.size(), diseaseSNPs.size() - linkedDiseaseSNPsMetaBrain.size(),
						linkedNonDiseaseSNPsMetaBrain.size(), nonDiseaseSNPs.size() - linkedNonDiseaseSNPsMetaBrain.size());
				double metabrainFisherExactZ = ZScores.pToZ(metabrainFisherExactP);
				double eqtlgenFisherExactP = fet.getFisherPValue(linkedDiseaseSNPsEQTLgen.size(), diseaseSNPs.size() - linkedDiseaseSNPsEQTLgen.size(),
						linkedNonDiseaseSNPsEQTLgen.size(), nonDiseaseSNPs.size() - linkedNonDiseaseSNPsEQTLgen.size());
				double eqtlgenFisherExactZ = ZScores.zToP(eqtlgenFisherExactP);

				double zdiff = metabrainFisherExactZ - eqtlgenFisherExactZ;
				double zdiffp = ZScores.zToP(zdiff);

				String outputln = elems[0]
						+ "\t" + diseaseSNPs.size()
						+ "\t" + nonDiseaseSNPs.size()
						+ "\t" + linkedDiseaseSNPsMetaBrain.size()
						+ "\t" + linkedNonDiseaseSNPsMetaBrain.size()
						+ "\t" + metabrainFisherExactP
						+ "\t" + metabrainFisherExactZ
						+ "\t" + linkedDiseaseSNPsEQTLgen.size()
						+ "\t" + linkedNonDiseaseSNPsEQTLgen.size()
						+ "\t" + eqtlgenFisherExactP
						+ "\t" + eqtlgenFisherExactZ
						+ "\t" + zdiff
						+ "\t" + zdiffp;
				System.out.println(outputln);
				output.writeln(outputln);
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		output.close();
	}

	private HashSet<String> asSet(HashMap<Chromosome, ArrayList<String>> diseasSNPsPerChr) {
		HashSet<String> output = new HashSet<String>();

		for (Chromosome chr : diseasSNPsPerChr.keySet()) {
			output.addAll(diseasSNPsPerChr.get(chr));
		}
		return output;
	}

	private HashSet<String> link(HashMap<Chromosome, ArrayList<String>> querySNPsPerChr, HashMap<Chromosome, ArrayList<EQTL>> eQTLSNPsPerChr,
								 int ldDistanceThreshold, double ldRSqThreshold,
								 String vcfTabixSampleLimitFile, String vcfTabixTemplate) throws IOException {

		HashSet<String> fullSetOfSnpsThatMayBeLinked = new HashSet<String>();

		Set<Chromosome> keyset = querySNPsPerChr.keySet();
		ArrayList<Chromosome> chrlist = new ArrayList<>();
		chrlist.addAll(keyset);

		IntStream.range(0, chrlist.size()).parallel().forEach(i -> {
			Chromosome chr = chrlist.get(i);
			ArrayList<String> querysnps = querySNPsPerChr.get(chr);
			ArrayList<EQTL> eqtlsnps = eQTLSNPsPerChr.get(chr);
			String filename = vcfTabixTemplate.replaceAll("CHR", "" + chr.getNumber());
			VCFTabix tabix = null;
			try {
				tabix = new VCFTabix(filename);

				HashSet<String> snpsThatMayBeLinked = new HashSet<String>();
				boolean[] sampleLimit = tabix.getSampleFilter(vcfTabixSampleLimitFile);
				HashSet<String> setThatMayBeLinked = new HashSet<>();

				if (eqtlsnps != null) {
					for (String querysnp : querysnps) {
						String[] querysnpelems = querysnp.split(":");
						String querysnpposstr = querysnpelems[1];
						Integer querysnppos = Integer.parseInt(querysnpposstr);

						for (EQTL e : eqtlsnps) {
							int distance = Math.abs(e.snppos - querysnppos);
							if (distance == 0) {
								snpsThatMayBeLinked.add(querysnp);
							} else if (distance < ldDistanceThreshold) {
								// test for LD
								VCFVariant diseaseVariant = getVariant(querysnp, tabix, sampleLimit); // tabix.getVariant(asFeature(querysnp), sampleLimit);
								VCFVariant nonDiseaseVariant = getVariant(e.snprs, tabix, sampleLimit); // tabix.getVariant(asFeature(e.snprs), sampleLimit);

								// remove snps that are not in the reference dataset.
								if (diseaseVariant != null && nonDiseaseVariant != null) {
									DetermineLD ld = new DetermineLD();
									Pair<Double, Double> ldobj = ld.getLD(diseaseVariant, nonDiseaseVariant);
									if (ldobj != null && ldobj.getRight() > ldRSqThreshold) {
										setThatMayBeLinked.add(querysnp);
									}
								}
							}
						}
					}
				}
				synchronized (fullSetOfSnpsThatMayBeLinked) {
					fullSetOfSnpsThatMayBeLinked.addAll(snpsThatMayBeLinked);
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		});
		return fullSetOfSnpsThatMayBeLinked;
	}

	private HashMap<Chromosome, ArrayList<String>> prune(HashMap<Chromosome, ArrayList<String>> diseasSNPsPerChr,
														 int ldDistanceThreshold, double ldRSqThreshold,
														 String vcfTabixSampleLimitFile, String vcfTabixTemplate) throws IOException {

		Set<Chromosome> keyset = diseasSNPsPerChr.keySet();
		ArrayList<Chromosome> chrlist = new ArrayList<>();
		chrlist.addAll(keyset);
		IntStream.range(0, chrlist.size()).parallel().forEach(i -> {
			Chromosome chr = chrlist.get(i);
			ArrayList<String> diseaseSNPs = diseasSNPsPerChr.get(chr);

			if (diseaseSNPs.size() > 1) {
				// open TABIX reference
				String filename = vcfTabixTemplate.replaceAll("CHR", "" + chr.getNumber());
				VCFTabix tabix = null;
				try {
					tabix = new VCFTabix(filename);

					boolean[] sampleLimit = tabix.getSampleFilter(vcfTabixSampleLimitFile);
					HashSet<String> setThatMayBeLinked = new HashSet<>();

					for (int s = 0; s < diseaseSNPs.size(); s++) {
						String diseaseSNP = diseaseSNPs.get(s);
						String[] diseaseSNPElems = diseaseSNP.split(":");
						String diseaseSNPPosStr = diseaseSNPElems[1];
						Integer diseaseSNPPos = Integer.parseInt(diseaseSNPPosStr);

						for (int s2 = s + 1; s < diseaseSNPs.size(); s++) {
							String diseaseSNP2 = diseaseSNPs.get(s2);
							if (!setThatMayBeLinked.contains(diseaseSNP2)) {
								String[] diseaseSNPElems2 = diseaseSNP2.split(":");
								String diseaseSNPPosStr2 = diseaseSNPElems2[1];
								Integer diseaseSNPPos2 = Integer.parseInt(diseaseSNPPosStr2);
								if (Math.abs(diseaseSNPPos - diseaseSNPPos2) < ldDistanceThreshold) {
									// check for linkage

									// test for LD
									VCFVariant diseaseVariant = getVariant(diseaseSNP, tabix, sampleLimit); //tabix.getVariant(asFeature(diseaseSNP), sampleLimit);
									VCFVariant nonDiseaseVariant = getVariant(diseaseSNP2, tabix, sampleLimit); // tabix.getVariant(asFeature(diseaseSNP2), sampleLimit);

									// remove snps that are not in the reference dataset.
									if (diseaseVariant == null) {
										setThatMayBeLinked.add(diseaseSNP);
									}
									if (nonDiseaseVariant == null) {
										setThatMayBeLinked.add(diseaseSNP2);
									}

									if (diseaseVariant != null && nonDiseaseVariant != null) {
										DetermineLD ld = new DetermineLD();
										Pair<Double, Double> ldobj = ld.getLD(diseaseVariant, nonDiseaseVariant);
										if (ldobj != null && ldobj.getRight() > ldRSqThreshold) {
											setThatMayBeLinked.add(diseaseSNP2);
										}
									}
								}
							}
						}
					}

					HashSet<String> diseaseSNPset = new HashSet<String>();
					diseaseSNPset.addAll(diseaseSNPs);
					diseaseSNPset = subtract(setThatMayBeLinked, diseaseSNPset);
					ArrayList<String> nonDiseasSNPList = new ArrayList<>();
					nonDiseasSNPList.addAll(diseaseSNPset);
					diseasSNPsPerChr.put(chr, nonDiseasSNPList);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}

		});

		return diseasSNPsPerChr;
	}

	private HashMap<Chromosome, ArrayList<String>> unlink(HashMap<Chromosome, ArrayList<String>> diseasSNPsPerChr,
														  HashMap<Chromosome, ArrayList<String>> nonDiseasSNPsPerChr,
														  int ldDistanceThreshold, double ldRSqThreshold,
														  String vcfTabixSampleLimitFile, String vcfTabixTemplate) throws IOException {
		// TODO: parallelize
		Set<Chromosome> keyset = diseasSNPsPerChr.keySet();
		ArrayList<Chromosome> chrlist = new ArrayList<>();
		chrlist.addAll(keyset);
		IntStream.range(0, chrlist.size()).parallel().forEach(i -> {
					Chromosome chr = chrlist.get(i);
					ArrayList<String> diseaseSNPs = diseasSNPsPerChr.get(chr);
					ArrayList<String> nonDiseaseSNPs = nonDiseasSNPsPerChr.get(chr);

					// open TABIX reference
					String filename = vcfTabixTemplate.replaceAll("CHR", "" + chr.getNumber());
					VCFTabix tabix = null;
					try {
						tabix = new VCFTabix(filename);
						boolean[] sampleLimit = tabix.getSampleFilter(vcfTabixSampleLimitFile);
						if (nonDiseaseSNPs != null) {

							// for each disease snp
							HashSet<String> setThatMayBeLinked = new HashSet<>();
							for (String diseaseSNP : diseaseSNPs) {
								String[] diseaseSNPElems = diseaseSNP.split(":");
								String diseaseSNPPosStr = diseaseSNPElems[1];
								Integer diseaseSNPPos = Integer.parseInt(diseaseSNPPosStr);
								// check whether there is something in the non-disease snpset that may be linked.
								for (String nonDiseaseSNP : nonDiseaseSNPs) {
									if (!setThatMayBeLinked.contains(nonDiseaseSNP)) {
										if (diseaseSNP.equals(nonDiseaseSNP)) {
											setThatMayBeLinked.add(diseaseSNP);
										} else {
											String[] nonDiseaseSNPElems = nonDiseaseSNP.split(":");
											String nonDiseaseSNPPosStr = nonDiseaseSNPElems[1];
											Integer nonDiseaseSNPPos = Integer.parseInt(nonDiseaseSNPPosStr);

											// if the SNPs are within ldDistanceThreshold
											if (Math.abs(diseaseSNPPos - nonDiseaseSNPPos) < ldDistanceThreshold) {
												// test for LD
												VCFVariant diseaseVariant = getVariant(diseaseSNP, tabix, sampleLimit);
												VCFVariant nonDiseaseVariant = getVariant(nonDiseaseSNP, tabix, sampleLimit); // tabix.getVariant(asFeature(nonDiseaseSNP), sampleLimit);

												if (diseaseVariant != null && nonDiseaseVariant != null) {
													DetermineLD ld = new DetermineLD();
													Pair<Double, Double> ldobj = ld.getLD(diseaseVariant, nonDiseaseVariant);
													if (ldobj != null && ldobj.getRight() > ldRSqThreshold) {
														setThatMayBeLinked.add(nonDiseaseSNP);
													}
												}
											}
										}
									}
								}
							}
							HashSet<String> nonDiseasSNPSet = new HashSet<String>();
							nonDiseasSNPSet.addAll(nonDiseaseSNPs);
							nonDiseasSNPSet = subtract(setThatMayBeLinked, nonDiseasSNPSet);
							ArrayList<String> nonDiseasSNPList = new ArrayList<>();
							nonDiseasSNPList.addAll(nonDiseasSNPSet);
							nonDiseasSNPsPerChr.put(chr, nonDiseasSNPList);
						}
					} catch (IOException e) {
						e.printStackTrace();
					}
				}

		);


		return nonDiseasSNPsPerChr;
	}

	private VCFVariant getVariant(String querySNP, VCFTabix tabix, boolean[] sampleLimit) throws IOException {
		if (nonExistingVariants.contains(querySNP)) {
			return null;
		}
		Chromosome chr = asChr(querySNP);
		VCFCache subcache = cache.get(chr);
		if (subcache == null) {
			subcache = new VCFCache();
		}
		synchronized (cache) {
			cache.put(chr, subcache);
		}

		VCFVariant var = subcache.get(querySNP);
		if (var != null) {
			return var;
		} else {
			var = tabix.getVariant(asFeature(querySNP), sampleLimit);
			if (var != null) {
				subcache.addVariant(querySNP, var);
				return var;
			} else {
				nonExistingVariants.add(querySNP);
			}
		}
		return null;
	}

	private Chromosome asChr(String querySNP) {

		String[] elems = querySNP.split(":");
		return Chromosome.parseChr(elems[0]);

	}

	private Feature asFeature(String diseaseSNP) {
		Feature f = new Feature();
		String[] elems = diseaseSNP.split(":");
		if (elems.length > 2) {
			f.setName(elems[2]);
		}
		f.setChromosome(Chromosome.parseChr(elems[0]));
		f.setStart(Integer.parseInt(elems[1]));
		f.setStop(Integer.parseInt(elems[1]));
		return f;
	}

	private HashMap<Chromosome, ArrayList<String>> groupSNPsPerChromosome(HashSet<String> snps) {

		HashMap<Chromosome, ArrayList<String>> output = new HashMap<>();
		for (String snp : snps) {
			String[] elems = snp.split(":");
			Chromosome chr = Chromosome.parseChr(elems[0]);
			ArrayList<String> list = output.get(chr);
			if (list == null) {
				list = new ArrayList<>();
			}
			list.add(snp);
			output.put(chr, list);
		}
		return output;

	}

	private HashSet<String> subtract(HashSet<String> diseaseSNPs, HashSet<String> uniqueDiseaseSNPs) {
		HashSet<String> output = new HashSet<String>();
		for (String snp : uniqueDiseaseSNPs) {
			if (!diseaseSNPs.contains(snp)) {
				output.add(snp);
			}
		}
		return output;
	}

	private HashSet<String> getUniqueDiseaseSNPs(String gwassnpfile, int minnrsnps) throws IOException {
		HashSet<String> output = new HashSet<>();
		TextFile tf2 = new TextFile(gwassnpfile, TextFile.R);
		tf2.readLine();
		String[] elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {

			String[] snps = elems[2].split(";");
			if (snps.length > minnrsnps) {
				// test snps for linkage in both sets
				for (String snp : snps) {
					String[] snpelems = snp.split(":");

					// remove non-autosomal variants
					Chromosome chr = Chromosome.parseChr(snpelems[0]);
					if (chr.isAutosome()) {
						output.add(snp);
					}
				}
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		return output;
	}

	private HashMap<Chromosome, ArrayList<EQTL>> readEQTLs(String eqtlfile, HashSet<String> allowedGenes) throws
			IOException {
		TextFile tf = new TextFile(eqtlfile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<Chromosome, ArrayList<EQTL>> output = new HashMap<>();
		int ctr = 0;
		while (elems != null) {

			String gene = elems[4];
			if (allowedGenes.contains(gene)) {
				Chromosome genechr = Chromosome.parseChr(elems[5]);

				String snprs = elems[1];
				String snpchr = elems[2];
				String snppos = elems[3];
				Chromosome chr = Chromosome.parseChr(snpchr);
				ArrayList<EQTL> chrEQTL = output.get(chr);
				if (chrEQTL == null) {
					chrEQTL = new ArrayList<>();
				}

				EQTL e = new EQTL();
				e.snpchr = snpchr;
				e.snppos = Integer.parseInt(snppos);
				if (snprs.startsWith("rs")) {
					snprs = snpchr + ":" + snppos + ":" + snprs;
				}
				e.snprs = snprs;
				e.gene = gene;

				chrEQTL.add(e);
				ctr++;

				output.put(chr, chrEQTL);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println(eqtlfile + " has " + ctr + " eqtls.");
		return output;
	}

	private void overlapGeneLists(String eqtlgengenelist, String metabraingenelist, String overlapout) throws
			IOException {
		HashSet<String> query = new HashSet<String>();
		TextFile tf = new TextFile(eqtlgengenelist, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String ensg = elems[0];
			query.add(ensg);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		TextFile tf2 = new TextFile(metabraingenelist, TextFile.R);
		TextFile outf = new TextFile(overlapout, TextFile.W);
		String ln = tf2.readLine();
		while (ln != null) {
			String ensg = ln.split("\\.")[0];
			if (query.contains(ensg)) {
				outf.writeln(ln);
			}
			ln = tf2.readLine();
		}
		tf2.close();
		outf.close();

	}

	public void mapGWASCatalogToSNPPositions(String catalog, String output, String snpmap, double threshold) throws
			IOException {
		GWASCatalog c = new GWASCatalog(catalog);

		HashSet<String> rsids = new HashSet<String>();

		GWASTrait[] traits = c.getTraits();

		for (GWASTrait t : traits) {
			GWASSNP[] snps = t.getSNPs(threshold);
			int inc = 0;
			int excl = 0;
			for (GWASSNP snp : snps) {
				if (snp.getName().contains("rs")) {
					rsids.add(snp.getName());
					inc++;
				} else {
					System.out.println("Chose to exclude: " + snp.getName());
					excl++;
				}
			}
			System.out.println(t.getName() + "\t" + inc + " included, " + excl + " excluded variants.");
		}

		System.out.println(rsids.size() + " rs ids found..");
		TextFile tf = new TextFile(snpmap, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, String> rsIdToPos = new HashMap<>();
		int ctr = 0;
		while (elems != null) {
			String rs = elems[2];
			if (rsids.contains(rs)) {
				rsIdToPos.put(rs, elems[0] + ":" + elems[1]);
			}
			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(ctr + " snp rsids parsed; " + rsIdToPos.size() + " / " + rsids.size() + " mapped.");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println("Could recover " + rsIdToPos.size() + " / " + rsids.size() + " rs ids");

		TextFile out = new TextFile(output, TextFile.W);

		out.writeln("Trait\tnrSNPs\tSNPs");
		for (GWASTrait t : traits) {
			GWASSNP[] snps = t.getSNPs(threshold);
			String ln = t.getName();
			ArrayList<String> positions = new ArrayList<String>();
			for (GWASSNP snp : snps) {
				if (snp.getName().contains("rs")) {
					String pos = rsIdToPos.get(snp.getName());
					if (pos != null) {
						pos = pos + ":" + snp.getName();
						positions.add(pos);
					}
				}
			}
			if (!positions.isEmpty()) {
				out.writeln(ln + "\t" + positions.size() + "\t" + Strings.concat(positions, Strings.semicolon));
			}
		}
		out.close();


	}


}
