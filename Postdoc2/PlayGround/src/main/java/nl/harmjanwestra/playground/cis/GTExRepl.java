package nl.harmjanwestra.playground.cis;

import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class GTExRepl {
	
	public static void main(String[] args) {
//		String fullcis = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\eQTLsFDR-ProbeLevel.txt.gz";
//		String gtextarball = "D:\\Sync\\SyncThing\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v7_eQTL\\";
//		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-06-05-cis-gtexrepl\\ConcordanceBetweenTopGTExAndEQTLFullCis.txt";
		boolean maponposition = true;
		GTExRepl r = new GTExRepl();

//		args = new String[]{
//				"D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz",
//				"D:\\Sync\\SyncThing\\Data\\eQTLs\\GTEx\\test\\",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-06-05-cis-gtexrepl\\test\\output.txt"
//		};
		
		try {
			if (args.length < 3) {
				System.out.println("Usage: fullcis gtexloc output");
			} else {
				r.determineConcordanceBetweenTopGTExAndFullCis(args[0], args[1], maponposition, args[2]);
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void determineConcordanceBetweenTopGTExAndFullCis(String fullcis, String gtextarball, boolean maponposition, String outputfile) throws IOException {
		ArrayList<Pair<String, ArrayList<EQTL>>> topgtexfx = getGTExTopFx(gtextarball,
				"egenes.txt.gz",
				null,
				false);
		
		// make a list of snp/gene combos to read from eqtlgen
		HashSet<String> toread = new HashSet<>();
		for (Pair<String, ArrayList<EQTL>> pair : topgtexfx) {
			for (EQTL e : pair.getRight()) {
				if (maponposition) {
					String query = Strings.cache(e.getRsChr() + "_" + e.getRsChrPos() + "_" + e.getProbe());
					toread.add(query);
				} else {
					String query = Strings.cache(e.getRsName() + "_" + e.getProbe());
					toread.add(query);
				}
				
			}
		}
		System.out.println(toread.size() + " total snp/gene combos to read in");
		HashMap<String, EQTL> eqtlmap = new HashMap<>();
		int ctr = 0;
		HashSet<String> refgenes = new HashSet<>();
		HashSet<String> refsnps = new HashSet<>();
		TextFile tf = new TextFile(fullcis, TextFile.R, 1048576);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = Strings.cache(elems[1]);
			String gene = Strings.cache(elems[4]);
			String query = null;
			if (maponposition) {
				query = Strings.cache(elems[2] + "_" + elems[3] + "_" + elems[4]);
			} else {
				query = elems[1] + "_" + elems[4];
			}
			refgenes.add(gene);
			refsnps.add(snp);
			if (toread.contains(query)) {
				EQTL subqtl = new EQTL();
				//
				if (elems.length >= 21) {
					subqtl.setRsChr(Byte.parseByte(elems[2]));
					subqtl.setRsChrPos(Integer.parseInt(elems[3]));
					subqtl.setRsName(snp);
					subqtl.setProbe(gene);
					subqtl.setProbeChr(Byte.parseByte(elems[5]));
					subqtl.setProbeChrPos(Integer.parseInt(elems[6]));
					subqtl.setAlleleAssessed(Strings.cache(elems[9]));
					subqtl.setAlleles(Strings.cache(elems[8]));
					subqtl.setZscore(Double.parseDouble(elems[10]));
					subqtl.setFDR(Double.parseDouble(elems[20]));
					eqtlmap.put(Strings.cache(query), subqtl);
				}
			}
			if (ctr % 100000 == 0) {
				System.out.print("\r" + ctr + "\tqtls processed. " + eqtlmap.size() + " combos found so far");
			}
			ctr++;
//			if (eqtlmap.size() >= 3000) {
//				break;
//			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println();
		System.out.println(eqtlmap.size() + " eqtls read from eqtlgen");
		
		TextFile notfound = new TextFile(outputfile + "-notfound.txt", TextFile.W);
		notfound.writeln("snp\tsnpchr\tsnppos\tgene\tsnppresent\tgenepresent");
		
		// now check concordance
		TextFile out = new TextFile(outputfile, TextFile.W);
		String header = "tissue" +
				"\tnrInGTEx" +
				"\tnrGTExSNPsInEQTLGen" +
				"\tnrGTExGenesInEQTLGen" +
				"\tnrGTExeQTLsPossibleInEQTLGen" +
				"\tnrSignificantGTEx" +
				"\tShared" +
				"\tSharedSignificantInGTEx" +
				"\tSharedSignificantInEQTLGen" +
				"\tSharedSignificantInBoth" +
				"\tConcordant" +
				"\tConcordantSignificantInGTEx" +
				"\tConcordantSignificantInEQTLGen" +
				"\tConcordantSignificantInBoth";
		
		out.writeln(header);
		
		for (Pair<String, ArrayList<EQTL>> pair : topgtexfx) {
			
			String tissue = pair.getLeft();
			int nringtex = pair.getRight().size();
			
			
			int nrsigingtex = 0;
			int shared = 0;
			int sharedsiggtex = 0;
			int sharedsigeqtlgen = 0;
			int sharedsigboth = 0;
			int concordant = 0;
			int concordantsiggtex = 0;
			int concordantsigeqtlgen = 0;
			int concordantsigboth = 0;
			int nrsharedgenes = 0;
			int nrsharedsnps = 0;
			
			int nrpossiblydetectedineqtlgen = 0;
			for (EQTL e : pair.getRight()) {
				if (refsnps.contains(e.getRsName()) || refgenes.contains(e.getProbe())) {
					nrpossiblydetectedineqtlgen++;
				}
				if (refsnps.contains(e.getRsName())) {
					nrsharedsnps++;
				}
				if (refgenes.contains(e.getProbe())) {
					nrsharedgenes++;
				}
			}
			TextFile nonc = new TextFile(outputfile + tissue + "-nonconcordant.txt.gz", TextFile.W);
			for (EQTL e : pair.getRight()) {
				String query = null;
				if (maponposition) {
					query = Strings.cache(e.getRsChr() + "_" + e.getRsChrPos() + "_" + e.getProbe());
				} else {
					query = e.getRsName() + "_" + e.getProbe();
				}
				
				
				EQTL ref = eqtlmap.get(query);
				if (e.getFDR() < 0.05) {
					nrsigingtex++;
				}
				if (ref == null) {
					String notfoundstr = tissue + "\t" + e.getRsName() + "\t" + e.getRsChr() + "\t" + e.getRsChrPos() + "\t" + e.getProbe() + "\t" + refsnps.contains(e.getRsName()) + "\t" + refgenes.contains(e.getProbe());
					notfound.writeln(notfoundstr);
				} else {
					
					
					shared++;
					if (e.getFDR() < 0.05) {
						sharedsiggtex++;
					}
					// check concordance
					if (ref.getFDR() < 0.05) {
						sharedsigeqtlgen++;
						if (e.getFDR() < 0.05) {
							sharedsigboth++;
						}
					}
					
					Boolean flipalleles = BaseAnnot.flipalleles(ref.getAlleles(), ref.getAlleleAssessed(), e.getAlleles(), e.getAlleleAssessed());
					boolean concordance = false;
					double z = e.getZscore();
					if (flipalleles != null) {
						if (flipalleles) {
							e.setZscore(-z);
						}
						
						if ((e.getZscore() >= 0 && ref.getZscore() >= 0) || (e.getZscore() < 0 && ref.getZscore() < 0)) {
							concordance = true;
						}
					}
					
					String outln = e.getProbe()
							+ "\t" + e.getRsName()
							+ "\t" + e.getRsChr()
							+ "\t" + e.getRsChrPos()
							+ "\t" + e.getAlleles()
							+ "\t" + e.getAlleleAssessed()
							+ "\t" + z
							+ "\t" + e.getFDR()
							+ "\t" + ref.getRsName()
							+ "\t" + ref.getRsChr()
							+ "\t" + ref.getRsChrPos()
							+ "\t" + ref.getAlleles()
							+ "\t" + ref.getAlleleAssessed()
							+ "\t" + ref.getZscore()
							+ "\t" + ref.getFDR();
					
					if (concordance) {
//						System.out.println(outln);
						concordant++;
						if (e.getFDR() < 0.05) {
							concordantsiggtex++;
							if (ref.getFDR() < 0.05) {
								concordantsigboth++;
							}
						}
						if (ref.getFDR() < 0.05) {
							concordantsigeqtlgen++;
						}
					} else {
						
						nonc.writeln(outln);
						
					}
					
				}
				
			}
			nonc.close();
			String ln = tissue
					+ "\t" + nringtex
					+ "\t" + nrsigingtex
					+ "\t" + shared
					+ "\t" + sharedsiggtex
					+ "\t" + sharedsigeqtlgen
					+ "\t" + sharedsigboth
					+ "\t" + concordant
					+ "\t" + concordantsiggtex
					+ "\t" + concordantsigeqtlgen
					+ "\t" + concordantsigboth;
			out.writeln(ln);
		}
		out.close();
		notfound.close();
	}
	
	public void determineLDBetweenTopCisFx(String cisfxfile, String gtextarball, String ensembl, String tabixfile, String samplefilter) throws IOException {
		
		// read top cis fx
		QTLTextFile tf = new QTLTextFile(cisfxfile, TextFile.R);
		
		HashMap<String, EQTL> topeqtlgencispergene = new HashMap<>();
		
		Iterator<EQTL> it = tf.getEQtlIterator();
		while (it.hasNext()) {
			EQTL e = it.next();
			EQTL e2 = topeqtlgencispergene.get(e.getProbe());
			if (e2 == null) {
				topeqtlgencispergene.put(e.getProbe(), e);
			} else {
				if (Math.abs(e2.getZscore()) < Math.abs(e.getZscore())) {
					topeqtlgencispergene.put(e.getProbe(), e);
				}
			}
		}
		tf.close();
		
		// sort the genes by chr (for easier lookup when calculating LD)
		HashMap<String, Integer> geneMap = new HashMap<>();
		EQTL[] reference = null;
		{
			GTFAnnotation a = new GTFAnnotation(ensembl);
			ArrayList<Gene> selectedGenes = new ArrayList<Gene>();
			for (String key : topeqtlgencispergene.keySet()) {
				Gene g = a.getStrToGene().get(key);
				if (g != null) {
					selectedGenes.add(g);
				}
			}
			Collections.sort(selectedGenes, new FeatureComparator(true));
			
			int ctr = 0;
			reference = new EQTL[selectedGenes.size()];
			for (Gene key : selectedGenes) {
				EQTL e = topeqtlgencispergene.get(key.getName());
				if (e != null) {
					geneMap.put(key.getName(), ctr);
					reference[ctr] = e;
				}
				ctr++;
			}
		}
		
		
		// order the cis-eQTL using the genemap
		ArrayList<Pair<String, ArrayList<EQTL>>> topgtexfx = getGTExTopFx(gtextarball,
				"egenes.txt.gz",
				topeqtlgencispergene.keySet(),
				false
		);
		
		// for convenience, plug into array of arrays
		EQTL[][] testqtls = new EQTL[topgtexfx.size()][];
		String[] tissues = new String[topgtexfx.size()];
		for (int j = 0; j < topgtexfx.size(); j++) {
			Pair<String, ArrayList<EQTL>> pair = topgtexfx.get(j);
			tissues[j] = pair.getLeft();
			ArrayList<EQTL> eqtls = pair.getRight();
			for (int i = 0; i < eqtls.size(); i++) {
				EQTL e = eqtls.get(i);
				Integer geneindex = geneMap.get(e.getProbe());
				if (geneindex != null) {
					testqtls[j][geneindex] = e;
				}
			}
		}
		
		// eQTLs loaded.. perform comparisons
		VCFTabix tabixref = null; // new VCFTabix(tabixfile);
		boolean[] filter = null; //t.getSampleFilter(samplefilter);
		Byte prevchr = -1;
		DetermineLD ldcalc = new DetermineLD();
		for (int i = 0; i < reference.length; i++) {
			EQTL refqtl = reference[i];
			Byte chr = refqtl.getRsChr();
			
			if (tabixref == null || chr > prevchr) {
				tabixref = new VCFTabix(tabixfile);
				filter = tabixref.getSampleFilter(samplefilter);
				prevchr = chr;
			}
			SNPFeature snpfeat1 = new SNPFeature(Chromosome.parseChr("" + chr), refqtl.getRsChrPos(), refqtl.getRsChrPos());
			VCFVariant var1 = tabixref.getVariant(snpfeat1, filter);
			
			for (int j = 0; j < testqtls.length; j++) {
				EQTL testqtl = testqtls[j][i];
				Double rsq = null;
				if (testqtl != null) {
					// check whether the variants are in ld
					if (refqtl.getRsName().equals(testqtl.getRsName()) || refqtl.getRsChrPos() == testqtl.getRsChrPos()) {
						// ld == 1
						rsq = 1d;
					} else {
						// check LD
						SNPFeature snpfeat2 = new SNPFeature(Chromosome.parseChr("" + chr), testqtl.getRsChrPos(), testqtl.getRsChrPos());
						VCFVariant var2 = tabixref.getVariant(snpfeat2, filter);
						nl.harmjanwestra.utilities.legacy.genetica.containers.Pair<Double, Double> ldvals = ldcalc.getLD(var1, var2);
						rsq = ldvals.getRight();
					}
				}
			}
		}
		
		
	}
	
	private ArrayList<Pair<String, ArrayList<EQTL>>> getGTExTopFx(String gtextarball,
																  String filefilter,
																  Set<String> genelimit,
																  boolean includesexchr) throws IOException {
		
		ArrayList<Pair<String, ArrayList<EQTL>>> eqtlspertissue = new ArrayList<>();
		
		// now iterate the GTEx file
//		TarArchiveInputStream tarInput = new TarArchiveInputStream(new GzipCompressorInputStream(new FileInputStream(gtextarball)));
//		TarArchiveEntry currentEntry = tarInput.getNextTarEntry();
//		BufferedReader br = null;
//		StringBuilder sb = new StringBuilder();
		File folder = new File(gtextarball);
		System.out.println(gtextarball + "\texists: " + folder.exists());
		if (!folder.exists()) {
			System.exit(-1);
		}
		File[] listOfFiles = folder.listFiles();
		
		for (File file : listOfFiles) {
			if (file.getName().endsWith(filefilter)) {
				TextFile br = new TextFile(file, TextFile.R); // new BufferedReader(new InputStreamReader(tarInput)); // Read directly from tarInput
				br.readLine();
//				System.out.println("For File = " + currentEntry.getName());
				String tissue = file.getName();
				HashMap<String, EQTL> topeqtlgencispergenegtextissue = new HashMap<>();
				String line;
				while ((line = br.readLine()) != null) {
					String[] elems = line.split("\t");
					/*
0	gene_id
2	gene_chr
3	gene_start
4	gene_end
10	variant_id
12	chr
13	snp_pos
14	ref
15	alt
16	rs_id_dbSNP142_GRCh37p13
22	pval_nominal
23	slope
24	slope_se
27	qval
28	pval_nominal_threshold
					 */
/* v7:
0	gene_id
1	gene_name
2	gene_chr
3	gene_start
4	gene_end
5	strand
6	num_var
7	beta_shape1
8	beta_shape2
9	true_df
10	pval_true_df
11	variant_id
12	tss_distance
13	chr
14	pos
15	ref
16	alt
17	num_alt_per_site
18	rs_id_dbSNP147_GRCh37p13
19	minor_allele_samples
20	minor_allele_count
21	maf
22	ref_factor
23	pval_nominal
24	slope
25	slope_se
26	pval_perm
27	pval_beta
28	qval
29	pval_nominal_threshold
30	log2_aFC
31	log2_aFC_lower
32	log2_aFC_upper

*/
					String geneid = Strings.cache(elems[0].split("\\.")[0]);
					if (genelimit == null || genelimit.contains(geneid)) {
						// check whether this is the top fx
						
						EQTL top = topeqtlgencispergenegtextissue.get(geneid);
						double beta = Double.parseDouble(elems[24]);
						double se = Double.parseDouble(elems[25]);
						double zscore = beta / se;
						if (top == null) {
							// if not use this as the top
							
							Chromosome chr = Chromosome.parseChr(elems[13]);
							Chromosome chr2 = Chromosome.parseChr(elems[2]);
							if (chr.isAutosome() || (!chr.isAutosome() && includesexchr)) {
								if (chr2.isAutosome() || (!chr2.isAutosome() && includesexchr)) {
									EQTL e = new EQTL();
									e.setRsName(Strings.cache(elems[18]));
									e.setRsChr((byte) chr.getNumber());
									e.setRsChrPos(Integer.parseInt(elems[14]));
									e.setProbe(geneid);
									e.setFDR(Double.parseDouble(elems[28]));
									e.setProbeChr((byte) chr2.getNumber());
									e.setAlleles(Strings.cache(elems[15] + "/" + elems[16]));
									e.setAlleleAssessed(Strings.cache(elems[16]));
									e.setZscore(zscore);
									topeqtlgencispergenegtextissue.put(geneid, e);
								}
							}
							
						} else {
							// compare absolute z-scores
							double zref = Math.abs(top.getZscore());
							double ztest = Math.abs(zscore);
							if (zref < ztest) {
								Chromosome chr = Chromosome.parseChr(elems[13]);
								Chromosome chr2 = Chromosome.parseChr(elems[2]);
								if (chr.isAutosome() || (!chr.isAutosome() && includesexchr)) {
									if (chr2.isAutosome() || (!chr2.isAutosome() && includesexchr)) {
										EQTL e = new EQTL();
										e.setRsName(Strings.cache(elems[18]));
										e.setRsChr((byte) chr.getNumber());
										e.setRsChrPos(Integer.parseInt(elems[14]));
										e.setProbe(geneid);
										e.setFDR(Double.parseDouble(elems[28]));
										e.setProbeChr((byte) chr2.getNumber());
										e.setAlleles(Strings.cache(elems[15] + "/" + elems[16]));
										e.setAlleleAssessed(Strings.cache(elems[16]));
										e.setZscore(zscore);
										topeqtlgencispergenegtextissue.put(geneid, e);
									}
								}
							}
						}
					}
				}
//				br.close();
				ArrayList<EQTL> eqtls = new ArrayList<>();
				for (String key : topeqtlgencispergenegtextissue.keySet()) {
					eqtls.add(topeqtlgencispergenegtextissue.get(key));
				}
				eqtlspertissue.add(new Pair<>(tissue, eqtls));
				System.out.println(tissue + "\t" + eqtls.size());
			}
		}
		return eqtlspertissue;
	}
}
