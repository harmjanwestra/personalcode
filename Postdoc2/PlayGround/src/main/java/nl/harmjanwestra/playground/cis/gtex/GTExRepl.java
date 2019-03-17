package nl.harmjanwestra.playground.cis.gtex;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import nl.harmjanwestra.playground.legacy.vcf.DetermineLD;
import nl.harmjanwestra.playground.legacy.vcf.VCFTabix;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import org.apache.commons.io.comparator.NameFileComparator;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.FeatureComparator;
import umcg.genetica.features.Gene;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
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
//				"D:\\Sync\\SyncThing\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v7_eQTL\\",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-06-05-cis-gtexrepl\\test\\output-withannot.txt",
//				"D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-06-12-biosstats\\BIOS_expression_summary_20180515.txt",
//				"D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz"
//		};
//
		try {
			if (args.length < 3) {
				System.out.println("Usage: fullcis gtexloc output [avgexp] [gtf]");
			} else {
				String avgexp = null;
				String gtf = null;
				boolean includeblood = true;
				boolean usebonferroni = true;
				if (args.length > 3) {
					avgexp = args[3];
					gtf = args[4];
				}
				r.determineConcordanceBetweenTopGTExAndFullCis(args[0], args[1], maponposition, args[2], includeblood, usebonferroni);
//				r.signficantGenesPerBin(args[0], args[1], maponposition, args[2], avgexp, gtf);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	
	public void signficantGenesPerBin(String fullcis, String gtextarball, boolean maponposition, String outputfile, String avgexpressionfile, String gtf, boolean includeblood) throws IOException {
		int nrbins = 10;
		boolean usebonferroni = true;
//		Double customthreshold = 1.8e-5;
		Double customthreshold = null;
		HashMap<String, Integer> rankedgenes = null;
		int[] nrgenesperbin = null;
		if (avgexpressionfile != null) {
			Pair<int[], HashMap<String, Integer>> generanks = loadAndRankGenes(avgexpressionfile, nrbins, true, gtf);
			rankedgenes = generanks.getRight();
			nrgenesperbin = generanks.getLeft();
			System.out.println(rankedgenes.size() + " genes ranked from " + avgexpressionfile);
		}
		
		ArrayList<Pair<String, ArrayList<EQTL>>> topgtexfx = getGTExTopFx(gtextarball,
				"egenes.txt.gz",
				null,
				false,
				includeblood);
		
		
		TextFile tf = new TextFile(fullcis, TextFile.R, 10485760);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> significantGenes = new HashSet<String>();
		
		int nrsig = 0;
		int ctr = 0;
		while (elems != null) {
			String snp = (elems[1]);
			String gene = (elems[4]);
			String query = null;
			
			double fdr = Double.parseDouble(elems[20]);
			if (fdr < 0.05) {
				nrsig++;
				significantGenes.add(gene);
			}
			
			
			if (ctr % 100000 == 0) {
				System.out.print("\r" + ctr + "\tqtls processed. " + significantGenes.size() + " significant genes found so far");
			}
			ctr++;
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println();
		
		HashMap<Integer, HashSet<String>> genesWithoutEQTLInEQTLGEN = new HashMap<>();
		for (String gene : rankedgenes.keySet()) {
			Integer bin = rankedgenes.get(gene);
			if (!significantGenes.contains(gene)) {
				HashSet<String> genes = genesWithoutEQTLInEQTLGEN.get(bin);
				if (genes == null) {
					genes = new HashSet<>();
				}
				genes.add(gene);
				genesWithoutEQTLInEQTLGEN.put(bin, genes);
			}
		}
		
		
		HashMap<Integer, HashSet<String>> genesWithEQTLInTissue = new HashMap<>();
		
		int tctr = 0;
		
		ArrayList<HashMap<String, EQTL>> geneToEQTL = new ArrayList<>();
		
		for (Pair<String, ArrayList<EQTL>> pair : topgtexfx) {
			String tissue = pair.getLeft();
			
			HashMap<String, EQTL> emap = new HashMap<>();
			HashSet<String> bloodGenesUniqueToTissue = new HashSet<>();
			for (EQTL e : pair.getRight()) {
				String gene = e.getProbe();
				Integer bin = rankedgenes.get(gene);
				boolean significant = false;
				if (customthreshold != null) {
					if (e.getPvalue() < customthreshold) {
						significant = true;
					}
				} else if (usebonferroni) {
					double threshold = 0.05 / pair.getRight().size();
					if (e.getPvalue() < threshold) {
						significant = true;
					}
				} else {
					if (e.getFDR() < 0.05) {
						significant = true;
					}
				}
				
				if (significant && bin != null && !significantGenes.contains(gene)) {
					HashSet<String> genes = genesWithEQTLInTissue.get(bin);
					if (genes == null) {
						genes = new HashSet<>();
					}
					genes.add(gene);
					bloodGenesUniqueToTissue.add(gene);
					genesWithEQTLInTissue.put(bin, genes);
					emap.put(gene, e);
					
				}
			}
			geneToEQTL.add(emap);
			tctr++;
			
			
		}
		
		for (int bin = 0; bin < nrbins; bin++) {
			HashSet<String> genes = genesWithoutEQTLInEQTLGEN.get(bin);
			HashSet<String> gtexgenes = genesWithEQTLInTissue.get(bin);
			
			HashSet<String> shared = new HashSet<>();
			for (String gene : genes) {
				if (gtexgenes.contains(gene)) {
					shared.add(gene);
				}
			}
//			System.out.println(bin + "\t" + genes.size() + "\t" + shared.size());
//
//			if (bin >= 5) {
//				for (String gene : genes) {
//					if (!shared.contains(gene)) {
//						System.out.println(bin + "\t" + gene);
//					}
//				}
//			}
			
			if (bin >= 5) {
//				System.out.println();
				// print the eQTLs in sign. tissues
				for (String gene : genes) {
					ArrayList<String> tissues = new ArrayList<>();
					boolean print = false;
					for (int t = 0; t < topgtexfx.size(); t++) {
						HashMap<String, EQTL> emap = geneToEQTL.get(t);
						EQTL e = emap.get(gene);
						
						if (e != null) {
							tissues.add(topgtexfx.get(t).getLeft());
							print = true;
//							System.out.println(gene + "\t" + topgtexfx.get(t).getLeft() + "\t" + e.getRsName() + "\t" + e.getPvalue() + "\t" + e.getFDR());
						}
					}
					if (print) {
						System.out.println(bin + "\t" + gene + "\t" + Strings.concat(tissues, Strings.semicolon));
					}
				}
			}
		}
	}
	
	public void determineConcordanceBetweenTopGTExAndFullCis(String fullcis, String gtextarball, boolean maponposition, String outputfile, boolean includeblood, boolean usebonferroni) throws Exception {
		ArrayList<Pair<String, ArrayList<EQTL>>> topgtexfx = getGTExTopFx(gtextarball,
				"egenes.txt.gz",
				null,
				false,
				includeblood);
		
		// make a list of snp/gene combos to read from eqtlgen
		HashSet<String> toread = new HashSet<>();
		for (Pair<String, ArrayList<EQTL>> pair : topgtexfx) {
			for (EQTL e : pair.getRight()) {
				if (maponposition) {
					String query = (e.getRsChr() + "_" + e.getRsChrPos() + "_" + e.getProbe());
					toread.add(query);
				} else {
					String query = (e.getRsName() + "_" + e.getProbe());
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
		HashSet<String> significantGenes = new HashSet<String>();
		
		int nrsig = 0;
		while (elems != null) {
			String snp = (elems[1]);
			String gene = (elems[4]);
			String query = null;
			if (maponposition) {
				query = (elems[2] + "_" + elems[3] + "_" + elems[4]);
			} else {
				query = elems[1] + "_" + elems[4];
			}
			
			refgenes.add(gene);
			refsnps.add(snp);
			
			double fdr = Double.parseDouble(elems[20]);
			if (fdr < 0.05) {
				nrsig++;
				significantGenes.add(gene);
			}
			
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
					subqtl.setAlleleAssessed((elems[9]));
					subqtl.setAlleles((elems[8]));
					subqtl.setZscore(Double.parseDouble(elems[10]));
					subqtl.setFDR(fdr);
					eqtlmap.put((query), subqtl);
				}
			}
			if (ctr % 1000000 == 0) {
				System.out.print("\r" + ctr + "\tqtls processed. " + eqtlmap.size() + " combos found so far");
			}
			ctr++;
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println();
		System.out.println(eqtlmap.size() + " eqtls read from eqtlgen");
		
		HashSet<String> gtexgenes = new HashSet<>();
		HashSet<String> significantgtexgenes = new HashSet<>();
		HashSet<String> gtexsnps = new HashSet<>();
		
		for (Pair<String, ArrayList<EQTL>> pair : topgtexfx) {
			double bonferroni = 0.05 / pair.getRight().size();
			for (EQTL e : pair.getRight()) {
				
				if ((usebonferroni && e.getPvalue() < bonferroni) || (!usebonferroni && e.getFDR() < 0.05)) {
					significantgtexgenes.add(e.getProbe());
				}
				gtexgenes.add(e.getProbe());
				gtexsnps.add(e.getRsName());
			}
		}
		
		TextFile tfo3 = new TextFile(outputfile + "-00-gtex-genesnotfound.txt", TextFile.W);
		for (String gene : gtexgenes) {
			if (!refgenes.contains(gene)) {
				tfo3.writeln(gene);
			}
		}
		tfo3.close();
		
		TextFile tfo4 = new TextFile(outputfile + "-00-eqtlgen-snpsnotfound.txt", TextFile.W);
		for (String snp : gtexsnps) {
			if (!refsnps.contains(snp)) {
				tfo4.writeln(snp);
			}
		}
		tfo4.close();
		
		
		TextFile notfound = new TextFile(outputfile + "-00-gtex-eQTLs-notfound.txt", TextFile.W);
		notfound.writeln("snp\tsnpchr\tsnppos\tgene\tsnppresent\tgenepresent");
		
		// now check concordance
		TextFile out = new TextFile(outputfile, TextFile.W);
		String header = "tissue" +
				"\tnrInGTEx" +
				"\tnrSignificantGTEx" +
				"\tnrGTExSNPsInEQTLGen" +
				"\tnrGTExGenesInEQTLGen" +
				"\tnrGTExeQTLsPossibleInEQTLGen" +
				"\tShared" +
				"\tSharedSignificantInGTEx" +
				"\tSharedSignificantInGTExAlsoSignificantInEQTLGen" +
				"\tSharedSignificantInEQTLGen" +
				"\tSharedSignificantInEQTLGenAlsoSignificantInGTEx" +
				"\tSharedSignificantInBoth" +
				"\tConcordant" +
				"\tConcordantSignificantInGTEx" +
				"\tConcordantSignificantInGTExAlsoSignificantInEQTLGen" +
				"\tConcordantSignificantInEQTLGen" +
				"\tConcordantSignificantInEQTLGenAlsoSignificantInGTEx" +
				"\tConcordantSignificantInBoth";
		
		out.writeln(header);
		
		int tissuectr = 0;
		
		String[] datasetnames = new String[topgtexfx.size() + 1];
		datasetnames[datasetnames.length - 1] = "eQTLGen";
		
		Range range = new Range(-100, -50, 100, 50);
		
		Grid grid = new Grid(300, 300, 10, 5, 100, 100);
		Grid gridsiggtex = new Grid(300, 300, 10, 5, 100, 100);
		Grid gridsigeqtlgen = new Grid(300, 300, 10, 5, 100, 100);
		Grid gridsigboth = new Grid(300, 300, 10, 5, 100, 100);
		
		
		tissuectr = 0;
		for (Pair<String, ArrayList<EQTL>> pair : topgtexfx) {
			
			String tissue = pair.getLeft();
			String nicetissue = tissue.replaceAll("\\.v7\\.egenes\\.txt\\.gz", "");
			datasetnames[tissuectr] = nicetissue;
			int nringtex = pair.getRight().size();
			
			int nrsigingtex = 0;
			int shared = 0;
			int sharedsiggtex = 0;
			int sharedsiggtexalsosigeqtlgen = 0;
			int sharedsigeqtlgen = 0;
			int sharedsigeqtlgenalsosiggtex = 0;
			int sharedsigboth = 0;
			int concordant = 0;
			int concordantsiggtex = 0;
			int concordantsiggtexalsosigeqtlgen = 0;
			int concordantsigeqtlgen = 0;
			int concordantsigeqtlgenalsosiggtex = 0;
			int concordantsigboth = 0;
			int nrsharedgenes = 0;
			int nrsharedsnps = 0;
			
			int nrpossiblydetectedineqtlgen = 0;
			for (EQTL e : pair.getRight()) {
				if (refsnps.contains(e.getRsName()) && refgenes.contains(e.getProbe())) {
					nrpossiblydetectedineqtlgen++;
				}
				if (refsnps.contains(e.getRsName())) {
					nrsharedsnps++;
				}
				if (refgenes.contains(e.getProbe())) {
					nrsharedgenes++;
				}
			}
			
			double bonferroni = 0.05 / pair.getRight().size();
			
			TextFile zscoreout = new TextFile(outputfile + "-" + tissue + "-zcomp.txt.gz", TextFile.W);
			String head = "Gene\tRsId\tRefAlleles\tRefAssessed\tRefZ\tRefFDR\tTestAlleles\tTestAssessed\tTestZ\tTestFDR\tConcordant";
			zscoreout.writeln(head);
			
			ArrayList<Double> xval = new ArrayList<>();
			ArrayList<Double> yval = new ArrayList<>();
			ArrayList<Double> xvalsiggtex = new ArrayList<>();
			ArrayList<Double> yvalsiggtex = new ArrayList<>();
			ArrayList<Double> xvalsigeqtlgen = new ArrayList<>();
			ArrayList<Double> yvalsigeqtlgen = new ArrayList<>();
			ArrayList<Double> xvalsigeboth = new ArrayList<>();
			ArrayList<Double> yvalsigeboth = new ArrayList<>();
			
			for (EQTL e : pair.getRight()) {
				String query = null;
				if (maponposition) {
					query = (e.getRsChr() + "_" + e.getRsChrPos() + "_" + e.getProbe());
				} else {
					query = e.getRsName() + "_" + e.getProbe();
				}
				
				
				EQTL ref = eqtlmap.get(query);
				boolean significantingtex = false;
				if ((usebonferroni && e.getPvalue() < bonferroni) || (!usebonferroni && e.getFDR() < 0.05)) {
					significantingtex = true;
					nrsigingtex++;
				}
				if (ref == null) {
					String notfoundstr = tissue + "\t" + e.getRsName() + "\t" + e.getRsChr() + "\t" + e.getRsChrPos() + "\t" + e.getProbe() + "\t" + refsnps.contains(e.getRsName()) + "\t" + refgenes.contains(e.getProbe());
					notfound.writeln(notfoundstr);
				} else {
					
					
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
					
					if (significantingtex) {
						if (ref.getFDR() < 0.05) {
							sharedsiggtexalsosigeqtlgen++;
							if (concordance) {
								concordantsiggtexalsosigeqtlgen++;
							}
						}
					}
					
					xval.add(ref.getZscore());
					yval.add(e.getZscore());
					
					shared++;
					if (significantingtex) {
						sharedsiggtex++;
						xvalsiggtex.add(ref.getZscore());
						yvalsiggtex.add(e.getZscore());
					}
					
					// check concordance
					boolean significantinref = false;
					if (ref.getFDR() < 0.05) {
						if (significantingtex) {
							sharedsigeqtlgenalsosiggtex++;
							if (concordance) {
								concordantsigeqtlgenalsosiggtex++;
							}
						}
						significantinref = true;
						xvalsigeqtlgen.add(ref.getZscore());
						yvalsigeqtlgen.add(e.getZscore());
						
						sharedsigeqtlgen++;
						if (significantingtex) {
							sharedsigboth++;
							xvalsigeboth.add(ref.getZscore());
							yvalsigeboth.add(e.getZscore());
						}
					}
					
					if (concordance) {
//						System.out.println(outln);
						concordant++;
						if (significantingtex) {
							concordantsiggtex++;
							if (significantinref) {
								concordantsigboth++;
							}
						}
						if (significantinref) {
							concordantsigeqtlgen++;
						}
					}
					
					String outln = e.getProbe()
							+ "\t" + e.getRsName()
							+ "\t" + ref.getRsName()
							+ "\t" + ref.getAlleles()
							+ "\t" + ref.getAlleleAssessed()
							+ "\t" + ref.getZscore()
							+ "\t" + ref.getFDR()
							+ "\t" + e.getAlleles()
							+ "\t" + e.getAlleleAssessed()
							+ "\t" + z
							+ "\t" + e.getFDR()
							+ "\t" + concordance;
					zscoreout.writeln(outln);
				}
			}
			zscoreout.close();
			
			String ln = tissue
					+ "\t" + nringtex
					+ "\t" + nrsigingtex
					+ "\t" + nrsharedsnps
					+ "\t" + nrsharedgenes
					+ "\t" + nrpossiblydetectedineqtlgen
					+ "\t" + shared
					+ "\t" + sharedsiggtex
					+ "\t" + sharedsiggtexalsosigeqtlgen
					+ "\t" + sharedsigeqtlgen
					+ "\t" + sharedsigeqtlgenalsosiggtex
					+ "\t" + sharedsigboth
					+ "\t" + concordant
					+ "\t" + concordantsiggtex
					+ "\t" + concordantsiggtexalsosigeqtlgen
					+ "\t" + concordantsigeqtlgen
					+ "\t" + concordantsigeqtlgenalsosiggtex
					+ "\t" + concordantsigboth;
			out.writeln(ln);
			tissuectr++;
			
			double xmax = Primitives.max(Primitives.toPrimitiveArr(xval));
			double xmin = Primitives.min(Primitives.toPrimitiveArr(xval));
			System.out.println("Max X: " + xmax);
			System.out.println("Min X: " + xmin);
			double ymax = Primitives.max(Primitives.toPrimitiveArr(yval));
			double ymin = Primitives.min(Primitives.toPrimitiveArr(yval));
			System.out.println("Max Y: " + ymax);
			System.out.println("Min Y: " + ymin);
			System.out.println();
			xval = clip(xval, range.getMinX(), range.getMaxX());
			yval = clip(yval, range.getMinY(), range.getMaxY());
			xvalsigeboth = clip(xvalsigeboth, range.getMinX(), range.getMaxX());
			yvalsigeboth = clip(yvalsigeboth, range.getMinY(), range.getMaxY());
			xvalsiggtex = clip(xvalsiggtex, range.getMinX(), range.getMaxX());
			yvalsiggtex = clip(yvalsiggtex, range.getMinY(), range.getMaxY());
			xvalsigeqtlgen = clip(xvalsigeqtlgen, range.getMinX(), range.getMaxX());
			yvalsigeqtlgen = clip(yvalsigeqtlgen, range.getMinY(), range.getMaxY());
			
			
			ScatterplotPanel p = new ScatterplotPanel(1, 1);
			p.setData(Primitives.toPrimitiveArr(xval), Primitives.toPrimitiveArr(yval));
			
			p.setDataRange(range);
			p.setPlotElems(true, false);
			p.setAlpha(0.2f);
			DecimalFormat df = new DecimalFormat("##.##");
			String percconcordant = df.format(((double) concordant / shared) * 100);
			String psigtitle = "Shared: " + shared + "; concordant: " + concordant + " (" + percconcordant + "%)";
			p.setTitle(psigtitle);
			p.setLabels("eQTLGen", nicetissue);
			
			grid.addPanel(p);
			
			ScatterplotPanel psigboth = new ScatterplotPanel(1, 1);
			percconcordant = df.format(((double) concordantsigboth / sharedsigboth) * 100);
			psigtitle = "Shared: " + sharedsigboth + "; concordant: " + concordantsigboth + " (" + percconcordant + "%)";
			psigboth.setTitle(psigtitle);
			psigboth.setData(Primitives.toPrimitiveArr(xvalsigeboth), Primitives.toPrimitiveArr(yvalsigeboth));
			psigboth.setDataRange(range);
			psigboth.setPlotElems(true, false);
			psigboth.setAlpha(0.2f);
			psigboth.setLabels("eQTLGen", nicetissue);
			gridsigboth.addPanel(psigboth);
			
			ScatterplotPanel psiggtex = new ScatterplotPanel(1, 1);
			percconcordant = df.format(((double) concordantsiggtex / sharedsiggtex) * 100);
			psigtitle = "Shared: " + sharedsigboth + "; concordant: " + concordantsigboth + " (" + percconcordant + "%)";
			psiggtex.setTitle(psigtitle);
			psiggtex.setData(Primitives.toPrimitiveArr(xvalsiggtex), Primitives.toPrimitiveArr(yvalsiggtex));
			psiggtex.setDataRange(range);
			psiggtex.setPlotElems(true, false);
			psiggtex.setAlpha(0.2f);
			psiggtex.setLabels("eQTLGen", nicetissue);
			gridsiggtex.addPanel(psiggtex);
			
			ScatterplotPanel psigeqtlgen = new ScatterplotPanel(1, 1);
			percconcordant = df.format(((double) concordantsigeqtlgen / sharedsigeqtlgen) * 100);
			psigtitle = "Shared: " + sharedsigboth + "; concordant: " + concordantsigboth + " (" + percconcordant + "%)";
			psigeqtlgen.setTitle(psigtitle);
			psigeqtlgen.setData(Primitives.toPrimitiveArr(xvalsigeqtlgen), Primitives.toPrimitiveArr(yvalsigeqtlgen));
			psigeqtlgen.setDataRange(range);
			psigeqtlgen.setPlotElems(true, false);
			psigeqtlgen.setAlpha(0.2f);
			psigeqtlgen.setLabels("eQTLGen", nicetissue);
			gridsigeqtlgen.addPanel(psigeqtlgen);
		}
		out.close();
		notfound.close();
		
		
		grid.draw(outputfile + "-zscoreplot.pdf");
		gridsigboth.draw(outputfile + "-zscoreplot-sig-both.pdf");
		gridsigeqtlgen.draw(outputfile + "-zscoreplot-sig-eqtlgen.pdf");
		gridsiggtex.draw(outputfile + "-zscoreplot-sig-gtex.pdf");
		
		System.out.println(nrsig + " eqtls significant; " + significantGenes.size() + " genes significant our of " + refgenes.size());
		
		
	}
	
	private ArrayList<Double> clip(ArrayList<Double> xval, double min, double max) {
		ArrayList<Double> output = new ArrayList<>();
		for (Double z : xval) {
			if (z < min) {
				output.add(min);
			} else if (z > max) {
				output.add(max);
			} else {
				output.add(z);
			}
		}
		return output;
		
	}
	
	public void determineLDBetweenTopCisFx(String cisfxfile, String gtextarball, String ensembl, String tabixfile, String samplefilter, boolean includeblood) throws IOException {
		
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
				false,
				includeblood
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
						Pair<Double, Double> ldvals = ldcalc.getLD(var1, var2);
						rsq = ldvals.getRight();
					}
				}
			}
		}
		
		
	}
	
	private ArrayList<Pair<String, ArrayList<EQTL>>> getGTExTopFx(String gtextarball,
																  String filefilter,
																  Set<String> genelimit,
																  boolean includesexchr,
																  boolean includeblood) throws IOException {
		
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
		
		// sort files
		Arrays.sort(listOfFiles, NameFileComparator.NAME_COMPARATOR);
		
		for (File file : listOfFiles) {
			if (file.getName().endsWith(filefilter)) {
				TextFile br = new TextFile(file, TextFile.R); // new BufferedReader(new InputStreamReader(tarInput)); // Read directly from tarInput
				br.readLine();
//				System.out.println("For File = " + currentEntry.getName());
				String tissue = file.getName();
				if (includeblood || (!includeblood && !tissue.contains("Whole_Blood.v7.egenes"))) {
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
						String geneid = (elems[0].split("\\.")[0]);
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
										e.setRsName((elems[18]));
										e.setRsChr((byte) chr.getNumber());
										e.setRsChrPos(Integer.parseInt(elems[14]));
										e.setProbe(geneid);
										e.setFDR(Double.parseDouble(elems[28]));
										e.setProbeChr((byte) chr2.getNumber());
										e.setAlleles((elems[15] + "/" + elems[16]));
										e.setAlleleAssessed((elems[16]));
										e.setZscore(zscore);
										
										Double pval = Double.parseDouble(elems[23]);
										e.setPvalue(pval);
										
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
											e.setRsName((elems[18]));
											e.setRsChr((byte) chr.getNumber());
											e.setRsChrPos(Integer.parseInt(elems[14]));
											e.setProbe(geneid);
											e.setFDR(Double.parseDouble(elems[28]));
											e.setProbeChr((byte) chr2.getNumber());
											e.setAlleles((elems[15] + "/" + elems[16]));
											e.setAlleleAssessed((elems[16]));
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
		}
		return eqtlspertissue;
	}
	
	class GenePair implements Comparable<GenePair> {
		String gene;
		Double exp;
		
		public GenePair(String gene, Double exp) {
			this.gene = gene;
			this.exp = exp;
		}
		
		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			
			GenePair genePair = (GenePair) o;
			
			return exp.equals(genePair.exp);
		}
		
		@Override
		public int hashCode() {
			return exp.hashCode();
		}
		
		
		@Override
		public int compareTo(GenePair o) {
			if (this.equals(o)) {
				return 0;
			} else if (this.exp > o.exp) {
				return 1;
			} else {
				return -1;
			}
		}
	}
	
	public Pair<int[], HashMap<String, Integer>> loadAndRankGenes(String expfile, int nrbins, boolean onlyautosomal, String gtf) throws IOException {
		GTFAnnotation f = new GTFAnnotation(gtf);
		
		HashMap<String, Integer> geneToBin = new HashMap<String, Integer>();
		TextFile tf = new TextFile(expfile, TextFile.R);
		tf.readLine(); // header;
		String[] elems = tf.readLineElems(TextFile.tab);
		ArrayList<GenePair> pairs = new ArrayList<>();
		while (elems != null) {
			String gene = elems[0];
			
			String tested = elems[4];
			if (!tested.equals("not_tested")) {
				Double exp = Double.parseDouble(elems[1]);
				
				Gene geneObj = f.getStrToGene().get(gene);
				
				if (geneObj == null) {
					System.out.println("Could not find gene: " + gene);
				} else {
					if (onlyautosomal && geneObj.getChromosome().isAutosome()) {
						pairs.add(new GenePair(gene, exp));
					}
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println("File has " + pairs.size() + " genes, all autosome: " + onlyautosomal);
		
		Collections.sort(pairs);
		
		int nrgenes = pairs.size();
		int genesperbin = pairs.size() / nrbins;
		int ctr = 0;
		int bin = 0;
		int[] nrgenesperbin = new int[nrbins];
		while (ctr < nrgenes) {
			String gene = pairs.get(ctr).gene;
			geneToBin.put(gene, bin);
			nrgenesperbin[bin]++;
			ctr++;
			if (ctr % genesperbin == 0) {
				if (bin < nrbins - 1) {
					bin++;
				}
			}
		}
		for (int i = 0; i < nrbins; i++) {
			System.out.println(i + "\t" + nrgenesperbin[i]);
		}
		
		System.out.println("max bin: " + bin);
		return new Pair<>(nrgenesperbin, geneToBin);
		
	}
}
