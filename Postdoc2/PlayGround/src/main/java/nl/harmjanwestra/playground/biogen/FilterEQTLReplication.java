package nl.harmjanwestra.playground.biogen;

import eqtlmappingpipeline.util.QTLFileSorter;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

public class FilterEQTLReplication {
	
	public static void main(String[] args) {
		FilterEQTLReplication f = new FilterEQTLReplication();
		String annotation = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\annotation\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";
		String snpannot = "D:\\Sync\\SyncThing\\Data\\Ref\\dbsnp\\SNPMappings-dbsnp151-b38.txt.gz";
//        {
////            // Brainseq
//            String input = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\brainseq\\allEqtls-rewrite.txt.gz";
//            String rewriteGenesOutput = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\brainseq\\allEqtls-rewritegenes.txt.gz";
//            String rewriteSNPOutput = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\brainseq\\output\\allEqtls.txt.gz";
//            try {
//                f.rewriteGeneNames(input, annotation, rewriteGenesOutput);
//                f.rewriteSNPPositions(rewriteGenesOutput, snpannot, rewriteSNPOutput);
//                String allsignificantout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\brainseq\\output\\allEqtls-FDR0.05.txt.gz";
//                String topfxsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\brainseq\\output\\allEqtls-topfxpergene-FDR0.05.txt.gz";
//                String topfxnonsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\brainseq\\output\\allEqtls-topfxpergene.txt.gz";
//                f.filterFDR(rewriteSNPOutput, 0.05, allsignificantout, topfxsig, topfxnonsig);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }
//        System.exit(0);
//        {
//            // CMC no SVA
//            String input = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedNoSVA-rewrite.txt.gz";
//            String rewriteGenesOutput = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedNoSVA-rewritegenes.txt.gz";
//            try {
//                f.rewriteGeneNames(input, annotation, rewriteGenesOutput);
//                String rewriteSNPOutput = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\output\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedNoSVA.txt.gz";
//                f.rewriteSNPPositions(rewriteGenesOutput, snpannot, rewriteSNPOutput);
//                String allsignificantout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\output\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedNoSVA-FDR0.05.txt.gz";
//                String topfxsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\output\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedNoSVA-topfxpergene-FDR0.05.txt.gz";
//                String topfxnonsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\output\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedNoSVA-topfxpergene.txt.gz";
//                f.filterFDR(rewriteSNPOutput, 0.05, allsignificantout, topfxsig, topfxnonsig);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }
//        {
//            // CMC SVA
//            String input = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA-rewrite.txt.gz";
//            String rewriteGenesOutput = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA-rewritegenes.txt.gz";
//
//
//            try {
//                f.rewriteGeneNames(input, annotation, rewriteGenesOutput);
//                String rewriteSNPOutput = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\output\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt.gz";
//                f.rewriteSNPPositions(rewriteGenesOutput, snpannot, rewriteSNPOutput);
//                String allsignificantout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\output\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA-FDR0.05.txt.gz";
//                String topfxsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\output\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA-topfxpergene-FDR0.05.txt.gz";
//                String topfxnonsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\cmc\\output\\CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA-topfxpergene.txt.gz";
//                f.filterFDR(rewriteSNPOutput, 0.05, allsignificantout, topfxsig, topfxnonsig);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }
		{
			// Mostavafi
//            try {
//
//                String input = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\mostavafi\\eQTLs_all-rewrite.txt.gz";
//
//                String rewriteGenesOutput = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\mostavafi\\eQTLs_all-rewritegenes.txt.gz";
//                f.rewriteGeneSymbols(input, annotation, rewriteGenesOutput);
//
//
//                String rewritesnpidout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\mostavafi\\eQTLs_all-rewritegenes-rewritesnpids.txt.gz";
//                String b37annot = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz";
//                 f.rewriteSNPIds(rewriteGenesOutput, b37annot, rewritesnpidout);
//
//                String rewritesnpout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\mostavafi\\eQTLs_all-rewritegenes-rewritesnppos-rewritesnpids.txt.gz";
//                 f.rewriteSNPPositions(rewritesnpidout, snpannot, rewritesnpout);
//
//                // calculate actual Z-scores?
//                // perform BH FDR
//                String outpustatsfixed = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\mostavafi\\output\\eQTLs_alltxt.gz";
//                f.recalculateSummaryStats(rewritesnpout, 494, outpustatsfixed);

//                String allsignificantout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\mostavafi\\output\\eQTLs_all-FDR0.05.txt.gz";
//                String topfxsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\mostavafi\\output\\eQTLs_all-topfxpergene-FDR0.05.txt.gz";
//                String topfxnonsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\mostavafi\\output\\eQTLs_all-topfxpergene.txt.gz";
//                f.filterFDR(outpustatsfixed, 0.05, allsignificantout, topfxsig, topfxnonsig);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
		
		}
		{
			// eqtlgen trans
//            String input = "D:\\trans\\2018-04-03-trans-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05-CohortInfoRemoved.txt.gz";
//            String rewritesnpoutput = "D:\\trans\\2018-04-03-trans-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05-CohortInfoRemoved-b38-snpupdate.txt.gz";
//            try {
//                f.rewriteSNPPositions(input, snpannot, rewritesnpoutput);
//                String rewritegeneoutput = "D:\\trans\\2018-04-03-trans-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05-CohortInfoRemoved-b38-snpupdate-geneupdate.txt.gz";
//                f.rewriteGeneNames(rewritesnpoutput, annotation, rewritegeneoutput);
//
////                String allsignificantout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\eqtlgen\\output\\2018-04-03-trans-eQTLsFDR.txt.gz";
////                String topfxsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\eqtlgen\\output\\2018-04-03-trans-eQTLsFDR-topfxpergene-FDR0.05.txt.gz";
////                String topfxnonsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\eqtlgen\\output\\2018-04-03-trans-eQTLsFDR-topfxpergene.txt.gz";
////                f.filterFDR(rewritegeneoutput, 0.05, allsignificantout, topfxsig, topfxnonsig);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
			
			// eqtlgen trans
//			String input = "D:\\trans\\2018-09-04-trans-eQTLsFDR-CohortInfoRemoved.txt.gz";
//			String rewritesnpoutput = "D:\\trans\\2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-b38-snpupdate.txt.gz";
//			try {
//				f.rewriteSNPPositions(input, snpannot, rewritesnpoutput);
//				String rewritegeneoutput = "D:\\trans\\2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-snpupdate-geneupdate.txt.gz";
//				f.rewriteGeneNames(rewritesnpoutput, annotation, rewritegeneoutput);
//
////                String allsignificantout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\eqtlgen\\output\\2018-04-03-trans-eQTLsFDR.txt.gz";
////                String topfxsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\eqtlgen\\output\\2018-04-03-trans-eQTLsFDR-topfxpergene-FDR0.05.txt.gz";
////                String topfxnonsig = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-11-07-replicationdata\\eqtlgen\\output\\2018-04-03-trans-eQTLsFDR-topfxpergene.txt.gz";
////                f.filterFDR(rewritegeneoutput, 0.05, allsignificantout, topfxsig, topfxnonsig);
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//
			
			// eqtlgen high confidence
//			String in = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2019-01-transrepl\\highconfidencecombos.txt";
//			String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2019-01-transrepl\\highconfidencecombos-b38.txt";
//			try {
//				f.rewriteCombos(in, annotation, out);
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//
//        }

//			String in = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.txt.gz";
//			String out = "D:\\biogen\\eqtlgenfx\\cis\\2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-b38-geneupdate.txt.gz";
//			try {
//				f.rewriteGeneNames(in, annotation, out);
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
			
			String in = "D:\\biogen\\metabrainfx\\cis\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
			String out = "D:\\biogen\\metabrainfx\\cis\\eQTLProbesFDR0.05-ProbeLevel-b37snppos.txt.gz";
			String annot = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz";
			try {
				f.rewriteSNPPositions(in, annot, out);
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
	}
	
	private void rewriteCombos(String in, String annotation, String out) throws IOException {
		GTFAnnotation annot = new GTFAnnotation(annotation);
		Collection<Gene> genes = annot.getGenes();
		HashMap<String, Gene> oldIdToNewID = new HashMap<String, Gene>();
		for (Gene gene : genes) {
			String geneName = gene.getName();
			String oldid = geneName.split("\\.")[0];
			oldIdToNewID.put(oldid, gene);
		}
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String old = elems[1];
			Gene newid = oldIdToNewID.get(old);
			if (newid != null) {
				tfo.writeln(elems[0] + "\t" + newid.getName());
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		
		
		tf.close();
		tfo.close();
		
	}
	
	public void recalculateSummaryStats(String in, int n, String out) throws IOException {
		Correlation.correlationToZScore(n);
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile outf = new TextFile(out + "-statsfixed.txt.gz", TextFile.W);
		String[] elems = tf.readLineElems(TextFile.tab);
		outf.writeln(Strings.concat(elems, Strings.tab, 0, elems.length - 2));
		
		
		elems = tf.readLineElems(TextFile.tab);
		System.out.println("Rewriting summary stats.");
		int ctr = 0;
		while (elems != null) {
			// 9.491118674287809531e-01        rs58552927      10      22752857        -       10      23729830        -       G/A     G       -       ROSMAP  -       494     -       -       OTUD1   -2.878779700448074255e-03       -2.878779700448074255e-03       -2.878779700448074255e-03
			double r = Double.parseDouble(elems[17]);
			if (!Double.isNaN(r) && Double.isFinite(r)) {
				double z = Correlation.convertCorrelationToZScore(n, r);
				elems[0] = "" + ZScores.zToP(z);
				elems[10] = "" + z;
				outf.writeln(Strings.concat(elems, Strings.tab, 0, elems.length - 2));
			}
			ctr++;
			if (ctr % 10000000 == 0) {
				System.out.println(ctr + " lines processed.");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		outf.close();
		
		System.out.println("Sorting eqTLs");
		QTLFileSorter sorter = new QTLFileSorter();
		sorter.run(out + "-statsfixed.txt.gz", out + "-sorted.txt.gz");
		
		Gpio.moveFile(out + "-sorted.txt.gz", out);
		
		System.out.println("Doing BHFDR");
		runBHFDR(out, 0.05, out + "sorted-FDR.txt.gz");
		
	}
	
	public void runBHFDR(String in, double significance, String out) throws IOException {
		System.out.println("Benjamini Hochberg correction");
		System.out.println("in: " + in);
		System.out.println("out: " + out);
		System.out.println("Threshold: " + significance);
		umcg.genetica.io.text.TextFile tf = new umcg.genetica.io.text.TextFile(in, umcg.genetica.io.text.TextFile.R);
		int nrlines = tf.countLines();
		tf.close();
		
		int nrtests = nrlines - 1;
		QTLTextFile tf2 = new QTLTextFile(in, umcg.genetica.io.text.TextFile.R);
		umcg.genetica.io.text.TextFile outf = new umcg.genetica.io.text.TextFile(out, umcg.genetica.io.text.TextFile.W);
		umcg.genetica.io.text.TextFile outfsig = new umcg.genetica.io.text.TextFile(out + "-FDR" + significance + ".txt.gz", umcg.genetica.io.text.TextFile.W);
		String header = tf2.readLine();
		outf.writeln(header);
		outfsig.writeln(header);
		Iterator<EQTL> it = tf2.getEQtlIterator();
		
		int ctr = 1;
		double prevp = -1;
		int nrsig = 0;
		
		while (it.hasNext()) {
			EQTL e = it.next();
			if (prevp == -1) {
				prevp = e.getPvalue();
			} else {
				if (e.getPvalue() < prevp) {
					System.out.println("QTL File not properly sorted! Faulty line: " + ctr + "\tSNP: " + e.getRsName() + "\tGene: " + e.getProbe());
					System.exit(-1);
				}
				prevp = e.getPvalue();
			}
			double qvalue = e.getPvalue() * (nrtests / ctr);
//			double qvalue = ((ctr + 1d) / nrlines) * e.getPvalue();
//			System.out.println(ctr + "\t" + qvalue);
//			if (ctr > 10) {
//				System.exit(-1);
//			}
			e.setFDR(qvalue);
			outf.writeln(e.toString());
			if (qvalue < significance) {
				outfsig.writeln(e.toString());
				nrsig++;
			}
			ctr++;
		}
		tf2.close();
		outf.close();
		outfsig.close();
		System.out.println(nrsig + " out of " + nrtests + " are signifciant.");
	}
	
	public void rewriteGeneSymbols(String in, String annot, String out) throws IOException {
		GTFAnnotation gtf = new GTFAnnotation(annot);
		
		HashMap<String, Gene> symboltoensg = new HashMap<>();
		for (Gene g : gtf.getGenes()) {
			String symbol = g.getGeneSymbol();
			symboltoensg.put(symbol, g);
		}
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		tfo.writeln(tf.readLine());
		String[] elems = tf.readLineElems(TextFile.tab);
		int lctr = 0;
		int written = 0;
		while (elems != null) {
			String symbol = elems[16];
			Gene newId = symboltoensg.get(symbol);
			if (newId != null) {
				int midpos = newId.getStart() + ((newId.getStop() - newId.getStart()) / 2);
				elems[4] = newId.getName();
				elems[5] = "" + newId.getChromosome().getNumber();
				elems[6] = "" + midpos;
				tfo.writeln(Strings.concat(elems, Strings.tab));
				written++;
			}
			elems = tf.readLineElems(TextFile.tab);
			lctr++;
			if (lctr % 1000000 == 0) {
				System.out.println(lctr + " lines, " + written + " written");
			}
		}
		tfo.close();
		tf.close();
		
	}
	
	public void rewriteSNPIds(String in, String annot, String out) throws IOException {
		
		System.out.println("About to rewrite snp ids");
		TextFile tf = new TextFile(in, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> toupdate = new HashSet<String>();
		int ctr = 0;
		while (elems != null) {
			
			String snp = elems[1];
			if (!snp.startsWith("rs")) {
				toupdate.add(snp);
			}
			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(ctr + " snps read. " + toupdate.size() + " snps to update");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(toupdate.size() + " snps to update");
		
		System.out.println("Reading snp map: " + annot);
		HashMap<String, String> map = new HashMap<String, String>();
		TextFile tf2 = new TextFile(annot, TextFile.R);
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String q = "chr" + elems[0] + ":" + elems[1];
			if (toupdate.contains(q)) {
				map.put(q, elems[2]);
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		System.out.println("Done reading map. " + map.size() + " snps will be re-IDed out of " + toupdate.size());
		
		tf.open();
		TextFile tfo = new TextFile(out, TextFile.W);
		tfo.writeln(tf.readLine());
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String snp = elems[1];
			if (!snp.startsWith("rs") && map.containsKey(snp)) {
				elems[1] = map.get(snp);
			}
			tfo.writeln(Strings.concat(elems, Strings.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		tfo.close();
		tf.close();
		System.out.println("Done re-IDing snps");
	}
	
	public void filtertrans(String in, String cisout, String transout2) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		String header = tf.readLine();
		TextFile tfc = new TextFile(cisout, TextFile.W);
		TextFile tft = new TextFile(transout2, TextFile.W);
		tfc.writeln(header);
		tft.writeln(header);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String chr1 = elems[2];
			String chr2 = elems[5];
			
			if (chr1.equals(chr2)) {
				// test distance?
				Integer pos1 = Integer.parseInt(elems[3]);
				Integer pos2 = Integer.parseInt(elems[6]);
				int dist = Math.abs(pos1 - pos2);
				if (dist <= 1000000) {
					tfc.writeln(Strings.concat(elems, Strings.tab));
				} else if (dist >= 5000000) {
					tft.writeln(Strings.concat(elems, Strings.tab));
				}
			} else {
				tft.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tfc.close();
		tf.close();
		tft.close();
		
		
	}
	
	public void rewriteSNPPositions(String in, String annot, String out) throws IOException {
		HashMap<String, Integer> snppos = new HashMap<String, Integer>();
		HashMap<String, Integer> snpchr = new HashMap<String, Integer>();
		
		TextFile tf = new TextFile(in, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			snppos.put(elems[1], null);
			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(ctr + " snppositions read");
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println("Done reading eqtl ");
		System.out.println(snppos.size() + " snps requested. ");
		int loaded = 0;
		int ctr2 = 0;
		TextFile tf2 = new TextFile(annot, TextFile.R);
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			
			String snp = elems[2];
			if (snppos.containsKey(snp)) {
				loaded++;
				
				Integer chr = Chromosome.parseChr(elems[0]).getNumber();
				Integer pos = Integer.parseInt(elems[1]);
				snppos.put(snp, pos);
				snpchr.put(snp, chr);
			}
			ctr2++;
			if (ctr2 % 1000000 == 0) {
				System.out.println(ctr2 + "\t" + loaded + " snps loaded / " + snppos.size());
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		System.out.println("Done reading SNPs");
		System.out.println(loaded + " out of " + snppos.size() + " loaded.");
		
		tf.open();
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln(tf.readLine());
		elems = tf.readLineElems(TextFile.tab);
		ctr = 0;
		
		while (elems != null) {
			
			String snp = elems[1];
			if (snppos.get(snp) != null) {
				elems[2] = "" + snpchr.get(snp);
				elems[3] = "" + snppos.get(snp);
				outf.writeln(Strings.concat(elems, Strings.tab));
				
			}
			
			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(ctr);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		outf.close();
		System.out.println("Done remapping snps");
	}
	
	public void filterFDR(String in, double threshold, String out, String outpairssig, String outpairsnonsig) throws
			IOException {
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		String header = tf.readLine();
		tfo.writeln(header);
		
		HashMap<String, Pair<Double, String>> topfxpergene = new HashMap<String, Pair<Double, String>>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			String snp = elems[1];
			if (snp.startsWith("rs")) {
				double fdr = Double.parseDouble(elems[elems.length - 1]);
				String gene = elems[4];
				double absZ = Math.abs(Double.parseDouble(elems[10]));
				if (elems[8].length() < 4) {
					if (!topfxpergene.containsKey(gene)) {
						topfxpergene.put(gene, new Pair<Double, String>(absZ, Strings.concat(elems, Strings.tab)));
					} else {
						Pair<Double, String> p = topfxpergene.get(gene);
						if (p.getLeft() < absZ) {
							topfxpergene.put(gene, new Pair<Double, String>(absZ, Strings.concat(elems, Strings.tab)));
						}
					}
					if (fdr < threshold) {
						tfo.writeln(Strings.concat(elems, Strings.tab));
					}
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		
		tf.close();
		tfo.close();
		
		
		// write the top fx
		TextFile tfo2 = new TextFile(outpairssig, TextFile.W);
		tfo2.writeln(header);
		TextFile tfo3 = new TextFile(outpairsnonsig, TextFile.W);
		tfo3.writeln(header);
		for (Pair<Double, String> p : topfxpergene.values()) {
			String e = p.getRight();
			elems = e.split("\t");
			double fdr = Double.parseDouble(elems[elems.length - 1]);
			if (fdr < threshold) {
				tfo2.writeln(e);
			}
			tfo3.writeln(e);
		}
		tfo2.close();
		tfo3.close();
	}
	
	public void rewriteGeneNames(String in, String geneannot, String out) throws IOException {
		
		GTFAnnotation annot = new GTFAnnotation(geneannot);
		Collection<Gene> genes = annot.getGenes();
		HashMap<String, Gene> oldIdToNewID = new HashMap<String, Gene>();
		for (Gene gene : genes) {
			String geneName = gene.getName();
			String oldid = geneName.split("\\.")[0];
			oldIdToNewID.put(oldid, gene);
		}
		
		TextFile tf = new TextFile(in, TextFile.R);
		TextFile tfo = new TextFile(out, TextFile.W);
		String header = tf.readLine();
		tfo.writeln(header);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		int written = 0;
		while (elems != null) {
			
			String gene = elems[4];
			Gene newId = oldIdToNewID.get(gene);
			if (newId != null) {
				int midpos = newId.getStart() + ((newId.getStop() - newId.getStart()) / 2);
				elems[4] = newId.getName();
				elems[5] = "" + newId.getChromosome().getNumber();
				elems[6] = "" + midpos;
				tfo.writeln(Strings.concat(elems, Strings.tab));
				written++;
			}
			ctr++;
			
			if (ctr % 10000 == 0) {
				System.out.print("\r" + ctr + " lines parsed, " + written + " written.");
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfo.close();
		System.out.println("Done");
		
	}
}
