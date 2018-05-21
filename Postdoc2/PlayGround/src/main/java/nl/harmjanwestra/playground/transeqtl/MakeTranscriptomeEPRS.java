package nl.harmjanwestra.playground.transeqtl;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Feature;
import nl.harmjanwestra.utilities.features.FeatureComparator;
import nl.harmjanwestra.utilities.features.SNPFeature;
import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import org.apache.tools.ant.Executor;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;

public class MakeTranscriptomeEPRS extends MakeTranscriptome {
	
	
	public static void main(String[] args) {
		
		String signprs = "D:\\Work\\eprsqc\\eQTLsFDR-Significant-0.05-snpprobecombos.txt.gz";
		String dbloc = "D:\\Work\\eprsqc\\db\\";
		String snplist = "D:\\Work\\eprsqc\\allprssnps.txt";
		String snpmap = "D:\\Work\\eprsqc\\GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz";
		String regiondeffolder = "D:\\Work\\eprsqc\\regiondef\\";
		String snpseqoutputdir = "D:\\Work\\eprsqc\\seq\\snps\\";
		String alndir = "D:\\Work\\eprsqc\\aln\\";
		
		int windowsize = 10000000;
		int clusterdistance = 1000000;
		
		String genomeFasta = "D:\\Sync\\SyncThing\\Data\\HumanGenome\\human_g1k_v37.fasta.gz";
		String geneoutputdir = "D:\\Work\\eprsqc\\seq\\genes\\";
		String prsfastaoutdir = "D:\\Work\\eprsqc\\seq\\prs\\";
		String ensemblannotation = "D:\\Sync\\SyncThing\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
		boolean exportindividualgenes = true;
		boolean exportfullexonsequences = false;
		int readlen = 35;
		int shift = 2;
		
		MakeTranscriptomeEPRS m = new MakeTranscriptomeEPRS();
		try {
			
			Gpio.createDir(regiondeffolder);
			Gpio.createDir(snpseqoutputdir);
			Gpio.createDir(geneoutputdir);

//			m.createGeneSequences(signprs, genomeFasta, geneoutputdir,
//					ensemblannotation, exportindividualgenes, exportfullexonsequences,
//					readlen, shift);
			int threads = 16;
//			m.createSNPList(signprs, dbloc, snplist);
			boolean onlyexportsnps = true;
			boolean overwriteexistingsnps = false;
//			m.createPRSRegions(signprs, dbloc, snpmap, snplist, clusterdistance, windowsize,
//					genomeFasta, snpseqoutputdir, regiondeffolder, threads, onlyexportsnps, overwriteexistingsnps);
			
			m.createAlignmentScripts(signprs, snpseqoutputdir,
					geneoutputdir, prsfastaoutdir,
					regiondeffolder, alndir);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void createAlignmentScripts(String significantEPRSFile, String snpfastadir,
										String genefastadir, String prsfastaoutdir,
										String regionDefFolder, String scriptOutputFolder) throws IOException {
		Gpio.createDir(scriptOutputFolder);
		Gpio.createDir(prsfastaoutdir);
		
		TextFile tf = new TextFile(significantEPRSFile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashMap<String, ArrayList<String>> prsgenes = new HashMap<>();
		while (elems != null) {
			if (elems.length >= 2) {
				String prs = elems[0];
				String gene = elems[1];
				ArrayList<String> genes = prsgenes.get(prs);
				if (genes == null) {
					genes = new ArrayList<>();
				}
				genes.add(gene);
				prsgenes.put(prs, genes);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		// translate prs into snps
		TextFile aln = new TextFile(scriptOutputFolder + "align.sh", TextFile.W);
		TextFile sam = new TextFile(scriptOutputFolder + "samse.sh", TextFile.W);
		aln.writeln("mkdir -p ./alignments/");
		int prsctr = 0;
		ExecutorService s = Executors.newWorkStealingPool(16);
		ArrayList<Future> outputs = new ArrayList<>();
		for (String prs : prsgenes.keySet()) {
			Runnable r = () -> {
				try {
					String[] prselems = prs.split("_");
					String pvalstr = prselems[prselems.length - 1];
					pvalstr = pvalstr.replaceAll("P", "");
					Double pvalthresh = Double.parseDouble(pvalstr);
					
					
					int lastunderscore = prs.lastIndexOf("_");
					String traitname = prs.substring(0, lastunderscore);
					// outputfolder + traitname + "_P" + pvalstr + ".txt.gz"
					String regiondef = regionDefFolder + traitname + "_P" + pvalstr + ".txt.gz";
					ArrayList<String> keysnps = new ArrayList<>();
					TextFile tf2 = new TextFile(regiondef, TextFile.R);
					String[] elems2 = tf2.readLineElems(TextFile.tab);
					while (elems2 != null) {
						if (elems2.length > 3) {
							if (Boolean.parseBoolean(elems2[3])) {
								String snpstr = elems2[1].split("_")[2];
								keysnps.add(snpstr);
							}
						}
						elems2 = tf2.readLineElems(TextFile.tab);
					}
					tf2.close();
					
					// write snps into single file
					// regiondef = regionDefFolder + traitname + "_P" + pvalstr + ".txt.gz";
					String fastaoutfile = prsfastaoutdir + traitname + "_P" + pvalstr + ".fa.gz";
					TextFile fastaout = new TextFile(fastaoutfile, TextFile.W);
					System.out.println("Processing: " + fastaoutfile + " for " + keysnps.size() + " snps");
					int ct = 0;
					for (String keysnp : keysnps) {
						String fastainfile = snpfastadir + keysnp + ".fa.gz";
						TextFile fa = new TextFile(fastainfile, TextFile.R);
						String ln = fa.readLine();
						while (ln != null) {
							fastaout.writeln(ln);
							ln = fa.readLine();
						}
						System.out.println(traitname + "\t" + fastainfile + "\t" + ct + "/" + keysnps.size());
						ct++;
						fa.close();
					}
					
					fastaout.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			};
			Future<?> output = s.submit(r);
			outputs.add(output);
		}
		
		boolean alldone = false;
		while (!alldone) {
			alldone = true;
			int nrdone = 0;
			
			for (Future f : outputs) {
				if (!f.isDone()) {
					alldone = false;
				} else {
					nrdone++;
				}
			}
			
			try {
				Thread.sleep(10000);
				System.out.println();
				System.out.println("-------------------------------------");
				System.out.println(nrdone + "/" + outputs.size() + " done");
				System.out.println("-------------------------------------");
				System.out.println();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		s.shutdown();
		
		for (String prs : prsgenes.keySet()) {
			ArrayList<String> genelist = prsgenes.get(prs);
			// write alignment script
			
			String[] prselems = prs.split("_");
			String pvalstr = prselems[prselems.length - 1];
			pvalstr = pvalstr.replaceAll("P", "");
			Double pvalthresh = Double.parseDouble(pvalstr);
			
			
			int lastunderscore = prs.lastIndexOf("_");
			String traitname = prs.substring(0, lastunderscore);
			
			String fastaoutfile = prsfastaoutdir + traitname + pvalstr + ".fa.gz";
			for (String gene : genelist) {
				
				String genefasta = genefastadir + gene + ".fa.gz";
				
				// bwa aln -t 1 -l 10 -0 ./snps/rs1131017.fa.gz ./genes/ENSG00000196656.fa.gz > ./alignments/rs1131017_ENSG00000196656.sai
				String saiout = "./alignments/" + traitname + "_" + gene + ".sai";
				String bamout = "./alignments/" + traitname + "_" + gene + ".bam";
				String script = "bwa aln -t 1 -l 10 -0 " + fastaoutfile + " " + genefasta + " > " + saiout;
				aln.writeln(script);
				// bwa samse ./snps/rs1131017.fa.gz ./alignments/rs1131017_ENSG00000196656.sai ./genes/ENSG00000196656.fa.gz > ./alignments/rs1131017_ENSG00000196656.sam
				String script2 = "bwa samse " + fastaoutfile + " " + saiout + " " + genefasta + " | samtools sort -O BAM -o " + bamout + " -";
				sam.writeln(script2);
			}
		}
		sam.close();
		aln.close();
	}
	
	private void createSNPList(String significantEPRSFile, String dbloc, String snplistoutputfile) throws IOException {
		ArrayList<String> uniqueprs = new ArrayList<>();
		{
			HashSet<String> significantPRS = new HashSet<String>();
			HashSet<String> significantPRSWoThresh = new HashSet<String>();
			
			TextFile tf = new TextFile(significantEPRSFile, TextFile.R);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 2) {
					String prs = elems[0];
					significantPRS.add(prs);
					
					int lastunderscore = prs.lastIndexOf("_");
					String traitname = prs.substring(0, lastunderscore);
					significantPRSWoThresh.add(traitname);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			System.out.println(significantPRS.size() + " significant PRS");
			uniqueprs.addAll(significantPRS);
		}
		
		HashSet<String> snps = new HashSet<>();
		ArrayList<Callable<ArrayList<String>>> tasks = new ArrayList<>();
		
		for (int i = 0; i < uniqueprs.size(); i++) {
			String prs = uniqueprs.get(i);
			String[] prselems = prs.split("_");
			String pvalstr = prselems[prselems.length - 1];
			pvalstr = pvalstr.replaceAll("P", "");
			Double pvalthresh = Double.parseDouble(pvalstr);
			
			
			int lastunderscore = prs.lastIndexOf("_");
			String traitname = prs.substring(0, lastunderscore);
			
			
			Callable<ArrayList<String>> c = () -> {
				String traitloc = dbloc + traitname;
				TextFile tf = new TextFile(traitloc, TextFile.R);
				tf.readLine();
				System.out.println("Parsing " + traitloc + " with threshold: " + pvalthresh);
				String[] elems = tf.readLineElems(TextFile.tab);
				ArrayList<String> snps1 = new ArrayList<>();
				while (elems != null) {
					if (elems.length >= 4) {
						String snp = new String(elems[0].getBytes(), "UTF-8");
						double p = Double.parseDouble(elems[3]);
						if (p < pvalthresh) {
							// includes
							
							snps1.add(snp);
							
						}
					}
					elems = tf.readLineElems(TextFile.tab);
				}
				
				System.out.println(snps1.size() + " significant snps.");
				tf.close();
				return snps1;
			};
			tasks.add(c);
		}
		{
			ExecutorService ex = Executors.newFixedThreadPool(16);
			ArrayList<Future<ArrayList<String>>> output = new ArrayList<>();
			for (Callable<ArrayList<String>> c : tasks) {
				output.add(ex.submit(c));
			}
			
			
			for (Future<ArrayList<String>> c : output) {
				try {
					while (!c.isDone()) {
						Thread.sleep(1000);
					}
					snps.addAll(c.get());
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
			}
			ex.shutdown();
		}
		
		System.out.println(snps.size() + " snps in total");
		TextFile out = new TextFile(snplistoutputfile, TextFile.W);
		for (String s : snps) {
			out.writeln(s);
		}
		out.close();
	}
	
	
	private void createGeneSequences(String significantEPRSFile, String genomeFasta, String geneoutputdir,
									 String ensemblannotation, boolean exportindividualgenes, boolean exportfullexonsequences,
									 int readlen, int shift) throws IOException {
		
		TextFile tf = new TextFile(significantEPRSFile, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		HashSet<String> significantPRSGenes = new HashSet<String>();
		while (elems != null) {
			if (elems.length >= 2) {
				String prs = elems[1];
				significantPRSGenes.add(prs);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		
		super.exportGenes(genomeFasta, geneoutputdir, significantPRSGenes, ensemblannotation,
				exportindividualgenes, exportfullexonsequences, readlen, shift);
		
		
	}
	
	
	public void createPRSRegions(String significantEPRSFile, String dbloc, String snpmap, String snplist, int clusterdistance, int windowsize,
								 String genomeFasta, String snpseqoutputdir,
								 String regionDefinitionOutputFolder, int threads, boolean exportonly, boolean overwriteexistingsnps) throws IOException {
		
		ArrayList<Triple<String, Chromosome, Integer>> snps = null;
		if (!exportonly) {
			ArrayList<String> uniqueprs = new ArrayList<>();
			
			HashSet<String> significantPRS = new HashSet<String>();
			HashSet<String> significantPRSWoThresh = new HashSet<String>();
			
			TextFile tf = new TextFile(significantEPRSFile, TextFile.R);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 2) {
					String prs = elems[0];
					significantPRS.add(prs);
					
					int lastunderscore = prs.lastIndexOf("_");
					String traitname = prs.substring(0, lastunderscore);
					significantPRSWoThresh.add(traitname);
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			System.out.println(significantPRS.size() + " significant PRS");
			uniqueprs.addAll(significantPRS);
			
			
			loadSNPMap(snpmap, snplist);
			ExecutorService ex = Executors.newFixedThreadPool(threads);
			ExecutorCompletionService<ArrayList<Feature>> ecs = new ExecutorCompletionService<>(ex);
			
			for (int i = 0; i < uniqueprs.size(); i++) {
				String prs = uniqueprs.get(i);
				String[] prselems = prs.split("_");
				String pvalstr = prselems[prselems.length - 1];
				pvalstr = pvalstr.replaceAll("P", "");
				Double pvalthresh = Double.parseDouble(pvalstr);
				
				
				int lastunderscore = prs.lastIndexOf("_");
				String traitname = prs.substring(0, lastunderscore);
				
				RegionCollectionTask task = new RegionCollectionTask(dbloc, traitname, pvalthresh, pvalstr,
						clusterdistance, windowsize, regionDefinitionOutputFolder, genomeFasta, snpseqoutputdir, i);
				ecs.submit(task);
				
			}
			int returned = 0;
			HashSet<Feature> tagsnps = new HashSet<>();
			while (returned < uniqueprs.size()) {
				try {
					ArrayList<Feature> taskout = ecs.take().get();
					tagsnps.addAll(taskout);
					returned++;
					System.out.println();
					System.out.println(returned + " / " + uniqueprs.size() + "\t" + tagsnps.size() + " tags so far");
					System.out.println();
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
				
			}
			
			ex.shutdown();
			
			System.out.println(tagsnps.size() + " total tags");
			System.out.println("Exporting to disk.");
			
			// create regions
			
			snps = new ArrayList<>();
			TextFile tfo = new TextFile(regionDefinitionOutputFolder + "tags.txt", TextFile.W);
			for (Feature tag : tagsnps) {
				Triple<String, Chromosome, Integer> tags = new Triple<>(tag.getName(), tag.getChromosome(), tag.getStart());
				tfo.writeln(tags.getMiddle().toString() + "\t" + tags.getRight() + "\t" + tags.getLeft());
				snps.add(tags);
			}
			tfo.close();
			Gpio.createDir(snpseqoutputdir);
		} else {
			snps = new ArrayList<>();
			TextFile tfo = new TextFile(regionDefinitionOutputFolder + "tags.txt", TextFile.R);
			String[] elems = tfo.readLineElems(TextFile.tab);
			while (elems != null) {
				if (elems.length >= 3) {
					Triple<String, Chromosome, Integer> t = new Triple<>(elems[2], Chromosome.parseChr(elems[0]), Integer.parseInt(elems[1]));
					
					String filename = snpseqoutputdir + elems[2] + ".fa.gz";
					if (!Gpio.exists(filename)) {
						snps.add(t);
					}
				}
				elems = tfo.readLineElems(TextFile.tab);
			}
		}
		MakeTranscriptome t = new MakeTranscriptome();
		t.exportSNPs(snps, genomeFasta, windowsize, snpseqoutputdir, overwriteexistingsnps);
	}
	
	
	HashMap<String, Feature> snpmap;
	
	private void loadSNPMap(String snpmapfile, String snplist) throws IOException {
		
		HashSet<String> query = null;
		if (snplist != null) {
			query = new HashSet<>();
			TextFile tf = new TextFile(snplist, TextFile.R);
			query.addAll(tf.readAsArrayList());
			tf.close();
			System.out.println(query.size() + " snps loaded from " + snplist);
		}
		
		System.out.println("Loading snp map: " + snpmapfile);
		TextFile tf = new TextFile(snpmapfile, TextFile.R);
		
		snpmap = new HashMap<String, Feature>();
		
		String[] elems = tf.readLineElems(TextFile.tab);
		int lnctr = 0;
		while (elems != null) {
			if (elems.length >= 3) {
				String snp = elems[2];
				
				if (query == null || query.contains(snp)) {
					snp = new String(elems[2].getBytes(), "UTF-8");
					Chromosome chr = Chromosome.parseChr(elems[0]);
					Integer pos = Integer.parseInt(elems[1]);
					Feature f = new Feature();
					f.setChromosome(chr);
					f.setStart(pos);
					f.setStop(pos);
					f.setName(snp);
					snpmap.put(snp, f);
				}
			}
			lnctr++;
			if (lnctr % 10000 == 0) {
				System.out.print("\r" + lnctr + " lines parsed");
				
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		System.out.println();
		
		System.out.println("SNP map has " + snpmap.size() + " elements.");
		
	}
	
	private class RegionCollectionTask implements Callable<ArrayList<Feature>> {
		private final String genomeFasta;
		private final String seqoutputfolder;
		int clusterdistance;
		String dbloc;
		String traitname;
		double pvalthresh;
		String pvalstr;
		String outputfolder;
		int windowsize;
		int id;
		
		public RegionCollectionTask(String dbloc, String traitname, Double pvalthresh, String pvalstr, int clusterdistance,
									int windowsize, String regionDefinitionOutputFolder, String genomeFasta, String seqoutputfolder, int i) throws IOException {
			this.dbloc = dbloc;
			this.traitname = traitname;
			this.pvalstr = pvalstr;
			this.pvalthresh = pvalthresh;
			this.clusterdistance = clusterdistance;
			this.windowsize = windowsize;
			this.outputfolder = regionDefinitionOutputFolder;
			this.id = i;
			this.genomeFasta = genomeFasta;
			this.seqoutputfolder = seqoutputfolder;
		}
		
		
		@Override
		public ArrayList<Feature> call() {
			ArrayList<Feature> tagsnps = new ArrayList<>();
			try {
				
				
				String traitloc = dbloc + traitname;
				TextFile tf = new TextFile(traitloc, TextFile.R);
				tf.readLine();
				System.out.println("Parsing " + traitloc + " with threshold: " + pvalthresh);
				String[] elems = tf.readLineElems(TextFile.tab);
				ArrayList<SNPFeature> tmp = new ArrayList<>();
				while (elems != null) {
					if (elems.length >= 4) {
						String snp = new String(elems[0].getBytes(), "UTF-8");
						double p = Double.parseDouble(elems[3]);
						if (p < pvalthresh) {
							// includes
							Feature snpf = snpmap.get(snp);
							if (snpf != null && snpf.getChromosome().isAutosome()) {
								SNPFeature f = new SNPFeature(snpf.getChromosome(), snpf.getStart(), snpf.getStop());
								f.setStart(f.getStart());
								if (f.getStart() < 0) {
									f.setStart(0);
								}
								f.setName(snp);
								f.setStop(f.getStop());
								f.setP(p);
								//
								f.useNameForComparison(false);
								f.useAllelesForComparison(false);
								
								tmp.add(f);
							}
						}
					}
					elems = tf.readLineElems(TextFile.tab);
				}
				tf.close();
				System.out.println("Found " + tmp.size() + " snps that match threshold: " + pvalthresh);
				try {
					Collections.sort(tmp, new FeatureComparator(true));
				} catch (IllegalArgumentException e) {
					System.out.println("Meh");
					FeatureComparator c = new FeatureComparator(true);
					for (Feature f1 : tmp) {
						for (Feature f2 : tmp) {
							for (Feature f3 : tmp) {
								if (c.compare(f1, f2) > 0 && c.compare(f2, f3) > 0 && c.compare(f1, f3) < 0) {
									System.out.println("This should break (bigger than)");
									System.out.println(f1.toString());
									System.out.println(f2.toString());
									System.out.println(f3.toString());
								} else if (c.compare(f1, f2) < 0 && c.compare(f2, f3) < 0 && c.compare(f1, f3) > 0) {
									System.out.println("This should also break (smaller than)");
									System.out.println(f1.toString());
									System.out.println(f2.toString());
									System.out.println(f3.toString());
								} else if (c.compare(f1, f2) == 0 && c.compare(f2, f3) == 0 && c.compare(f1, f3) != 0) {
									System.out.println("This should also break (equals)");
									System.out.println(f1.toString());
									System.out.println(f2.toString());
									System.out.println(f3.toString());
								}
							}
						}
					}
					
					e.printStackTrace();
					System.exit(-1);
				}
				
				
				// make lists of overlapping features
				TextFile outf = new TextFile(outputfolder + traitname + "_P" + pvalstr + ".txt.gz", TextFile.W);
				if (tmp.size() == 1) {
					SNPFeature hit = tmp.get(0);
					Feature f = snpmap.get(hit.getName());
					Feature snpfeature = new Feature(f);
					Feature region = new Feature(f);
					region.setStart(region.getStart() - windowsize);
					region.setStop(region.getStop() + windowsize);
					tagsnps.add(snpfeature);
					outf.writeln(region.toString() + "\t" + hit.toString() + "\t" + hit.getP() + "\tTRUE");
				} else if (tmp.size() > 1) {
					
					// add a certain window around each snp
					for (SNPFeature f : tmp) {
						f.setStart(f.getStart() - clusterdistance);
						if (f.getStart() < 0) {
							f.setStart(0);
						}
						f.setStop(f.getStop() + clusterdistance);
					}
					
					ArrayList<ArrayList<SNPFeature>> clusters = new ArrayList<>();
					
					SNPFeature current = tmp.get(0);
					ArrayList<SNPFeature> currentCluster = new ArrayList<>();
					currentCluster.add(current);
					int index = 1;
					while (index < tmp.size()) {
						SNPFeature fi = tmp.get(index);
						if (!fi.overlaps(current)) {
							clusters.add(currentCluster);
							currentCluster = new ArrayList<>();
						}
						currentCluster.add(fi);
						current = fi;
						index++;
					}
					clusters.add(currentCluster);
					System.out.println(clusters.size() + " clusters of snps detected");
					
					// select lowest p in each cluster
					for (ArrayList<SNPFeature> fa : clusters) {
						double lowestP = 2;
						SNPFeature tophit = null;
						for (SNPFeature f : fa) {
							if (f.getP() < lowestP) {
								lowestP = f.getP();
								tophit = f;
							}
						}
						
						// define region around top hit
						Feature f = snpmap.get(tophit.getName());
						Feature snpfeature = new Feature(f);
						Feature region = new Feature(f);
						region.setName(f.getName());
						region.setStart(region.getStart() - windowsize);
						region.setStop(region.getStop() + windowsize);
						tagsnps.add(snpfeature);
						for (SNPFeature hit : fa) {
							outf.writeln(region.toString() + "\t" + hit.toString() + "\t" + hit.getP() + "\t" + tophit.getName().equals(hit.getName()));
						}
						
					}
					
				}
				outf.close();
				
			} catch (Exception e) {
				e.printStackTrace();
			}
			System.out.println("Defined " + tagsnps.size() + " tag snps for trait " + traitname);
			return tagsnps;
		}
	}
}
