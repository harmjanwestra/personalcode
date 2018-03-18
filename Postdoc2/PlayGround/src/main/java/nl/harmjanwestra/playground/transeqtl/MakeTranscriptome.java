package nl.harmjanwestra.playground.transeqtl;

import htsjdk.samtools.SAMRecord;
import javafx.scene.control.ProgressBar;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;
import nl.harmjanwestra.utilities.annotation.ensembl.EnsemblStructures;
import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.*;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class MakeTranscriptome {
	
	public static void main(String[] args) {
		MakeTranscriptome mt = new MakeTranscriptome();
		if (args.length == 0) {
			System.out.println("Usage: makeseq || align || quantify || combine ");
		} else if (args[0].equals("makeseq")) {
			
			if (args.length < 8) {
				System.out.println("Usage: makeseq eqtlfile genome ensembl.gtf.gz outdir windowsize readlen shift");
				System.out.println("Defaults: windowsize=2000000 readlen=25 shift=2");
				System.exit(-1);
			} else {
				String eqtlfile = args[1]; // "D:\\Work\\eqtlgeneremap\\eQTLs-allsnpprobecombos.txt.gz";
				String genomefasta = args[2]; // "D:\\Data\\HumanGenome\\human_g1k_v37.fasta.gz";
				String ensemblannotation = args[3]; // "D:\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
				String outdir = args[4]; // "D:\\Work\\eqtlgeneremap\\";
				int windowsize = Integer.parseInt(args[5]); // 2000000;
				int readlen = Integer.parseInt(args[6]); // 25;
				int shift = Integer.parseInt(args[7]); // 2;
				boolean exportindividualgenes = Boolean.parseBoolean(args[8]);
				boolean exportexonseqs = Boolean.parseBoolean(args[9]);
				try {
					mt.createSequences(eqtlfile, genomefasta, ensemblannotation, outdir, windowsize, readlen, shift, exportindividualgenes, exportexonseqs);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}

//			String eqtlfile = "D:\\Work\\eqtlgeneremap\\eQTLs-allsnpprobecombos.txt.gz";
////		String eqtlfile = "D:\\Work\\eqtlgeneremap\\test.txt";
//			String genomefasta = "D:\\Data\\HumanGenome\\human_g1k_v37.fasta.gz";
//			String ensemblannotation = "D:\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
//			String outdir = "D:\\Work\\eqtlgeneremap\\";
//			int windowsize = 2000000;
//			int readlen = 25;
//			int shift = 2;
		
		
		} else if (args[0].equals("align")) {
			
			if (args.length < 3) {
				System.out.println("Usage: align indir eqtlfile");
			} else {
				String indir = args[1];
				String eqtlfile = args[2];
				try {
					mt.align(indir, eqtlfile);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			
		} else if (args[0].equals("quantify")) {
			if (args.length < 4) {
				System.out.println("Usage quantify eqtlfile indir out");
			} else {
				String eqtlfile = args[1];
				String indir = args[2];
				String out = args[3];
				try {
					mt.quantify(eqtlfile, indir, out);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		} else if (args[0].equals("combine")) {
			
			if (args.length < 4) {
				System.out.println("Usage: combine alignout eqtl output");
			} else {
				
				String alignout = args[1];
				String eqtl = args[2];
				String output = args[3];
				try {
					mt.combine(alignout, eqtl, output);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
		}
	}
	
	public void combine(String alignmentoutput, String eqtlfile, String output) throws IOException {
		
		HashMap<String, String> snpprobeToOutput = new HashMap<String, String>();
		TextFile tf = new TextFile(alignmentoutput, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			snpprobeToOutput.put(elems[0] + "_" + elems[1], Strings.concat(elems, Strings.tab, 2, elems.length));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile out = new TextFile(output, TextFile.W);
		QTLTextFile qtf = new QTLTextFile(eqtlfile, QTLTextFile.R);
		out.writeln("SNP\tSNPChr\tSNPChrPos\tGene\tGeneChr\tGeneChrPos\tZ\tAbsZ\tNrReads\t%BasesMapped\t%ReadsMapped");
		
		Iterator<EQTL> it = qtf.getEQtlIterator();
		while (it.hasNext()) {
			EQTL e = it.next();
			String query = e.getRsName() + "_" + e.getProbe();
			String append = snpprobeToOutput.get(query);
			
			String lnout = e.getRsName()
					+ "\t" + e.getRsChr()
					+ "\t" + e.getRsChrPos()
					+ "\t" + e.getProbe()
					+ "\t" + e.getProbeChr()
					+ "\t" + e.getProbeChrPos()
					+ "\t" + e.getZscore()
					+ "\t" + e.getZscoreAbs()
					+ "\t" + append;
			out.writeln(lnout);
		}
		tf.close();
		out.close();
		
	}
	
	public void quantify(String eqtlfile, String indir, String out) throws IOException {
		
		// inventory for genes.
		TextFile tf2 = new TextFile(eqtlfile, TextFile.R);
		HashSet<String> genes = new HashSet<>();
		tf2.readLine();
		String[] tf2elems = tf2.readLineElems(TextFile.tab);
		while (tf2elems != null) {
			String genename = tf2elems[3];
			genes.add(genename);
			tf2elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		HashMap<String, Integer> geneToNrReads = new HashMap<String, Integer>();
		for (String genename : genes) {
			String fasta = indir + "/genes/" + genename + ".fa.gz";
			TextFile faIN = new TextFile(fasta, TextFile.R);
			int nrlines = faIN.countLines();
			nrlines /= 2;
			geneToNrReads.put(genename, nrlines);
		}
		
		// iterate the QTL file once more
		TextFile tf = new TextFile(eqtlfile, TextFile.R);
		tf.readLine(); // skip header
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		
		TextFile outtf = new TextFile(out, TextFile.W);
		outtf.writeln("SNP\tGene\tNrReads\t%BasesMapped\t%ReadsMapped");
		
		while (elems != null) {
			Chromosome chr = Chromosome.parseChr(elems[1]);
			Integer snppos = Integer.parseInt(elems[2]);
			String snpname = elems[0];
			String genename = elems[3];
			
			// get the alignment
			String alignment = indir + "/alignments/" + snpname + "_" + genename + ".sam";
			if (Gpio.exists(alignment)) {
//				System.out.println(alignment);
				SAMFileReader reader = new SAMFileReader(new File(alignment));
				SAMRecordIterator it = reader.iterator();
				double totalperc = 0;
				
				// iterate all alignments, take the ones that are properly mapping
				while (it.hasNext()) {
					net.sf.samtools.SAMRecord record = it.next();
					List<AlignmentBlock> blocks = record.getAlignmentBlocks();
					int readlen = record.getReadLength();
					int mappedlen = 0;
					for (AlignmentBlock b : blocks) {
						mappedlen += b.getLength();
					}
					
					double perc = (double) mappedlen / readlen;
					totalperc += perc;
				}
				
				int nrReadsTotal = geneToNrReads.get(genename);
				
				outtf.writeln(snpname + "\t" + genename + "\t" + nrReadsTotal + "\t" + totalperc + "\t" + (totalperc / nrReadsTotal));
				
				reader.close();
			} else {
				System.out.println("Could not find alignment: " + alignment);
			}
			
			
			elems = tf.readLineElems(TextFile.tab);
			ctr++;
			if (ctr % 100 == 0) {
				System.out.print("\r" + ctr + " lines parsed..");
			}
		}
		tf.close();
		outtf.close();
	}
	
	
	public void align(String indir, String eqtlfile) throws IOException {
		
		ExecutorService ex = Executors.newFixedThreadPool(16);
		
		TextFile tf = new TextFile(eqtlfile, TextFile.R);
		
		tf.readLine();
		int nrAlignments = tf.countLines();
		tf.close();
		tf.open();
		
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		int nrsubmitted = 0;
		AtomicInteger nrdone = new AtomicInteger();
		nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar pb = new nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar(nrAlignments);
		while (elems != null) {
			
			String snpname = elems[0];
			Gpio.createDir(indir + "/alignments/" + snpname + "/");
			String genename = elems[3];
			
			AlignmentTask t = new AlignmentTask();
			t.gene = genename;
			t.snp = snpname;
			t.dir = indir;
			t.nrdone = nrdone;
			
			ex.submit(t);
			nrsubmitted++;
			
			// don't let the buffer run too full
//			if (nrdone.get() - nrsubmitted > 128) {
//				try {
//					Thread.sleep(1000);
//				} catch (InterruptedException e) {
//					e.printStackTrace();
//				}
//			}

//			if (nrsubmitted == 10000) {
//				while (nrdone.get() < nrsubmitted) {
//					try {
//
//						pb.set(nrdone.get());
//						Thread.sleep(1000);
//					} catch (InterruptedException e) {
//						e.printStackTrace();
//					}
//				}
//				System.exit(-1);
//			}
			
			if (nrsubmitted % 1000 == 0) {
				pb.set(nrdone.get());
			}
			
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		while (nrdone.get() < nrsubmitted) {
			try {
				pb.set(nrdone.get());
//				System.out.print("\r" + nrdone.get() + "/" + nrsubmitted);
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		ex.shutdown();
		pb.close();
		
	}
	
	
	public class AlignmentTask implements Runnable {
		
		String snp;
		String gene;
		String dir;
		public AtomicInteger nrdone;
		
		@Override
		public void run() {


//			String ln = "bwa aln -l 10 -0 " + dir + "/snps/" + snp + ".fa.gz " + dir + "/genes/" + gene + ".fa.gz > " + dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sai";
			String[] cmd1 = new String[]{
					"bwa",
					"aln",
					"-l 10",
					"-0",
					dir + "/snps/" + snp + ".fa.gz",
					dir + "/genes/" + gene + ".fa.gz"
				
			};


//			String ln2 = "bwa samse " + dir + "/snps/" + snp + ".fa.gz " + dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sai " + dir + "/genes/" + gene + ".fa.gz > " + dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sam";
			String[] cmd2 = new String[]{
					"bwa",
					"samse",
					dir + "/snps/" + snp + ".fa.gz",
					dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sai",
					dir + "/genes/" + gene + ".fa.gz"
				
			};
			try {
				ProcessBuilder pb = new ProcessBuilder(cmd1);
				pb.redirectOutput(new File(dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sai"));
				Process process = pb.start();
				
				int errCode = process.waitFor();
//				System.out.println("Echo command executed, any errors? " + (errCode == 0 ? "No" : "Yes"));
				ProcessBuilder pb2 = new ProcessBuilder(cmd2);
				pb2.redirectOutput(new File(dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sam"));
				Process process2 = pb2.start();
				int errCode2 = process2.waitFor();
//				System.out.println("Echo command executed, any errors? " + (errCode2 == 0 ? "No" : "Yes"));
			} catch (IOException e) {
				e.printStackTrace();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			nrdone.getAndIncrement();
			
		}
	}
	
	
	public void createSequences(String eqtlfile, String genomeFasta, String ensemblannotation, String outputdir, int windowsize, int readlen, int shift, boolean exportindividualgenes, boolean exportfullexonsequences) throws IOException {
		
		
		String geneoutputdir = outputdir + "/genes/";
		String snpoutputdir = outputdir + "/snps/";
		Gpio.createDir(geneoutputdir);
		Gpio.createDir(snpoutputdir);
		
		System.out.println("Output of genes: " + geneoutputdir);
		System.out.println("Output of snps: " + snpoutputdir);
		
		// create reference genome from snps
		TextFile tf = new TextFile(eqtlfile, TextFile.R);
		tf.readLine(); // skip header
		
		
		ArrayList<Triple<String, Chromosome, Integer>> snps = new ArrayList<>();
		HashSet<String> snpids = new HashSet<String>();
		HashSet<String> geneIds = new HashSet<String>();
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		while (elems != null) {
			Chromosome chr = Chromosome.parseChr(elems[1]);
			Integer snppos = Integer.parseInt(elems[2]);
			String snpname = elems[0];
			String genename = elems[3];
			if (!snpids.contains(snpname)) {
				Triple<String, Chromosome, Integer> snp = new Triple<>(snpname, chr, snppos);
				snps.add(snp);
				snpids.add(snpname);
			}
			
			geneIds.add(genename);
			elems = tf.readLineElems(TextFile.tab);
			ctr++;
			if (ctr % 100000 == 0) {
				System.out.print("\r" + ctr + " lines parsed..");
			}
		}
		tf.close();
		System.out.println();
		System.out.println(snps.size() + " snps to export");
		System.out.println(geneIds.size() + " genes to export");
		
		GTFAnnotation s = new GTFAnnotation(ensemblannotation);
		
		// convert gene strings to ensembl gene ids
		HashMap<String, Gene> geneHash = s.getStrToGene();
		ArrayList<Gene> genesToExport = new ArrayList<>();
		for (String g : geneIds) {
			Gene gobj = geneHash.get(g);
			if (gobj != null) {
				genesToExport.add(gobj);
			}
		}
		
		
		System.out.println(genesToExport.size() + " genes found in annotation.");
		
		if (genesToExport.isEmpty()) {
			System.err.println("No genes to export..?");
			System.exit(-1);
		}
		
		// create 1mb regions around snps
		System.out.println("Loading genome from: " + genomeFasta);
		htsjdk.samtools.reference.FastaSequenceFile fastaFile = new htsjdk.samtools.reference.FastaSequenceFile(new File(genomeFasta), false);
		htsjdk.samtools.reference.ReferenceSequence seq = fastaFile.nextSequence();
//		System.out.println("available sequences");
//		while (seq != null) {
//			String seqname = seq.getName();
//			byte[] bases = seq.getBases();
//			byte[] base10 = new byte[10];
//			System.arraycopy(bases, 0, base10, 0, 10);
//			System.out.println(seqname + "\t" + bases.length + "\t" + byteToStr(base10) + "\t" + seq.length());
//			seq = fastaFile.nextSequence();
//
//		}
//		System.exit(-1);
		
		
		int nrsnpsexported = 0;
		int nrgenesexported = 0;
		TextFile allgenefa = new TextFile(geneoutputdir + "allgenes.fa.gz", TextFile.W);
		
		while (seq != null) {
			
			byte[] bases = seq.getBases();
			
			String seqname = seq.getName();
			while (seqname.contains("  ")) {
				seqname = seqname.replaceAll("  ", " ");
			}
			String[] seqnamelems = seqname.split(" ");
			seqname = seqnamelems[0];
			System.out.println("Getting sequence: " + seq.getName() + "\t" + seq.length() + "\tParsed seq: " + seqname);
			
			Chromosome seqchr = Chromosome.parseChr(seqname);
			for (Triple<String, Chromosome, Integer> snp : snps) {
				if (snp.getMiddle().equals(seqchr)) {
					
					int snppos = snp.getRight();
					int windowleft = snppos - (windowsize / 2);
					int windowright = snppos + (windowsize / 2);
					
					if (windowleft < 0) {
						windowleft = 0;
					}
					if (windowright > bases.length - 1) {
						windowright = bases.length - 1;
					}
					
					
					int nrbases = windowright - windowleft;
					if (nrbases <= 0) {
						System.err.println("Error: <=0 bases..");
						System.err.println(nrsnpsexported + "/" + snps.size() + " | SNP: " + snp + " pos: " + snppos + " chr " + seqchr + " left: " + windowleft + " right: " + windowright + " size: " + (windowright - windowleft));
						System.exit(-1);
					}
					byte[] snpwindow = new byte[nrbases];
					System.arraycopy(bases, windowleft, snpwindow, 0, nrbases);
					
					// make into fasta
					TextFile snpfastaout = new TextFile(snpoutputdir + snp.getLeft() + ".fa.gz", TextFile.W);
					String fastaStr = ">" + snp.getMiddle().toString() + "_" + snp.getRight() + "_" + snp.getLeft() + "_" + windowleft + "_" + windowright
							+ "\n" + byteToStr(snpwindow);
					snpfastaout.writeln(fastaStr);
					snpfastaout.close();
					System.out.println(nrsnpsexported + "/" + snps.size() + " | SNP: " + snp + " chr " + seqchr + " left: " + windowleft + " right: " + windowright + " size: " + (windowright - windowleft));
					nrsnpsexported++;
				}
			}
			
			System.out.println("Done exporting snps for chr " + seqchr);
			
			System.out.println("Now exporting genes");
			// now parse the ensemblgenes on this chr
			for (Gene g : genesToExport) {
				if (g.getChromosome().equals(seqchr)) {
					TextFile genesfastaout = null;
					if (exportindividualgenes) {
						genesfastaout = new TextFile(geneoutputdir + g.getName() + ".fa.gz", TextFile.W);
					}
					ArrayList<Transcript> transcripts = g.getTranscripts();
					HashSet<Exon> allexons = new HashSet<Exon>();
					for (Transcript t : transcripts) {
						ArrayList<Exon> exons = t.getExons();
						allexons.addAll(exons);
					}
//					System.out.println("Gene: " + g.getName() + "\tTotal exons: " + allexons.size());
					
					if (exportfullexonsequences) {
						TextFile genesfastaoutFull = new TextFile(geneoutputdir + g.getName() + "_exons.fa.gz", TextFile.W);
						for (Exon currentexon : allexons) {
							int len = currentexon.getStop() - currentexon.getStart();
							byte[] b = new byte[len];
							System.arraycopy(bases, currentexon.getStart(), b, 0, len);
							String fastaStr = ">" + g.getName() + "_" + currentexon.getName() + "_" + currentexon.getChromosome().toString() + "_" + currentexon.getStart() + "_" + currentexon.getStop()
									+ "\n" + byteToStr(b);
							genesfastaoutFull.writeln(fastaStr);
						}
						genesfastaoutFull.close();
					}
					// split up the transcript in readlen reads
					int nrseqs = 0;
					for (Transcript t : transcripts) {
						ArrayList<Exon> exons = t.getExons();
						
						Collections.sort(exons, new FeatureComparator(true)); // sort the exons
						// output exon sequences
						int exonctr = 0;
						int tstart = exons.get(0).getStart();

//						if (exons.size() > 1) {
//							System.out.println("Bench this");
//						}
						
						int tstop = exons.get(exons.size() - 1).getStop();
						int currentPos = tstart;
						Exon currentexon = exons.get(0);
						while (currentPos < tstop) {
							
							// iterate current exon
							while (currentPos + readlen <= currentexon.getStop()) {
								byte[] b = new byte[readlen];
								System.arraycopy(bases, currentPos, b, 0, readlen);
								String fastaStr = ">" + g.getName() + "_" + t.getName() + "_" + currentexon.getName() + "_" + currentexon.getChromosome().toString() + "_" + currentPos + "_" + (currentPos + readlen)
										+ "\n" + byteToStr(b);
								if (exportindividualgenes) {
									genesfastaout.writeln(fastaStr);
								}
								allgenefa.writeln(fastaStr);
// System.out.println(nrseqs + "\texon: " + exonctr + "\tsta: " + currentPos + "\tsto: " + (currentPos + readlen) + "\texonstop: " + currentexon.getStop());
								nrseqs++;
								currentPos += shift;
								
							}
							
							// process exon exon boundary
							
							// if there are exons remaining
							if (exonctr + 1 < exons.size()) {
								while (currentPos < currentexon.getStop()) { // iterate until we run out of currentexon
									int basesFromCurrentExon = currentexon.getStop() - currentPos;
									byte[] b = new byte[readlen];
									System.arraycopy(bases, currentPos, b, 0, basesFromCurrentExon);
									
									int basesRemaining = readlen - basesFromCurrentExon;
									// now we should check whether the next exon actually has enough bases
									int tmpexonctr = exonctr + 1;
									while (basesRemaining > 0 && tmpexonctr < exons.size()) {
										// get remaining bases from next exons
										Exon tmpexon = exons.get(tmpexonctr);
										int tmpexonstart = tmpexon.getStart();
										int tmpexonstop = tmpexon.getStop();
										
										int basesToTakeFromExon = basesRemaining;
										if (tmpexonstart + basesRemaining > tmpexonstop) {
											basesToTakeFromExon = tmpexonstop - tmpexonstart;
										}
										
										
										System.arraycopy(bases, tmpexonstart, b, basesFromCurrentExon, basesToTakeFromExon);
										basesRemaining -= basesToTakeFromExon;
//										System.out.println("Boundary: " + nrseqs + "\texon: " + exonctr + "\ttmpexon: " + tmpexonctr + "\tsta: " + currentPos + "\texonstop: " + currentexon.getStop() + "\t" + basesRemaining);
										tmpexonctr++;
									}
									if (basesRemaining == 0) {
										// write the resulting freakshow of a sequence
										String fastaStr = ">" + g.getName() + "_" + t.getName() + "_" + currentexon.getName() + "_" + currentPos + "_" + (currentPos + readlen)
												+ "\n" + byteToStr(b);
										if (exportindividualgenes) {
											genesfastaout.writeln(fastaStr);
										}
										allgenefa.writeln(fastaStr);
										nrseqs++;
									}
									
									currentPos += shift;
								}
							}
							
							// start at next exon
							if (exonctr < exons.size()) {
//								System.out.println("Investigating exon: " + exonctr);
								currentexon = exons.get(exonctr);
								currentPos = currentexon.getStart();
							} else {
//								System.out.println("Reached end of transcript");
								currentPos = currentexon.getStop();
							}
							exonctr++;
						}
						
					}
					if (exportindividualgenes) {
						genesfastaout.close();
					}
					System.out.println(nrgenesexported + "/" + genesToExport.size() + " | Gene: " + g.getName() + " Chr: " + g.getChromosome().toString() + " start: " + g.getStart() + " stop: " + g.getStop() + " Strand: " + g.getStrand() + " nr sequences: " + nrseqs);
					nrgenesexported++;
				}
			}
			System.out.println("Done exporting genes");
			seq = fastaFile.nextSequence();
		}
		allgenefa.close();
		
		// make index script
		TextFile tfs1 = new TextFile(outputdir + "index.sh", TextFile.W);
		for (String snp : snpids) {
			System.out.println(snpoutputdir + snp + ".fa.gz");
			if (Gpio.exists(snpoutputdir + snp + ".fa.gz")) {
				tfs1.writeln("bwa index ./snps/" + snp + ".fa.gz");
			}
		}
		tfs1.close();
		
		TextFile tfsa1 = new TextFile(outputdir + "alignall.sh", TextFile.W);
		TextFile tfsa2 = new TextFile(outputdir + "samseall.sh", TextFile.W);
		tfsa1.writeln("mkdir -p alignments");
		for (String snp : snpids) {
			System.out.println(snpoutputdir + snp + ".fa.gz");
			if (Gpio.exists(snpoutputdir + snp + ".fa.gz")) {
				String ln = "bwa aln -t 16 -l 10 -0 ./snps/" + snp + ".fa.gz ./genes/allgenes.fa.gz > ./alignments/" + snp + "_allgenes.sai";
				String ln2 = "bwa samse ./snps/" + snp + ".fa.gz ./alignments/" + snp + "_allgenes.sai ./genes/allgenes.fa.gz > ./alignments/" + snp + "_allgenes.sam";
				tfsa1.writeln(ln);
				tfsa2.writeln(ln2);
			}
		}
		tfsa1.close();
		tfsa2.close();
		
		
		// make alignment script
		tf.open();
		tf.readLine();
		elems = tf.readLineElems(TextFile.tab);
		
		TextFile tfs2 = new TextFile(outputdir + "align.sh", TextFile.W);
		TextFile tfs3 = new TextFile(outputdir + "samse.sh", TextFile.W);
		tfs2.writeln("mkdir -p alignments");
		while (elems != null) {
			
			String snp = elems[0];
			String gene = elems[3];
			if (Gpio.exists(snpoutputdir + snp + ".fa.gz") && Gpio.exists(geneoutputdir + gene + ".fa.gz")) {
				String ln = "bwa aln -t 1 -l 10 -0 ./snps/" + snp + ".fa.gz ./genes/" + gene + ".fa.gz > ./alignments/" + snp + "_" + gene + ".sai";
				String ln2 = "bwa samse ./snps/" + snp + ".fa.gz ./alignments/" + snp + "_" + gene + ".sai ./genes/" + gene + ".fa.gz > ./alignments/" + snp + "_" + gene + ".sam";
				tfs2.writeln(ln);
				tfs3.writeln(ln2);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		tfs2.close();
		tfs3.close();
	}
	
	private String byteToStr(byte[] bytes) throws UnsupportedEncodingException {
		return new String(bytes, "UTF-8");
	}
}
