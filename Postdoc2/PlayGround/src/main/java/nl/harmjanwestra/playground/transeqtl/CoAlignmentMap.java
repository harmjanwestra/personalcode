package nl.harmjanwestra.playground.transeqtl;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import nl.harmjanwestra.playground.legacy.GTFAnnotation;

import org.apache.commons.cli.*;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Exon;
import umcg.genetica.features.FeatureType;
import umcg.genetica.features.Gene;
import umcg.genetica.features.Transcript;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;

public class CoAlignmentMap {

	public static void main(String[] args) {

		Options options = new Options();
		Option option = Option.builder("g")
				.longOpt("gtf")
				.argName("annotation.gtf.gz")
				.hasArg()
				.required()
				.desc("GTF file annotation")
				.build();
		options.addOption(option);
		option = Option.builder("o")
				.longOpt("out")
				.argName("outdir")
				.hasArg()
				.required()
				.desc("Output directory")
				.build();
		options.addOption(option);
		option = Option.builder("f")
				.longOpt("fasta")
				.argName("")
				.hasArg()
				.required()
				.desc("Fasta file")
				.build();
		options.addOption(option);
		option = Option.builder("w")
				.longOpt("window")
				.argName("int")
				.hasArg()
				.desc("Cis window size [default 5mb]")
				.build();
		options.addOption(option);
		option = Option.builder("r")
				.longOpt("readlen")
				.argName("int")
				.hasArg()
				.desc("Read length [default 50]")
				.build();
		options.addOption(option);
		option = Option.builder("s")
				.longOpt("readslide")
				.argName("int")
				.hasArg()
				.desc("Read sliding window size [default 10]")
				.build();
		options.addOption(option);


		CommandLineParser parser = new DefaultParser();

		try {
			CommandLine cmd = parser.parse(options, args);
			String gtfFile = cmd.getOptionValue("g");
			String genomeFasta = cmd.getOptionValue("f");
			String outdir = cmd.getOptionValue("o");

			int windowsize = 5000000;
			int readLength = 50;
			int readOverlap = 10;
			if (cmd.hasOption("w")) {
				windowsize = Integer.parseInt(cmd.getOptionValue('w'));
			}
			if (cmd.hasOption("r")) {
				readLength = Integer.parseInt(cmd.getOptionValue('r'));
			}
			if (cmd.hasOption("s")) {
				readOverlap = Integer.parseInt(cmd.getOptionValue('s'));
			}

			if (!outdir.endsWith("/")) {
				outdir = outdir + "/";
			}

			String refOutdirPrefix = outdir + "cisgenes/";
			String readOutputDir = outdir + "genereads/";
			String jobOutputDir = outdir + "jobs/";

			CoAlignmentMap c = new CoAlignmentMap();
			c.createReferenceSequences(gtfFile, genomeFasta, windowsize, refOutdirPrefix);
			c.createFakeReads(gtfFile, genomeFasta, readLength, readOverlap, readOutputDir);
			c.createAlignmentJobs(gtfFile, refOutdirPrefix, readOutputDir, jobOutputDir);

		} catch (ParseException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("ant", options);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public void createAlignmentJobs(String gtfFile, String refOutdirPrefix, String readOutputDir, String jobOutputDir) throws IOException {
		GTFAnnotation annotation = new GTFAnnotation(gtfFile);

		HashSet<Chromosome> chrs = new HashSet<>();
		for (Gene g : annotation.getGenes()) {
			Chromosome chr = g.getChromosome();
			chrs.add(chr);
		}

		TextFile tfA = new TextFile(jobOutputDir + "align.sh", TextFile.W);
		TextFile tfB = new TextFile(jobOutputDir + "samse.sh", TextFile.W);
		for (Gene g : annotation.getGenes()) {
			String refOutputDir = refOutdirPrefix + "/" + g.getName() + "/";
			String reffile = refOutputDir + "ref-" + g.getName() + ".fq.gz";

			for (Chromosome chr : chrs) {
				String chrFile = readOutputDir + "genesequences-chr" + chr.getNumber() + ".fq.gz";
				String alignmentOutSaiFile = refOutdirPrefix + "alignments-chr" + chr.getNumber() + ".sai";
				String alignmentOutSamFile = refOutdirPrefix + "alignments-chr" + chr.getNumber() + ".sam";
				String alingmentStr = "bwa aln -t 12 -l 10 -0 " + reffile + " " + chrFile + " > " + alignmentOutSaiFile;
				String samseStr = "bwa samse " + reffile + " " + alignmentOutSaiFile + " " + chrFile + " > " + alignmentOutSamFile;
				tfA.writeln(alingmentStr);
				tfB.writeln(samseStr);
			}
		}
		tfA.close();
		tfB.close();
	}

	public void createFakeReads(String gtfFile,
								String genomeFasta,
								int readLength,
								int readOverlap,
								String outdir) throws IOException {
		System.out.println("Creating fake reads with:");
		System.out.println("Read length: " + readLength);
		System.out.println("Read sliding window: " + readOverlap);
		GTFAnnotation annotation = new GTFAnnotation(gtfFile);

		Gpio.createDir(outdir);

		System.out.println("Loading genome from: " + genomeFasta);
		FastaSequenceFile fastaFile = new FastaSequenceFile(new File(genomeFasta), false);
		ReferenceSequence seq = fastaFile.nextSequence();

		HashMap<Chromosome, ArrayList<Gene>> genesPerChr = new HashMap<>();
		for (Gene g : annotation.getGenes()) {
			Chromosome chr = g.getChromosome();
			ArrayList<Gene> genesForChr = genesPerChr.get(chr);
			if (genesForChr == null) {
				genesForChr = new ArrayList<>();
			}
			genesForChr.add(g);
			genesPerChr.put(chr, genesForChr);
		}


		ArrayList<Pair<Chromosome, byte[]>> referenceGenome = new ArrayList<Pair<Chromosome, byte[]>>();
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
			referenceGenome.add(new Pair<>(seqchr, bases));
			seq = fastaFile.nextSequence();
		}

		AtomicInteger genesWritten = new AtomicInteger(0);
		referenceGenome.parallelStream().forEach(refseq -> {

			Chromosome seqchr = refseq.getLeft();
			byte[] bases = refseq.getRight();
			System.out.println("Now exporting genes");
			ArrayList<Gene> genesForChr = genesPerChr.get(seqchr);

			try {

				TextFile chrFastaOut = new TextFile(outdir + "genesequences-chr" + seqchr.getNumber() + ".fq.gz", TextFile.W);

				for (Gene g : genesForChr) {

					ArrayList<Transcript> transcripts = g.getTranscripts();
					for (Transcript t : transcripts) {
						int transcriptLength = 0;
						Exon[] exons = t.getExonsRanked();
						for (Exon e : exons) {
							if (e.getType().equals(FeatureType.EXON)) {
								transcriptLength += e.getSize();
							}
						}
						byte[] transcriptBases = new byte[transcriptLength];
						int basesCopied = 0;
						for (Exon e : exons) {
							if (e.getType().equals(FeatureType.EXON)) {
								int start = e.getStart();
								System.arraycopy(bases, start, transcriptBases, basesCopied, e.getSize());
								basesCopied += e.getSize();
							}
						}

						// generate reads
						int curBasePos = 0;
						byte[] read = new byte[readLength];
						int readctr = 0;
						while (curBasePos + readLength < transcriptBases.length) {
							System.arraycopy(transcriptBases, curBasePos, read, 0, readLength);

							String outFastaLn = ">" + g.getName() + "-" + t.getName() + "-" + readctr + "\n" + byteToStr(read);
							readctr += 1;
							chrFastaOut.writeln(outFastaLn);
							curBasePos += readOverlap;
						}
					}
					int gw = genesWritten.getAndIncrement();
					if (gw % 100 == 0) {
						System.out.print(gw + " genes written sofar out of " + annotation.getGenes().size() + "\n");
					}
				}
				chrFastaOut.close();
			} catch (Exception e) {
				e.printStackTrace();
			}


		});
		System.out.print(genesWritten + " genes written sofar out of " + annotation.getGenes().size() + "\n");

	}

	public void createReferenceSequences(String gtfFile, String genomeFasta, int windowsize, String outdir) throws
			IOException {
		GTFAnnotation annotation = new GTFAnnotation(gtfFile);

		Gpio.createDir(outdir);

		System.out.println("Loading genome from: " + genomeFasta);
		FastaSequenceFile fastaFile = new FastaSequenceFile(new File(genomeFasta), false);
		ReferenceSequence seq = fastaFile.nextSequence();

		HashMap<Chromosome, ArrayList<Gene>> genesPerChr = new HashMap<>();
		for (Gene g : annotation.getGenes()) {
			Chromosome chr = g.getChromosome();
			ArrayList<Gene> genesForChr = genesPerChr.get(chr);
			if (genesForChr == null) {
				genesForChr = new ArrayList<>();
			}
			genesForChr.add(g);
			genesPerChr.put(chr, genesForChr);
		}

		TextFile jobIndexOut = new TextFile(outdir + "index.sh", TextFile.W);
		AtomicInteger genesWrittenAtomic = new AtomicInteger(0);

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

			System.out.println("Now exporting genes");
			ArrayList<Gene> genesForChr = genesPerChr.get(seqchr);
			ReferenceSequence finalSeq = seq;
			genesForChr.parallelStream().forEach(g -> {
				int start = g.getStart();
				int stop = g.getStop();
				start = start - windowsize;
				stop = stop + windowsize;
				if (start < 0) {
					start = 0;
				}
				if (stop > finalSeq.length()) {
					stop = finalSeq.length() - 1;
				}
				int len = stop - start;
				if (len > 0) {
					byte[] b = new byte[len];
					System.arraycopy(bases, start, b, 0, len);
					String baseoutput = new String(b, StandardCharsets.UTF_8);
					String fastaStr = ">" + g.getName() + "_"
							+ g.getChromosome().toString() + ":" + g.getStart() + "-" + g.getStop()
							+ "+w" + windowsize
							+ "\n" + baseoutput;


					String gOutputDir = outdir + "/" + g.getName() + "/";
					try {
						createDir(gOutputDir);

						String gOutputFile = gOutputDir + "ref-" + g.getName() + ".fq.gz";
						TextFile tf = new TextFile(gOutputFile, TextFile.W);
						tf.writeln(fastaStr);
						tf.close();
						jobIndexOut.writelnsynced("bwa index " + gOutputFile);
						int genesWritten = genesWrittenAtomic.getAndIncrement();
						if (genesWritten % 10 == 0) {
							System.out.print(genesWritten + "/" + genesForChr.size() + " genes written sofar \r");
						}
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
				}
			});

			System.out.print(genesWrittenAtomic.get() + "/" + genesForChr.size() + " genes written sofar \n");
			seq = fastaFile.nextSequence();
		}
		jobIndexOut.close();
	}

	private String byteToStr(byte[] bytes) throws UnsupportedEncodingException {
		return new String(bytes, StandardCharsets.UTF_8);
	}

	public static void createDir(String dirName) throws IOException {
		if (!dirName.endsWith(getFileSeparator())) {
			dirName = dirName + getFileSeparator();
		}

		if (!exists(dirName)) {
			boolean success = (new File(dirName)).mkdirs();
		}

	}

	public static String getFileSeparator() {
		return System.getProperty("file.separator");
	}

	public static boolean exists(String dir) {
		return existsAndReadable(new File(dir));
	}

	public static boolean existsAndReadable(File file) {
		return file.exists() && file.canRead();
	}

}
