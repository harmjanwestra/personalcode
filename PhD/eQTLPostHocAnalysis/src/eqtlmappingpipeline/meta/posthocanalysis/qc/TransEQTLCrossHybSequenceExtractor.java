/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.qc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.*;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class TransEQTLCrossHybSequenceExtractor {

    int windowsize = 2500000; // 2.5mb window

    public static void main(String[] args) {
	try {
	    if (args.length < 2) {
		System.out.println("Usage: file.jar eqtlfilename outputdir");
	    } else {
		TransEQTLCrossHybSequenceExtractor t = new TransEQTLCrossHybSequenceExtractor();
		t.run(args[0], args[1]);
	    }
	} catch (IOException e) {
	    e.printStackTrace();
	}

    }
    private HashMap<String, String> probeSequences;
    String outfilename = "";
    private HashSet<String> eSNPs;
    private HashMap<String, Pair<Integer, Pair<Integer, Integer>>> eSNPPositions;
    private HashMap<String, ArrayList<String>> eQTLProbes;
    private String outdir;
    private String snp;

    public void run(String eQTLFileName, String dir) throws IOException {

	if (!dir.endsWith("/")) {
	    dir += "/";
	}

        this.outdir = dir;
	Gpio.createDir(dir);
	Gpio.createDir(this.outdir + "/probeseq/");
	Gpio.createDir(this.outdir + "/alignments/");
	Gpio.createDir(this.outdir + "/snpseq/");


	loadProbeSequences("/home/users/westra/MetaQC/2011-09-14-ProbeTranslationTable.log",
		this.outdir + "/probeseq/");


	eSNPs = new HashSet<String>();
	eSNPPositions = new HashMap<String, Pair<Integer, Pair<Integer, Integer>>>();


	eQTLProbes = new HashMap<String, ArrayList<String>>();

	TextFile eQTLFile = new TextFile(eQTLFileName, TextFile.R);

	// read all trans-eQTL snps
	// determine regions to extract from genome
	eQTLFile.readLineElemsReturnObjects(TextFile.tab); // headerline
	String[] elems = eQTLFile.readLineElemsReturnObjects(TextFile.tab);
	while (elems != null) {
	    String snp = elems[1];
	    eSNPs.add(snp);
	    Integer chrPos = Integer.parseInt(elems[3]);
	    String probe = elems[4];

	    ArrayList<String> probes = eQTLProbes.get(snp);
	    if (probes == null) {
		probes = new ArrayList<String>();
	    }

	    probes.add(probe);
	    eQTLProbes.put(snp, probes);

	    Pair<Integer, Integer> p = new Pair<Integer, Integer>(chrPos - windowsize, chrPos + windowsize);
	    Pair<Integer, Pair<Integer, Integer>> p2 = new Pair<Integer, Pair<Integer, Integer>>(Integer.parseInt(elems[2]), p);
	    eSNPPositions.put(snp, p2);
	    elems = eQTLFile.readLineElemsReturnObjects(TextFile.tab);
	}
	eQTLFile.close();

	// extract regions from genome by iterating over reference sequence..

	System.out.println(eSNPPositions.size() + " eSNPs loaded!");

	int nrSNPs = eSNPs.size();
	int snpsCtr = 0;

//	ExecutorService threadPool = Executors.newFixedThreadPool(16);
//	CompletionService pool = new ExecutorCompletionService(threadPool);

	TextFile outcommands = new TextFile(outdir + "runAlignments.sh", TextFile.W);
	int ctr = 0;
	for (String esnp : eSNPs) {
	    this.snp = esnp;

	    Pair<Integer, Pair<Integer, Integer>> positionAndRegion = eSNPPositions.get(esnp);
	    Integer chr = positionAndRegion.getLeft();
	    Pair<Integer, Integer> window = positionAndRegion.getRight();

	    String sequence = getPositionFromGenome(chr, window.getLeft(), window.getRight());


	    System.out.println("Writing for SNP\t" + snp + "\t" + snpsCtr + "/" + nrSNPs + "\t" + sequence.length());

	    TextFile out = new TextFile(outdir + "/snpseq/" + esnp + ".fa", TextFile.W);
	    out.writeln(">" + esnp);
	    out.writeln(sequence);
	    out.close();

//	    ArrayList<String> probes = eQTLProbes.get(callablesnp);
//	    Callable callable = new Callable() {
//
//		String callablesnp = snp;
//
//		@Override
//		public Object call() throws Exception {
//		    try {
//			callablesnp = snp;
//			Pair<Integer, Pair<Integer, Integer>> positionAndRegion = eSNPPositions.get(callablesnp);
//			Integer chr = positionAndRegion.getLeft();
//			Pair<Integer, Integer> window = positionAndRegion.getRight();
//
//			String sequence = getPositionFromGenome(chr, window.getLeft(), window.getRight());
//
//
////			System.out.println("Writing for SNP\t" + snp + "\t" + snpsCtr + "/" + nrSNPs + "\t" + sequence.length());
//			ArrayList<String> probes = eQTLProbes.get(callablesnp);

//			TextFile out = new TextFile(outdir + "/snpseq/" + callablesnp + ".fa", TextFile.W);
//			out.writeln(">" + callablesnp);
//			out.writeln(sequence);
//			out.close();

	    //
//			//			    out.writeln(">" + p);
////			    out.writeln(probeSequences.get(p));
////			for (String p : probes) {
////			    outfilename = outdir + "/snpseq/" + callablesnp + "-Probe" + p + ".fa";
////
////
////
////			}
//
////			for (String p : probes) {
////			    System.out.println("Command to run:\n"+"/home/users/westra/MetaQC/shrimp/SHRiMP_2_2_2/bin/gmapper " + outdir + "/probeseq/" + p + ".fa " + outdir + "/snpseq/" + callablesnp + ".fa -N 4 -o 5 -h 80% >"+outdir+"/alignments/"+callablesnp+"-"+p+".out 2>"+outdir+"/alignments/"+callablesnp+"-"+p+".log");
////
////			    Process proc = Runtime.getRuntime().exec("/home/users/westra/MetaQC/shrimp/SHRiMP_2_2_2/bin/gmapper " + outdir + "/probeseq/" + p + ".fa " + outdir + "/snpseq/" + callablesnp + ".fa -N 4 -o 5 -h 80% >"+outdir+"/alignments/"+callablesnp+"-"+p+".out 2>"+outdir+"/alignments/"+callablesnp+"-"+p+".log");
////			    proc.waitFor();
////			}
//
//		    } catch (IOException e) {
//			e.printStackTrace();
//		    }
//		    return null;
//		}
//
//		private String getPositionFromGenome(Integer chr, Integer left, Integer right) throws IOException {
//		    TextFile chromosomeSequenceFile = new TextFile("/home/users/westra/MetaQC/genomesequence/Homo_sapiens.NCBI36.54.dna.chromosome." + chr + ".fa.gz", TextFile.R);
//
//		    String line = chromosomeSequenceFile.readLine();
//		    int basesRead = 0;
//
//		    StringBuilder seq = new StringBuilder();
//		    while (line != null) {
//			if (line.startsWith(">")) {
//			    // fasta header line
//			} else {
//
//
//			    for (int i = 0; i < line.length(); i++) {
//
//				char c = line.charAt(i);
//				if (c == 'A' || c == 'T' || c == 'G' || c == 'C' || c == 'N') {
//				    if (basesRead >= left && basesRead < right) {
//					seq.append(c);
//				    } else if (basesRead >= left && basesRead >= right) {
//					break;
//				    }
//
//				} else {
//				    System.err.println("Could not parse char: " + c);
//				}
//				basesRead++;
//			    }
//			}
//
//			line = chromosomeSequenceFile.readLine();
//		    }
//		    chromosomeSequenceFile.close();
//		    return seq.toString();
//		}
//	    };

//	    pool.submit(callable);



	    snpsCtr++;
	    ArrayList<String> probes = eQTLProbes.get(esnp);
	    for (String p : probes) {
		String callablesnp = esnp;
		outcommands.writeln("/home/users/westra/MetaQC/shrimp/SHRiMP_2_2_2/bin/gmapper " + outdir + "/probeseq/" + p + ".fa " + outdir + "/snpseq/" + callablesnp + ".fa -N 4 -m 10 -i 0 -q -250 -f -100 -h 30% --shrimp-format >" + outdir + "/alignments/" + callablesnp + "-" + p + ".out 2>" + outdir + "/alignments/" + callablesnp + "-" + p + ".log &");

		ctr++;
		if (ctr % 16 == 0) {
		    outcommands.writeln("wait");
		}
	    }


	    System.out.println(snpsCtr);

	}

	outcommands.close();
//	threadPool.shutdown();


    }

    private String getPositionFromGenome(Integer chr, Integer left, Integer right) throws IOException {
	TextFile chromosomeSequenceFile = new TextFile("/home/users/westra/MetaQC/genomesequence/Homo_sapiens.NCBI36.54.dna.chromosome." + chr + ".fa.gz", TextFile.R);

	String line = chromosomeSequenceFile.readLine();
	int basesRead = 0;

	StringBuilder seq = new StringBuilder();
	while (line != null) {
	    if (line.startsWith(">")) {
		// fasta header line
	    } else {


		for (int i = 0; i < line.length(); i++) {

		    char c = line.charAt(i);
		    if (c == 'A' || c == 'T' || c == 'G' || c == 'C' || c == 'N' || c == 'M' || c == 'R' ) {
			if (basesRead >= left && basesRead < right) {
			    seq.append(c);
			} else if (basesRead >= left && basesRead >= right) {
			    break;
			}

		    } else {
			System.err.println("Could not parse char: " + c);
		    }
		    basesRead++;
		}
	    }

	    line = chromosomeSequenceFile.readLine();
	}
	chromosomeSequenceFile.close();
	return seq.toString();
    }

    private void loadProbeSequences(String probefile, String directoryToWriteTo) throws IOException {
	System.out.println("Loading probe sequences from " + probefile);
	probeSequences = new HashMap<String, String>();
	TextFile tf = new TextFile(probefile, TextFile.R);

	String[] elems = tf.readLineElemsReturnObjects(TextFile.tab);

	ArrayList<String> probeIds = new ArrayList<String>();
	while (elems != null) {
	    probeIds.add(elems[0]);
	    probeSequences.put(elems[0], elems[1]);
	    elems = tf.readLineElemsReturnObjects(TextFile.tab);
	}
	System.out.println(probeSequences.size() + " probes read");
	tf.close();

	for (int i = 0; i < probeIds.size(); i++) {
	    String probe = probeIds.get(i);
	    TextFile fasta = new TextFile(directoryToWriteTo + "/" + probe + ".fa", TextFile.W);

	    fasta.writeln(">" + probe);
	    fasta.writeln(probeSequences.get(probe));
	    fasta.close();
	}

    }
}
