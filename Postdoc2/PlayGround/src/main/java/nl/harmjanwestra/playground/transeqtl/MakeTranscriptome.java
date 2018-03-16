package nl.harmjanwestra.playground.transeqtl;

import nl.harmjanwestra.utilities.annotation.ensembl.EnsemblStructures;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.features.Exon;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.features.Transcript;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class MakeTranscriptome {
	
	public static void main(String[] args) {
		
		String eqtlfile = "";
		String genomefasta = "";
		String ensemblannotation = "";
		String snpoutputdir = "";
		int windowsize = 2000000;
		int readlen = 25;
		int shift = 2;
		
	}
	
	public void run(String eqtlfile, String genomeFasta, String ensemblannotation, String geneoutputdir, String snpoutputdir, int windowsize, int readlen, int shift) throws IOException {
		
		// create reference genome from snps
		TextFile tf = new TextFile(eqtlfile, TextFile.R);
		tf.readLine();
		
		
		ArrayList<Triple<String, Chromosome, Integer>> snps = new ArrayList<>();
		HashSet<String> snpids = new HashSet<String>();
		HashSet<String> geneIds = new HashSet<String>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			Chromosome chr = Chromosome.parseChr(elems[1]);
			Integer snppos = Integer.parseInt(elems[2]);
			String snpname = elems[3];
			String genename = elems[4];
			if (!snpids.contains(snpname)) {
				Triple<String, Chromosome, Integer> snp = new Triple<>(snpname, chr, snppos);
				snps.add(snp);
				snpids.add(snpname);
			}
			
			geneIds.add(genename);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		EnsemblStructures s = new EnsemblStructures(ensemblannotation);
		
		// convert gene strings to ensembl gene ids
		HashMap<String, Gene> geneHash = s.getStrToGene();
		ArrayList<Gene> genesToExport = new ArrayList<>();
		for (String g : geneIds) {
			Gene gobj = geneHash.get(g);
			if (gobj != null) {
				genesToExport.add(gobj);
			}
		}
		
		// create 1mb regions around snps
		htsjdk.samtools.reference.FastaSequenceFile fastaFile = new htsjdk.samtools.reference.FastaSequenceFile(new File(genomeFasta), false);
		htsjdk.samtools.reference.ReferenceSequence seq = fastaFile.nextSequence();
		
		while (seq != null) {
			System.out.println(seq.getName() + "\t" + seq.length());
			byte[] bases = seq.getBases();
			
			Chromosome seqchr = Chromosome.parseChr(seq.getName());
			
			for (Triple<String, Chromosome, Integer> snp : snps) {
				if (snp.getMiddle().equals(seqchr)) {
					TextFile snpfastaout = new TextFile(snpoutputdir + "snpseqs.fa.gz", TextFile.W);
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
					byte[] snpwindow = new byte[nrbases];
					System.arraycopy(bases, windowleft, snpwindow, 0, nrbases);
					
					// make into fasta
					String fastaStr = ">" + snp.getMiddle().toString() + "_" + snp.getRight() + "_" + snp.getLeft() + "_" + windowleft + "_" + windowright
							+ "\n" + byteToStr(snpwindow);
					snpfastaout.writeln(fastaStr);
					snpfastaout.close();
				}
			}
			
			// now parse the ensemblgenes on this chr
			for (Gene g : genesToExport) {
				if (g.getChromosome().equals(seqchr)) {
					TextFile genesfastaout = new TextFile(geneoutputdir + g.getName() + ".fa.gz", TextFile.W);
					ArrayList<Transcript> transcripts = g.getTranscripts();
					// split up the transcript in readlen reads
					for (Transcript t : transcripts) {
						Exon[] exons = t.getExonsRanked();
						// output exon sequences
						int exonctr = 0;
						int tstart = exons[0].getStart();
						int tstop = exons[exons.length - 1].getStop();
						int currentPos = tstart;
						Exon currentexon = exons[0];
						while (currentPos < tstop) {
							
							// iterate current exon
							while (currentPos + readlen < currentexon.getStop()) {
								byte[] b = new byte[readlen];
								System.arraycopy(bases, currentPos, b, 0, readlen);
								String fastaStr = ">" + g.getName() + "_" + t.getName() + "_" + currentexon.getName() + "_" + currentPos + "_" + (currentPos + readlen)
										+ "\n" + byteToStr(b);
								genesfastaout.writeln(fastaStr);
								currentPos += shift;
							}
							
							// process exon exon boundary
							
							// if there are exons remaining
							if (exonctr + 1 < exons.length) {
								while (currentPos < currentexon.getStop()) { // iterate until we run out of currentexon
									int basesFromCurrentExon = currentexon.getStop() - currentPos;
									byte[] b = new byte[readlen];
									System.arraycopy(bases, currentPos, b, 0, basesFromCurrentExon);
									
									int basesRemaining = readlen - basesFromCurrentExon;
									// now we should check whether the next exon actually has enough bases
									int tmpexonctr = exonctr;
									while (basesRemaining > 0 && tmpexonctr < exons.length - 1) {
										// get remaining bases from next exons
										Exon tmpexon = exons[tmpexonctr];
										int tmpexonstart = tmpexon.getStart();
										int tmpexonstop = tmpexon.getStop();
										
										int basesToTakeFromExon = basesRemaining;
										if (tmpexonstart + basesRemaining > tmpexonstop) {
											basesToTakeFromExon = tmpexonstop - tmpexonstart;
										}
										
										System.arraycopy(bases, tmpexonstart, b, basesFromCurrentExon, basesToTakeFromExon);
										basesRemaining -= basesToTakeFromExon;
										
										tmpexonctr++;
									}
									
									currentPos += shift;
								}
							}
							
							// start at next exon
							currentexon = exons[exonctr];
							currentPos = currentexon.getStart();
							exonctr++;
						}
						
					}
					genesfastaout.close();
				}
			}
			seq = fastaFile.nextSequence();
		}
		
	}
	
	private String byteToStr(byte[] bytes) throws UnsupportedEncodingException {
		return new String(bytes, "UTF-8");
	}
}
