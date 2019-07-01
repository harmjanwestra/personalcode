package nl.harmjanwestra.playground.biogen;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Exon;
import umcg.genetica.features.Gene;
import umcg.genetica.features.Transcript;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

public class MakeSNPAndProbeAnnotationFiles {

	public static void main(String[] args) {
		String gtfin = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\annotation\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\GeneAnnotationFile-Gencodev24-Ensembl83-b38.txt.gz";
		String platformstr = "Gencodev24-Ensembl83-b38";

		MakeSNPAndProbeAnnotationFiles p = new MakeSNPAndProbeAnnotationFiles();
		try {
			p.gtfToProbeAnnotationFile(gtfin, platformstr, out);

			out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\GeneAnnotationFile-Gencodev24-Ensembl83-b38-transcripts.txt.gz";
			p.gtfToProbeAnnotationTranscriptFile(gtfin, platformstr, out);
			out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\GeneAnnotationFile-Gencodev24-Ensembl83-b38-exons.txt.gz";
			p.gtfToProbeAnnotationExonFile(gtfin, platformstr, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
//
//		String dbsnpvcf = args[0];
//		String output = args[1];
//		try {
//			p.dbsnpVCFToSNPMappings(dbsnpvcf, output);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

	}

	public void gtfToProbeAnnotationFile(String in, String platformstr, String out) throws IOException {
		GTFAnnotation gtfAnnotation = new GTFAnnotation(in);

		Collection<Gene> genes = gtfAnnotation.getGenes();
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("Platform\tArrayAddress\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tSeq");
		for (Gene g : genes) {
			String outln = platformstr + "\t" + g.getName() + "\t" + g.getGeneSymbol() + "\t" + g.getChromosome().getNumber() + "\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
			if (g.getChromosome().equals(Chromosome.X)) {
				outln = platformstr + "\t" + g.getName() + "\t" + g.getGeneSymbol() + "\tX\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
			} else if (g.getChromosome().equals(Chromosome.Y)) {
				outln = platformstr + "\t" + g.getName() + "\t" + g.getGeneSymbol() + "\tY\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
			}

			outf.writeln(outln);
		}
		outf.close();
	}

	public void gtfToProbeAnnotationTranscriptFile(String in, String platformstr, String out) throws IOException {
		GTFAnnotation gtfAnnotation = new GTFAnnotation(in);


		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("Platform\tArrayAddress\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tSeq");
		for (Gene gene : gtfAnnotation.getGenes()) {
			for (Transcript g : gene.getTranscripts()) {
				String outln = platformstr + "\t" + g.getName() + "\t" + g.getGene().getName() + "\t" + g.getChromosome().getNumber() + "\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
				if (g.getChromosome().equals(Chromosome.X)) {
					outln = platformstr + "\t" + g.getName() + "\t" + g.getGene().getName() + "\tX\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
				} else if (g.getChromosome().equals(Chromosome.Y)) {
					outln = platformstr + "\t" + g.getName() + "\t" + g.getGene().getName() + "\tY\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
				}
				outf.writeln(outln);
			}


		}
		outf.close();
	}

	public void gtfToProbeAnnotationExonFile(String in, String platformstr, String out) throws IOException {
		GTFAnnotation gtfAnnotation = new GTFAnnotation(in);


		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln("Platform\tArrayAddress\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tSeq");
		for (Gene gene : gtfAnnotation.getGenes()) {
			for (Transcript t : gene.getTranscripts()) {
				for (Exon g : t.getExons()) {
					String outln = platformstr + "\t" + g.getName() + "\t" + gene.getName() + ";" + t.getName() + "\t" + g.getChromosome().getNumber() + "\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
					if (g.getChromosome().equals(Chromosome.X)) {
						outln = platformstr + "\t" + g.getName() + "\t" + gene.getName() + ";" + t.getName() + "\tX\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
					} else if (g.getChromosome().equals(Chromosome.Y)) {
						outln = platformstr + "\t" + g.getName() + "\t" + gene.getName() + ";" + t.getName() + "\tY\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
					}

					outf.writeln(outln);
				}
			}

		}
		outf.close();
	}

	public void dbsnpVCFToSNPMappings(String vcf, String out) throws IOException {
		TextFile outf = new TextFile(out, TextFile.W, 32 * 1024);
		TextFile tfin = new TextFile(vcf, TextFile.R, 320 * 1024);
		String ln = tfin.readLine();
		HashSet<String> dups = new HashSet<String>(1000000);
		int ctr = 0;
		String prevchr = null;
		while (ln != null) {
			if (!ln.startsWith("#")) {
				VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.HEADER);
				if (!var.isBiallelic() && !var.isIndel()) {
					String chr = var.getChr();
					if (prevchr == null || !chr.equals(prevchr)) {
						dups = new HashSet<>();
						System.out.println();
						System.out.println(chr);
					}
					prevchr = chr;
					if (Chromosome.parseChr(chr).getNumber() < 25) {
						Integer pos = var.getPos();
						String rsid = var.getId();
						if (!dups.contains(rsid)) {
							outf.writeln(chr + "\t" + pos + "\t" + rsid);
							dups.add(rsid);
						}
					}

				}
			}
			ln = tfin.readLine();
			ctr++;
			if (ctr % 10000 == 0) {
				System.out.print("\r" + ctr);
			}
		}

		tfin.close();
		outf.close();
	}
}
