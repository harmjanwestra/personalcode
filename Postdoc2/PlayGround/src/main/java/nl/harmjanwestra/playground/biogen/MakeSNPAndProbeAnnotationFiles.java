package nl.harmjanwestra.playground.biogen;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;
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
//			p.gtfToProbeAnnotationFile(gtfin, platformstr, out);

            out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\GeneAnnotationFile-Gencodev24-Ensembl83-b38-transcripts.txt.gz";
//			p.gtfToProbeAnnotationTranscriptFile(gtfin, platformstr, out);
            out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\GeneAnnotationFile-Gencodev24-Ensembl83-b38-exons.txt.gz";
//			p.gtfToProbeAnnotationExonFile(gtfin, platformstr, out);

//			p.gtfToProbeAnnotationFile("D:\\Sync\\SyncThing\\Data\\Ref\\GTEX\\gencode.v19.genes.v7.patched_contigs.gtf.gz",
//					"gencode.v19.genes.v7",
//					"D:\\Sync\\SyncThing\\Data\\Ref\\GTEX\\GTEx-gencode.v19.genes.v7.patched_contigs.ProbeAnnotation.txt.gz");

//            p.gtfToProbeAnnotationFile("D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz",
//                    "gencode.v32",
//                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz",
//                    true);
//            p.writeListOfProteinCodingGenes("D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz",
//                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.proteincoding.txt.gz");

//            p.gtfToProbeAnnotationFile("D:\\tmp\\Homo_sapiens.GRCh37.75.gtf",
//                    "ensemble.75",
//                    "D:\\tmp\\Homo_sapiens.GRCh37.75.ProbeAnnotation.txt.gz", true);

//            p.gtfToProbeAnnotationFile("D:\\tmp\\genes.gtf",
//                    "refdata-gex-GRCh38-2020-A",
//                    "D:\\tmp\\refdata-gex-GRCh38-2020-A-GeneAnnotation.txt.gz", true);
//            p.writeListOfProteinCodingGenes("D:\\tmp\\genes.gtf",
//                    "D:\\tmp\\refdata-gex-GRCh38-2020-A-proteincoding.txt");
//                    "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.proteincoding.txt.gz");

            p.gtfToProbeAnnotationFile("/Users/harm-jan/SyncData/TMP/orfeas/transtest/data/gencode.v45.primary_assembly.annotation-copy.gtf.gz","GencodeV45PrimaryAssembly", "/Users/harm-jan/SyncData/TMP/orfeas/transtest/data/ProbeAnnotation-gencode.v45.primary_assembly.annotation.txt.gz",true, true);
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

    private void writeListOfProteinCodingGenes(String in, String out) throws IOException {
        GTFAnnotation gtfAnnotation = new GTFAnnotation(in);

        Collection<Gene> genes = gtfAnnotation.getGenes();
        HashSet<String> proteincoding = new HashSet<String>();
        for (Gene g : genes) {
            String type = g.getType();
            if (type != null && type.contains("protein")) {
//                proteincoding.add(g.getName());
                proteincoding.add(g.getGeneSymbol());
            }
        }

        TextFile outf = new TextFile(out, TextFile.W);
        for (String s : proteincoding) {
            outf.writeln(s);
        }
        outf.close();
    }

    public void gtfToProbeAnnotationFile(String in, String platformstr, String out, boolean printTSS, boolean stripGeneVersion) throws IOException {
        GTFAnnotation gtfAnnotation = new GTFAnnotation(in);

        Collection<Gene> genes = gtfAnnotation.getGenes();
        TextFile outf = new TextFile(out, TextFile.W);
        outf.writeln("Platform\tArrayAddress\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tStrand");
        for (Gene g : genes) {

            if (printTSS) {
                int tss = 0;
                if (g.getStrand().equals(Strand.NEG)) {
                    tss = g.getStop();
                } else {
                    tss = g.getStart();
                }
                String genename = g.getName();
                if(stripGeneVersion){
                    genename = genename.split("\\.")[0];
                }
                String outln = platformstr + "\t" + genename + "\t" + g.getGeneSymbol() + "\t" + g.getChromosome().getNumber() + "\t" + tss + "\t" + tss + "\t" + g.getName() + "\t" + g.getStrand();
                if (g.getChromosome().equals(Chromosome.X)) {
                    outln = platformstr + "\t" + genename + "\t" + g.getGeneSymbol() + "\tX\t" + tss + "\t" + tss + "\t" + g.getName() + "\t" + g.getStrand();
                } else if (g.getChromosome().equals(Chromosome.Y)) {
                    outln = platformstr + "\t" + genename + "\t" + g.getGeneSymbol() + "\tY\t" + tss + "\t" + tss + "\t" + g.getName() + "\t" + g.getStrand();
                }
                outf.writeln(outln);
            } else {
                String genename = g.getName();
                if(stripGeneVersion){
                    genename = genename.split("\\.")[0];
                }
                String outln = platformstr + "\t" + genename + "\t" + g.getGeneSymbol() + "\t" + g.getChromosome().getNumber() + "\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
                if (g.getChromosome().equals(Chromosome.X)) {
                    outln = platformstr + "\t" + genename + "\t" + g.getGeneSymbol() + "\tX\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
                } else if (g.getChromosome().equals(Chromosome.Y)) {
                    outln = platformstr + "\t" + genename + "\t" + g.getGeneSymbol() + "\tY\t" + g.getStart() + "\t" + g.getStop() + "\t" + g.getName() + "\tN";
                }
                outf.writeln(outln);
            }

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
