package nl.harmjanwestra.playground;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.enums.Strand;
import umcg.genetica.features.Exon;
import umcg.genetica.features.Gene;
import umcg.genetica.features.Transcript;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

public class FilterGTF {

    public static void main(String[] args) {
        String gtf = "U:\\2020-05-MetaBrainExonquant\\gencode.v32.primary_assembly.annotation.gtf.gz";
        String out = "U:\\2020-05-MetaBrainExonquant\\probeannotations\\ProbeAnnotation-gencode.v32.primary_assembly";
        FilterGTF g = new FilterGTF();
        try {
            g.run(gtf, out, "gencode.v32");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void run(String gtffile, String outf, String platformstr) throws IOException {


        GTFAnnotation gtf = new GTFAnnotation(gtffile);

        TextFile outputt = new TextFile(outf + "-transcripts.txt.gz", TextFile.W);
        TextFile outputt2 = new TextFile(outf + "-transcripts2.txt.gz", TextFile.W);
        TextFile outputg = new TextFile(outf + "-genes.txt.gz", TextFile.W);
        TextFile outpute = new TextFile(outf + "-exons.txt.gz", TextFile.W);
        TextFile outpute2 = new TextFile(outf + "-exons2.txt.gz", TextFile.W);
        String header = "Platform\tArrayAddress\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tStrand";
        outputt.writeln(header);
        outputt2.writeln(header);
        outputg.writeln(header);
        outpute.writeln(header);
        outpute2.writeln(header);

        Collection<Gene> genes = gtf.getGenes();
        for (Gene g : genes) {
            Strand strand = g.getStrand();
            int tss = g.getStart();
            if (strand.equals(Strand.NEG)) {
                tss = g.getStop();
            }
            outputg.writeln(platformstr + "\t" + g.getName() + "\t" + g.getGeneSymbol() + "\t" + g.getChromosome().getNumber() + "\t" + tss + "\t" + tss + "\t" + g.getName() + "\t" + strand);
            ArrayList<Transcript> transcripts = g.getTranscripts();
            for (Transcript t : transcripts) {
                tss = t.getStart();
                if (strand.equals(Strand.NEG)) {
                    tss = t.getStop();
                }
                outputt.writeln(platformstr + "\t" + t.getName() + "\t" + g.getGeneSymbol() + "\t" + t.getChromosome().getNumber() + "\t" + tss + "\t" + tss + "\t" + g.getName() + "\t" + strand);
                String transcriptid = t.getName() + "_chr" + t.getChromosome().getNumber() + "_" + t.getStart() + "_" + t.getStop();
                outputt2.writeln(platformstr + "\t" + t.getName() + "\t" + g.getGeneSymbol() + "\t" + t.getChromosome().getNumber() + "\t" + tss + "\t" + tss + "\t" + g.getName() + "\t" + strand);
                ArrayList<Exon> exons = t.getExons();
                for (Exon e : exons) {
                    tss = e.getStart();
                    if (strand.equals(Strand.NEG)) {
                        tss = e.getStop();
                    }
                    outpute.writeln(platformstr + "\t" + e.getName() + "\t" + g.getGeneSymbol() + "\t" + e.getChromosome().getNumber() + "\t" + tss + "\t" + tss + "\t" + g.getName() + "\t" + strand);
                    String exonid = e.getName() + "_chr" + e.getChromosome().getNumber() + "_" + e.getStart() + "_" + e.getStop();
                    outpute2.writeln(platformstr + "\t" + exonid + "\t" + g.getGeneSymbol() + "\t" + e.getChromosome().getNumber() + "\t" + tss + "\t" + tss + "\t" + g.getName() + "\t" + strand);
                }
            }
        }
        outputt.close();
        outputt2.close();
        outputg.close();
        outpute.close();
        outpute2.close();
    }
}
