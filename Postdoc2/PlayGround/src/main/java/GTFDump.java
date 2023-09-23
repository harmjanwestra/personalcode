import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;
import umcg.genetica.features.Exon;
import umcg.genetica.features.Gene;
import umcg.genetica.features.Transcript;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class GTFDump {


    public static void main(String[] args) {
//		String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz";
        String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.gtf.gz";
        try {
            GTFAnnotation g = new GTFAnnotation(gtf);
            TextFile tf = new TextFile("d:\\tmp\\GencodeV32-positions.txt.gz", TextFile.W);
            tf.writeln("Gene\tSymbol\tStrand\tGeneStart\tGeneEnd\tTranscript\tTranscriptStart\tTranscriptEnd\tExon\tExonStart\tExonEnd\tExonRankInTranscript");
            for (Gene gene : g.getGenes()) {
                for (Transcript t : gene.getTranscripts()) {
                    for (Exon e : t.getExons()) {
                        tf.writeln(gene.getName() + "\t" + gene.getGeneSymbol() + "\t" + gene.getStrand() + "\t" + gene.getStart() + "\t" + gene.getStop() + "\t" + t.getName() + "\t" + t.getStart() + "\t" + t.getStop() + "\t" + e.getName() + "\t" + e.getStart() + "\t" + e.getStop() + "\t" + t.getExonRank(e));
                    }
                }
            }
            tf.close();
            Gene
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
