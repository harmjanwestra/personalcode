package nl.harmjanwestra.playground;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Exon;
import umcg.genetica.features.FeatureComparator;
import umcg.genetica.features.Gene;
import umcg.genetica.features.Transcript;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

public class ListExons {
    public static void main(String[] args) {
        try {
            String query = "ENSG00000138398.16";

            GTFAnnotation gtf = new GTFAnnotation("D:\\Sync\\SyncThing\\Data\\Ref\\gencode\\gencode.v32.primary_assembly.annotation.gtf.gz");
            for (Gene g : gtf.getGenes()) {
                if (g.getChromosome().isAutosome() || g.getChromosome().equals(Chromosome.X) || g.getChromosome().equals(Chromosome.Y) || g.getChromosome().equals(Chromosome.MT)) {
                    if (g.getName().equals(query)) {
                        ArrayList<Transcript> transcripts = g.getTranscripts();
                        for (Transcript t : transcripts) {
                            ArrayList<Exon> exons = t.getExons();

                            for (Exon e : exons) {
                                System.out.println(g.getName() + "\t" + t.getName() + "\t" + e.getName() + "\t" + e.getStart() + "\t" + e.getStop() + "\t" + e.getStrand() + "\t" + t.getExonRank(e));
                            }
                        }
                    }
                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
