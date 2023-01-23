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

public class GeneLength {

    public static void main(String[] args) {
        try {
            TextFile tf = new TextFile("D:\\Sync\\SyncThing\\Data\\Ref\\gencode\\gencode.v32.primary_assembly.annotation-genelengths.txt.gz", TextFile.W);
            tf.writeln("Gene\tSymbol\tNrTranscripts\tNrExons\tNrMergedExons\tTotalGeneLength\tMergedExonLength");
            GTFAnnotation gtf = new GTFAnnotation("D:\\Sync\\SyncThing\\Data\\Ref\\gencode\\gencode.v32.primary_assembly.annotation.gtf.gz");
            for (Gene g : gtf.getGenes()) {
                if (g.getChromosome().isAutosome() || g.getChromosome().equals(Chromosome.X) || g.getChromosome().equals(Chromosome.Y) || g.getChromosome().equals(Chromosome.MT) ) {
                    ArrayList<Transcript> transcripts = g.getTranscripts();
                    HashSet<Exon> allexons = new HashSet<Exon>();
                    HashSet<String> exonNames = new HashSet<>();
                    for (Transcript t : transcripts) {
                        if (t.getChromosome().isAutosome() || t.getChromosome().equals(Chromosome.X) || t.getChromosome().equals(Chromosome.Y) || t.getChromosome().equals(Chromosome.MT)) {
                            ArrayList<Exon> exons = t.getExons();
                            for (Exon e : exons) {
                                String ename = e.getName();
                                if (!exonNames.contains(ename)) {
                                    allexons.add(e);
                                    exonNames.add(ename);
                                }
                            }
                        }
                    }

                    // overlap exons
                    ArrayList<Exon> exonArr = new ArrayList<>();
                    exonArr.addAll(allexons);
                    FeatureComparator c = new FeatureComparator();
                    Collections.sort(exonArr, c);

                    HashSet<Exon> merged = new HashSet<>();
                    ArrayList<Exon> combinedExons = new ArrayList<>();
                    for (int i = 0; i < exonArr.size(); i++) {
                        Exon a = exonArr.get(i);
                        if (!merged.contains(a)) {
                            Exon comb = new Exon("", a.getChromosome(), a.getStrand(), a.getGene(), a.getStart(), a.getStop());
                            combinedExons.add(comb);
                            for (int j = i + 1; j < exonArr.size(); j++) {
                                Exon b = exonArr.get(j);
                                if (!merged.contains(b)) {
                                    if (b.overlaps(comb)) {
                                        if (b.getStart() < comb.getStart()) {
                                            comb.setStart(b.getStart());
                                        }
                                        if (b.getStop() > comb.getStop()) {
                                            comb.setStop(b.getStop());
                                        }
                                        merged.add(b);
                                    }
                                }
                            }
                        }
                    }
                    int nrbases = 0;
                    for (Exon e : combinedExons) {
                        int bases = e.getStop() - e.getStart();
                        nrbases += bases;
                    }

                    tf.writeln(g.getName() + "\t" + g.getGeneSymbol() + "\t" + transcripts.size() + "\t" + exonArr.size() + "\t" + combinedExons.size() + "\t" + (g.getStop() - g.getStart()) + "\t" + nrbases);
                    System.out.println(g.getName() + "\t" + g.getGeneSymbol() + "\t" + transcripts.size() + "\t" + exonArr.size() + "\t" + combinedExons.size() + "\t" + (g.getStop() - g.getStart()) + "\t" + nrbases);
                }

            }
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
