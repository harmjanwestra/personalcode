import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class GTFToHGNC {

    public static void main(String[] args) {
        String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz";
        try {
            GTFAnnotation g = new GTFAnnotation(gtf);
            TextFile tf = new TextFile("d:\\tmp\\ENSGToHGNC-GencodeV32.txt", TextFile.W);
            for (Gene gene : g.getGenes()) {
                String[] geneelems = gene.getName().split("\\.");
                String genename = geneelems[0];
                tf.writeln(genename + "\t" + gene.getGeneSymbol());
            }
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
