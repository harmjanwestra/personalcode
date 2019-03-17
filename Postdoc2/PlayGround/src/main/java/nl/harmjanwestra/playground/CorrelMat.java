package nl.harmjanwestra.playground;

import nl.harmjanwestra.utilities.annotation.gtf.GTFAnnotation;
import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.enums.Strand;
import nl.harmjanwestra.utilities.features.Gene;
import nl.harmjanwestra.utilities.legacy.genetica.console.ProgressBar;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;

import java.io.IOException;

public class CorrelMat {

    public static void main(String[] args) {

        String mat = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\mergednorm\\TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.10PCAsOverSamplesRemoved.txt.gz";
        String query = "ENSG00000147894.14";

        String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-03-06-ENA\\c9orfgenes.txt";

        String annot = "D:\\Sync\\SyncThing\\Data\\Ref\\gencode\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";

        try {
            GTFAnnotation g = new GTFAnnotation(annot);
            TextFile tf = new TextFile(out, TextFile.W);
            tf.writeln("Ensembl\tPosition(b38)\tGeneSymbol\tSpearmanCorrelationCoefficient\tZ\tP-val");

            DoubleMatrixDataset<String, String> data = DoubleMatrixDataset.loadDoubleData(mat);
            Integer id = data.getRowIndex(query);
            if (id > -1) {

                ProgressBar pb = new ProgressBar(data.rows());
                for (int i = 0; i < data.rows(); i++) {
//                    if (i != id) {
                    String gname = data.getRowObjects().get(i);

                    SpearmansCorrelation c = new SpearmansCorrelation();
                    Gene gene = g.getStrToGene().get(gname);
                    String symbol = gname;
                    Chromosome chr = Chromosome.NA;
                    Strand strand = Strand.NA;
                    int start = 0;
                    int stop = 0;
                    if (gene != null) {
                        symbol = gene.getGeneSymbol();
                        chr = gene.getChromosome();
                        start = gene.getStart();
                        stop = gene.getStop();
                        strand = gene.getStrand();
                    }
                    double corr = c.correlation(data.getRow(id).toArray(), data.getRow(i).toArray());

                    Correlation.correlationToZScore(data.columns());
                    double z = Correlation.convertCorrelationToZScore(data.columns(), corr);
                    double p = ZScores.zToP(z);

                    tf.writeln(gname + "\t" + chr.toString() + ":" + start + "-" + stop + ":" + strand + "\t" + symbol + "\t" + corr + "\t" + z + "\t" + p);
//                    }
                    pb.iterate();
                }
                pb.close();
            }
            tf.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
