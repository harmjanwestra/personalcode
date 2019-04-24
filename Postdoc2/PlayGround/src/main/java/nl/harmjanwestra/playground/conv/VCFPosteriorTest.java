package nl.harmjanwestra.playground.conv;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class VCFPosteriorTest {
    public static void main(String[] args) {
        String vcf = "D:\\ENA\\ENA_chr21.gg-posterior-b37.vcf.gz";

        TextFile tf = null;
        try {
            tf = new TextFile(vcf, TextFile.R);

            String ln = tf.readLine();
            while (ln != null) {

                if (!ln.startsWith("#")) {
                    VCFVariant var = new VCFVariant(ln, VCFVariant.PARSE.ALL);

                    DoubleMatrix2D genotypes = var.getGenotypeAllelesAsMatrix2D();
                    DoubleMatrix2D probs = var.getPosteriorProbabilities();

                    DoubleMatrix2D dosages = var.calculateDosageFromProbabilities(probs);
                    for (int i = 0; i < var.getNrSamples(); i++) {
                        double a1 = genotypes.getQuick(i, 0);
                        double a2 = genotypes.getQuick(i, 1);

                        double probAA = probs.getQuick(i, 0);
                        double probAB = probs.getQuick(i, 1);
                        double probBB = probs.getQuick(i, 2);
                        double dosAA = dosages.getQuick(i, 0);


                        double probsum = probAA + (probAB) + (probBB);

                        System.out.println(i + "\t" + a1 + "/" + a2 + "\t" + dosAA + "\t" + probAA + "\t" + probAB + "\t" + probBB + "\t" + probsum);
                    }
                    System.exit(-1);
                }
                ln = tf.readLine();
            }
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }
}
