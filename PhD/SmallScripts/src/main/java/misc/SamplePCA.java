/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import eqtlmappingpipeline.normalization.Normalizer;
import java.io.IOException;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.concurrent.ConcurrentCorrelation;

/**
 *
 * @author harmjan
 */
public class SamplePCA {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String in = "";
        String out= "";
        try{
            SamplePCA pca = new SamplePCA();
            pca.run(in, out);
        } catch (IOException e){
            
        }
    }

    public void run(String in, String out) throws IOException {
        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(in);
        eqtlmappingpipeline.normalization.Normalizer normalizer = new Normalizer();

        ConcurrentCorrelation c = new ConcurrentCorrelation(2);
        double[][] correlationMatrix = c.pairwiseCorrelation(dataset.getRawDataTransposed());

        //outputFileNamePrefix = outdir + expressionFileName;
        Pair<DoubleMatrixDataset<String, String>, DoubleMatrixDataset<String, String>> PCAResults = normalizer.calculatePCA(dataset, correlationMatrix, out, correlationMatrix.length-1);
        
    }
}
