/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package correlator2;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class Correlator2 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        if (args.length > 1) {
            String ds = args[0];
            String cov = args[1];

            String output = ds + "-CorrelationWithCovariates.txt";
            try {
                Correlator2 covc = new Correlator2(ds, cov, output);
            } catch (IOException e) {
                e.printStackTrace();
            }
        } else {
            System.out.println("Usage: correlate dataset1 dataset2");
        }
    }

    private Correlator2(String dsloc, String covloc, String output) throws IOException {
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(dsloc);
        DoubleMatrixDataset<String, String> cov = new DoubleMatrixDataset<String, String>(covloc);
        cov.transposeDataset();

        int[] toSample = new int[ds.nrCols];
        int nrSharedSamples = 0;
        for (int col = 0; col < ds.nrCols; col++) {
            Integer id = cov.hashCols.get(ds.colObjects.get(col));
            if (id != null) {
                toSample[col] = id;
                nrSharedSamples++;
            } else {
                toSample[col] = -1;
            }
        }

        if (nrSharedSamples == 0) {
            ds.transposeDataset();
        }
        
        toSample = new int[ds.nrCols];
        nrSharedSamples = 0;
        for (int col = 0; col < ds.nrCols; col++) {
            Integer id = cov.hashCols.get(ds.colObjects.get(col));
            if (id != null) {
                toSample[col] = id;
                nrSharedSamples++;
            } else {
                toSample[col] = -1;
            }
        }
        if (nrSharedSamples == 0) {
            System.err.println("No shared samples!");
            System.exit(-1);
        }

        System.out.println("Shared samples: " + nrSharedSamples);
        TextFile outFile = new TextFile(output, TextFile.W);
        for (int row = 0; row < ds.nrRows; row++) {
            for (int covrow = 0; covrow < cov.nrRows; covrow++) {

                double[] valsX = new double[nrSharedSamples];
                double[] valsY = new double[nrSharedSamples];
                int ctr = 0;
                for (int col = 0; col < ds.nrCols; col++) {
                    int sample2Id = toSample[col];
                    if (sample2Id != -1) {
                        valsX[ctr] = ds.rawData[row][col];
                        valsY[ctr] = cov.rawData[covrow][sample2Id];
                        ctr++;
                    }
                }
                double corr = JSci.maths.ArrayMath.correlation(valsX, valsY);
                outFile.writeln(ds.rowObjects.get(row) + "\t" + cov.rowObjects.get(covrow) + "\t" + corr);
            }
        }
        outFile.close();
    }
}
