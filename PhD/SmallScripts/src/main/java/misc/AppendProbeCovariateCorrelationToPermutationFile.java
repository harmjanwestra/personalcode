/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class AppendProbeCovariateCorrelationToPermutationFile {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String permutationfile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/CisEffectsAndArrayQuality/2013-11-07-BloodHT12v3-WithCovariate-Neutroproxy-maf0.05-Parametric-YHatSpearman-MultiThread/PermutedEQTLsPermutationRound1.txt.gz";
        String expressionfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt";
        String covariatefile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-Groningen/Normalization/CellTypeProxyFile.txt";
        try {
            DoubleMatrixDataset<String, String> dscov = new DoubleMatrixDataset<String, String>(covariatefile);
            dscov.transposeDataset();
            DoubleMatrixDataset<String, String> dsexp = new DoubleMatrixDataset<String, String>(expressionfile);

            HashMap<String, Double> correlations = new HashMap<String, Double>();

            for (int row = 0; row < dsexp.nrRows; row++) {
                ArrayList<Double> x = new ArrayList<Double>();
                ArrayList<Double> y = new ArrayList<Double>();
                for (int sample = 0; sample < dsexp.nrCols; sample++) {
                    String sampleName = dsexp.colObjects.get(sample);
                    Integer sampleId = dscov.hashCols.get(sampleName);
                    if (sampleId != null) {
                        x.add(dsexp.rawData[row][sample]);
                        y.add(dscov.rawData[0][sampleId]);
                    }
                }

                double[] xarr = Primitives.toPrimitiveArr(x.toArray(new Double[0]));
                double[] yarr = Primitives.toPrimitiveArr(y.toArray(new Double[0]));

                double corr = JSci.maths.ArrayMath.correlation(xarr, yarr);
                correlations.put(dsexp.rowObjects.get(row), corr);
            }

            TextFile tfIn = new TextFile(permutationfile, TextFile.R);
            TextFile tfOut = new TextFile(permutationfile + "-CorrWCovariate.txt", TextFile.W);
            tfOut.writeln(tfIn.readLine() + "\tcorr");
            String[] elems = tfIn.readLineElems(TextFile.tab);

            while (elems != null) {
                tfOut.writeln(Strings.concat(elems, Strings.tab) + "\t" + correlations.get(elems[2]));
                elems = tfIn.readLineElems(TextFile.tab);
            }
            tfIn.close();
            tfOut.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
