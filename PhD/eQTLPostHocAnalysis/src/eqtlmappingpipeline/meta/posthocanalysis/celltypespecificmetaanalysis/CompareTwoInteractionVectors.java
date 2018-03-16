/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class CompareTwoInteractionVectors {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            CompareTwoInteractionVectors v = new CompareTwoInteractionVectors();
            DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-MetaAnalysis/MetaAnalysisZScoreMatrix.txt");

            Integer rowIdForInteractionTerm = metaMatrix.hashRows.get("CellTypeInteractionZScore");
            Integer rowIdForMainEffect = metaMatrix.hashRows.get("MainEffectZScore");

            ArrayList<Double> valsX = new ArrayList<Double>();
            ArrayList<Double> valsY = new ArrayList<Double>();
            for (int col = 0; col < metaMatrix.nrCols; col++) {
                valsX.add(metaMatrix.rawData[rowIdForInteractionTerm][col]);
                valsY.add(metaMatrix.rawData[rowIdForMainEffect][col]);
            }

            double[] xarr = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
            double[] yarr = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));

            double r = JSci.maths.ArrayMath.correlation(xarr, yarr);
            ScatterPlot plot = new ScatterPlot(500, 500, xarr, yarr, ScatterPlot.OUTPUTFORMAT.PDF, "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-MetaAnalysis/ComparisonBetweenInteractionAndMainEffect.pdf");
            System.out.println(r + "\t" + (r * r));

//            metaMatrix.save("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutputUsingRealNeutrophilCounts/CellTypeSpecificityMatrix.txt");
//            Integer rowIdForInteractionTerm = metaMatrix.hashRows.get("CellTypeInteractionZScore");
//            double[] interaction = metaMatrix.getRawData()[rowIdForInteractionTerm];
//            String[] matrixes = new String[]{
//                "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-EGCUT/2013-07-08-LymphocytesAsCovariate/CellTypeSpecificityMatrix.binary",
//                "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-EGCUT/2013-07-08-NeutrophilsAsCovariate/CellTypeSpecificityMatrix.binary"};
//            String output = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-EGCUT/ComparisonBetweenLymphocytesAndNeutrophils.pdf";
//            v.run(matrixes, output);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String[] matrixes, String output) throws IOException {
        System.out.println("Running..");
        ArrayList<double[]> vectors = new ArrayList<double[]>();
        ArrayList<List<String>> eQTLColNames = new ArrayList<List<String>>();
        ArrayList<Map<String, Integer>> eQTLIndexes = new ArrayList<Map<String, Integer>>();

        for (int m = 0; m < matrixes.length; m++) {
            DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(matrixes[m]);
            Integer rowIdForInteractionTerm = metaMatrix.hashRows.get("CellTypeInteractionZScore");
            Integer rowIdForMainEffect = metaMatrix.hashRows.get("MainEffectBeta");
            if (rowIdForMainEffect == null) {
                rowIdForMainEffect = metaMatrix.hashRows.get("MainEffectZScore");
            }

            double[] interaction = metaMatrix.getRawData()[rowIdForInteractionTerm];
            double[] beta = metaMatrix.getRawData()[rowIdForMainEffect];
            for (int i = 0; i < beta.length; i++) {
                if (beta[i] < 0) {
                    interaction[i] *= -1;
                }
            }

            eQTLColNames.add(metaMatrix.colObjects);
            vectors.add(interaction);
            eQTLIndexes.add(metaMatrix.hashCols);

        }

        for (int ds1Nr = 0; ds1Nr < matrixes.length; ds1Nr++) {
            double[] vector1 = vectors.get(ds1Nr);
            List<String> colNames = eQTLColNames.get(ds1Nr);
            for (int ds2Nr = ds1Nr + 1; ds2Nr < matrixes.length; ds2Nr++) {
                double[] vector2 = vectors.get(ds2Nr);

                ArrayList<Double> valsX = new ArrayList<Double>();
                ArrayList<Double> valsY = new ArrayList<Double>();

                Map<String, Integer> correspondingEQTL = eQTLIndexes.get(ds2Nr);
                for (int colNr = 0; colNr < vector1.length; colNr++) {
                    String eQTL = colNames.get(colNr);
                    Integer id = correspondingEQTL.get(eQTL);
                    if (id != null) {
                        valsX.add(vector1[colNr]);
                        valsY.add(vector2[id]);

                        System.out.println(eQTL + "\t" + vector1[colNr] + "\t" + vector2[id]);
                    }
                }

                double[] xarr = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                double[] yarr = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));

                double r = JSci.maths.ArrayMath.correlation(xarr, yarr);
                ScatterPlot plot = new ScatterPlot(500, 500, xarr, yarr, ScatterPlot.OUTPUTFORMAT.PDF, output);
                System.out.println(r + "\t" + (r * r));
            }
        }


    }
}
