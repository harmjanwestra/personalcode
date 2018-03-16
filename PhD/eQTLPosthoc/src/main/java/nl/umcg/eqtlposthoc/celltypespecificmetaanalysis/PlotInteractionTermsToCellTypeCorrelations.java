/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class PlotInteractionTermsToCellTypeCorrelations {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String metamatrix = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/2013-06-26-MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
            String endomatrix = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/RawEndoPhenotypes/EGCUTEndophenotypesValidSamples-NoSex.txt-asLoadedByNormalizerEndophenotypeVsExpressionCorrelationMatrix.txt";
            String output = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/DannyCorrelationAgainstMetaanalysis/";
            PlotInteractionTermsToCellTypeCorrelations c = new PlotInteractionTermsToCellTypeCorrelations();
            c.run(metamatrix, endomatrix, output);
        } catch (IOException e) {
        }
    }

    public void run(String infile, String endophenomatrix, String output) throws IOException {
        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(infile);
        DoubleMatrixDataset<String, String> endophenoCorrelationMatrix = new DoubleMatrixDataset<String, String>(endophenomatrix);

        Integer rowIdForInteractionTerm = metaMatrix.hashRows.get("CellTypeInteractionZScore");
        Integer rowIdForMainEffect = metaMatrix.hashRows.get("MainEffectBeta");

        Integer colIdforNeutrophils = endophenoCorrelationMatrix.hashCols.get("Neutrophils");

        ArrayList<Double> interactionTerms = new ArrayList<Double>();
        ArrayList<Double> interactionTermsUncorrected = new ArrayList<Double>();
        ArrayList<Double> probeCorrelationsToEndophenotypes = new ArrayList<Double>();


//        double[] interactionTerms = new double[metaMatrix.nrCols];
//        double[] probeCorrelationsToEndophenotypes = new double[metaMatrix.nrCols];
        for (int c = 0; c < metaMatrix.nrCols; c++) {
            String probe = metaMatrix.colObjects.get(c).split("-")[1];
            Integer probeId = endophenoCorrelationMatrix.hashRows.get(probe);

            if (probeId != null) {

                double maineffect = metaMatrix.rawData[rowIdForMainEffect][c];
                double direction = 1;
                if (maineffect < 0) {
                    direction = -1;
                }

                interactionTermsUncorrected.add(metaMatrix.rawData[rowIdForInteractionTerm][c]);
                interactionTerms.add(direction * metaMatrix.rawData[rowIdForInteractionTerm][c]);
                probeCorrelationsToEndophenotypes.add(endophenoCorrelationMatrix.rawData[probeId][colIdforNeutrophils]);
            } else {
                System.out.println("PROBE IS NULL: " + probe);
            }


        }

        double r = JSci.maths.ArrayMath.correlation(toPrimitiveArr(interactionTerms.toArray(new Double[0])), toPrimitiveArr(probeCorrelationsToEndophenotypes.toArray(new Double[0])));
        System.out.println(r + "\t" + (r * r));

        double rUnCorr = JSci.maths.ArrayMath.correlation(toPrimitiveArr(interactionTermsUncorrected.toArray(new Double[0])), toPrimitiveArr(probeCorrelationsToEndophenotypes.toArray(new Double[0])));
        System.out.println(rUnCorr + "\t" + (rUnCorr * rUnCorr));

        ScatterPlot plot = new ScatterPlot(500, 500, toPrimitiveArr(probeCorrelationsToEndophenotypes.toArray(new Double[0])), toPrimitiveArr(interactionTerms.toArray(new Double[0])), ScatterPlot.OUTPUTFORMAT.PDF, output + "EndophenotypeCorrelationVsCellTypeInteractionZScore.pdf");
        ScatterPlot plot2 = new ScatterPlot(500, 500, toPrimitiveArr(probeCorrelationsToEndophenotypes.toArray(new Double[0])), toPrimitiveArr(interactionTermsUncorrected.toArray(new Double[0])), ScatterPlot.OUTPUTFORMAT.PDF, output + "EndophenotypeCorrelationVsCellTypeInteractionZScoreUncorrected.pdf");


    }

    private double[] toPrimitiveArr(Double[] toArray) {
        double[] arr = new double[toArray.length];
        for (int i = 0; i < toArray.length; i++) {
            arr[i] = toArray[i];
        }
        return arr;
    }
}
