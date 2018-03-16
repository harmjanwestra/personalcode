/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import com.lowagie.text.DocumentException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import umcg.genetica.graphics.Heatmap;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class CorrelateCellCountsWithGeneExpression {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String expfile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-EGCUT/Normalization/ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.txt.gz";
            String cellfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/RawEndoPhenotypes/EGCUTEndophenotypesValidSamples-NoSex.txt-asLoadedByNormalizer.txt";
            String out = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-15-EGCUTCellCountCorrelations/CellCountCorrelationMatrix.txt";

            CorrelateCellCountsWithGeneExpression c = new CorrelateCellCountsWithGeneExpression();
            c.run(expfile, cellfile, out);

        } catch (IOException e) {
            e.printStackTrace();
        } catch (DocumentException e) {
            e.printStackTrace();
        }

    }

    public void run(String inexpFileName, String cellCountFileName, String output) throws IOException, DocumentException {
        DoubleMatrixDataset<String, String> expression = new DoubleMatrixDataset<String, String>(inexpFileName); // samples on the columns
        DoubleMatrixDataset<String, String> cellcounts = new DoubleMatrixDataset<String, String>(cellCountFileName); // samples somewhere else..

        // double[][] values, String[] rowHeaders, String[] colHeaders, int width, int height, String filename, Heatmap.Output output
        String[] cellTypes = cellcounts.rowObjects.toArray(new String[0]);
        String[] probes = expression.rowObjects.toArray(new String[0]);

        double[][] correlationMatrix = new double[probes.length][cellTypes.length];
        double[][] correlationMatrix2 = new double[cellTypes.length][probes.length];

        for (int cell = 0; cell < cellTypes.length; cell++) {
            for (int probe = 0; probe < probes.length; probe++) {
                ArrayList<Double> x = new ArrayList<Double>();
                ArrayList<Double> y = new ArrayList<Double>();
                for (int sample = 0; sample < cellcounts.nrCols; sample++) {
                    Integer sampleIdInOtherDataset = expression.hashCols.get(cellcounts.colObjects.get(sample));
                    if (sampleIdInOtherDataset != null) {
                        x.add(cellcounts.rawData[cell][sample]);
                        y.add(expression.rawData[probe][sampleIdInOtherDataset]);
                    }
                }
                
                double[] xarr = toPrimitiveArr(x.toArray(new Double[0]));
                double[] yarr = toPrimitiveArr(y.toArray(new Double[0]));
                double corr = JSci.maths.ArrayMath.correlation(xarr, yarr);
                correlationMatrix[probe][cell] = corr;
                correlationMatrix2[cell][probe] = corr;
            }
        }

        double[][] celltypecorrelationMatrix = new double[cellTypes.length][cellTypes.length];
        for (int cell = 0; cell < cellTypes.length; cell++) {
            for (int cell2 = 0; cell2 < cellTypes.length; cell2++) {
                celltypecorrelationMatrix[cell][cell2] = JSci.maths.ArrayMath.correlation(correlationMatrix2[cell], correlationMatrix2[cell2]);
            }
        }


        Heatmap.drawCorrelationHeatmap(celltypecorrelationMatrix, cellTypes, cellTypes, (cellTypes.length * 100) + 200, (cellTypes.length * 100) + 200, output + "CorrelationWithExpressionDataCorrelationsbetweencelltypesHeatMap.pdf", Heatmap.Output.PDF);
        
        double[][] celltypecorrelationMatrix2 = new double[cellTypes.length][cellTypes.length];
        for (int cell = 0; cell < cellTypes.length; cell++) {
            for (int cell2 = 0; cell2 < cellTypes.length; cell2++) {
                celltypecorrelationMatrix2[cell][cell2] = JSci.maths.ArrayMath.correlation(cellcounts.rawData[cell], cellcounts.rawData[cell2]);
            }
        }

        Heatmap.drawCorrelationHeatmap(celltypecorrelationMatrix2, cellTypes, cellTypes, (cellTypes.length * 100) + 200, (cellTypes.length * 100) + 200, output + "CorrelationBetweenCellTypes.pdf", Heatmap.Output.PDF);




        System.out.println(probes.length);
        System.out.println(cellTypes.length);
        System.out.println("");
        System.out.println((probes.length * 10));
        System.out.println((cellTypes.length * 10));
        
        DoubleMatrixDataset<String, String> outputds = new DoubleMatrixDataset<String, String>();
        outputds.rawData = correlationMatrix;
        outputds.colObjects = Arrays.asList(cellTypes);
        outputds.rowObjects = Arrays.asList(probes);
        outputds.recalculateHashMaps();
        outputds.save(output+"EndophenotypeVsExpressionCorrelationMatrix.txt");
        
        Heatmap.drawCorrelationHeatmap(correlationMatrix, probes, cellTypes, (cellTypes.length * 10) + 200, (probes.length * 10) + 400, output + "CorrelationWithExpressionDataHeatMap.pdf", Heatmap.Output.PDF);


    }

    private double[] toPrimitiveArr(Double[] toArray) {
        double[] arr = new double[toArray.length];
        for (int i = 0; i < toArray.length; i++) {
            arr[i] = toArray[i];
        }
        return arr;
    }
}
