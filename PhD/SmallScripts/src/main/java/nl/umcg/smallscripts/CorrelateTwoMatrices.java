/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class CorrelateTwoMatrices {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here
//            DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/ExpressionDataCovCorrected/ExpressionData.txt.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz.CovariatesRemoved.txt.gz");
//            DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/ExpressionDataCovCorrected/ExpressionData.txt.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz-SampleSizeCorrectedForCovariates.txt");
//
//            double[][] rawData1 = ds1.rawData;
//            for(int row =0 ; row < rawData1.length; row++){
//                double[] data1 = rawData1[row];
//                // get corresponding row in data2
//                Integer row2 = ds2.hashRows.get(ds1.rowObjects.get(row));
//                if(row2!=null){
//                    double[] data2 = ds2.rawData[row2];
//                    double[] newdata2 = new double[data2.length];
//                    for(int i=0 ;i<newdata2.length; i++){
//                        int col1 = ds1.hashCols.get(ds2.colObjects.get(i));
//                        newdata2[col1] = data2[i];
//                    }
//                    double r = JSci.maths.ArrayMath.correlation(data1, newdata2);
//                    System.out.println(ds1.rowObjects.get(row)+"\t"+r);
//                }
//            }

            DoubleMatrixDataset<String, String> uncorrected = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/ExpressionDataCovCorrected/ExpressionData.txt.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz-SampleSizeCorrectedForCovariates.txt");
            DoubleMatrixDataset<String, String> corrected = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/ExpressionDataCovCorrected/ExpressionData.txt.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz.CovariatesRemoved.txt.gz");
            DoubleMatrixDataset<String, String> endopheno = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-Endophenotypecorrelations/EGCUT/RawEndoPhenotypes/EGCUTEndophenotypesValidSamples.txt");

            endopheno.transposeDataset();
            double[][] uncorrectedRawData = uncorrected.rawData;
            double[][] correctedRawData = corrected.rawData;
            double[][] endoPhenoRawData = endopheno.rawData;

            for (int row = 0; row < uncorrectedRawData.length; row++) {
                double[] uncorrectedData = uncorrectedRawData[row];
                // get corresponding row in data2
                int correctedRow = corrected.hashRows.get(uncorrected.rowObjects.get(row));
                double[] correctedDataUnsorted = correctedRawData[correctedRow];
                double[] correctedDataSorted = new double[correctedDataUnsorted.length];

                for (int col = 0; col < correctedDataSorted.length; col++) {
                    int uncorrectedRowId = uncorrected.hashCols.get(corrected.colObjects.get(col));
                    correctedDataSorted[uncorrectedRowId] = correctedDataUnsorted[col];
                }

                for (int traitRow = 0; traitRow < endoPhenoRawData.length; traitRow++) {
                    double[] tmpPhenoRow = endopheno.rawData[traitRow];
                    double[] phenoRow = new double[tmpPhenoRow.length];
                    for (int i = 0; i < phenoRow.length; i++) {
                        int uncorrectedRowId = uncorrected.hashCols.get(endopheno.colObjects.get(i));
                        phenoRow[uncorrectedRowId] = tmpPhenoRow[i];
                    }
                    double rVsUnCorrected = JSci.maths.ArrayMath.correlation(uncorrectedData, phenoRow);
                    double rVsCorrected = JSci.maths.ArrayMath.correlation(correctedDataSorted, phenoRow);
                    double expressionComparison = JSci.maths.ArrayMath.correlation(correctedDataSorted, uncorrectedData);

                    if (rVsUnCorrected > 0.1) {
//                        System.out.println(uncorrected.rowObjects.get(row) + "\t" + endopheno.rowObjects.get(traitRow) + "\t" + rVsUnCorrected + "\t" + rVsCorrected + "\t" + expressionComparison);
                    }
                    if (rVsUnCorrected > 0.7) {
                        System.out.println(uncorrected.rowObjects.get(row) + "\t" + endopheno.rowObjects.get(traitRow));

                        for (int col = 0; col < uncorrectedData.length; col++) {
                            System.out.println(col + "\t" + uncorrected.colObjects.get(col) + "\t" + uncorrectedData[col] + "\t" + correctedDataSorted[col] + "\t" + phenoRow[col]);
                        }
                        System.out.println("");
//                        ScatterPlot p = new ScatterPlot(500,500,correctedDataSorted, uncorrectedData, 
//                                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-EGCUTTransCorrectedForCovariates/"+uncorrected.rowObjects.get(row)+"-"+endopheno.rowObjects.get(traitRow)+"-r"+expressionComparison+"-expressionComparison.png");
//                        
//                        ScatterPlot p2 = new ScatterPlot(500,500,correctedDataSorted, phenoRow,
//                                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-EGCUTTransCorrectedForCovariates/"+uncorrected.rowObjects.get(row)+"-"+endopheno.rowObjects.get(traitRow)+"-r"+rVsCorrected+"-covariateVsCorrected.png");
//                        
//                        ScatterPlot p3 = new ScatterPlot(500,500,uncorrectedData, phenoRow, 
//                                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-EGCUTTransCorrectedForCovariates/"+uncorrected.rowObjects.get(row)+"-"+endopheno.rowObjects.get(traitRow)+"-r"+rVsUnCorrected+"-covariateVsUnCorrected.png");
                        System.exit(0);
                    }
                }
            }

        } catch (IOException ex) {
            Logger.getLogger(CorrelateTwoMatrices.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
