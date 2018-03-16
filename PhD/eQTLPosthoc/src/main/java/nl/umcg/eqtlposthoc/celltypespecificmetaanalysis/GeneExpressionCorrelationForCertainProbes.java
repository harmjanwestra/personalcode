/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class GeneExpressionCorrelationForCertainProbes {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String[] datasetLocations = new String[]{
            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/RawExpressionData/ExpressionData.QuantileNormalized.Log2Transformed.txt.gz",
            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz"};
        String[] datasetNames = new String[]{"EGCUT","Groningen"};
        String[] probeList = new String[]{"7330546","2650750","4920202"};
        String outdir = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-10-09-MED8CovariateOverlap/";
        try {
            GeneExpressionCorrelationForCertainProbes c = new GeneExpressionCorrelationForCertainProbes();
            c.run(datasetLocations, datasetNames, probeList, outdir);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void run(String[] datasetLocations, String[] datasetNames, String[] probeList, String outdir) throws IOException {


        Gpio.createDir(outdir);
        HashSet<String> queryProbes = new HashSet<String>();
        queryProbes.addAll(Arrays.asList(probeList));

        SpearmansCorrelation correlation = new SpearmansCorrelation();

        for (int i = 0; i < datasetLocations.length; i++) {
            String loc = datasetLocations[i];
            String dsName = datasetNames[i];
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(loc, queryProbes);
            TextFile outputfile = new TextFile(outdir + dsName + "-CorrelationMatrix.txt", TextFile.W);
            outputfile.append("-");
            for (int row = 0; row < ds.nrRows; row++) {
                outputfile.append("\t" + ds.rowObjects.get(row));
            }
            outputfile.append("\n");
            for (int row = 0; row < ds.nrRows; row++) {
                double[] data = ds.rawData[row];
                outputfile.append(ds.rowObjects.get(row));
                for (int row2 = 0; row2 < ds.nrRows; row2++) {
                    double[] data2 = ds.rawData[row2];
                    double corr = correlation.correlation(data, data2);
                    outputfile.append("\t" + corr);
                }

                outputfile.append("\n");
            }
            outputfile.close();

        }





    }
}
