/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.WilcoxonMannWhitney;

/**
 *
 * @author harmjan
 */
public class CalculateMeanExpression {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//         TODO code application logic here
//        String[] matrixList = new String[]{
//            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataCHB/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
//            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataGIH/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
//            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataJPT/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
//            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataLWK/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
//            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataMEX/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
//            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataMKK/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
//            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataYRI/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz"
//        };
        String[] matrixList = new String[]{
            "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/Fairfax/BCells/ExpressionData.QuantileNormalized.Log2Transformed-HT12v3.txt.gz"
        };
      
        String outputFile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/Fairfax/BCells/MeanExpresion.txt";

//        
//        String[] matrixList = new String[]{"/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodH8v2/CovariateCorrectedData/ExpressionData.QuantileNormalized.Log2Transformed.CovariatesRemoved.txt.gz"};
//        String outputFile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-06-BloodH8v2/MeanExpresion.txt";
        String probeTranslationFile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";

//        String file1 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-04-CD4AndCD8Cells/CD4-ExpressionMeanExpressionVariance.txt";
//        String file1out = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-04-CD4AndCD8Cells/CD4-ExpressionMeanExpressionVariance-HT12v3.txt";
//        String file2 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-04-CD4AndCD8Cells/CD8-ExpressionMeanExpressionVariance.txt";
//        String file2out = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-04-CD4AndCD8Cells/CD8-ExpressionMeanExpressionVariance-HT12v3.txt";
        try {
//            Gpio.createDir("/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-29-FairfaxBCells");
            CalculateMeanExpression r = new CalculateMeanExpression();
            r.run(matrixList, outputFile, probeTranslationFile);
//            r.probeRewrite(file1, file1out, probeTranslationFile, "HT12v4.txt", "HT12v3.txt");
//            r.probeRewrite(file2, file2out, probeTranslationFile, "HT12v4.txt", "HT12v3.txt");
        } catch (IOException e) {
        }

    }

    public void run(String[] matrixlist, String output, String probeTranslation) throws IOException {

        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> wg6ToHt12v3 = pb.getProbeTranslation(probeTranslation, "HT12v3.txt", "HT12v3.txt");
//        HashMap<String, String> wg6ToHt12v3 = pb.getProbeTranslation(probeTranslation, "HT12v3.txt", "HT12v3.txt");

        HashSet<String> probes = new HashSet<String>();
        HashMap<String, Double> sumsAverage = new HashMap<String, Double>();
        HashMap<String, Double> sumsVariance = new HashMap<String, Double>();
        HashMap<String, Integer> count = new HashMap<String, Integer>();

        for (String s : matrixlist) {
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(s);
            for (int row = 0; row < ds.nrRows; row++) {
                String probe = ds.rowObjects.get(row);

                probes.add(probe);
                double median = JSci.maths.ArrayMath.mean(ds.rawData[row]);
                double variance = JSci.maths.ArrayMath.variance(ds.rawData[row]);
                Double runningAverage = sumsAverage.get(probe);
                Double runningVariance = sumsVariance.get(probe);
                Integer ct = count.get(probe);
                if (runningAverage == null) {
                    runningAverage = 0d;
                    runningVariance = 0d;
                    ct = 0;
                }
                ct++;

                runningAverage += median;
                runningVariance += variance;
                count.put(probe, ct);
                sumsAverage.put(probe, runningAverage);
                sumsVariance.put(probe, runningVariance);
            }
        }

        TextFile out = new TextFile(output, TextFile.W);
        for (String probe : probes) {
            Integer ct = count.get(probe);
            Double average = sumsAverage.get(probe);
            Double var = sumsVariance.get(probe);

            double avg = average / ct;
            double avgvar = var / ct;
            String ht12v43 = wg6ToHt12v3.get(probe);
            if (ht12v43 != null && !ht12v43.equals("-")) {
                out.writeln(ht12v43 + "\t" + avg + "\t" + avgvar);
            }

        }
        out.close();
    }

    public void probeRewrite(String file, String output, String probeTranslation, String src, String dst) throws IOException {
        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> probeToProbe = pb.getProbeTranslation(probeTranslation, src, dst);

        TextFile in = new TextFile(file, TextFile.R);
        in.readLine();
        TextFile out = new TextFile(output, TextFile.W);
        String[] elems = in.readLineElems(TextFile.tab);
        HashSet<String> probesOutput = new HashSet<String>();
        while (elems != null) {
            String probe = elems[0].split("_")[1];
            String val = elems[1];
            String val2 = elems[2];

            if (!val.equals("NA")) {
                if (!probesOutput.contains(probe)) {
                    out.writeln(probeToProbe.get(probe) + "\t" + val + "\t" + val2);
                    probesOutput.add(probe);
                }
            }

            elems = in.readLineElems(TextFile.tab);
        }
        out.close();
        in.close();


    }
}
