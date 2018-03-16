/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class KarinCorrelationWithDiseaseRisk {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String probeAnnotation = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-HT12v3.txt";
            String expressionfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz";
            String riskloadfile = "/Volumes/iSnackHD/Data/Projects/KarinFransen/2013-11-19-riskmodelrebuttalonlyriskscore.txt";
            String outputfile = "/Volumes/iSnackHD/Data/Projects/KarinFransen/2013-11-19-riskmodelrebuttalonlyriskscore-Correlations-WithPValues.txt";
            TextFile tf = new TextFile(riskloadfile, TextFile.R);
            HashMap<String, Double> sampleToRisk = new HashMap<String, Double>();
            HashMap<String, Double> sampleToRiskWeighted = new HashMap<String, Double>();

            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String sample = elems[1];
                double risk = Double.parseDouble(elems[6]);
                double riskweighed = Double.parseDouble(elems[7]);
                sampleToRisk.put(sample, risk);
                sampleToRiskWeighted.put(sample, riskweighed);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            HashMap<String, String> probeToAnnotation = new HashMap<String, String>();
            TextFile tf2 = new TextFile(probeAnnotation, TextFile.R);
            tf2.readLine();
            String[] elems2 = tf2.readLineElems(TextFile.tab);

            while (elems2 != null) {
                String probe = elems2[1];
                //HT12v3	2190220	TTF2	1	117644958	117645007	3	CAGCCATCTCTGCAGTTCTCTCAGTGCAGGCAGTTCTTCCTCTCAGGCTG
                String annotation = elems2[2] + "\t" + elems2[3] + "\t" + elems2[4] + "\t" + elems2[5];

                probeToAnnotation.put(probe, annotation);
                elems2 = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(expressionfile, null, sampleToRisk.keySet());

            ArrayList<Double> riskValues = new ArrayList<Double>();
            ArrayList<Double> riskValuesWeighted = new ArrayList<Double>();
            for (int col = 0; col < ds.nrCols; col++) {
                String sample = ds.colObjects.get(col);
                Double risk = sampleToRisk.get(sample);
                if (risk != null) {
                    riskValues.add(risk);
                    riskValuesWeighted.add(sampleToRiskWeighted.get(sample));
                }
            }
            System.out.println("Overlap: " + riskValuesWeighted.size());
            double[] riskValueArray = Primitives.toPrimitiveArr(riskValues.toArray(new Double[0]));
            double[] riskValueWeightedArray = Primitives.toPrimitiveArr(riskValuesWeighted.toArray(new Double[0]));
            SpearmansCorrelation correlator = new SpearmansCorrelation();

            TextFile output = new TextFile(outputfile, TextFile.W);
            output.writeln("Probe\tP-Unweighted\tZ-Unweighted\tSpearmanCorrelationRisk\tP-Weighted\tZ-Weighted\tSpearmanCorrelationRiskWeighted\tN\tGene\tChr\tChrStartPos\tChrStopPos");
            for (int row = 0; row < ds.nrRows; row++) {
                ArrayList<Double> expressionValues = new ArrayList<Double>();
                for (int col = 0; col < ds.nrCols; col++) {
                    String sample = ds.colObjects.get(col);
                    Double risk = sampleToRisk.get(sample);
                    if (risk != null) {
                        expressionValues.add(ds.rawData[row][col]);
                    }
                }

                double[] expressionValueArray = Primitives.toPrimitiveArr(expressionValues.toArray(new Double[0]));
                double correlation = correlator.correlation(riskValueArray, expressionValueArray);
                double correlationWeighed = correlator.correlation(riskValueWeightedArray, expressionValueArray);

                String probe = ds.rowObjects.get(row);
                Correlation.correlationToZScore(riskValueWeightedArray.length);
                double zunweighted = Correlation.convertCorrelationToZScore(riskValueWeightedArray.length, correlation);
                double zweighted = Correlation.convertCorrelationToZScore(riskValueWeightedArray.length, correlationWeighed);
                double punweighted = ZScores.zToP(zunweighted);
                double pweighted = ZScores.zToP(zweighted);
                output.writeln(ds.rowObjects.get(row) + "\t" + punweighted + "\t" + zunweighted + "\t" + correlation + "\t" + pweighted + "\t" + zweighted + "\t" + correlationWeighed + "\t" + expressionValueArray.length + "\t" + probeToAnnotation.get(probe));
            }
            output.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
