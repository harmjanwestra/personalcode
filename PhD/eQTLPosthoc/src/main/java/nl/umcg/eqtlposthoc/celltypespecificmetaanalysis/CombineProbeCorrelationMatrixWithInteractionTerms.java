/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class CombineProbeCorrelationMatrixWithInteractionTerms {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here

            String metaMatrix = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-12-GroningenHT12v3/CellTypeSpecificityMatrix.binary";
            String rawExpressionFile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-Groningen/Normalization/ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.CovariatesRemoved.txt.gz";
            String pcCorrectedExpressionFile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.40PCAsOverSamplesRemoved.txt.gz";
            String outputdir = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-08-02-CorrelateGeneExpressionCorrelationWithInteractionTermsGroningen/";

            CombineProbeCorrelationMatrixWithInteractionTerms q = new CombineProbeCorrelationMatrixWithInteractionTerms();
//            q.run(rawExpressionFile, pcCorrectedExpressionFile, metaMatrix, outputdir);
            String correlationOutput = outputdir + "CorrelationOfInteractionTermsWithProbeProbeCorrelation.txt";
            String correlationMatrix = outputdir + "ProbeCorrelationMatrix.txt";
            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            q.makeplots(metaMatrix, correlationMatrix, correlationOutput, outputdir, probetranslationfile);
        } catch (IOException ex) {
            Logger.getLogger(CombineProbeCorrelationMatrixWithInteractionTerms.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void run(String rawExpressionFile, String pcCorrectedExpressionFile, String metaMatrix, String outputdir) throws IOException {
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(metaMatrix); // eQTLs on columns
        DoubleMatrixDataset<String, String> dsPCCExp = new DoubleMatrixDataset<String, String>(pcCorrectedExpressionFile); // eQTLs on columns
        DoubleMatrixDataset<String, String> dsRawExp = new DoubleMatrixDataset<String, String>(rawExpressionFile); // eQTLs on columns

        HashSet<String> probesOfInterest = new HashSet<String>();
        for (int i = 0; i < ds.nrCols; i++) {
            String probe = ds.colObjects.get(i).split("-")[1];
            probesOfInterest.add(probe);
        }

        ArrayList<String> allProbes = new ArrayList<String>(probesOfInterest);
        SpearmansCorrelation corr = new SpearmansCorrelation();
        double[][] probeCorrelationMatrix = new double[allProbes.size()][allProbes.size()];
        System.out.println("Starting correlations");
        ProgressBar pb = new ProgressBar(allProbes.size());
        for (int i = 0; i < allProbes.size(); i++) {
            String probe1 = allProbes.get(i);
            Integer probeId1 = dsPCCExp.hashRows.get(probe1);
            if (probeId1 != null) {
                for (int j = 0; j < allProbes.size(); j++) {
                    String probe2 = allProbes.get(j);
                    Integer probeId2 = dsRawExp.hashRows.get(probe2);
                    if (probeId2 != null) {
                        ArrayList<Double> valsX = new ArrayList<Double>();
                        ArrayList<Double> valsY = new ArrayList<Double>();
                        for (int p1 = 0; p1 < dsPCCExp.nrCols; p1++) {
                            String person1 = dsPCCExp.colObjects.get(p1);
                            Integer person2Id = dsRawExp.hashCols.get(person1);
                            if (person2Id != null) {
                                valsX.add(dsPCCExp.rawData[probeId1][p1]);
                                valsY.add(dsRawExp.rawData[probeId2][person2Id]);
                            }
                        }

                        double[] xArr = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                        double[] yArr = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));
                        double correlation = corr.correlation(xArr, yArr);
                        probeCorrelationMatrix[i][j] = correlation;
                    }
                }
            }
            pb.set(i);
        }
        pb.close();

        DoubleMatrixDataset<String, String> corrMat = new DoubleMatrixDataset<String, String>();
        corrMat.rawData = probeCorrelationMatrix;
        corrMat.colObjects = allProbes;
        corrMat.rowObjects = allProbes;
        corrMat.recalculateHashMaps();
        corrMat.save(outputdir + "ProbeCorrelationMatrix.txt");

        // now iterate over all eQTLs, and correlate the individual covariate 
        // interaction term scores with the gene-gene correlations.

        System.out.println("Performing final calculations..");
        TextFile out = new TextFile(outputdir + "CorrelationOfInteractionTermsWithProbeProbeCorrelation.txt", TextFile.W);
        for (int i = 0; i < ds.nrCols; i++) {
            String eQTL = ds.colObjects.get(i);
            String probe = eQTL.split("-")[1];
            Integer correlationMatrixeQTLProbeId = corrMat.hashRows.get(probe);
            if (correlationMatrixeQTLProbeId != null) {

                ArrayList<Double> valsX = new ArrayList<Double>();
                ArrayList<Double> valsY = new ArrayList<Double>();

                ArrayList<Double> valsXAbs = new ArrayList<Double>();
                ArrayList<Double> valsYAbs = new ArrayList<Double>();

                for (int j = 0; j < ds.nrRows; j++) {
                    String covariateProbeName = ds.rowObjects.get(j);
                    Integer correlationMatrixCovariateProbeId = corrMat.hashRows.get(covariateProbeName);
                    if (correlationMatrixCovariateProbeId != null) {

                        valsX.add(ds.rawData[j][i]);
                        valsY.add(corrMat.rawData[correlationMatrixeQTLProbeId][correlationMatrixCovariateProbeId]);

                        valsXAbs.add(Math.abs(ds.rawData[j][i]));
                        valsYAbs.add(Math.abs(corrMat.rawData[correlationMatrixeQTLProbeId][correlationMatrixCovariateProbeId]));
                    }
                }

                double[] xArr = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                double[] yArr = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));
                double correlation = corr.correlation(xArr, yArr);

                xArr = Primitives.toPrimitiveArr(valsXAbs.toArray(new Double[0]));
                yArr = Primitives.toPrimitiveArr(valsYAbs.toArray(new Double[0]));
                double abscorrelation = corr.correlation(xArr, yArr);

                out.writeln(eQTL + "\t" + correlation + "\t" + abscorrelation + "\t" + xArr.length);

                System.out.println(eQTL + "\t" + correlation + "\t" + abscorrelation + "\t" + xArr.length);

            }

        }

        out.close();

    }

    private void makeplots(String metaMatrix, String probeCorrelationMatrix, String correlationOutput, String outputdir, String probetranslationfile) throws IOException {
        DoubleMatrixDataset<String, String> interactionTermMatrix = new DoubleMatrixDataset<String, String>(metaMatrix); // eQTLs on columns
        DoubleMatrixDataset<String, String> correlationMatrix = new DoubleMatrixDataset<String, String>(probeCorrelationMatrix);

        TextFile tf = new TextFile(correlationOutput, TextFile.R);
        SpearmansCorrelation corr = new SpearmansCorrelation();

        String[] elems = tf.readLineElems(TextFile.tab);

        String ht12v3String = "HumanHT-12_V3_0_R2_11283641_A.txt";
        String ht12v4String = "HumanHT-12_V4_0_R1_15002873_B.txt";

        String hugoStr = "HUGO";

        //       PROBE LOOKUP TABLES..
        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> ht12v3ToHugo = new HashMap<String, String>();
        ht12v3ToHugo = pbt.getProbeTranslation(probetranslationfile, ht12v3String, hugoStr);

        while (elems != null) {
            String eqtl = elems[0];
            Integer eqtlColId = interactionTermMatrix.hashCols.get(eqtl);
            double correlation = Double.parseDouble(elems[1]);
            if (Math.abs(correlation) > 0.75) {
                String probe1 = eqtl.split("-")[1];

                Integer probe1RowId = correlationMatrix.hashRows.get(probe1);
                if (probe1RowId != null) {
//                Integer probe2ColId = ds2.hashCols.get(probe1);
                    TextFile tfOut = new TextFile(outputdir + eqtl + "-" + ht12v3ToHugo.get(probe1) + ".txt", TextFile.W);

                    tfOut.writeln("EQTL\tInteractionTerm\tHUGO\tInteractionZscore\tCorrelationWithEQTLProbe");
                    ArrayList<Double> valsX = new ArrayList<Double>();
                    ArrayList<Double> valsY = new ArrayList<Double>();
                    for (int interactionTerm = 0; interactionTerm < interactionTermMatrix.nrRows; interactionTerm++) {
                        String interactionTermProbeName = interactionTermMatrix.rowObjects.get(interactionTerm);
                        Integer probe2ColId = correlationMatrix.hashCols.get(interactionTermProbeName);
                        if (probe2ColId != null) {
                            valsX.add(interactionTermMatrix.rawData[interactionTerm][eqtlColId]);
                            valsY.add(correlationMatrix.rawData[probe1RowId][probe2ColId]);
                            tfOut.writeln(eqtl + "\t" + interactionTermProbeName + "\t" + ht12v3ToHugo.get(interactionTermProbeName) + "\t" + interactionTermMatrix.rawData[interactionTerm][eqtlColId] + "\t" + correlationMatrix.rawData[probe1RowId][probe2ColId]);
                        }
                    }

                    tfOut.close();
                    double[] xArr = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                    double[] yArr = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));
                    double spearmancorr = corr.correlation(xArr, yArr);

                    new ScatterPlot(500, 500, xArr, yArr, ScatterPlot.OUTPUTFORMAT.PDF, outputdir + spearmancorr + "-" + eqtl + ".pdf");

                    System.out.println(correlation + "\t" + eqtl + "\t" + spearmancorr);
                }



            }




            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


    }
}
