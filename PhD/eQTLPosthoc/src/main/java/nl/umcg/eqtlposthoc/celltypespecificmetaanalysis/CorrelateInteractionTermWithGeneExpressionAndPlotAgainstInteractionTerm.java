/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.containers.Pair;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class CorrelateInteractionTermWithGeneExpressionAndPlotAgainstInteractionTerm {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            String outputdir = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-08-14-MetaAnalysis5DatasetsEffectsFlipped/CorrelationWithGeneExpressionForCellTypeProxies/";
            Gpio.createDir(outputdir);
            String[] datasetnames = new String[]{"EGCUT", "BloodHT12v3"};
            String[] expressionDatasets = new String[]{"/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-EGCUT/Normalization/ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.txt.gz",
                "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-Groningen/Normalization/ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.CovariatesRemoved.txt.gz"};
            String[] cellproxyfiles = new String[]{"/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-EGCUT/Normalization/CellTypeProxyFile.txt",
                "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-Groningen/Normalization/CellTypeProxyFile.txt"};
            String interactionterms = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-08-14-MetaAnalysis5DatasetsEffectsFlipped/MetaAnalysis/Vector-CellTypeInteractionZScore.txt";

            TextFile tf = new TextFile(interactionterms, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            ArrayList<Pair<String, Double>> snpprobes = new ArrayList<Pair<String, Double>>();
            HashSet<String> probes = new HashSet<String>();
            while (elems != null) {

                String snpprobe = elems[0];
                String[] snpprobeelems = snpprobe.split("-");
                probes.add(snpprobeelems[1]);

                Double d = Double.parseDouble(elems[1]);

                snpprobes.add(new Pair<String, Double>(snpprobe, d));

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            SpearmansCorrelation spearmancorrelation = new SpearmansCorrelation();

            HashSet<String> allProbesToMetaAnalyze = new HashSet<String>();
            ArrayList<HashMap<String, Double>> zScoresPerDataset = new ArrayList<HashMap<String, Double>>();
            ArrayList<HashMap<String, Double>> correlationsPerDataset = new ArrayList<HashMap<String, Double>>();
            Integer[] sampleSizes = new Integer[datasetnames.length];
            for (int dsNr = 0; dsNr < expressionDatasets.length; dsNr++) {
                String datasetname = datasetnames[dsNr];
                String expressionDataset = expressionDatasets[dsNr];
                String cellproxyfile = cellproxyfiles[dsNr];

                HashMap<String, Double> individualZScores = new HashMap<String, Double>();
                HashMap<String, Double> individualCorrelations = new HashMap<String, Double>();
                zScoresPerDataset.add(individualZScores);
                correlationsPerDataset.add(individualCorrelations);

                DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(expressionDataset);


                // load cell count proxy values
                TextFile tf2 = new TextFile(cellproxyfile, TextFile.R);
                tf2.readLine();
                String[] elems2 = tf2.readLineElems(TextFile.tab);
                HashMap<String, Double> cellcountproxyvalues = new HashMap<String, Double>();
                while (elems2 != null) {
                    String individual = elems2[0];
                    Double d = Double.parseDouble(elems2[1]);
                    cellcountproxyvalues.put(individual, d);
                    elems2 = tf2.readLineElems(TextFile.tab);
                }
                tf2.close();




                // now correlate..
                for (int expressionProbeId = 0; expressionProbeId < ds.nrRows; expressionProbeId++) {
                    String expressionProbe = ds.rowObjects.get(expressionProbeId);
                    if (probes.contains(expressionProbe)) {
                        ArrayList<Double> valsX = new ArrayList<Double>();
                        ArrayList<Double> valsY = new ArrayList<Double>();
                        for (int col = 0; col < ds.nrCols; col++) {
                            String name2 = ds.colObjects.get(col);

                            Double d = cellcountproxyvalues.get(name2);
                            if (d != null) {
                                valsX.add(ds.rawData[expressionProbeId][col]);
                                valsY.add(d);
                            }

                        }

                        double[] xarrd = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                        double[] yarrd = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));
                        sampleSizes[dsNr] = xarrd.length;
                        double spearman = spearmancorrelation.correlation(xarrd, yarrd);
                        Correlation.correlationToZScore(sampleSizes[dsNr]);

                        double z = Correlation.convertCorrelationToZScore(sampleSizes[dsNr], spearman);
                        individualZScores.put(expressionProbe, z);
                        individualCorrelations.put(expressionProbe, spearman);
                        allProbesToMetaAnalyze.add(expressionProbe);

                    }
                }

                // plot
                ArrayList<Double> valsX = new ArrayList<Double>();
                ArrayList<Double> valsY = new ArrayList<Double>();
                ArrayList<Double> valsY2 = new ArrayList<Double>();

                for (Pair<String, Double> p : snpprobes) {
                    String snpprobe = p.getLeft();
                    String probe = snpprobe.split("-")[1];
                    Double d = p.getRight();
                    Double z = individualZScores.get(probe);
                    Double c = individualCorrelations.get(probe);
                    if (z != null) {
                        valsX.add(d);
                        valsY.add(z);
                        valsY2.add(c);
                    }
                }
                double[] xarrd = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                double[] yarrd = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));
                double[] yarrd2 = Primitives.toPrimitiveArr(valsY2.toArray(new Double[0]));

                System.out.println("Zscore correlation " + datasetname + "\t" + spearmancorrelation.correlation(xarrd, yarrd));

                new ScatterPlot(500, 500, xarrd, yarrd, ScatterPlot.OUTPUTFORMAT.PDF, outputdir + datasetname + ".pdf");
                new ScatterPlot(500, 500, xarrd, yarrd2, ScatterPlot.OUTPUTFORMAT.PDF, outputdir + datasetname + "-correlations.pdf");
            }

            for (int d = 0; d < datasetnames.length; d++) {
                HashMap<String, Double> correlations1 = correlationsPerDataset.get(d);
                for (int d2 = d + 1; d2 < datasetnames.length; d2++) {
                    HashMap<String, Double> correlations2 = correlationsPerDataset.get(d2);
                    ArrayList<Double> valsX = new ArrayList<Double>();
                    ArrayList<Double> valsY = new ArrayList<Double>();
                    for (String probe : probes) {
                        Double val1 = correlations1.get(probe);
                        Double val2 = correlations2.get(probe);
                        if (val1 != null && val2 != null) {
                            valsX.add(val1);
                            valsY.add(val2);
                        }
                    }
                    double[] xarrd = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
                    double[] yarrd = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));
                    new ScatterPlot(500, 500, xarrd, yarrd, ScatterPlot.OUTPUTFORMAT.PDF, outputdir + datasetnames[d] + "-" + datasetnames[d2] + ".pdf");
                    System.out.println("Correlation correlation " + datasetnames[d] + "-" + datasetnames[d2] + "\t" + spearmancorrelation.correlation(xarrd, yarrd));
                }
            }

            // meta-analyze..
            HashMap<String, Double> metaZPerProbe = new HashMap<String, Double>();
            for (String probe : allProbesToMetaAnalyze) {
                double[] zarr = new double[datasetnames.length];
                int[] narr = new int[datasetnames.length];
                for (int d = 0; d < zarr.length; d++) {
                    Integer n = sampleSizes[d];
                    HashMap<String, Double> zs = zScoresPerDataset.get(d);
                    Double z = zs.get(probe);
                    if (z != null) {

                        zarr[d] = z;
                        narr[d] = n;
                    } else {
                        zarr[d] = Double.NaN;
                        narr[d] = 0;
                    }
                }

                double metaZ = ZScores.getWeightedZ(zarr, narr);
                metaZPerProbe.put(probe, metaZ);
            }

            ArrayList<Double> valsX = new ArrayList<Double>();
            ArrayList<Double> valsY = new ArrayList<Double>();

            for (Pair<String, Double> p : snpprobes) {
                String snpprobe = p.getLeft();
                String probe = snpprobe.split("-")[1];
                Double d = p.getRight();
                Double z = metaZPerProbe.get(probe);
                if (z != null) {
                    valsX.add(d);
                    valsY.add(z);
                }
            }
            double[] xarrd = Primitives.toPrimitiveArr(valsX.toArray(new Double[0]));
            double[] yarrd = Primitives.toPrimitiveArr(valsY.toArray(new Double[0]));

            System.out.println();
            new ScatterPlot(500, 500, xarrd, yarrd, ScatterPlot.OUTPUTFORMAT.PDF, outputdir + "MetaAnalysis.pdf");

            System.out.println("InteractionTerm MetaCorrelation\t" + spearmancorrelation.correlation(xarrd, yarrd));



        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
