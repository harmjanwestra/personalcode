/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import umcg.genetica.containers.Pair;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author harmjan
 */
public class CorrelateMetaAnalysisInteractionTermsWithCellTypeSpecificVectors {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String dannyfile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/DannyCorrelationAgainstMetaanalysis/COR_MetaAnalysisZScore_CellTypeTScore.txt";
        String metafile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/2013-06-26-MetaAnalysis/MetaAnalysisZScoreMatrix.txt";
        String output = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/DannyCorrelationAgainstMetaanalysis/";

        try {
            CorrelateMetaAnalysisInteractionTermsWithCellTypeSpecificVectors c = new CorrelateMetaAnalysisInteractionTermsWithCellTypeSpecificVectors();
            c.run(metafile, dannyfile, output);
        } catch (IOException e) {
        }
    }

    public void run(String metaFile, String vectorFile, String output) throws IOException {

        output = Gpio.formatAsDirectory(output);
        Gpio.createDir(output);

        DoubleMatrixDataset<String, String> metaMatrix = new DoubleMatrixDataset<String, String>(metaFile);

        DoubleMatrixDataset<String, String> correlationMatrix = new DoubleMatrixDataset<String, String>(vectorFile, ",");

        List<String> rows = correlationMatrix.rowObjects;
        ArrayList<String> newRows = new ArrayList<String>();
        for (String row : rows) {
            newRows.add(row.replace(".", "-"));
        }

        correlationMatrix.rowObjects = newRows;
        correlationMatrix.recalculateHashMaps();

        Integer rowIdForInteractionTerm = metaMatrix.hashRows.get("CellTypeInteractionZScore");
        Integer rowIdForMainEffect = metaMatrix.hashRows.get("MainEffectBeta");

        TextFile tfOut = new TextFile(output + "Correlations.txt", TextFile.W);
        tfOut.writeln("CellType\tCorrelation\tRSquared\tN");

        HashSet<Pair<Double, String>> eqtls = new HashSet<Pair<Double, String>>();

        double[] interactionTerms = metaMatrix.rawData[rowIdForInteractionTerm];
        double[] mainEffects = metaMatrix.rawData[rowIdForMainEffect];
        for (int ct = 0; ct < correlationMatrix.nrCols; ct++) {
            String cellType = correlationMatrix.colObjects.get(ct);
            ArrayList<Double> valsX = new ArrayList<Double>();
            ArrayList<Double> valsY = new ArrayList<Double>();
            for (int r = 0; r < correlationMatrix.nrRows; r++) {
                String eQTL = correlationMatrix.rowObjects.get(r);
                Integer eQTLIdInMetaMatrix = metaMatrix.hashCols.get(eQTL);
                if (eQTLIdInMetaMatrix != null) {
                    valsX.add(correlationMatrix.rawData[r][ct]);

                    double maineffect = mainEffects[eQTLIdInMetaMatrix];
                    double direction = 1;
                    if (maineffect < 0) {
                        direction = -1;
                    }

                    // force the interaction term into the same direction as the main effect

                    valsY.add(direction * interactionTerms[eQTLIdInMetaMatrix]);
                    eqtls.add(new Pair<Double, String>(direction * interactionTerms[eQTLIdInMetaMatrix], eQTL, "\t", Pair.SORTBY.LEFT));
                    if (Math.abs(interactionTerms[eQTLIdInMetaMatrix]) > 20) {
                        System.out.println(eQTL);
                    }
                }

            }
            if (!valsX.isEmpty()) {
                double[] xArr = toPrimitiveArr(valsX.toArray(new Double[0]));
                double[] yArr = toPrimitiveArr(valsY.toArray(new Double[0]));

                double corr = JSci.maths.ArrayMath.correlation(xArr, yArr);

                tfOut.writeln(cellType + "\t" + corr + "\t" + (corr * corr) + "\t" + xArr.length);
                ScatterPlot plot = new ScatterPlot(500, 500, xArr, yArr, ScatterPlot.OUTPUTFORMAT.PDF, output + cellType + "-InteractionTerm.pdf");
            }
        }

        ArrayList<Pair<Double, String>> eqtlarr = new ArrayList<Pair<Double, String>>();
        eqtlarr.addAll(eqtls);
        Collections.sort(eqtlarr);


//        HashMap<String, Pair<Double, Double>> ludeEQTLs = new HashMap<String, Pair<Double, Double>>();
//        TextFile tftmp = new TextFile("/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/DannyCorrelationAgainstMetaanalysis/NeutrophilSpecificCis-eQTLsSubsetColumns.txt", TextFile.R);
//        String[] header = tftmp.readLineElems(TextFile.tab);
//        String[] lnelems = tftmp.readLineElems(TextFile.tab);
//        while (lnelems != null) {
//            double p = Double.parseDouble(lnelems[3]);
//            double z = ZScores.pToZ(p);
//            String probe = lnelems[0];
//            String snp = lnelems[1];
//            String eQTL = snp + "-" + probe;
//            ludeEQTLs.put(eQTL, new Pair<Double, Double>(p, -z, "\t"));
//            lnelems = tftmp.readLineElems(TextFile.tab);;
//        }
//        tftmp.close();
//
//        ProbeTranslation pbt = new ProbeTranslation();
//        HashMap<String, String> ht12ToGene = pbt.getProbeTranslation("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", "HumanHT-12_V3_0_R2_11283641_A.txt", "HUGO");
//        
//        TextFile tfeQTLOut = new TextFile(output + "AllEQTLsSortedByInteractionTerm.txt", TextFile.W);
//        tfeQTLOut.writeln("Pvalue\tZScore\teQTL\tLudeP\tLudeZ\tHUGO");
//        for (int i = 0; i < eqtlarr.size(); i++) {
//            Pair<Double, Double> p = ludeEQTLs.get(eqtlarr.get(i).getRight());
//            String ludeqtl = null;
//            if (p != null) {
//                ludeqtl = p.toString();
//            } else {
//                ludeqtl = "-\t-";
//            }
//            tfeQTLOut.writeln(ZScores.zToP(eqtlarr.get(i).getLeft()) + "\t" + eqtlarr.get(i).toString() + "\t" + ludeqtl+"\t"+ht12ToGene.get(eqtlarr.get(i).getRight().split("-")[1]));
//        }
//        tfeQTLOut.close();



        correlationMatrix.transposeDataset();

        for (int r = 0; r < correlationMatrix.nrRows; r++) {
            for (int r2 = r + 1; r2 < correlationMatrix.nrRows; r2++) {
                String cellType1 = correlationMatrix.rowObjects.get(r);
                String cellType2 = correlationMatrix.rowObjects.get(r2);
                double corr = JSci.maths.ArrayMath.correlation(correlationMatrix.rawData[r], correlationMatrix.rawData[r2]);
                ScatterPlot plot = new ScatterPlot(500, 500, correlationMatrix.rawData[r], correlationMatrix.rawData[r2], ScatterPlot.OUTPUTFORMAT.PDF, output + cellType1 + "-" + cellType2 + ".pdf");

                tfOut.writeln(cellType1 + "-" + cellType2 + "\t" + corr + "\t" + (corr * corr) + "\t" + correlationMatrix.rawData[r2].length);


            }
        }
        tfOut.close();

    }

    private double[] toPrimitiveArr(Double[] toArray) {
        double[] arr = new double[toArray.length];
        for (int i = 0; i < toArray.length; i++) {
            arr[i] = toArray[i];
        }
        return arr;
    }
}
