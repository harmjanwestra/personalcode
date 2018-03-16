/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import eqtlmappingpipeline.binarymeta.meta.graphics.ZScorePlot;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.CellTypeSpecificMetaAnalsysisDatasetContainer;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class ComparisonBetweenDatasets {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here

            String[] datasetFileDirs = new String[]{
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-12-GroningenHT12v3/CellTypeSpecificityMatrix.binary",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-23_SHIP-TREND/CellTypeSpecificityMatrix.binary",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-24_ROTTERDAM-STUDY/CellTypeSpecificityMatrix.binary",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-25-EGCUT/CellTypeSpecificityMatrix.binary",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-07-29-DILGOM/CellTypeSpecificityMatrix.binary",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-08-19-KORA-F4/CellTypeSpecificityMatrix.binary",
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/CellTypeSpecificEQTLMapping/Results/2013-08-19-INCHIANTI/CellTypeSpecificityMatrix.binary",
                "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysisZScoreMatrix.binary"
            };

            String[] datasetNames = new String[]{
                "GroningenHT12v3",
                "SHIP-TREND",
                "Rotterdam",
                "EGCUT",
                "DILGOM",
                "KORA-F4",
                "INCHIANTI",
                "METAANALYSIS"
            };

//            boolean[] platformIsHT12v3 = new boolean[]{true, true, false};

            boolean[] platformIsHT12v3 = new boolean[]{
                true,
                true,
                false,
                true,
                true, true, true, true
            };

            String output = "Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-10-03-ComparisonsBetweenDatasets/";
            String pbt = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            ComparisonBetweenDatasets c = new ComparisonBetweenDatasets();

            String query = "";
            c.run(datasetFileDirs, platformIsHT12v3, output, pbt, datasetNames, query);
        } catch (IOException ex) {
            Logger.getLogger(ComparisonBetweenDatasets.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void run(String[] datasetFileDirs, boolean[] platformIsHT12v3, String output, String probetranslationfile, String[] datasetNames, String query) throws IOException {
        Gpio.createDir(output);

        String ht12v3String = "HumanHT-12_V3_0_R2_11283641_A.txt";
        String ht12v4String = "HumanHT-12_V4_0_R1_15002873_B.txt";

        String hugoStr = "HUGO";

        //       PROBE LOOKUP TABLES..
        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> ht12v4ToHT12v3 = new HashMap<String, String>();
        HashMap<String, String> ht12v3ToHugo = new HashMap<String, String>();
        ht12v4ToHT12v3 = pbt.getProbeTranslation(probetranslationfile, ht12v4String, ht12v3String);
        ht12v3ToHugo = pbt.getProbeTranslation(probetranslationfile, ht12v3String, hugoStr);

        // load the datasets..
        CellTypeSpecificMetaAnalsysisDatasetContainer container = new CellTypeSpecificMetaAnalsysisDatasetContainer(datasetFileDirs, platformIsHT12v3, ht12v4ToHT12v3);
        SNP[][] snps = container.getSnps();

        ArrayList<DoubleMatrixDataset<String, String>> datasets = container.getDatasets();
        HashSet<String> allSNPs = container.getAllSNPs();
        HashSet<String> allEQTL = container.getAllEQTL();


        // hash the SNPs for easy lookup..
        Integer[][] snpIndex = container.getSnpIndex(); // for each SNP in the meta analysis, this index points to a SNP in a dataset (in the jagged allSNPs array)
        HashMap<String, Integer> snpToId = container.getSnpToId();

        // hash the eQTLs as well..
        Integer[][] eQTLIndex = container.geteQTLIndex();
        HashMap<String, Integer> eQTLToId = container.geteQTLToId();
        String[] allEQTLArr = allEQTL.toArray(new String[0]);

        // hash the probes / covariates...
        HashSet<String> allCovariates = container.getAllCovariates();
        String[] allCovariateArr = container.getAllCovariateArr();
        Integer[][] covariateIndex = container.getCovariateIndex();
        HashMap<String, Integer> covariateToId = container.getCovariateToId();


        // determine which eQTLs to FLIP
        Boolean[][] flipFxPerEQTL = new Boolean[datasetFileDirs.length][allEQTLArr.length];
        for (int e = 0; e < allEQTLArr.length; e++) {
            String eQTL = allEQTLArr[e];
            String snp = eQTL.split("-")[0];
            String probe = eQTL.split("-")[1];
            Integer snpId = snpToId.get(snp);

            if (snpId == null) {
                System.err.println("ERROR snpId " + snp + " is null for eQTL: " + eQTL);
            }

            String firstSNPAlleles = null;
            String firstSNPAssessedAllele = null;
            String firstSNPName = null;
            Boolean[] flipEffect = new Boolean[datasets.size()];

            String mafStr = "";
            String crStr = "";
            String nrCalledStr = "";
            String hwepStr = "";

            String dsNames = "";

            // determine which effects should be flipped
            for (int d = 0; d < datasets.size(); d++) {
                Integer snpIdInDataset = snpIndex[d][snpId];
                if (snpIdInDataset != null) {
                    SNP snpObj = snps[d][snpIdInDataset];
                    String SNPAlleles = BaseAnnot.getAllelesDescription(snpObj.getAlleles());
                    String SNPAssessedAllele = BaseAnnot.toString(snpObj.getMinorAllele());
                    if (firstSNPAlleles == null) {
                        firstSNPAlleles = SNPAlleles;
                        firstSNPAssessedAllele = SNPAssessedAllele;
                        flipEffect[d] = false;
                        flipFxPerEQTL[d][e] = false;
                        firstSNPName = snpObj.getName();
                    } else {
                        flipFxPerEQTL[d][e] = BaseAnnot.flipalleles(firstSNPAlleles, firstSNPAssessedAllele, SNPAlleles, SNPAssessedAllele);
                        String secondSNPName = snpObj.getName();
                        if (!firstSNPName.equals(secondSNPName)) {
                            System.err.println("SNPs not linked together properly");
                        }
                    }
                } else {

                    flipFxPerEQTL[d][e] = null;
                }
            }
        }

        SpearmansCorrelation spearman = new SpearmansCorrelation();
        int nrDs = 0;
        for (int d = 0; d < datasets.size(); d++) {
            for (int d2 = d + 1; d2 < datasets.size(); d2++) {
                nrDs++;
            }
        }

        DoubleMatrixDataset<String, String> outputData = new DoubleMatrixDataset<String, String>(allCovariateArr.length, nrDs);
        ArrayList<String> colNames = new ArrayList<String>();
        DecimalFormat formatter = new DecimalFormat("#.###");
        
        ZScorePlot zs = new ZScorePlot();
        zs.init(nrDs, datasetNames, true, output+query+"-ComparisonChart.pdf");
        
        for (int d = 0; d < datasets.size(); d++) {
            for (int d2 = d + 1; d2 < datasets.size(); d2++) {

                ProgressBar pb = new ProgressBar(allCovariateArr.length, datasetNames[d] + "-" + datasetNames[d2]);

                for (int cov = 0; cov < allCovariateArr.length; cov++) {
                    Integer covIndex1 = covariateIndex[d][cov];
                    Integer covIndex2 = covariateIndex[d2][cov];
                    if (covIndex1 != null && covIndex2 != null) {
                        // get all the values over all shared eQTLs..
                        ArrayList<Double> xArr = new ArrayList<Double>();
                        ArrayList<Double> yArr = new ArrayList<Double>();
                        for (int e = 0; e < allEQTLArr.length; e++) {
                            Integer eQTLIndex1 = eQTLIndex[d][e];
                            Integer eQTLIndex2 = eQTLIndex[d2][e];
                            if (eQTLIndex1 != null && eQTLIndex2 != null) {
                                double v1 = datasets.get(d).get(covIndex1, eQTLIndex1);
                                double v2 = datasets.get(d2).get(covIndex2, eQTLIndex2);
                                if (!Double.isNaN(v1) && !Double.isNaN(v2)) {
                                    if (flipFxPerEQTL[d][e]) {
                                        v1 *= -1;
                                    }
                                    if (flipFxPerEQTL[d2][e]) {
                                        v2 *= -1;
                                    }
                                    xArr.add(v1);
                                    yArr.add(v2);
                                }
                            }
                        }

                        double[] x = Primitives.toPrimitiveArr(xArr.toArray(new Double[0]));
                        double[] y = Primitives.toPrimitiveArr(yArr.toArray(new Double[0]));
                        //equals("CellTypeInteractionZScore") || ds.rowObjects.get(r).equals("MainEffectZScore") || ds.rowObjects.get(r).equals("CellTypeZScore") || ds.rowObjects.get(r).equals("CellTypeSNPZScore")
                        double corr = spearman.correlation(x, y);
                        if (allCovariateArr[cov].equals(query)) {
                            for(int i=0;i<x.length;i++){
                                zs.draw(x[i], y[i], d, d2);
                            }
//                            ScatterPlot p = new ScatterPlot(500, 500, x, y, ScatterPlot.OUTPUTFORMAT.PDF, output + allCovariateArr[cov] + "-" + datasetNames[d] + "-" + datasetNames[d2] + "-" + formatter.format(corr) + ".pdf");
                        }


                        outputData.rawData[cov][colNames.size()] = corr;
                    } else {
                        outputData.rawData[cov][colNames.size()] = Double.NaN;
                    }
                    pb.iterate();
                }
                pb.close();
                colNames.add(datasetNames[d] + "-" + datasetNames[d2]);
            }

        }
        try {
            zs.write(output+query+"-ComparisonChart.pdf");
        } catch (Exception ex) {
            Logger.getLogger(ComparisonBetweenDatasets.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        outputData.colObjects = colNames;
        ArrayList<String> covariatesplus = new ArrayList<String>();
        for (int cov = 0; cov < allCovariateArr.length; cov++) {
            String hugo = ht12v3ToHugo.get(allCovariateArr[cov]);
            if (hugo != null && !hugo.equals("-")) {
                covariatesplus.add(allCovariateArr[cov] + ";" + hugo);
            } else {
                covariatesplus.add(allCovariateArr[cov]);
            }
        }
        outputData.rowObjects = covariatesplus;
        outputData.recalculateHashMaps();
        outputData.save(output + "CorrelationOfCovariatesBetweenDatasets.txt");

    }
}
