/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ConvertMatrixProbeIdentifiersToEnsembl {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            String dsLoc = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/RawEndoPhenotypes/EGCUTEndophenotypesValidSamples-NoSex.txt-asLoadedByNormalizerEndophenotypeVsExpressionCorrelationMatrix.txt";
            String dsLocOut = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/RawEndoPhenotypes/EGCUTEndophenotypesValidSamples-NoSex.txt-asLoadedByNormalizerEndophenotypeVsExpressionCorrelationMatrix-Ensembl.txt";

            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            HashMap<String, String> ht12ToProbeId = new HashMap<String, String>();

            ProbeTranslation pbt = new ProbeTranslation();
            String ht12v3String = "HumanHT-12_V3_0_R2_11283641_A.txt";
            ht12ToProbeId = pbt.getProbeTranslation(probetranslationfile, ht12v3String, "Probe");

            TextFile tf = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-OnlyEnsemblAnnotation.txt", TextFile.R);
            HashMap<String, String> probeToEnsembl = (HashMap<String, String>) tf.readAsHashMap(0, 4);

            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(dsLoc);
            boolean[] hasId = new boolean[ds.nrRows];
            int nrWithEns = 0;

            HashMap<String, Integer> ensemblToNewRow = new HashMap<String, Integer>();
            ArrayList<String> newRowIds = new ArrayList<String>();
            for (int r = 0; r < ds.nrRows; r++) {
                String cov = ds.rowObjects.get(r);
                String probe = ht12ToProbeId.get(cov);
                String ens = probeToEnsembl.get(probe);
                for (int c = 0; c < ds.nrCols; c++) {
                    if (Double.isNaN(ds.rawData[r][c])) {
                        System.out.println(r + "\t" + c + "\t" + ds.rowObjects.get(r) + "\t" + ds.colObjects.get(c) + "\t is NaN\t" + ds.rawData[r][c]);
                    }
                }
                if (ens != null && ens.trim().length() > 0) {
                    hasId[r] = true;
                    String[] enselems = ens.split(" ");
                    ens = enselems[0];
                    if (!ensemblToNewRow.containsKey(ens)) {
                        ensemblToNewRow.put(ens, nrWithEns);
                        newRowIds.add(ens);
                        nrWithEns++;
                    }


                }
            }

            System.out.println(nrWithEns + " ensembl identifiers remain.");

            double[][] newData = new double[nrWithEns][ds.nrCols];
            int[] numMerged = new int[nrWithEns];

            int ctr = 0;
            // sum
            for (int r = 0; r < ds.nrRows; r++) {
                if (hasId[r]) {
                    String cov = ds.rowObjects.get(r);
                    String probe = ht12ToProbeId.get(cov);
                    String ens = probeToEnsembl.get(probe);
                    String[] enselems = ens.split(" ");
                    ens = enselems[0];
                    Integer rowId = ensemblToNewRow.get(ens);
                    for (int c = 0; c < ds.nrCols; c++) {

                        newData[rowId][c] += ds.rawData[r][c];
                        if (Double.isNaN(newData[rowId][c])) {
//                            System.out.println(r + "\t" + c + "\t is NaN\t" + ds.rawData[r][c]);
                        }
                    }
                    numMerged[rowId]++;
                }
            }

            // average 
            for (int r = 0; r < newData.length; r++) {
                for (int c = 0; c < newData[r].length; c++) {
                    newData[r][c] /= numMerged[r];

                }
            }

            DoubleMatrixDataset<String, String> dsOut = new DoubleMatrixDataset<String, String>();
            dsOut.rawData = newData;
            dsOut.colObjects = ds.colObjects;
            dsOut.rowObjects = newRowIds;
            dsOut.recalculateHashMaps();
            dsOut.save(dsLocOut);

            
            // check for NaNs
            CheckForNaNs.check(dsOut);
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
