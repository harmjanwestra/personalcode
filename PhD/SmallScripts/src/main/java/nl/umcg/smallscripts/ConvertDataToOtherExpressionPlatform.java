/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ConvertDataToOtherExpressionPlatform {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {
            String probet = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt";
            String fileIn = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataYRI/ExpressionData.txt.gz.10PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt.gz";
            String fileOut = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap3r2/ExpDataYRI/ExpressionData.txt.gz.10PCAsOverSamplesRemoved-GeneticVectorsNotRemovedHT12v3.txt.gz";
            String source = "HumanWG-6_V2_0_R4_11223189_A.txt";
            String destin = "HumanHT-12_V3_0_R2_11283641_A.txt";


            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(fileIn);

            ProbeTranslation pb = new ProbeTranslation();
            HashMap<String, String> ht12v4toHT12v3 = pb.getProbeTranslation(probet, source, destin);

            ArrayList<String> newRows = new ArrayList<String>();
            double[][] rawData = ds.rawData;
            boolean[] includerow = new boolean[rawData.length];

            for (int row = 0; row < rawData.length; row++) {
                String probe = ds.rowObjects.get(row);



                String ht12v3 = ht12v4toHT12v3.get(probe);
                if (ht12v3 == null) {
                    // try complement?
                    String complement = BaseAnnot.getComplement(probe);
                    ht12v3 = ht12v4toHT12v3.get(complement);
                    if (ht12v3 == null) {
                        complement = BaseAnnot.getReverseComplement(probe);
                        ht12v3 = ht12v4toHT12v3.get(complement);
                    }
                }
                if (ht12v3 == null || ht12v3.equals("-")) {
                    includerow[row] = false;
                } else {
                    newRows.add(ht12v3);
                    includerow[row] = true;

                }
            }

            double[][] newRawData = new double[newRows.size()][0];
            int ctr = 0;
            for (int row = 0; row < rawData.length; row++) {
                if (includerow[row]) {
                    newRawData[ctr] = rawData[row];
                    ctr++;
                }
            }

            System.out.println(ctr+" new rows. Was: "+includerow.length);

            DoubleMatrixDataset<String, String> dsOut = new DoubleMatrixDataset<String, String>(newRawData, newRows, ds.colObjects);
            dsOut.save(fileOut);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
