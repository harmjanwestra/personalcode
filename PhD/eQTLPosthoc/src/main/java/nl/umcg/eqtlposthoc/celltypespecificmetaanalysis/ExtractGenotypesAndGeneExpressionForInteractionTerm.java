/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ExtractGenotypesAndGeneExpressionForInteractionTerm {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String snp = "rs7135577";
        String probe1 = "3450180";
        String probe2 = "7040035";

        String genotypeData = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
        String rawExpData = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-08-Groningen/Normalization/ExpressionData/ExpressionDataSamplePCQC-QNormLog2Transform.CovariatesRemoved.txt.gz";
        String pcCorrectedData = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.40PCAsOverSamplesRemoved.txt.gz";

        ExtractGenotypesAndGeneExpressionForInteractionTerm f = new ExtractGenotypesAndGeneExpressionForInteractionTerm();
        try {
            f.run(snp, probe1, probe2, genotypeData, rawExpData, pcCorrectedData);
        } catch (IOException ex) {
            Logger.getLogger(ExtractGenotypesAndGeneExpressionForInteractionTerm.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void run(String snp, String probe1, String probe2, String genotypeDataLoc, String expressionData1Loc, String expressionData2Loc) throws IOException {
        TriTyperGenotypeData genotypeData = new TriTyperGenotypeData(genotypeDataLoc);
        SNPLoader loader = genotypeData.createSNPLoader();

        Integer snpId = genotypeData.getSnpToSNPId().get(snp);
        SNP snpObj = genotypeData.getSNPObject(snpId);
        loader.loadGenotypes(snpObj);
        loader.loadDosage(snpObj);

        HashSet<String> probes = new HashSet<String>();
        probes.add(probe1);
        probes.add(probe2);

        DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>(expressionData1Loc, probes);
        DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(expressionData2Loc, probes);

        System.out.println("Individual\tGenotype\tDosage\t" + ds1.rowObjects.get(0) + "\t" + ds1.rowObjects.get(1) + "\t" + ds2.rowObjects.get(0) + "\t" + ds2.rowObjects.get(1));
        for (int i = 0; i < genotypeData.getIndividuals().length; i++) {
            String individual = genotypeData.getIndividuals()[i];
            if (genotypeData.getIsIncluded()[i]) {
                Integer sampleIdExp1 = ds1.hashCols.get(individual);
                Integer sampleIdExp2 = ds2.hashCols.get(individual);

                if (sampleIdExp1 != null && sampleIdExp2 != null) {
                    System.out.println(genotypeData.getIndividuals()[i] + "\t" + snpObj.getGenotypes()[i] + "\t" + snpObj.getDosageValues()[i] + "\t" + ds1.rawData[0][sampleIdExp1] + "\t" + ds1.rawData[1][sampleIdExp1] + "\t" + ds2.rawData[0][sampleIdExp2] + "\t" + ds2.rawData[1][sampleIdExp2]);
                }

            }
        }


        snpObj.clearGenotypes();
    }
}
