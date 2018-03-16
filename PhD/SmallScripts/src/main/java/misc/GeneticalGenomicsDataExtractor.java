/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class GeneticalGenomicsDataExtractor {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String[] snps = {"rs10876864"};
        String[] probes = {""};

        String expression = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt";
        String genotype = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
        String covariate = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/2013-10-10-BloodHT12DataForInteractionTerms/Covariate.txt-DirectionCorrected.txt";

        try {
            umcg.genetica.io.trityper.TriTyperGenotypeData ds = new umcg.genetica.io.trityper.TriTyperGenotypeData(genotype);
            DoubleMatrixDataset<String, String> dsexp = new DoubleMatrixDataset<String, String>(expression);
            DoubleMatrixDataset<String, String> dscov = new DoubleMatrixDataset<String, String>(covariate);
            dscov.transposeDataset();

            for (int p = 0; p < snps.length; p++) {
                String snp = snps[p];
                String probe = probes[p];
                Integer snpid = ds.getSnpToSNPId().get(snp);
                if (snpid != null) {
                    String outputFile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/CisEffectsAndArrayQuality/" + snp + "-" + probe + "-data.txt";
                    TextFile outtf = new TextFile(outputFile, TextFile.W);
                    outtf.writeln("individual\tgenotype\texpression\tcovariate");
                    SNPLoader loader = ds.createSNPLoader();
                    Integer expProbeId = dsexp.hashRows.get(probe);
                    SNP snpObj = ds.getSNPObject(snpid);
                    loader.loadGenotypes(snpObj);
                    loader.loadDosage(snpObj);
                    double[] dosages = snpObj.getDosageValues();
                    String[] genoIndividuals = ds.getIndividuals();
                    for (int i = 0; i < genoIndividuals.length; i++) {
                        String individual = genoIndividuals[i];
                        String output = individual;

                        output += "\t" + dosages[i];

                        Integer expId = dsexp.hashCols.get(individual);
                        Integer covId = dscov.hashCols.get(individual);

                        if (expId != null && covId != null) {
                            output += "\t" + dsexp.rawData[expProbeId][expId];
                            output += "\t" + dscov.rawData[0][covId];
                            outtf.writeln(output);
                        }
                    }
                    snpObj.clearGenotypes();
                    outtf.close();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
