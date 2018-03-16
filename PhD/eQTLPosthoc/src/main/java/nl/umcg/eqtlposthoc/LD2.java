/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc;

import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class LD2 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            TriTyperGenotypeData ds = new TriTyperGenotypeData();
            ds.load("/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/");

            SNPLoader loader = ds.createSNPLoader();
            DetermineLD ld = new DetermineLD();

            Integer id1 = ds.getSnpToSNPId().get("rs8192916");
            Integer id2 = ds.getSnpToSNPId().get("rs2236337");

            SNP snp1 = ds.getSNPObject(id1);
            SNP snp2 = ds.getSNPObject(id2);
            loader.loadGenotypes(snp1);
            loader.loadGenotypes(snp2);
            System.out.println(ld.getRSquared(snp1, snp2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false));
            System.out.println(ld.getRSquared(snp1, snp2, ds, DetermineLD.RETURN_D_PRIME, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false));



        } catch (Exception e) {
        }
    }
}
