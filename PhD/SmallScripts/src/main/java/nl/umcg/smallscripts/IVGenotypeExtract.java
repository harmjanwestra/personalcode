/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.trityper.converters.TriTyperToMachImputedTransposed;

/**
 *
 * @author harmjan
 */
public class IVGenotypeExtract {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            TriTyperToMachImputedTransposed.convert("/Data/GeneticalGenomicsDatasets/BloodHT12ImputeTriTyper/", "/Data/GeneticalGenomicsDatasets/BloodHT12ImputeTriTyper/2013-01-09-IVAnalysisGenotypes.txt", "/Volumes/iSnackHD/Data/Projects/2013-IVAnalysis/ProbeFiles/cis_eQTLs_08012013_dummy_snpsOnly.txt");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
