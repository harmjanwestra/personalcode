/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.trityper.converters.TriTyperToVCF;

/**
 *
 * @author harmjan
 */
public class ConvertToVCF {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            TriTyperToVCF v = new TriTyperToVCF();
            v.convert("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap2r24-CEU/", "/Volumes/iSnackHD/tmp/VCFOUT/", null);
        } catch (IOException e){
            e.printStackTrace();
        }
    }
}
