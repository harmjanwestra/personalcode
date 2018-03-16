/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ConvertBinaryToText {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix.binary");
            
            ds.save("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/CellTypeSpecificityMatrix.txt");
        } catch (IOException e){
            e.printStackTrace();
        }
    }
}
