/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ConvertCellTypeMatrixToText {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            String infile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-09-Groningen/CellTypeSpecificEQTLWoInteractionEffect/CellTypeSpecificityMatrix.binary";
            String outfile = "/Volumes/iSnackHD/Data/Projects/2013-06-CellTypeSpecificEQTLs/Datasets/2013-07-09-Groningen/CellTypeSpecificEQTLWoInteractionEffect/CellTypeSpecificityMatrix.txt";
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(infile);
            ds.save(outfile);
                    
        } catch (IOException e){
            
        }
    }
}
