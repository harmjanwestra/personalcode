/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class DoubleMatrixDatasetSplitter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String in = "/Volumes/iSnackHD/DataForYangLi/Choy/AffymetrixExpressionDataOriginal.txt.TriTyperFormat.txt";
        String out = "/Volumes/iSnackHD/DataForYangLi/Choy/AffymetrixExpressionDataOriginal.txt.TriTyperFormat-CEU.txt";
        String f = "/Volumes/iSnackHD/DataForYangLi/Choy/GenotypeToExpressionIdOriginal-CEUWoDuplicates.txt";


        HashSet<String> samplesToInclude = new HashSet<String>();
        try {
            TextFile tf = new TextFile(f, TextFile.R);

            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                samplesToInclude.add(elems[1]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(in, null, samplesToInclude);
            ds.save(out);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
