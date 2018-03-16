/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author harmjan
 */
public class JosephPowellDataValidator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String infile = "";
            TextFile tf = new TextFile(infile, TextFile.R);
            int nrLines = tf.countLines();
            tf.close();
            tf.open();
            int nrCols = tf.countCols(TextFile.tab);
            tf.close();
            tf.open();
            tf.readLine(); // header
            String[] elems = tf.readLineElems(TextFile.tab);

            double[][] data = new double[nrCols - 2][nrLines - 1];

            int lnctr = 0;
            while (elems != null) {
                // probes are on the columns.. read in transposed
                for (int e = 2; e < elems.length; e++) {
                    data[e - 2][lnctr] = Double.parseDouble(elems[e]);
                }
                lnctr++;
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            tf.close();

            for (int l = 0; l < data.length; l++) {
                double mean = Descriptives.mean(data[l]);
                if (mean < 0) {
                    System.out.println(l + "\t" + mean);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
