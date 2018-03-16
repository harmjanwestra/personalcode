/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ColumnCounter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String f = "/Volumes/iSnackHD/tmp/chr22_top10Line.dose";
            TextFile tf = new TextFile(f, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                System.out.println(elems.length);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
