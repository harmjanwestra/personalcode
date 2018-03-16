/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class BatchiFier {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            TextFile tf = new TextFile("/Volumes/iSnackHD/Batches.txt", TextFile.R);
            TextFile tfout = new TextFile("/Volumes/iSnackHD/BatchesUnique.txt", TextFile.W);
            String[] elems = tf.readLineElems(TextFile.space);
            int sum = 0;
            while (elems != null) {

                HashSet<String> samples = new HashSet<String>();
                for (int i = 1; i < elems.length; i++) {
                    samples.add(elems[i]);
                }
                sum += samples.size();
                System.out.println(elems[0] + "\t" + samples.size());
                tfout.writeln(elems[0] + "\t" + Strings.concat(samples.toArray(new String[0]), Strings.tab));
                elems = tf.readLineElems(TextFile.space);
            }
            tf.close();

            tfout.close();
            System.out.println(sum);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
