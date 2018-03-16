/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class SuperDuperCombiner {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String f1 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/CorrelationsWithPC1-MetaZ.txt";
            String f2 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MED8Validation/2013-10-09-MED8CovariateOverlap-PC1OnlyMerged/PC1DirectionCorrectedCellTypeInteractionZScore-CellTypeInteractionZScore.txt";
            TextFile tf = new TextFile(f1, TextFile.R);
            HashMap<String, String> probeToCorr = new HashMap<String, String>();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                probeToCorr.put(elems[0], elems[elems.length - 2]+"\t"+elems[elems.length - 1]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile tf2 = new TextFile(f2, TextFile.R);
            TextFile tf2w = new TextFile(f2 + "-Appended.txt", TextFile.W);
            tf2w.writeln(tf2.readLine() + "\tPC1MetaZ\tPC1IndividualZScores");
            String line = tf2.readLine();
            while (line != null) {

                String[] lineElems = line.split("\t");
                String probe = lineElems[0];
                String corr = probeToCorr.get(probe);
                line+="\t"+corr;
                tf2w.writeln(line);
                
                line = tf2.readLine();
            }
            tf2w.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
