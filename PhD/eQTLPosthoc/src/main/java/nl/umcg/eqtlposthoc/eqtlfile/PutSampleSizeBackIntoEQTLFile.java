/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class PutSampleSizeBackIntoEQTLFile {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        
        // EGCUT;SHIP_TREND;Groningen-HT12;Groningen-H8v2;Rotterdam;DILGOM;INCHIANTI;HVH-v3;HVH-v4
        String[] datasetnames = new String[]{"EGCUT", "SHIP-Trend", "Fehrmann HT12", "Fehrmann H8v2", "Rotterdam Study", "DILGOM", "InCHIANTI", "HVH HT12v3", "HVH HT12v4", "Meta-analysis", "B-Cells", "Monocytes", "Kora", "BSGS"};
        int[] datasetweights = new int[]{891, 963, 1240, 229, 762, 509, 611, 43, 63, 5311, 10, 10, 10, 10};
        try {
            TextFile tf    = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/DataForRebuttal/2013-03-Trans40PCsCisFXNotRegressedOut/FDR0.5SNPProbes/eQTLs.txt", TextFile.R);
            TextFile tfout = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/DataForRebuttal/2013-03-Trans40PCsCisFXNotRegressedOut/FDR0.5SNPProbes/eQTLs.txt-WSampleSize.txt.gz", TextFile.W);
            
            tfout.writeln(tf.readLine());
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                String dsNames = elems[12];
                String[] dsElems = dsNames.split(";");

                String[] sampleSElems = new String[dsElems.length];
                for (int i = 0; i < dsElems.length; i++) {
                    if (dsElems[i].equals("-")) {
                        sampleSElems[i] = "-";
                    } else {
                        sampleSElems[i] = "" + datasetweights[i];
                    }
                }

                String sampleSizes = Strings.concat(sampleSElems, Strings.semicolon);

                elems[13] = sampleSizes;

                tfout.writeln(Strings.concat(elems, Strings.tab));
                elems = tf.readLineElems(TextFile.tab);
            }
            tfout.close();
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
