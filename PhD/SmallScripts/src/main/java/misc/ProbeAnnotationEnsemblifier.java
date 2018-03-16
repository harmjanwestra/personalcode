/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ProbeAnnotationEnsemblifier {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String ensemblAnnotation = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/All.0.96%Identity-Merged-PerProbe-UniqueMappings-Ensembl70.txt";
        String ht12v3reannotated = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-HT12v3.txt";
        String outfile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-HT12v3-Ensembl70.txt";
        try {
            TextFile tf = new TextFile(ensemblAnnotation, TextFile.R);
            HashMap<String, String> probeToEns = new HashMap<String, String>();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                String probe = elems[0];
                String ens = elems[6];
                probeToEns.put(probe, ens);

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile tfout = new TextFile(outfile, TextFile.W);
            TextFile tfin = new TextFile(ht12v3reannotated, TextFile.R);
            tfout.writeln(tfin.readLine());
            elems = tfin.readLineElems(TextFile.tab);
            while (elems != null) {
                String probe = elems[6];
                String ens = probeToEns.get(probe);
                if (ens == null) {
                    ens = "-";
                }
                elems[2] = ens;
                tfout.writeln(Strings.concat(elems, Strings.tab));
                elems = tfin.readLineElems(TextFile.tab);
            }

            tfin.close();
            tfout.close();
// 6
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
