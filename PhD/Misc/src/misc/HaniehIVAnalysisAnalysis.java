/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class HaniehIVAnalysisAnalysis {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {


            HashSet<String> eProbes = new HashSet<String>();
            TextFile efile = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/PostQC/eQTLProbesFDR0.05.txt-HT12v3Ids.txt", TextFile.R);
            String[] efileElems = efile.readLineElems(TextFile.tab);
            efileElems = efile.readLineElems(TextFile.tab);
            while (efileElems != null) {
                String probe = efileElems[4];
                eProbes.add(probe);
                efileElems = efile.readLineElems(TextFile.tab);
            }

            efile.close();

            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/IVGeneral/2012-07-31-IV/trans_ma_results_31072012.txt", TextFile.R);

            String[] elems = tf.readLineElems(TextFile.tab);
            elems = tf.readLineElems(TextFile.tab);
            HashMap<String, HashSet<String>> probesPerClass = new HashMap<String, HashSet<String>>();
            HashSet<String> uniqueProbes = new HashSet<String>();
            ArrayList<String> classes = new ArrayList<String>();
            while (elems != null) {

                String classOfIV = elems[elems.length - 1];
                String probe = elems[1];
                uniqueProbes.add(probe);
                HashSet<String> probesForClass = probesPerClass.get(classOfIV);
                if (probesForClass == null) {
                    probesForClass = new HashSet<String>();
                    classes.add(classOfIV);
                }

                probesForClass.add(probe);

                probesPerClass.put(classOfIV, probesForClass);

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            System.out.println(uniqueProbes.size() + " unique probes");
            String cls1 = classes.get(0);
            HashSet<String> probes1 = probesPerClass.get(cls1);
            for (String cls2 : classes) {
                HashSet<String> probes2 = probesPerClass.get(cls2);
                int ctr2 = 0;
                for (String p : probes1) {
                    if (probes2.contains(p)) {
                        ctr2++;
                    }
                }
                System.out.println(cls1 + "\t" + cls2 + "\t" + ctr2);
            }

            System.out.println("");

            for (String cls : classes) {
                HashSet<String> probes = probesPerClass.get(cls);
                int ctr = 0;
                for (String p : probes) {
                    if (eProbes.contains(p)) {
                        ctr++;
                    }
                }

                double perc = (double) ctr / probes.size();
                System.out.println(cls + "\t" + probes.size() + "\t" + ctr + "\t" + perc);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
