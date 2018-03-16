/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class RewriteProbeAnnotation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String probetranslationfile = "/Volumes/iSnackHD/Data/ProbeMappingB37E70/tmp2/2013-07-18-ProbeAnnotationFile.txt";
            String source = "HT12v3.txt";
            String destination = "Gene";
            int probecolumn = 1;
            int annotationcolumn = 2;
           
            String infile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-12-04-2-NoGWASPValueFilter-BinomialDiseaseTestFDR0.05/ListOfEffectsForCrohnsDisease.txt";
            String outfile = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/MetaAnalysisResults/2013-12-04-2-NoGWASPValueFilter-BinomialDiseaseTestFDR0.05/ListOfEffectsForCrohnsDisease.txt-UpdatedGeneAnnotation.txt";

            ProbeTranslation pbt = new ProbeTranslation();
            HashMap<String, String> probeToAnnotation = pbt.getProbeTranslation(probetranslationfile, source, destination);

            TextFile tf = new TextFile(infile, TextFile.R);
            TextFile tf2 = new TextFile(outfile, TextFile.W);
            tf2.writeln(tf.readLine());

            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                String probe = elems[probecolumn];
                String annotation = probeToAnnotation.get(probe);
                if (annotation == null) {
                    annotation = "-";
                }
                elems[annotationcolumn] = annotation;

                tf2.writeln(Strings.concat(elems, Strings.tab));

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            tf2.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
