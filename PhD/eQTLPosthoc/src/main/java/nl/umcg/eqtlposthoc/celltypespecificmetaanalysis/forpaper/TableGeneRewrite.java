/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class TableGeneRewrite {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String probeTranslationFile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-18-ProbeAnnotationFile.txt";
        String tableToConvert = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/2013-12-NatureMethods-Draft2/Supplementary/Supplementary Table 2 - Results of cell type specificity analysis.txt-WithR.txt";
        String tableOut = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/2013-12-NatureMethods-Draft2/Supplementary/Supplementary Table 2 - Results of cell type specificity analysis.txt-WithR-UpdatedGeneNames.txt";
        int probeCol = 4;
        int geneCol = 7;

        try {

            ProbeTranslation pb = new ProbeTranslation();
            HashMap<String, String> probeToHUGO = pb.getProbeTranslation(probeTranslationFile, "HT12v3.txt", "Gene");

            TextFile tf = new TextFile(tableToConvert, TextFile.R);
            TextFile tfout = new TextFile(tableOut, TextFile.W);

            tfout.writeln(tf.readLine());
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String probe = elems[probeCol];
                String gene = elems[geneCol];
                String assocGene = probeToHUGO.get(probe);
                elems[geneCol] = assocGene;
                tfout.writeln(Strings.concat(elems, Strings.tab));
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            tfout.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
