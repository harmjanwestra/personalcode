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
public class AddGeneName {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String pbtfile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-18-ProbeAnnotationFile.txt";
        String infile = "/Volumes/iSnackHD/Data/Projects/LudeFranke/2014-02-12-TransInteractionTerms/Meta/MetaAnalysisZScoreMatrix.txt";
        String outfile = "/Volumes/iSnackHD/Data/Projects/LudeFranke/2014-02-12-TransInteractionTerms/Meta/EtLeMatrixZScore-AvecLeGenes.txt";
        
        HashMap<String, String> probeToGene = new HashMap<String, String>();
        ProbeTranslation pbt = new ProbeTranslation();
        try {
            probeToGene = pbt.getProbeTranslation(pbtfile, "HT12v3.txt", "Gene");
            TextFile tf = new TextFile(infile, TextFile.R);
            TextFile tfout = new TextFile(outfile, TextFile.W);
            tfout.writeln(tf.readLine()+ "\tGene");
            String[] elems = tf.readLineElems(TextFile.tab);
            while(elems!=null){
                String probe = elems[0];
                String gene = probeToGene.get(probe);
                tfout.writeln(Strings.concat(elems, Strings.tab)+ "\t"+gene);
                elems = tf.readLineElems(TextFile.tab);
            }
            tfout.close();
            tf.close();
        } catch (IOException e) {

        }
    }

}
