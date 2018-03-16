/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CellTypeSpecificVectorProbeRewrite {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            String probetranslationfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            String src = "HumanHT-12_V3_0_R2_11283641_A.txt";
            String dst = "Probe";
            String infile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/2013-06-21-EGCUT-Vector-rs12057769-2000128.txt";
            String outfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/CellTypeSpecificTestEQTLOutput/2013-06-21-EGCUT-Vector-rs12057769-2000128-Ensembl.txt";
            String ensemblfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-OnlyEnsemblAnnotation.txt";

            ProbeTranslation pbt = new ProbeTranslation();

            HashMap<String, String> ht12toprobe = new HashMap<String, String>();
            HashMap<String, String> probetoEns = new HashMap<String, String>();
            ht12toprobe = pbt.getProbeTranslation(probetranslationfile, src, dst);

            TextFile tfens = new TextFile(ensemblfile, TextFile.R);
            String[] elemsens = tfens.readLineElems(TextFile.tab);
            while (elemsens != null) {
                if (elemsens.length >= 5) {
                    String probe = elemsens[0];
                    String ens = elemsens[4];
                    probetoEns.put(probe, ens);
                }
                elemsens = tfens.readLineElems(TextFile.tab);
            }

            tfens.close();

            TextFile tf = new TextFile(infile, TextFile.R);
            TextFile tfout = new TextFile(outfile, TextFile.W);
            tfout.writeln(tf.readLine());
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                String ht12probe = elems[0];
                String probe = ht12toprobe.get(ht12probe);
                String ens = probetoEns.get(probe);

                if (ens != null) {
                    elems[0] = ens;
                    String output = Strings.concat(elems, Strings.tab);
                    tfout.writeln(output);
                }



                elems = tf.readLineElems(TextFile.tab);
            }

            tfout.close();
            tf.close();




        } catch (IOException e) {
        }

    }
}
