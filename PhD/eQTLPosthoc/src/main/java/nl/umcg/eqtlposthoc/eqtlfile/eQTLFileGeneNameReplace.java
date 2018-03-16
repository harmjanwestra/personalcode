/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class eQTLFileGeneNameReplace {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            // cis
//            eQTLFileGeneNameReplace.run(
//                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/",
//                    "/Volumes/iSnackHD/SkyDrive/MetaAnalysisAnnotationFiles/2012-08-08-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt");
//            // trans
            // /Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR0.05.txt
            eQTLFileGeneNameReplace.run(
                    "/Volumes/iSnackHD/SkyDrive/latesteQTLs/transFile/",
                    "/Volumes/iSnackHD/SkyDrive/MetaAnalysisAnnotationFiles/2012-08-08-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt", true);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void run(String eQTLFileLoc, String ensemblLoc, boolean useEnsIDs) throws IOException {


        HashMap<String, String> probeToGene = new HashMap<String, String>();
        HashMap<String, String> probeToEnsembl = new HashMap<String, String>();
        TextFile tf = new TextFile(ensemblLoc, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String probe = elems[0];
            String hugo = elems[5];
            String ens = elems[4];
            probeToEnsembl.put(probe, ens);
            probeToGene.put(probe, hugo);
            elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();

        String[] files = Gpio.getListOfFiles(eQTLFileLoc, "txt");

        for (String f : files) {
            System.out.println("Parsing file: " + f);
            String efileIn = f;

            TextFile eQTLFileIn = new TextFile(efileIn, TextFile.R);

            String header = eQTLFileIn.readLine();

            String[] eqtlElems = eQTLFileIn.readLineElems(TextFile.tab);
            if (eqtlElems != null && eqtlElems.length > 16) {
                TextFile eQTLFileOut = new TextFile(efileIn + "-UpdatedEnsemblAnnotation.txt", TextFile.W);
                eQTLFileOut.writeln(header);
                while (eqtlElems != null) {
                    String probe = eqtlElems[4];
                    if (!useEnsIDs) {
                        String hugo = probeToGene.get(probe);
                        eqtlElems[eQTLTextFile.HUGO] = hugo;
                    } else {
                        String hugo = probeToEnsembl.get(probe);
                        eqtlElems[eQTLTextFile.HUGO] = hugo;
                    }
                    eQTLFileOut.writeln(Strings.concat(eqtlElems, Strings.tab));
                    eqtlElems = eQTLFileIn.readLineElems(TextFile.tab);
                }
                eQTLFileOut.close();

            }
            eQTLFileIn.close();
        }

    }
}
