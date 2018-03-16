/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.replication;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class FilterForSignificantProbes {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        FilterForSignificantProbes p = new FilterForSignificantProbes();
        try {

            p.run(
                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLProbesFDR0.05.txt",
                    "/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/MonocyteReplication/SubsetOfProbesUsedDuringNormalization/",
                    10,
                    "/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/MonocyteReplication/SubsetOfProbesUsedDuringNormalization/FilteredForSignificantProbes/",
                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt",
                    "HumanHT-12_V4_0_R1_15002873_B.txt",
                    "Probe");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String significantEQTLFile, String indir, int perm, String outdir, String probeTranslation, String source, String dest) throws IOException {

        // read probe annotation
        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> probeToMetaId = pb.getProbeTranslation(probeTranslation, source, dest);

        System.out.println(probeToMetaId.size() + " Probe translations read");
        HashSet<String> significantProbes = new HashSet<String>();
        TextFile tf = new TextFile(significantEQTLFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            significantProbes.add(elems[4]);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        System.out.println(significantProbes.size() + " significant probes");

        Gpio.createDir(outdir);
        for (int p = 0; p < perm + 1; p++) {
            String infile = "eQTLs.txt";
            if (p > 0) {
                infile = "PermutedEQTLsPermutationRound" + p + ".txt.gz";
                //PermutedEQTLsPermutationRound1.txt.gz
            }

            TextFile tf2 = new TextFile(indir + infile, TextFile.R);

            TextFile tfout = new TextFile(outdir + infile, TextFile.W);

            System.out.println("Processing " + indir + infile);
            System.out.println("Output: " + outdir + infile);

            tfout.writeln(tf2.readLine());

            elems = tf2.readLineElems(TextFile.tab);
            while (elems != null) {
                String probe = elems[4];
                String metaId = probeToMetaId.get(probe);
                if (significantProbes.contains(metaId)) {
                    elems[4] = metaId;
                    String output = Strings.concat(elems, Strings.tab);
                    tfout.writeln(output);
                }

                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();
            tfout.close();


        }
    }
}
