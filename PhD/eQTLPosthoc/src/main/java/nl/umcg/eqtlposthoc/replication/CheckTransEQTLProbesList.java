/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.replication;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.lang.String;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class CheckTransEQTLProbesList {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        /*
         * "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLProbesFDR0.05.txt",
         "/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/MonocyteReplication/SubsetOfProbesUsedDuringNormalization/",
         10,
         "/Volumes/iSnackHD/eQTLMetaAnalysis/BenFairFax/MonocyteReplication/SubsetOfProbesUsedDuringNormalization/FilteredForSignificantProbes/",
         "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt",
         "HumanHT-12_V4_0_R1_15002873_B.txt",
         ""
         */
        try{
        CheckTransEQTLProbesList.run(
                "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt",
                "Probe", 
                "HumanHT-12_V4_0_R1_15002873_B.txt", 
                "/Volumes/iSnackHD/SkyDrive/MissingProbes.txt", 
                "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/", 
                10);
        } catch (IOException e){
            e.printStackTrace();
        }
    }

    public static void run(String probetranslation, String source, String dest, String probesthataremissing, String indir, int perm) throws IOException {
        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> pt = pb.getProbeTranslation(probetranslation, source, dest);

        TextFile tf = new TextFile(probesthataremissing, TextFile.R);
        HashSet<String> probesToQuery = new HashSet<String>();
        probesToQuery.addAll(tf.readAsArrayList());
        tf.close();
        HashSet<String> probesThatAreNotMissing = new HashSet<String>();
        for (int p = 0; p < perm + 1; p++) {
            String infile = "eQTLs.txt";
            if (p > 0) {
                infile = "PermutedEQTLsPermutationRound" + p + ".txt.gz";
                //PermutedEQTLsPermutationRound1.txt.gz
            }
            
            System.out.println(indir + infile);

            TextFile tfe = new TextFile(indir + infile, TextFile.R);
            String[] elems = tfe.readLineElems(TextFile.tab);
            elems = tfe.readLineElems(TextFile.tab);

            TextFile out = null;
            
            if(p==0){
                out = new TextFile("/Volumes/iSnackHD/SkyDrive/ProbesThatWereMissingButHaveTrans.txt", TextFile.W);
            }
            
            while (elems != null) {
                String meta = elems[4];
                String probe = pt.get(meta);
                if (probesToQuery.contains(probe)) {
                    if (p == 0) {
                        out.writeln(Strings.concat(elems, Strings.tab));
                    }
                }
                probesThatAreNotMissing.add(probe);
                elems = tfe.readLineElems(TextFile.tab);
            }
            if(p==0){
                out.close();
            }
            tfe.close();
        }

        System.out.println("");
        System.out.println("Probes that are not in the trans-eQTL data");
        for (String p : probesToQuery) {
            if (!probesThatAreNotMissing.contains(p)) {
                System.out.println(p);
            }
        }
    }
}
