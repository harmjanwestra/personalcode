/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import eqtlmappingpipeline.meta.posthocanalysis.eqtlfile.eQTLFileProbeNameReplace;
import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class eQTLToEffectSize {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        String filterfor = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt";
        String infile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/DataForRebuttal/2013-03-Trans40PCsCisFXNotRegressedOut/FDR0.5SNPProbes/eQTLs.txt-WSampleSize.txt.gz";
        String outfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/DataForRebuttal/2013-03-Trans40PCsCisFXNotRegressedOut/FDR0.5SNPProbes/eQTLs.txt-WSampleSize.txt.gzFilteredForMetaCisFXRemovedFDR0.05.txt";

        try {

            HashSet<Pair<String, String>> eqtls = new HashSet<Pair<String, String>>();
            TextFile tf = new TextFile(filterfor, TextFile.R);

            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[1];
                String probe = elems[4];

                eqtls.add(new Pair<String, String>(snp, probe));
                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();

            TextFile tfout = new TextFile(outfile, TextFile.W);
            TextFile tf2 = new TextFile(infile, TextFile.R);
            tfout.writeln(tf2.readLine()+"\tr2");
            elems = tf2.readLineElems(TextFile.tab);
            while (elems != null) {

                String snp = elems[1];
                String probe = elems[4];
                Pair<String, String> p = new Pair<String, String>(snp, probe);
                if(eqtls.contains(p)){
                    
                    int nrSamples = 0; 
                    
                    String[] sampleSizeStr = elems[13].split(";");
                    for(int i=0; i<sampleSizeStr.length; i++){
                        if(!sampleSizeStr[i].equals("-")){
                            int ct = Integer.parseInt(sampleSizeStr[i]);
                            nrSamples+=ct;
                        }
                    }
                    Double obsZ = Double.parseDouble(elems[10]);
                    double r = ZScores.zScoreToCorrelation(obsZ, nrSamples);
                    double r2 = r*r;
                    tfout.writeln(Strings.concat(elems, Strings.tab) + "\t"+r2);
                }
                

                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();
            tfout.close();


            eQTLFileProbeNameReplace.rewriteProbeIds(
                    outfile,
                    outfile + "ht12v3.txt",
                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",
                    "Probe",
                    "HumanHT-12_V3_0_R2_11283641_A.txt", true);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
