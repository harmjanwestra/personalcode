/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class GetEQTLSNPs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {





            if (1 == 1) {
                HashSet<String> probes = new HashSet<String>();


                TextFile in = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/2012-09-24-GroningenProbesFDR0.05vsAllButGroningenFDR0.05.txt-OppositeEQTLs.txt", TextFile.R);
                String[] elems = in.readLineElems(TextFile.tab);
                while (elems != null) {
                    probes.add(elems[1]);
                    elems = in.readLineElems(TextFile.tab);
                }

                in.close();

                TextFile tf2 = new TextFile("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", TextFile.R);
                String header = tf2.readLine();
                String[] elems2 = tf2.readLineElems(TextFile.tab);
                while (elems2 != null) {
                    String probe = elems2[0];
                    String hugo = elems2[4];
                    if (probes.contains(probe)) {
                        System.out.println(probe + "\t" + hugo);
                    }
                    elems2 = tf2.readLineElems(TextFile.tab);
                }
                tf2.close();
                System.exit(0);
            }







            // /Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/2012-08-08-AllButGroningen/Sorted/eQTLsFDR0.05.txt
//            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/2012-08-08-GroningenOnly/Sorted/eQTLs.txt.gz", TextFile.R);
            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/2012-08-08-AllButGroningen/Sorted/eQTLs.txt.gz", TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab); // header
            elems = tf.readLineElems(TextFile.tab);
            HashSet<String> snpsInFile = new HashSet<String>();
            HashSet<String> probesInFile = new HashSet<String>();
            HashSet<String> snpProbePairsInFile = new HashSet<String>();
            while (elems != null) {
                String snp = elems[1];
                String probe = elems[4];
                probesInFile.add(probe);
                snpsInFile.add(snp);
                snpProbePairsInFile.add(snp + "\t" + probe);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            System.out.println("Nr unique SNPs tested in datasets excluding Groningen:\t" + snpsInFile.size());
            System.out.println("Nr unique probes tested in datasets excluding Groningen:\t" + probesInFile.size());

//            String[] list = snpsInFile.toArray(new String[0]);

//            TextFile out = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/cisSNPsInGroningenData.txt", TextFile.W);
//            TextFile out = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/cisSNPsInGroningenData.txt", TextFile.W);
//            TextFile out = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/cisSNPsInGroningenData.txt", TextFile.W);
//            out.writeList(Arrays.asList(list));
//            out.close();

            HashSet<String> snpProbePairsInGroningenFile = new HashSet<String>();
            if (1 == 1) {
                TextFile tf2 = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/2012-08-08-GroningenOnly/Sorted/eQTLProbesFDR0.05.txt", TextFile.R);
                TextFile tfout = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/2012-08-08-GroningenOnly/Sorted/eQTLProbesFDR0.05-SNPsAndProbePairsTestedInNonGroningenDatasets.txt", TextFile.W);
                String header = tf2.readLine();
                tfout.writeln(header);
                String[] tf2elems = tf2.readLineElems(TextFile.tab);
                while (tf2elems != null) {
                    String snp = tf2elems[1];
                    String probe = tf2elems[4];
                    //if(snpsInFile.contains(snp) && probesInFile.contains(probe)){
                    if (snpProbePairsInFile.contains(snp + "\t" + probe)) {
                        snpProbePairsInGroningenFile.add(snp + "\t" + probe);
                        tfout.writeln(Strings.concat(tf2elems, Strings.tab));
                    }

                    tf2elems = tf2.readLineElems(TextFile.tab);
                }
                tf2.close();
                tfout.close();
            }

            if (1 == 1) {
                TextFile tf2 = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/2012-08-08-AllButGroningen/Sorted/eQTLs.txt.gz", TextFile.R);
                TextFile tfout = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisOnSubsets/cis/2012-08-08-AllButGroningen/Sorted/eQTLs.txt.gz-OnlySNPProbesPairsSignificantEQTLProbesFDR0.05Groningen.txt", TextFile.W);
                String header = tf2.readLine();
                tfout.writeln(header);
                String[] tf2elems = tf2.readLineElems(TextFile.tab);
                while (tf2elems != null) {
                    String snp = tf2elems[1];
                    String probe = tf2elems[4];
                    if (snpProbePairsInGroningenFile.contains(snp + "\t" + probe)) {
                        tfout.writeln(Strings.concat(tf2elems, Strings.tab));
                    }

                    tf2elems = tf2.readLineElems(TextFile.tab);
                }
                tf2.close();
                tfout.close();
            }


        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
