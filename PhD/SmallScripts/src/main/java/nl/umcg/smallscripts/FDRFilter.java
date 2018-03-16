/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class FDRFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
//            for (int perm = 0; perm < 11; perm++) {
//                String indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/";
//                String outdir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/";
//                String infile = "eQTLs.txt.gz";
//                if (perm > 0) {
//                    infile = "PermutedEQTLsPermutationRound" + perm + ".txt.gz";
//                }
//
//                TextFile tfout = new TextFile(outdir + infile, TextFile.W);
//                TextFile tf = new TextFile(indir + infile, TextFile.R);
//                HashSet<String> probesWritten = new HashSet<String>();
//                String header = tf.readLine();
//                tfout.writeln(header);
//                String[] elems = tf.readLineElems(TextFile.tab);
//                while (elems != null) {
//
//                    String probe = elems[4];
//                    if (!probesWritten.contains(probe)) {
//                        tfout.writeln(Strings.concat(elems, Strings.tab));
//                        probesWritten.add(probe);
//                    }
//
//                    elems = tf.readLineElems(TextFile.tab);
//                }
//                tf.close();
//                tfout.close();
//            }

            for (int perm = 0; perm < 101; perm++) {
                String indir = "/Volumes/iSnackHD/Data/Projects/Vinod/eQTLResults/2012-10-16-LincRNARerun-GRNG+EGCUT-AllSNPProbePairs/";
                String outdir = "/Volumes/iSnackHD/Data/Projects/Vinod/eQTLResults/2012-10-16-LincRNARerun-GRNG+EGCUT-AllSNPProbePairs/FDRForProbesOnly/";
                Gpio.createDir(outdir);
                String infile = "eQTLs.txt";
                if (perm > 0) {
                    infile = "PermutedEQTLsPermutationRound" + perm + ".txt.gz";
                }

                TextFile tfout = new TextFile(outdir + infile, TextFile.W);
                TextFile tf = new TextFile(indir + infile, TextFile.R);
                HashSet<String> probesWritten = new HashSet<String>();
                String header = tf.readLine();
                tfout.writeln(header);
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {

                    String probe = elems[4];
                    if (!probesWritten.contains(probe)) {
                        tfout.writeln(Strings.concat(elems, Strings.tab));
                        probesWritten.add(probe);
                    }

                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
                tfout.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
