/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class VinodTopSNPFinder {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            HashSet<String> probesToSelect = new HashSet<String>();
            TextFile tf = new TextFile("/Volumes/iSnackHD/Data/Projects/Vinod/2012-10-03-LincRNA_ProbesFDR0.05.txt", TextFile.R);
            probesToSelect.addAll(tf.readAsArrayList());
            tf.close();
            String outdir = "/Volumes/iSnackHD/Data/Projects/Vinod/lincRNAPermutedTopSNPs/";
            Gpio.createDir(outdir);
            String indir = "/Volumes/iSnackHD/Users/vinod/LincRNA_probes_PermutedFiles/";
            for (int perm = 1; perm < 101; perm++) {
                System.out.println("Processing: "+perm);
                String fin = indir + "PermutedEQTLsPermutationRound"+perm+".txt";

                TextFile tfin = new TextFile(fin, TextFile.R);
                TextFile tfout = new TextFile(outdir+"PermutationRound"+perm+"-TopSNPs.txt", TextFile.W);
                HashSet<String> visitedProbes = new HashSet<String>();

                tfin.readLine();
                String[] elems = tfin.readLineElems(TextFile.tab);
                while (elems != null) {
                    String probe = elems[4];
                    if(!visitedProbes.contains(probe) && probesToSelect.contains(probe)){
                        tfout.writeln(elems[1]+"\t"+probe);
                        visitedProbes.add(probe);
                    }
                    elems = tfin.readLineElems(TextFile.tab);
                }
                tfin.close();
                tfout.close();

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
