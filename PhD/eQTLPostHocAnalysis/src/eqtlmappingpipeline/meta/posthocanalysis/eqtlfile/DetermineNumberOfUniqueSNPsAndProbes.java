/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class DetermineNumberOfUniqueSNPsAndProbes {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            TextFile in = new TextFile("/Volumes/iSnackHD/Data/Projects/Vinod/eQTLResults/2012-10-16-LincRNARerun/eQTLsFDR0.05.txt", TextFile.R);
            String[] elems = in.readLineElems(TextFile.tab);
            HashSet<String> snps = new HashSet<String>();
            HashSet<String> probes = new HashSet<String>();
            HashSet<String> genes = new HashSet<String>();
            while (elems != null) {

                snps.add(elems[1]);
                probes.add(elems[4]);
                
                
                elems = in.readLineElems(TextFile.tab);
            }

            in.close();
            
            System.out.println("snps: "+snps.size());
            System.out.println("prob: "+probes.size());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
