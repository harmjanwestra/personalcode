/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.trans;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harmjan
 */
public class GetAssociatedSNPs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        GetAssociatedSNPs g = new GetAssociatedSNPs();
        int nrSelected = 0;
        try {
            String indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/";
            String outdir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/SNPSelectionsForFuncAnalysis/";
            Gpio.createDir(outdir);
            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLSNPsFDR0.05.txt", TextFile.R);

            HashSet<String> selectedProbes = new HashSet<String>();
            String[] elems = tf.readLineElems(TextFile.tab);// skip the header.....
            elems = tf.readLineElems(TextFile.tab);
            HashSet<String> snps = new HashSet<String>();
            while (elems != null) {
                snps.add(elems[1]);
                selectedProbes.add(elems[eQTLTextFile.PROBE]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile tfout = new TextFile(outdir + "RealData.txt", TextFile.W);
            tfout.writeList(Arrays.asList(snps.toArray(new String[0])));
            tfout.close();

            for (int perm = 1; perm < 11; perm++) {
                String fileIn = indir + "PermutedEQTLsPermutationRound" + perm + ".txt.gz";
                String out = outdir + "PermutationRound" + perm + ".txt";
                int nr = g.collect(fileIn, out, 346);
                System.out.println(perm + "\t" + nr);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public int collect(String fileIn, String fileout, int nrToSelect) throws IOException {
        TextFile in = new TextFile(fileIn, TextFile.R);
        HashSet<String> snps = new HashSet<String>();
        String[] elems = in.readLineElems(TextFile.tab);
        elems = in.readLineElems(TextFile.tab);
        if (nrToSelect == -1) {
            nrToSelect = Integer.MAX_VALUE;
        }
        while (elems != null) {

            if (snps.size() >= nrToSelect) {
                break;
            } else {
                snps.add(elems[1]);
            }
            elems = in.readLineElems(TextFile.tab);
        }

        in.close();

        TextFile out = new TextFile(fileout, TextFile.W);
        out.writeList(Arrays.asList(snps.toArray(new String[0])));
        out.close();
        in.close();
        return snps.size();
    }
}
