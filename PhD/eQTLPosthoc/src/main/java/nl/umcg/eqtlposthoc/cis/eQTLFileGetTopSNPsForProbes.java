/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.cis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.regex.Pattern;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class eQTLFileGetTopSNPsForProbes {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {///Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/PreQC
//            String indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/";
//
//            String outdir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/SNPSelectionsForFuncAnalysis/";
//
//            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLProbesFDR0.05.txt", TextFile.R);
            String indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/";

            String outdir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC-SNPsForFuncAnnot/";

            TextFile tf = new TextFile("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLProbesFDR0.05-ProbeLevel.txt", TextFile.R);
                        
            HashSet<String> selectedProbes = new HashSet<String>();

            String[] elems = tf.readLineElems(TextFile.tab);// skip the header.....
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                selectedProbes.add(elems[eQTLTextFile.PROBE]);
                elems = tf.readLineElems(TextFile.tab);
            }
            System.out.println("Detected: " + selectedProbes.size() + " probes");
            tf.close();

            eQTLFileGetTopSNPsForProbes n = new eQTLFileGetTopSNPsForProbes();

            n.collect(tf.getFileName(),
                    selectedProbes, 
                    outdir+"RealData");
//            n.collect("/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PostQC/eQTLProbesFDR0.05.txt",
//                    selectedProbes, 
//                    outdir+"RealData");
            
            for (int perm = 1; perm < 11; perm++) {
                System.out.println("Running perm " + perm);
                String fileIn = indir + "PermutedEQTLsPermutationRound" + perm + ".txt.gz";

                n.collect(fileIn, selectedProbes, outdir + "PermutationRound" + perm);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void collect(String infile, HashSet<String> selectedProbes, String out) throws IOException {
        TextFile tf = new TextFile(infile, TextFile.R);


        String[] elems = tf.readLineElems(TextFile.tab);
        elems = tf.readLineElems(TextFile.tab);
        HashSet<String> probesVisited = new HashSet<String>();

        ArrayList<String> selectedSNPs = new ArrayList<String>();
        while (elems != null) {
            String probe = elems[eQTLTextFile.PROBE];
            if (selectedProbes.contains(probe) && !probesVisited.contains(probe)) {
                selectedSNPs.add(elems[1]);
                probesVisited.add(probe);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        System.out.println("Done reading data..");
        String[] snpArr = selectedSNPs.toArray(new String[0]);
        Pattern p = Pattern.compile("\n");
        Gpio.createDir(out);
        for (int i = 0; i < selectedSNPs.size(); i += 1000) {
            TextFile tfout = new TextFile(out + "/Bin" + i + "-" + (i + 1000) + ".txt", TextFile.W);
            int max = i + 1001;
            if (max > selectedSNPs.size()) {
                max = i + (selectedSNPs.size() - i);
            }
            System.out.println(i + "\t" + max);
            String output = Strings.concat(snpArr, p, i, max);
            tfout.write(output);
            tfout.close();
        }
        
        TextFile allSNPs = new TextFile(out+"/AllSNPs.txt", TextFile.W);
        for(String snp: snpArr){
            allSNPs.writeln(snp);
        }
        allSNPs.close();


        tf.close();
    }
}
