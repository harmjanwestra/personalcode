/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harmjan
 */
public class Bonferronifier {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        double nrTests = 11172453;
        String filename = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt-UpdatedEnsemblAnnotation.txt.gz";
        HashSet<String> allowedSNPs = new HashSet<String>();
        String snpfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-TestedGWASSNPs.txt";
//        double nrTests = 153134630;
//        String filename = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt.gz";

        //"/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt.gz"
        try {
            
            TextFile tfsnps = new TextFile(snpfile, TextFile.R);
            allowedSNPs.addAll(tfsnps.readAsArrayList());
            tfsnps.close();
            
            double threshold = 0.05 / nrTests;
            TextFile tf = new TextFile(filename, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            HashSet<String> snps = new HashSet<String>();
            HashSet<String> probe = new HashSet<String>();
            HashSet<String> gene = new HashSet<String>();
            HashSet<String> probewogene = new HashSet<String>();
            int ctr = 0;
            double minZ = Double.MAX_VALUE;
            while (elems != null) {
                double p = Double.parseDouble(elems[0]);
                if (p < threshold && allowedSNPs.contains(elems[1])) {
                    ctr++;
                    double z = Math.abs(Double.parseDouble(elems[eQTLTextFile.METAZ]));
                    if(z < minZ){
                        minZ = z;
                    }
                    snps.add(elems[1]);
                    probe.add(elems[4]);
                    String geneStr = elems[eQTLTextFile.HUGO];
                    if (geneStr.equals("-")) {
                        probewogene.add(elems[4]);
                    } else {
                        String[] genesArr = geneStr.split(",");
                        for(String g:genesArr){
                            gene.add(g);
                        }
                    }

                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            
            System.out.println("Threshold\t"+threshold);
            System.out.println("Z\t"+minZ);
            System.out.println("Pairs\t"+ctr);
            System.out.println("snps\t"+snps.size());
            System.out.println("probes\t"+probe.size());
            System.out.println("genes\t"+gene.size());
            System.out.println("probeswogenes\t"+probewogene.size());
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
