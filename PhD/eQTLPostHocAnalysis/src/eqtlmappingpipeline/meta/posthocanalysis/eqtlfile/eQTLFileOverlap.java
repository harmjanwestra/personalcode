/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harmjan
 */
public class eQTLFileOverlap {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            eQTLFileOverlap.determineOverlap(
                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PostQC/eQTLsFDR0.05.txt",
                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/trans/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt");

            eQTLFileOverlap.determineOverlapSNPs(
                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-GWAS-SNPs.txt", 
                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PostQC/eQTLsFDR0.5.txt");
            
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void determineOverlap(String file1, String file2) throws IOException {
        HashSet<String> snpsInFile1 = new HashSet<String>();
        HashSet<String> snpsInFile2 = new HashSet<String>();

        HashSet<String> genesInFile1 = new HashSet<String>();
        HashSet<String> genesInFile2 = new HashSet<String>();

        HashSet<String> probesInFile1 = new HashSet<String>();
        HashSet<String> probesInFile2 = new HashSet<String>();


        int nrUniqueSNPsFile1 = 0;
        int nrUniqueProbesFile1 = 0;
        int nrUniqueGenesFile1 = 0;

        int nrUniqueSNPsFile2 = 0;
        int nrUniqueProbesFile2 = 0;
        int nrUniqueGenesFile2 = 0;

        int nrSharedSNPs = 0;
        int nrSharedProbes = 0;
        int nrSharedGenes = 0;


        TextFile tf = new TextFile(file1, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        elems = tf.readLineElemsReturnObjects(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];
            String gene = elems[eQTLTextFile.HUGO];

            if (!gene.equals("-")) {
                genesInFile1.add(gene);
            }

            probesInFile1.add(probe);
            snpsInFile1.add(snp);
            elems = tf.readLineElemsReturnObjects(TextFile.tab);
        }


        tf.close();

        TextFile tf2 = new TextFile(file2, TextFile.R);
        elems = tf2.readLineElems(TextFile.tab);
        elems = tf2.readLineElemsReturnObjects(TextFile.tab);




        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];


            String gene = elems[eQTLTextFile.HUGO];

            if (!gene.equals("-")) {
                genesInFile2.add(gene);
            }
            probesInFile2.add(probe);
            snpsInFile2.add(snp);

            elems = tf2.readLineElemsReturnObjects(TextFile.tab);
        }
        tf2.close();

        String[] genesInFile1Arr = genesInFile1.toArray(new String[0]);
        String[] snpsInFile1Arr = snpsInFile1.toArray(new String[0]);
        String[] probesInFile1Arr = probesInFile1.toArray(new String[0]);


        String[] genesInFile2Arr = genesInFile2.toArray(new String[0]);
        String[] snpsInFile2Arr = snpsInFile2.toArray(new String[0]);
        String[] probesInFile2Arr = probesInFile2.toArray(new String[0]);


        for (String s : genesInFile1Arr) {
            if (!genesInFile2.contains(s)) {
                nrUniqueGenesFile1++;
            } else {
                nrSharedGenes++;
            }
        }

        for (String s : snpsInFile1Arr) {
            if (!snpsInFile2.contains(s)) {
                nrUniqueSNPsFile1++;
            } else {
                nrSharedSNPs++;
            }
        }

        for (String s : probesInFile1Arr) {
            if (!probesInFile2.contains(s)) {
                nrUniqueProbesFile1++;
            } else {
                nrSharedProbes++;
            }
        }

        /// again for file2 
        for (String s : genesInFile2Arr) {
            if (!genesInFile1.contains(s)) {
                nrUniqueGenesFile2++;
            }
        }

        for (String s : snpsInFile2Arr) {
            if (!snpsInFile1.contains(s)) {
                nrUniqueSNPsFile2++;
            }
        }

        for (String s : probesInFile2Arr) {
            if (!probesInFile1.contains(s)) {
                nrUniqueProbesFile2++;
            }

        }
        System.out.println("");
        System.out.println("Statistics:");
        System.out.println("-----------");
        System.out.println("Type\tFile1\tUniqueToFile1\tFile2\tUniqueToFile2\tShared");
        System.out.println("SNPs\t" + snpsInFile1.size() + "\t" + nrUniqueSNPsFile1 + "\t" + snpsInFile2.size() + "\t" + nrUniqueSNPsFile2 + "\t" + nrSharedSNPs);
        System.out.println("Probes\t" + probesInFile1.size() + "\t" + nrUniqueProbesFile1 + "\t" + probesInFile2.size() + "\t" + nrUniqueProbesFile2 + "\t" + nrSharedProbes);
        System.out.println("Genes\t" + genesInFile1.size() + "\t" + nrUniqueGenesFile1 + "\t" + genesInFile2.size() + "\t" + nrUniqueGenesFile2 + "\t" + nrSharedGenes);
        System.out.println("");
    }

    private static void determineOverlapSNPs(String snpfile, String eqtlfile) throws IOException {
        TextFile tf = new TextFile(snpfile, TextFile.R);
        HashSet<String> snps = new HashSet<String>();
        snps.addAll(tf.readAsArrayList());
        tf.close();

        TextFile tf2 = new TextFile(eqtlfile, TextFile.R);
        String[] elems = tf2.readLineElems(TextFile.tab);

        HashSet<String> snpsInFile = new HashSet<String>();
        HashSet<String> probesInFile = new HashSet<String>();
        HashSet<String> genesInFile = new HashSet<String>();


        elems = tf2.readLineElems(TextFile.tab);
        while (elems != null) {

            String snp = elems[1];
            if (snps.contains(snp)) {
                String gene = elems[eQTLTextFile.HUGO];
                snpsInFile.add(snp);
                if (!gene.equals("-")) {
                    genesInFile.add(gene);
                }
                String probe = elems[4];
                probesInFile.add(probe);
            }

            elems = tf2.readLineElems(TextFile.tab);
        }

        tf2.close();
        System.out.println("");
        System.out.println("Shared with SNP list:");
        System.out.println("---------------------");
        System.out.println("SNPs\t"+snpsInFile.size());
        System.out.println("Probes\t"+probesInFile.size());
        System.out.println("Genes\t"+genesInFile.size());
        System.out.println("");
    }
}
