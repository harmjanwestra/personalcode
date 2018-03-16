/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harmjan
 */
public class eQTLFileSummary {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // full cis
        String cisDir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/";
        // cis on trait assoc snps only.
//        cisDir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-08-22-CisOnTraitAssocSNPs/Sorted/";
        // /Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz-UpdatedEnsemblAnnotation.txt.gz"
        String transDir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/";
        try {
            // allcis
//            eQTLFileSummary.countSNPsProbesAndGenes(cisDir + "PreQC/eQTLs.txt.gz-UpdatedEnsemblAnnotation.txt.gz");
//            eQTLFileSummary.countSNPsProbesAndGenes(cisDir + "PreQC/eQTLs.txt-UpdatedEnsemblAnnotation.txt");

            // cisFDR
            eQTLFileSummary.countSNPsProbesAndGenes(cisDir + "PostQC/eQTLsFDR0.05.txt-UpdatedEnsemblAnnotation.txt");

            // trans FDR
            eQTLFileSummary.countSNPsProbesAndGenes(transDir + "PostQC/eQTLsFDR0.05.txt-UpdatedEnsemblAnnotation.txt");


            // intersect
            eQTLFileSummary.intersect(
                    cisDir + "PostQC/eQTLsFDR0.05.txt-UpdatedEnsemblAnnotation.txt",
//                    cisDir + "PostQC/eQTLsFDR0.5.txt-UpdatedEnsemblAnnotation.txt",
                    transDir + "PostQC/eQTLsFDR0.05.txt-UpdatedEnsemblAnnotation.txt",
                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-GWAS-SNPs.txt");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void countSNPsProbesAndGenes(String file) throws IOException {

        System.out.println("Summary of file: " + file);
        HashSet<String> genes = new HashSet<String>();
        HashSet<String> probes = new HashSet<String>();
        HashSet<String> snps = new HashSet<String>();

        TextFile t = new TextFile(file, TextFile.R);
        String[] elems = t.readLineElems(TextFile.tab);
        elems = t.readLineElems(TextFile.tab);
        int nrProbesNotMappingToGenes = 0;
        HashSet<String> probesNotMappingToGenes = new HashSet<String>();
        while (elems != null) {
            String geneStr = elems[eQTLTextFile.HUGO];

            probes.add(elems[4]);
            snps.add(elems[1]);
            if (geneStr.equals("-")) {
                probesNotMappingToGenes.add(elems[4]);
            } else {
                String[] geneStrArr = geneStr.split(",");
                genes.addAll(Arrays.asList(geneStrArr));
            }
            elems = t.readLineElems(TextFile.tab);
        }
        t.close();

        System.out.println("Genes: " + genes.size());
        System.out.println("Probes: " + probes.size());
        System.out.println("SNPs: " + snps.size());
        System.out.println("Probes not mapping to genes: " + probesNotMappingToGenes.size());

    }

    public static void intersect(String file1, String file2, String snpsToConfineTo) throws IOException {
        HashSet<String> snpsHashToConfineTo = null;
        if (snpsToConfineTo != null) {
            snpsHashToConfineTo = new HashSet<String>();
            TextFile snpfile = new TextFile(snpsToConfineTo, TextFile.R);
            snpsHashToConfineTo.addAll(snpfile.readAsArrayList());
            snpfile.close();
        }

        System.out.println("file1: " + file1);
        System.out.println("file2: " + file2);


        HashSet<String> probesInFile1 = new HashSet<String>();
        HashSet<String> snpsInFile1 = new HashSet<String>();
        HashSet<String> genesInFile1 = new HashSet<String>();
        HashSet<String> probesWithoutGenesInFile1 = new HashSet<String>();

        TextFile t = new TextFile(file1, TextFile.R);
        String[] elems = t.readLineElems(TextFile.tab);//header
        elems = t.readLineElems(TextFile.tab);
        int nrProbesNotMappingToGenes = 0;
        HashSet<String> probesNotMappingToGenes = new HashSet<String>();
        while (elems != null) {

            String snp = elems[1];
            if (snpsHashToConfineTo == null || snpsHashToConfineTo.contains(snp)) {
                String geneStr = elems[eQTLTextFile.HUGO];
                probesInFile1.add(elems[4]);
                snpsInFile1.add(elems[1]);
                if (geneStr.equals("-")) {
                    probesWithoutGenesInFile1.add(elems[4]);
                } else {
                    String[] geneStrArr = geneStr.split(",");
                    genesInFile1.addAll(Arrays.asList(geneStrArr));
                }
            }
            elems = t.readLineElems(TextFile.tab);
        }
        t.close();

        HashSet<String> probesInFile2 = new HashSet<String>();
        HashSet<String> snpsInFile2 = new HashSet<String>();
        HashSet<String> genesInFile2 = new HashSet<String>();
        HashSet<String> probesWithoutGenesInFile2 = new HashSet<String>();

        TextFile t2 = new TextFile(file2, TextFile.R);
        String[] elems2 = t2.readLineElems(TextFile.tab); //header
        elems2 = t2.readLineElems(TextFile.tab);

        while (elems2 != null) {
            String snp = elems2[1];
            if (snpsHashToConfineTo == null || snpsHashToConfineTo.contains(snp)) {
                String geneStr = elems2[eQTLTextFile.HUGO];

                probesInFile2.add(elems2[4]);
                snpsInFile2.add(elems2[1]);
                if (geneStr.equals("-")) {
                    probesWithoutGenesInFile2.add(elems2[4]);
                } else {
                    String[] geneStrArr = geneStr.split(",");
                    genesInFile2.addAll(Arrays.asList(geneStrArr));

                }
            }
            elems2 = t2.readLineElems(TextFile.tab);

        }
        t2.close();

        int probesShared = 0;
        int snpssShared = 0;
        int genesShared = 0;
        int probesWithoutGenesShared = 0;

        int probesUniqueToFile1 = 0;
        int snpsUniqueToFile1 = 0;
        int genesUniqueToFile1 = 0;
        int probesWithoutGenesUniqueToFile1 = 0;

        int probesUniqueToFile2 = 0;
        int snpsUniqueToFile2 = 0;
        int genesUniqueToFile2 = 0;
        int probesWithoutGenesUniqueToFile2 = 0;


        for (String s : probesInFile1) {
            if (probesInFile2.contains(s)) {
                probesShared++;
            } else {
                probesUniqueToFile1++;
            }
        }
        for (String s : snpsInFile1) {
            if (snpsInFile2.contains(s)) {
                snpssShared++;
            } else {
                snpsUniqueToFile1++;
            }
        }
        for (String s : genesInFile1) {
            if (genesInFile2.contains(s)) {
//                System.out.println(s+"\tshared");
                genesShared++;
            } else {
                genesUniqueToFile1++;
            }
        }
        for (String s : probesWithoutGenesInFile1) {
            
            if (probesWithoutGenesInFile2.contains(s)) {
                probesWithoutGenesShared++;
            } else {
                probesWithoutGenesUniqueToFile1++;
            }
        }

// difference file 2
        for (String s : probesInFile2) {
            if (!probesInFile1.contains(s)) {
                
                probesUniqueToFile2++;
            }
        }
        for (String s : snpsInFile2) {
            if (!snpsInFile1.contains(s)) {
                
                snpsUniqueToFile2++;
            }
        }
        for (String s : genesInFile2) {
            if (!genesInFile1.contains(s)) {
                
                genesUniqueToFile2++;
            }
        }
        for (String s : probesWithoutGenesInFile2) {
            if (!probesWithoutGenesInFile1.contains(s)) {
                
                probesWithoutGenesUniqueToFile2++;
            }
        }

        HashSet<String> snpsInBothFiles = new HashSet<String>();
        snpsInBothFiles.addAll(snpsInFile1);
        snpsInBothFiles.addAll(snpsInFile2);
        
        HashSet<String> genesInBothFiles = new HashSet<String>();
        genesInBothFiles.addAll(genesInFile1);
        genesInBothFiles.addAll(genesInFile2);
        
        HashSet<String> probesInBothFiles = new HashSet<String>();
        probesInBothFiles.addAll(probesInFile1);
        probesInBothFiles.addAll(probesInFile2);
        
        HashSet<String> probesWithoutGeneNameInBothFiles = new HashSet<String>();
        probesWithoutGeneNameInBothFiles.addAll(probesWithoutGenesInFile1);
        probesWithoutGeneNameInBothFiles.addAll(probesWithoutGenesInFile2);
        
        
        System.out.println("Summary of intersect: ");
        System.out.println("Type\t\t\tNrInBothFiles\tNrInFile1\tNrInFile2\tShared\tUniqueToFile1\tUniqueToFile2");

        System.out.println("SNPs\t\t\t"+snpsInBothFiles.size()+"\t\t"+snpsInFile1.size()+"\t\t"+snpsInFile2.size()+"\t\t"+snpssShared+"\t"+snpsUniqueToFile1+"\t\t"+snpsUniqueToFile2);
        System.out.println("Probes\t\t\t"+probesInBothFiles.size()+"\t\t"+probesInFile1.size()+"\t\t"+probesInFile2.size()+"\t\t"+probesShared+"\t"+probesUniqueToFile1+"\t\t"+probesUniqueToFile2);
        System.out.println("Genes\t\t\t"+genesInBothFiles.size()+"\t\t"+genesInFile1.size()+"\t\t"+genesInFile2.size()+"\t\t"+genesShared+"\t"+genesUniqueToFile1+"\t\t"+genesUniqueToFile2);
        System.out.println("Probes w/o genes\t"+probesWithoutGeneNameInBothFiles.size()+"\t\t"+probesWithoutGenesInFile1.size()+"\t\t"+probesWithoutGenesInFile2.size()+"\t\t"+probesWithoutGenesShared+"\t"+probesWithoutGenesUniqueToFile1+"\t\t"+probesWithoutGenesUniqueToFile2);
        
    }
}
