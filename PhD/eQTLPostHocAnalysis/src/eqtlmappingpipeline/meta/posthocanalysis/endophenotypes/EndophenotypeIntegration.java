/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.endophenotypes;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class EndophenotypeIntegration {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            EndophenotypeIntegration.run(
                    "/Volumes/iSnackHD/Dropbox/eQTLMeta/EndophenotypeCorrelations/2012-09-27-EGCUT-spearman-LOGp-TransFDR-FADS1-Exp40PCA-genVecNotRem_109012.txt",
                    "/Volumes/iSnackHD/Dropbox/eQTLMeta/EndophenotypeCorrelations/RS-RESULTS-endophenotypes-and-probes-correlations-P-VALUES-AND-RHO.txt",
                    "/Volumes/iSnackHD/Dropbox/eQTLMeta/EndophenotypeCorrelations/EndoPhenotypeTranslationTable2.txt");
        } catch (IOException e) {
        }
    }

    private static void run(String f1, String f2, String translation) throws IOException {

        // load all pvals in f1

        // inventorize the set of genes available.
        TextFile tf1 = new TextFile(f1, TextFile.R);
        HashMap<String, String> probeToGene = new HashMap<String, String>();
        String[] traits1 = tf1.readLineElems(TextFile.tab);
        HashMap<String, Integer> trait1Map = new HashMap<String, Integer>();
        for (int i = 0; i < traits1.length; i++) {
            trait1Map.put(traits1[i], i);
        }

        String[] elems1 = tf1.readLineElems(TextFile.tab);
        HashMap<String, Integer> geneMap1 = new HashMap<String, Integer>();
        ArrayList<String> genes1 = new ArrayList<String>();
        int geneCTR = 0;
        while (elems1 != null) {
            String gene = elems1[1];
            Integer geneId = geneMap1.get(gene);
            if (geneId == null) {
                geneId = geneCTR;
                probeToGene.put(gene, elems1[5]);
                geneMap1.put(gene, geneId);
                genes1.add(gene);
                geneCTR++;
            }
            elems1 = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();

        double[][] pvals1 = new double[geneMap1.size()][traits1.length];
        for (int i = 0; i < pvals1.length; i++) {
            for (int j = 0; j < pvals1[i].length; j++) {
                pvals1[i][j] = 2;
            }

        }

        tf1.open();

        // read in the pvalues
        elems1 = tf1.readLineElems(TextFile.tab); // header
        elems1 = tf1.readLineElems(TextFile.tab);
        while (elems1 != null) {
            String gene = elems1[1];
            Integer id = geneMap1.get(gene);
            for (int i = 7; i < elems1.length; i++) {

                Double d = 2d;
                try {
                    d = Double.parseDouble(elems1[i]);
                    // convert to pval using Power(10,-d)
                    d = Math.pow(10, -d);
                } catch (NumberFormatException e) {
//                    System.out.println(e.getMessage());
                }

                pvals1[id][i] = d;
            }

            elems1 = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();


        // load all pvals in f2

        TextFile tf2 = new TextFile(f2, TextFile.R);
        HashMap<String, Integer> geneMap2 = new HashMap<String, Integer>();
        String[] traits2 = tf2.readLineElems(TextFile.tab);
        HashMap<String, Integer> trait2Map = new HashMap<String, Integer>();
        for (int i = 0; i < traits2.length; i++) {
            trait2Map.put(traits2[i], i);
        }

        ArrayList<String> genes2 = new ArrayList<String>();
        String[] elems2 = tf2.readLineElems(TextFile.tab);
        geneCTR = 0;
        while (elems2 != null) {
            if (elems2.length > 0) {
                String gene = elems2[0];
                Integer geneId = geneMap2.get(gene);
                if (geneId == null) {
                    geneId = geneCTR;
                    geneMap2.put(gene, geneId);
                    genes2.add(gene);
                    geneCTR++;
//                System.out.println(geneCTR);
                }

            }
//            String ln = Strings.concat(elems2, Strings.tab);
//            System.out.println(ln);
            elems2 = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();

        double[][] pvals2 = new double[geneMap2.size()][traits2.length];
        for (int i = 0; i < pvals2.length; i++) {
            for (int j = 0; j < pvals2[i].length; j++) {
                pvals2[i][j] = 2;
            }

        }
        tf2.open();
        elems2 = tf2.readLineElems(TextFile.tab); // header
        elems2 = tf2.readLineElems(TextFile.tab);
        while (elems2 != null) {

            if (elems2.length > 0) {
                String gene = elems2[0];
                Integer id = geneMap2.get(gene);

                for (int i = 0; i < elems2.length; i++) {
                    Double d = 2d;
                    double prevd = pvals2[id][i];

                    try {
                        d = Double.parseDouble(elems2[i]);
                        // convert to pval using Power(10,-d)
                    } catch (NumberFormatException e) {
//                    System.out.println(e.getMessage());
                    }

                    if (prevd < 2 && d >= 2) {
                        // do not replace
                    } else if (prevd < 2 && d < 2 && prevd != d) {
//                        System.out.println("WARNING: same trait / genec combo, but different pvalue?: " + prevd + "\t" + d + "\t" + traits2[i] + "\t" + gene);
                    } else {
                        pvals2[id][i] = d;
                    }
                }
            }

            elems2 = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();

        // load endophenotype translation table
        int[] EGToRsFeno = new int[traits1.length];
        int[] rsToEGFeno = new int[traits2.length];

        // format is one (RS) --> many (EGCUT)
        TextFile endop = new TextFile(translation, TextFile.R);
        String[] endoelems = endop.readLineElems(TextFile.tab);
        while (endoelems != null) {
            if (endoelems.length > 1) {
                String egFeno = endoelems[0];
                String rsFeno = endoelems[1];

                Integer egFenoId = trait1Map.get(egFeno);
                Integer rsFenoId = trait2Map.get(rsFeno);

                if (rsFenoId != null && egFenoId != null) {
                    EGToRsFeno[egFenoId] = rsFenoId;
                    rsToEGFeno[rsFenoId] = egFenoId;
                }
            }
            endoelems = endop.readLineElems(TextFile.tab);
        }
        endop.close();


        // pick correlations where p < 0.05

        for (int gene1Id = 0; gene1Id < pvals1.length; gene1Id++) {
            for (int trait1 = 0; trait1 < pvals1[gene1Id].length; trait1++) {
                double egpval = pvals1[gene1Id][trait1];
                if (egpval < 0.05) {
                    // we have a hit. print fenotype and get rs accompanying data, when present
                    String trait = traits1[trait1];
                    String gene = genes1.get(gene1Id);
                    int rsFeno = EGToRsFeno[trait1];
                    Integer rsGene = geneMap2.get(gene);

                    String output = gene + "\t" + probeToGene.get(gene) + "\t" + trait + "\t" + egpval;

                    if (rsFeno > 0 && rsGene != null) {

                        double pvalRS = pvals2[rsGene][rsFeno];
                        if (pvalRS < 2) {
                            output += "\t" + traits2[rsFeno] + "\t" + pvalRS;
                        } else {
                            output += "\tPheno not tested in RS for gene\t-";
                        }
                    } else if (rsGene == null) {
                        output += "\tGene not tested in RS\t-";
                    } else {
                        output += "\tPheno not present in RS\t-";
                    }
                    System.out.println(output);
                }
            }
        }
        System.out.println("");
        System.out.println("---");
        System.out.println("");
        for (int geneId2 = 0; geneId2 < pvals2.length; geneId2++) {
            for (int trait2 = 3; trait2 < pvals2[geneId2].length; trait2++) {
                double rspval = pvals2[geneId2][trait2];
                if (rspval < 0.05) {
                    // see if there are other 
                    String gene = genes2.get(geneId2);
                    String trait = traits2[trait2];
                    String output = gene + "\t" + probeToGene.get(gene) + "\t" + trait + "\t" + rspval;


                    // get trait in EGCUT
                    int egtraitid = rsToEGFeno[trait2];
                    Integer egGene = geneMap1.get(gene);
                    if (egGene != null) {
                        if (pvals1[egGene][egtraitid] < 2) {
                            output += "\t" + traits1[egtraitid] + "\t" + pvals1[egGene][egtraitid];
                        }
                    } else {
                        output += "\tgene not tested in EGCUT\t-";
                    }
                    System.out.println(output);
                }
            }
        }


    }
}
