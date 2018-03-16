/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.cis;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class GetAllProxiesForSNPsFromOlderProxyFile {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        // load in proxies from older file

        // load snps in new file
        String pf1F = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/SNPFunctionalAnnotation/Cis-SNPs/AllSNPsPerm.txt-WithProxies.txt";
        String pf2F = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/SNPFunctionalAnnotation/Cis-SNPs/RealData/AllSNPs.txt-WithProxies.txt";
        String indir = "";

        
        try {

            for (int perm = 0; perm < 11; perm++) {
                System.out.println("Perm: "+perm);
                if (perm == 0) {
                    indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC-SNPsForFuncAnnot/RealData/";
                } else {
                    indir = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC-SNPsForFuncAnnot/PermutationRound" + perm + "/";
                }

                for (int i = 0; i < 9000; i += 1000) {
                    HashSet<String> query = new HashSet<String>();
                    String queryfile = indir + "Bin" + i + "-" + (i + 1000) + ".txt";
                    System.out.println(queryfile);

                    TextFile in = new TextFile(queryfile, TextFile.R);
                    query.addAll(in.readAsArrayList());
                    in.close();

                    HashSet<String> proxiesToOutput = new HashSet<String>();
                    TextFile pf1 = new TextFile(pf1F, TextFile.R);
                    String[] elems = pf1.readLineElems(TextFile.tab);
                    while (elems != null) {
                        if (query.contains(elems[0])) {
                            proxiesToOutput.add(elems[0]);
                            proxiesToOutput.add(elems[1]);
                        }
                        elems = pf1.readLineElems(TextFile.tab);
                    }
                    pf1.close();

                    pf1 = new TextFile(pf2F, TextFile.R);
                    elems = pf1.readLineElems(TextFile.tab);
                    while (elems != null) {
                        if (query.contains(elems[0])) {
                            proxiesToOutput.add(elems[0]);
                            proxiesToOutput.add(elems[1]);
                        }
                        elems = pf1.readLineElems(TextFile.tab);
                    }
                    pf1.close();
                    
                    TextFile out = new TextFile(indir + "Bin" + i + "-" + (i + 1000) + "-WithProxies.txt", TextFile.W);
                    out.writeList(Arrays.asList(proxiesToOutput.toArray(new String[0])));
                    out.close();
                }



            }

        } catch (IOException e) {
            e.printStackTrace();
        }

//        HashSet<Pair<String, String>> proxy = new HashMap<String, HashSet<String>>();

    }
}
