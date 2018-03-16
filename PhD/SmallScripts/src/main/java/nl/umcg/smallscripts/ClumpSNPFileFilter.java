/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ClumpSNPFileFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String queryFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/TestedSNPs.txt";
            HashSet<String> querySNPs = new HashSet<String>();
            TextFile tf = new TextFile(queryFile, TextFile.R);
            querySNPs.addAll(tf.readAsArrayList());
            tf.close();

            String gwasSNPFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt";
            TextFile tf2 = new TextFile(gwasSNPFile, TextFile.R);
            HashMap<String, Double> lowestPForGWASSNP = new HashMap<String, Double>();

            String[] elems = tf2.readLineElems(TextFile.tab);
            while (elems != null) {
                if (elems.length > 27) {
                    String snpStr = elems[21];
                    snpStr = snpStr.replaceAll(" ", "");
                    String[] multisnp = snpStr.split(",");
                    for (String snp : multisnp) {
                        String pStr = elems[27];
                        System.out.println(snp + "\t" + pStr);
                        try {
                            Double p = Double.parseDouble(pStr);

                            Double d = lowestPForGWASSNP.get(snp);
                            if (d == null) {
                                lowestPForGWASSNP.put(snp, p);
                            } else {
                                if (d > p) {
                                    lowestPForGWASSNP.put(snp, p);
                                }
                            }
                        } catch (NumberFormatException e) {
                        }
                    }

                }
                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            System.out.println(lowestPForGWASSNP.size());

            System.out.println("Now printing to: " + gwasSNPFile + "FilteredForTestedSNPs.txt");
            TextFile tfout = new TextFile(gwasSNPFile + "FilteredForTestedSNPs.txt", TextFile.W);
            tfout.writeln("P\tSNP");
            int ctr = 0;
            for (String snp : querySNPs) {
                Double p = lowestPForGWASSNP.get(snp);
                if(p == null){
                    p = 0.9;
                }
                System.out.println(p + "\t" + snp);
                tfout.writeln(p + "\t" + snp);
                ctr++;
            }
            System.out.println(ctr);
            tfout.close();


//            tfout.writeln("P\tSNP");
//            String[] elems = tf2.readLineElems(TextFile.tab);
//            elems = tf2.readLineElems(TextFile.tab);
//            while(elems!=null){
//                String p = elems[0];
//                String snps = elems[1];
//                String[] snpelems = snps.split(",");
//                for(String s: snpelems){
//                    if(querySNPs.contains(s)){
//                        tfout.writeln(p+"\t"+s);
//                    }
//                }
//                        
//                elems = tf2.readLineElems(TextFile.tab);
//            }
//            tfout.close();
//            tf2.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
