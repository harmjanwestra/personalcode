/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.smallscripts;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class GWASCatalogTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String loc = "/Users/harmjan/Downloads/tmp.txt";
            GWASCatalog c = new GWASCatalog(loc);

            GWASTrait t = c.getTraitToObj().get("Celiac disease");
            GWASSNP[] data = t.getSNPs(5E-8);
            System.out.println(data.length);
            
//            TextFile tf = new TextFile(loc, TextFile.R);
//            tf.readLine();
//            String[] elems = tf.readLineElems(TextFile.tab);
//            HashMap<String, HashSet<String>> gwasloci = new HashMap<String, HashSet<String>>();
//            HashSet<String> traits = new HashSet<String>();
//            while (elems != null) {
//                if (elems.length > 20) {
//                    String trait = elems[7];
//                    while(trait.endsWith(" ")){
//                        trait = trait.substring(0, trait.length()-1);
//                    }
//        
//                    HashSet<String> loci = gwasloci.get(trait);
//                    if (loci == null) {
//                        loci = new HashSet<String>();
//                    }
//
//                    String snps = elems[20];
//                    snps = snps.replaceAll(" ", "");
//                    String[] snpElems = snps.split("-");
//                    for (String s : snpElems) {
//                        if (s.startsWith("rs")) {
//                            loci.add(s);
//                        }
//                    }
//
//                    snps = elems[21];
//                    snps = snps.replaceAll(" ", "");
//                    snpElems = snps.split(",");
//
//                    for (String s : snpElems) {
//                        if (s.startsWith("rs")) {
//                            loci.add(s);
//                        }
//                    }
//                    traits.add(trait);
//                    gwasloci.put(trait, loci);
//
//                }
//                elems = tf.readLineElems(TextFile.tab);
//            }
//            tf.close();
//
//            for (String trait : traits) {
//                boolean outputtrait = false;
//                GWASTrait catalogTrait = c.getTraitToObj().get(trait);
//                HashSet<String> loci = gwasloci.get(trait);
//                if (catalogTrait == null) {
//                    System.out.println(trait + "\t" + 0 + "\t" + loci.size() + "\tTrait not in catalog object..");
//                } else {
//                    GWASSNP[] catalogSNPs = catalogTrait.getSNPs();
//
//                    HashSet<String> snpsInCatalog = new HashSet<String>();
//                    for (GWASSNP snp : catalogSNPs) {
//                        snpsInCatalog.add(snp.getName());
//                    }
//
//                    if (snpsInCatalog.size() != loci.size()) {
//                        outputtrait = true;
//                        for (String s : loci) {
//                            if (!snpsInCatalog.contains(s)) {
//                                System.out.println(trait + "\t" + s + "\tnot found in catalog object");
//                            }
//                        }
//
//                        for (String s : snpsInCatalog) {
//                            if (!loci.contains(s)) {
//                                System.out.println(trait + "\t" + s + "\tnot found in test");
//                            }
//                        }
//                    }
//                    if (outputtrait) {
//                        System.out.println("");
//                    }
////                    System.out.println(trait + "\t" + snpsInCatalog.size() + "\t" + loci.size());
//                }
//
//
//
//
//            }


        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
