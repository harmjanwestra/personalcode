/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.smallscripts;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class CrohnsDiseaseTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String gwasCatalog = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-09-23-gwascatalog.txt";
        String eqtlfile = "/Volumes/iSnackHD/AeroFS/table.txt";
        try {
            TextFile tf = new TextFile(eqtlfile, TextFile.R);
            HashSet<String> eqtlSNPs = new HashSet<String>();
            HashSet<String> neutroeqtlSNPs = new HashSet<String>();
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[1];
                String neutro = elems[elems.length - 1];
                if (neutro.equals("Neutrophils")) {
                    neutroeqtlSNPs.add(snp);
                }
                eqtlSNPs.add(snp);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            GWASCatalog catalog = new GWASCatalog(gwasCatalog);
            GWASTrait trait = catalog.getTraitToObj().get("Crohn's disease");
            GWASSNP[] snps = trait.getSNPs(5E-8);
            HashSet<String> uniqueSNPs = new HashSet<String>();
            for (GWASSNP snp : snps) {
                uniqueSNPs.add(snp.getName());
            }

            int neutro = 0;
            int eqtl = 0;
            for (String snp : uniqueSNPs) {
                if (neutroeqtlSNPs.contains(snp)) {
                    neutro++;
                }
                if (eqtlSNPs.contains(snp)) {
                    eqtl++;
                }
            }

            System.out.println(uniqueSNPs.size() + "\t" + eqtl + "\t" + neutro);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
