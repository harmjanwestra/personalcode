/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class MergeProxySNPLists {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {// 
            String list1 = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/IBDSNPsWithHapMapProxies-WithPositions.txt";
            String list2 = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/IBDSNPsWith1KgProxies-WithPositions.txt";
            String output = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-09-23-gwascatalog-WithAllIBDLoci-Proxies.txt";

            String gwascatalog = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-09-23-gwascatalog-WithAllIBDLoci.txt";
            HashMap<String, HashSet<String>> proxiesPerSNP = new HashMap<String, HashSet<String>>();

            TextFile tf = new TextFile(list1, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[0];
                String snp2 = elems[1];
                if (snp2.startsWith("rs")) {
                    HashSet<String> q = proxiesPerSNP.get(snp);
                    if (q == null) {
                        q = new HashSet<String>();
                    }
                    q.add(snp2);
                    proxiesPerSNP.put(snp, q);
                }

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            tf = new TextFile(list2, TextFile.R);
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[0];
                String snp2 = elems[1];
                if (snp2.startsWith("rs")) {
                    HashSet<String> q = proxiesPerSNP.get(snp);
                    if (q == null) {
                        q = new HashSet<String>();
                    }
                    q.add(snp2);
                    proxiesPerSNP.put(snp, q);
                }

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile tf2 = new TextFile(gwascatalog, TextFile.R);
            TextFile out = new TextFile(output, TextFile.W);
            String[] elems2 = tf2.readLineElems(TextFile.tab);
            int traitcol = 7;
            int snpcol = 20;
            while (elems2 != null) {

                String trait = elems2[traitcol];
                if (trait.equals("Nature2012-IBD+Crohn-NoUCSpecific")) {
                    String snp = elems2[snpcol];
                    elems2[traitcol] = trait+"WithPerfectProxies";
                    HashSet<String> proxies = proxiesPerSNP.get(snp);
                    if (proxies == null) {
                        out.writeln(Strings.concat(elems2, Strings.tab));
                    } else {
                        for(String proxy: proxies){
                            elems2[snpcol] = proxy;
                            out.writeln(Strings.concat(elems2, Strings.tab));
                        }
                    }
                }

                elems2 = tf2.readLineElems(TextFile.tab);
            }

            out.close();
            tf2.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
