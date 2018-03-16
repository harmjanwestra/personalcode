/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.iv;

import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harm-jan
 */
public class DetermineNumberOFIVWithMultipleCis {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            TextFile tf = new TextFile("D:\\SkyDrive\\latesteQTLs\\trans_ma_results_31072012.txt-Filtered2.txt", TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            elems = tf.readLineElems(TextFile.tab);

            HashMap<String, HashSet<String>> snpprobes = new HashMap<String, HashSet<String>>();

            HashSet<Triple<String, String, String>> triples = new HashSet<Triple<String, String, String>>();
            HashSet<Triple<String, String, String>> signtriples = new HashSet<Triple<String, String, String>>();
            HashSet<String> snps = new HashSet<String>();
            while (elems != null) {

                String snp = elems[0];
                String probe = elems[1];
                String transprobe = elems[2];
                Double ivP = Double.parseDouble(elems[7]);

                HashSet<String> probes = snpprobes.get(snp);
                if (probes == null) {
                    probes = new HashSet<String>();
                }

                snps.add(snp);
                probes.add(probe);
                snpprobes.put(snp, probes);
                triples.add(new Triple<String, String, String>(snp, probe, transprobe));
                if (ivP < 0.05) {
                    signtriples.add(new Triple<String, String, String>(snp, probe, transprobe));
                }
                elems = tf.readLineElems(TextFile.tab);
            }

            int nrSig = 0;
            int nrMultiCis = 0;
            for(String snp: snps){
                if(snpprobes.get(snp).size() > 1){
                    // multi cis
                    
                    // get number of IV effects
                    
                    for(Triple<String, String, String> t: triples){
                        if(t.getLeft().equals(snp)){
                            nrMultiCis++;
                            if(signtriples.contains(t)){
                                nrSig++;
                            }
                        }
                    }
                }
            }
            
            System.out.println(nrMultiCis);
            System.out.println(nrSig);
            
            tf.close();
        } catch (Exception e) {
        }
    }
}
