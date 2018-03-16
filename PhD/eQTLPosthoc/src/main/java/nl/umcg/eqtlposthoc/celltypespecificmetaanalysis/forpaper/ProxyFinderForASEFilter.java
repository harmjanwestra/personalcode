/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ProxyFinderForASEFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String[] proxyfiles = new String[]{"/Volumes/iSnackHD/AeroFS/2013-12-09-ProxiesForASE/100GenomesCEU-Proxies-0.99.txt", "/Volumes/iSnackHD/AeroFS/2013-12-09-ProxiesForASE/HapMap2CEU-Proxies-0.99.txt"};
        String snpFile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap-Altschuler/SNPs.txt";
        String outputfile = "/Volumes/iSnackHD/AeroFS/2013-12-09-ProxiesForASE/Hap300Proxies.txt";

        try {
            ProxyFinderForASEFilter filter = new ProxyFinderForASEFilter();
            filter.filter(snpFile, proxyfiles, outputfile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void filter(String snpFile, String[] proxyfiles, String outputfile) throws IOException {
        HashMap<String, HashSet<String>> proxiesPerSNP = new HashMap<String, HashSet<String>>();
        for (String proxyfile : proxyfiles) {
            TextFile tf = new TextFile(proxyfile, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String inputSNP = elems[0];
                String proxies = elems[3];
                String[] proxyElems = proxies.split(";");
                HashSet<String> proxieSet = proxiesPerSNP.get(inputSNP);
                if (proxieSet == null) {
                    proxieSet = new HashSet<String>();
                }
                proxieSet.addAll(Arrays.asList(proxyElems));
                proxiesPerSNP.put(inputSNP, proxieSet);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
        }

        TextFile tf2 = new TextFile(snpFile, TextFile.R);
        Set<String> snpSelection = tf2.readAsSet(0, TextFile.tab);
        tf2.close();
        
        TextFile out = new TextFile(outputfile, TextFile.W);
        out.writeln("InputSNP\tIsOnHap300\tProxiesOnHap300");
        Set<String> keys = proxiesPerSNP.keySet();
        for(String s: keys){
            HashSet<String> proxies = proxiesPerSNP.get(s);
            HashSet<String> proxiesOnHap300 = new HashSet<String>();
            for(String p: proxies){
                if(snpSelection.contains(p)){
                    proxiesOnHap300.add(p);
                }
            }
            String outputStr = s+"\t"+snpSelection.contains(s)+"\t"+Strings.concat(proxiesOnHap300.toArray(new String[0]), Strings.semicolon);
            out.writeln(outputStr);
        }
        out.close();

    }

}
