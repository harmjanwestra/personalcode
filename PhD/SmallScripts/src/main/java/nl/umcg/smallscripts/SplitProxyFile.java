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
public class SplitProxyFile {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            SplitProxyFile.convert(
                        "/Volumes/iSnackHD/SkyDrive/Cis-SNPs/RealData/AllSNPs.txt-WithProxies.txt",
                        "/Volumes/iSnackHD/SkyDrive/Cis-SNPs/RealData/AllSNPs.txt");
//            for (int perm = 1; perm < 11; perm++) {
//                SplitProxyFile.convert(
//                        "/Volumes/iSnackHD/SkyDrive/Cis-SNPs/AllSNPsPerm.txt-WithProxies.txt",
//                        "/Volumes/iSnackHD/SkyDrive/Cis-SNPs/PermutationRound" + perm + "/AllSNPs.txt");
//            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void convert(String proxyFile, String fileToConvert) throws IOException {
        TextFile proxy = new TextFile(proxyFile, TextFile.R);
        HashMap<String, HashSet<String>> proxiesPerSNP = new HashMap<String, HashSet<String>>();
        String[] elems = proxy.readLineElems(TextFile.tab);
        while (elems != null) {
            String input = elems[0];
            String proxySNP = elems[1];
            HashSet<String> proxies = proxiesPerSNP.get(input);
            if (proxies == null) {
                proxies = new HashSet<String>();
            }
            proxies.add(proxySNP);
            elems = proxy.readLineElems(TextFile.tab);
        }
        proxy.close();

        TextFile tfOut = new TextFile(fileToConvert + "-SNPNexusFormat.txt", TextFile.W);
        TextFile tf = new TextFile(fileToConvert, TextFile.R);

        String snp = tf.readLine();
        HashSet<String> snpVisited = new HashSet<String>();
        while (snp != null) {
            if (!snpVisited.contains(snp)) {
                HashSet<String> proxies = proxiesPerSNP.get(snp);
                if (proxies != null) {
                    String[] proxiesArr = proxies.toArray(new String[0]);
                    for (String prox : proxiesArr) {
                        if (!snpVisited.contains(prox)) {
                            tfOut.writeln("dbsnp\t" + prox);
                            snpVisited.add(prox);
                        }
                    }

                }
                if (!snpVisited.contains(snp)) {
                    tfOut.writeln("dbsnp\t" + snp);
                }
                snpVisited.add(snp);
            }
            snp = tf.readLine();
        }

        tf.close();
        tfOut.close();


    }
}
