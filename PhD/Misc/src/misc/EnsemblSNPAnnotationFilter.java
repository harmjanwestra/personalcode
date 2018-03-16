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
public class EnsemblSNPAnnotationFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            int ctr = 0;

            TextFile snpFileOut = new TextFile("/Volumes/iSnackHD/Data/GenomeSequences/Ensembl70_HG19/snps/2013-07-18-SNPMappings.txt", TextFile.W);
            for (int chr = 1; chr < 26; chr++) {
                String file = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl70_HG19/snps/SNPs_chr" + chr + ".txt.gz";
                if (chr == 23) {
                    file = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl70_HG19/snps/SNPs_chrX.txt.gz";
                }
                if (chr == 24) {
                    file = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl70_HG19/snps/SNPs_chrY.txt.gz";
                }
                if (chr == 25) {
                    file = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl70_HG19/snps/SNPs_MT.txt.gz";
                }

                HashMap<String, Integer> snpToPosition = new HashMap<String, Integer>();
                HashSet<String> excludeSNPs = new HashSet<String>();
                HashSet<String> snps = new HashSet<String>();

                TextFile tf = new TextFile(file, TextFile.R);

                tf.readLine();//skip header;
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {

                    String snp = elems[0];
                    String pos = elems[3];
                    Integer posI = Integer.parseInt(pos);
                    Integer posOld = snpToPosition.get(snp);
                    if (posOld == null) {
                        snpToPosition.put(snp, posI);
                    } else {
                        if (posI.equals(posOld)) {
                            excludeSNPs.add(snp);
                        }
                    }
                    snps.add(snp);

                    if (snps.size() % 100000 == 0) {
                        System.out.println("Chr " + chr + " \t " + snps.size() + "\t" + excludeSNPs.size());
                    }

                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();


                String chrStr = "" + chr;
                if (chr == 23) {
                    chrStr = "X";
                }
                if (chr == 24) {
                    chrStr = "Y";
                }
                if (chr == 25) {
                    chrStr = "MT";
                }
                for (String snp : snps) {
                    if (!excludeSNPs.contains(snp)) {
                        snpFileOut.writeln(chr + "\t" + snpToPosition.get(snp) + "\t" + snp);
                        ctr++;
                    }
                }

                System.out.println("RUNNING TOTAL: " + ctr);

            }

            snpFileOut.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
