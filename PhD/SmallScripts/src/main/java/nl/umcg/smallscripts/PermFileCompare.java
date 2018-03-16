/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class PermFileCompare {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            HashSet<Pair<String, String>> snpProbesInFirstFile = new HashSet<Pair<String, String>>();
            HashSet<String> snpsInFile1 = new HashSet<String>();
            HashSet<String> probesInFile1 = new HashSet<String>();

            for (int perm = 1; perm < 101; perm++) {
                String fileIn = "/Volumes/iSnackHD/Data/Projects/VaezBarzani/eqtlresults/2012-09-28-eQTLs-trans/PermutedEQTLsPermutationRound" + perm + ".txt.gz";


                HashSet<Pair<String, String>> snpProbesInThisFile = new HashSet<Pair<String, String>>();

                HashSet<String> snpsInFile2 = new HashSet<String>();
                HashSet<String> probesInFile2 = new HashSet<String>();


                TextFile tf = new TextFile(fileIn, TextFile.R);

                tf.readLine();
                String[] elems = tf.readLineElems(TextFile.tab);
                while (elems != null) {
                    String snp = elems[1];
                    String probe = elems[4];
                    Pair<String, String> p = new Pair<String, String>(snp, probe);
                    if (perm == 1) {
                        snpsInFile1.add(snp);
                        probesInFile1.add(probe);
                        
                        if(snpProbesInFirstFile.contains(p)){
                            System.out.println("Double pair: "+p);
                        } else {
                            snpProbesInFirstFile.add(p);
                        }
                    } else {
                        snpsInFile2.add(snp);
                        probesInFile2.add(probe);
                        snpProbesInThisFile.add(p);
                    }
                    elems = tf.readLineElems(TextFile.tab);
                }

                tf.close();

                // now compare both files
                if (perm > 1) {
                    int shared = 0;
                    int newPairs = 0;
                    for (Pair<String, String> p : snpProbesInThisFile) {
                        if (snpProbesInFirstFile.contains(p)) {
                            shared++;
                        } else {
                            newPairs++;
//                            System.out.println(p.toString() + "\tis new in file " + perm);
                        }
                    }
                    System.out.println("New: " + newPairs + "\tShared: " + shared);
                    int probesShared = 0;
                    int snpsShared = 0;

                    int newProbesF1 = 0;
                    int newProbesF2 = 0;

                    int newSNPsF1 = 0;
                    int newSNPsF2 = 0;

                    for (String p : probesInFile1) {
                        if (probesInFile2.contains(p)) {
                            probesShared++;
                        } else {
                            System.out.println(p);
                        }
                    }
                    for (String p : probesInFile2) {
                        if (!probesInFile1.contains(p)) {
                            System.out.println(p);
                        }
                    }

                    for (String p : snpsInFile1) {
                        if (snpsInFile2.contains(p)) {
                            snpsShared++;
                        } else {
                            System.out.println("1\t" + p);
                        }
                    }
                    for (String p : snpsInFile2) {
                        if (!snpsInFile1.contains(p)) {
                            System.out.println("2\t" + p);
                        }
                    }
                    System.out.println(perm + "\t" + snpsShared + "\t" + snpsInFile1.size() + "\t" + snpsInFile2.size() + "\t" + probesShared + "\t" + probesInFile1.size() + "\t" + probesInFile2.size() + "\t" + shared + "\t" + snpProbesInFirstFile.size() + "\t" +snpProbesInThisFile.size());
                } else {
                    System.out.println("Initial file contains: " + snpProbesInFirstFile.size());
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
