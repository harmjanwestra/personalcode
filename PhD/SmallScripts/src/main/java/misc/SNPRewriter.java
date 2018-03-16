/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class SNPRewriter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String snpmapfile = "/Volumes/iSnackHD/AeroFS/AnnotationFiles/2013-08-05-SNPMappings_Filtered.txt.gz";
            String snpfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/EGCUT1000GenomesImputed/SNPs-Original2.txt.gz";
            String out = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/EGCUT1000GenomesImputed/";
            String snpsplitz0r = ":";
            String newsnpsplitz0r = ":";

            TextFile tf = new TextFile(snpfile, TextFile.R);
            Set<String> snps = tf.readAsSet(0, TextFile.tab);
            tf.close();
            System.out.println(snps.size() + "SNPs loaded.");

            TextFile tf2 = new TextFile(snpmapfile, TextFile.R);
            String[] mapElems = tf2.readLineElems(TextFile.tab);
            HashMap<String, String> posToRs = new HashMap<String, String>();
            while (mapElems != null) {
                String snpAsInDs = new String(mapElems[0] + snpsplitz0r + mapElems[1]);
                String rs = new String(mapElems[2]);
                if (snps.contains(snpAsInDs)) {
                    posToRs.put(snpAsInDs, mapElems[2]);
                }
                if (snps.contains(rs)) {
                    posToRs.put(rs, snpAsInDs);
                }
                mapElems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            System.out.println(posToRs.size() + " positions read.");

            HashSet<String> visitedRSNumbers = new HashSet<String>();
            TextFile out1 = new TextFile(out + "SNPs.txt.gz", TextFile.W);
            TextFile out2 = new TextFile(out + "SNPMappings.txt.gz", TextFile.W);
            tf.open();

            String snp = tf.readLine();
            while (snp != null) {

                if (visitedRSNumbers.contains(snp)) {

                    while (visitedRSNumbers.contains(snp)) {
                        snp += "_dup";
                    }
                    out1.writeln(snp);
                    visitedRSNumbers.add(snp);
                    out2.writeln(0 + "\t" + 0 + "\t" + snp);

                } else {
                    String[] snpElems = snp.split(snpsplitz0r);
                    if (snpElems.length == 2) {
                        String rs = posToRs.get(snp);
                        if (rs == null || visitedRSNumbers.contains(rs)) {
                            out1.writeln(snpElems[0] + newsnpsplitz0r + snpElems[1]);
                            visitedRSNumbers.add(snp);
                        } else {
                            out1.writeln(rs);
                            visitedRSNumbers.add(rs);
                        }
                        out2.writeln(snpElems[0] + "\t" + snpElems[1] + "\t" + snpElems[0] + newsnpsplitz0r + snpElems[1]);
                    } else {
                        System.out.println("ERROR could not parse snp: " + snp);

                        if (visitedRSNumbers.contains(snp)) {
                            while (visitedRSNumbers.contains(snp)) {
                                snp += "_dup";
                            }
                        }

                        String rs = posToRs.get(snp);
                        if (rs != null && snp.startsWith("rs")) {
                            snpElems = rs.split(snpsplitz0r);
                            out1.writeln(snp);
                            System.out.println("Resolved: " + snpElems[0] + "\t" + snpElems[1] + "\t" + snp);
                            out2.writeln(snpElems[0] + "\t" + snpElems[1] + "\t" + snp);
                        } else {
                            out1.writeln(snp);
                            out2.writeln(0 + "\t" + 0 + "\t" + snp);
                        }
                        visitedRSNumbers.add(snp);
                    }
                }
                snp = tf.readLine();
            }
            tf.close();

            out1.close();
            out2.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
