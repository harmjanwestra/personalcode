/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.iv;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.Triple;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harm-jan
 */
public class IVSNPProbePairFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String cisFile = "D:\\Skydrive\\latesteQTLs\\cisFDR0.05.txt.gz";
        String transFile = "D:\\Skydrive\\latesteQTLs\\transFDR0.05.txt.gz";
        String ivfile = "D:\\Skydrive\\latesteQTLs\\trans_ma_results_31072012.txt";
        String probetrans = "D:\\Skydrive\\MetaAnalysisAnnotationFiles\\2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
        String probetransSource = "HumanHT-12_V3_0_R2_11283641_A.txt";
        String probetransDest = "Probe";
        try {
            HashMap<String, String> probeMap = new ProbeTranslation().getProbeTranslation(probetrans, probetransSource, probetransDest);

            HashSet<Pair<String, String>> cisPairs = new HashSet<Pair<String, String>>();
            HashSet<Pair<String, String>> transPairs = new HashSet<Pair<String, String>>();

            HashSet<String> transSNPs = new HashSet<String>();
            TextFile trans = new TextFile(transFile, TextFile.R);
            String[] data = trans.readLineElems(TextFile.tab);
            data = trans.readLineElems(TextFile.tab);
            while (data != null) {
                transPairs.add(new Pair<String, String>(data[1], data[4]));
                transSNPs.add(data[1]);
                data = trans.readLineElems(TextFile.tab);
            }
            trans.close();

            System.out.println(transPairs.size() + " trans pairs");
            TextFile cis = new TextFile(cisFile, TextFile.R);
            data = cis.readLineElems(TextFile.tab);
            data = cis.readLineElems(TextFile.tab);
            HashSet<String> cisSNPs = new HashSet<String>();
            while (data != null) {
                double pval = Double.parseDouble(data[0]);
                if (transSNPs.contains(data[1])) {
                    if (pval <= 1.3141052405861965E-4) {
                        cisPairs.add(new Pair<String, String>(data[1], data[4]));
                        cisSNPs.add(data[1]);
                    }
                }
                data = cis.readLineElems(TextFile.tab);
            }
            cis.close();
            System.out.println(cisPairs.size() + " cis pairs");

            HashSet<Triple<String, String, String>> snpCisTransProbePairs = new HashSet<Triple<String, String, String>>();
            HashSet<String> cisTransSNPs = new HashSet<String>();
            for (String transSNP : transSNPs) {
                ArrayList<Pair<String, String>> transPairsForSNP = new ArrayList<Pair<String, String>>();
                ArrayList<Pair<String, String>> cisPairsForSNP = new ArrayList<Pair<String, String>>();
                for (Pair<String, String> pair : transPairs) {
                    if (pair.getLeft().equals(transSNP)) {
                        transPairsForSNP.add(pair);
                    }
                }

                for (Pair<String, String> pair : cisPairs) {
                    if (pair.getLeft().equals(transSNP)) {
                        cisPairsForSNP.add(pair);
                    }
                }
                for (Pair<String, String> pair1 : cisPairsForSNP) {
                    for (Pair<String, String> pair2 : transPairsForSNP) {
                        snpCisTransProbePairs.add(new Triple<String, String, String>(transSNP, pair1.getRight(), pair2.getRight()));
                        cisTransSNPs.add(transSNP);
                    }
                }
            }

            System.out.println(snpCisTransProbePairs.size() + " snp - cis - trans combinations");
            System.out.println(cisTransSNPs.size() + " cistrans SNPs");

            TextFile iv = new TextFile(ivfile, TextFile.R);
            TextFile ivout = new TextFile(ivfile + "-Filtered2.txt", TextFile.W);
            data = iv.readLineElems(TextFile.tab);
            ivout.writeln(Strings.concat(data, Strings.tab));
            data = iv.readLineElems(TextFile.tab);
            while (data != null) {
                String snp = data[0];

                if (cisTransSNPs.contains(snp)) {
                    String cisProbe = probeMap.get(data[1]);
                    String transProbe = probeMap.get(data[2]);
                    boolean comboFound = false;
                    for (Triple<String, String, String> t : snpCisTransProbePairs) {
                        if (t.getLeft().equals(snp) && t.getMiddle().equals(cisProbe) && t.getRight().equals(transProbe)) {
                            comboFound = true;
                        }
                    }
                    if (comboFound) {
                        // output
                        ivout.writeln(Strings.concat(data, Strings.tab));
                    }
                } else {
                    // don't output
                }
                data = iv.readLineElems(TextFile.tab);
            }
            ivout.close();
            iv.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
