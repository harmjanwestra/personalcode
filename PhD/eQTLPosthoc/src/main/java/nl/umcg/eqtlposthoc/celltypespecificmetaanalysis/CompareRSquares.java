/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class CompareRSquares {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String file1 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/CisEffectsAndArrayQuality/2013-11-04-BloodHT12v3-WithCovariate-Parametric/eQTLProbesFDR0.05-ProbeLevel.txt";
            String file2 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/CisEffectsAndArrayQuality/2013-11-04-BloodHT12v3-WithRandomCovariate-Parametric/eQTLProbesFDR0.05-ProbeLevel.txt";
            String file3 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/CisEffectsAndArrayQuality/2013-11-04-BloodHT12v3-WithOutCovariate-Parametric/eQTLProbesFDR0.05-ProbeLevel.txt";
            String file4 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/CisEffectsAndArrayQuality/2013-11-04-BloodHT12v3-WithOutCovariateNormalModel-Parametric/eQTLs.txt.gz";
            String outfile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/CisEffectsAndArrayQuality/Comparison-ProbeLevel-RandomCovariate-ComparedWithNormalModel.txt";
            HashMap<Pair<String, String>, Double> eqtls1 = new HashMap<Pair<String, String>, Double>();
            HashMap<Pair<String, String>, Double> eqtls2 = new HashMap<Pair<String, String>, Double>();
            HashMap<Pair<String, String>, Double> eqtls3 = new HashMap<Pair<String, String>, Double>();

            TextFile tf = new TextFile(file1, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                Double r = Double.parseDouble(elems[17]);
                String snp = elems[1];
                String probe = elems[4];
                eqtls1.put(new Pair<String, String>(snp, probe), r);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            TextFile tf3 = new TextFile(file3, TextFile.R);
            tf3.readLine();
            String[] elems3 = tf3.readLineElems(TextFile.tab);
            while (elems3 != null) {
                Double r = Double.parseDouble(elems3[17]);
                String snp = elems3[1];
                String probe = elems3[4];
                eqtls2.put(new Pair<String, String>(snp, probe), r);
                elems3 = tf3.readLineElems(TextFile.tab);
            }
            tf3.close();

            TextFile tf4 = new TextFile(file4, TextFile.R);
            tf4.readLine();
            String[] elems4 = tf4.readLineElems(TextFile.tab);
            while (elems4 != null) {
                Double r = Double.parseDouble(elems4[17]);
                String snp = elems4[1];
                String probe = elems4[4];
                eqtls3.put(new Pair<String, String>(snp, probe), r);
                elems4 = tf4.readLineElems(TextFile.tab);
            }
            tf4.close();

            TextFile tf2 = new TextFile(file2, TextFile.R);
            tf2.readLine();
            String[] elems2 = tf2.readLineElems(TextFile.tab);
            int notshared = 0;
            int shared = 0;

            TextFile out = new TextFile(outfile, TextFile.W);

            out.writeln("snp\tprobe\tRandomCovariate\tCovariate\tNoCovariate");
            while (elems2 != null) {
                Double r = Double.parseDouble(elems2[17]);
                String snp = elems2[1];
                String probe = elems2[4];
                Double otherfx = eqtls1.get(new Pair<String, String>(snp, probe));
                Double otherfx2 = eqtls2.get(new Pair<String, String>(snp, probe));
                Double otherfx3 = eqtls3.get(new Pair<String, String>(snp, probe));
                if (otherfx != null && otherfx2 != null && otherfx3 != null) {
                    out.writeln(snp + "\t" + probe + "\t" + Math.abs(r) + "\t" + Math.abs(otherfx) + "\t" + Math.abs(otherfx2) + "\t" + (otherfx3 * otherfx3));
                    shared++;
                } else {
                    notshared++;
                }
                elems2 = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();

            out.close();
            System.out.println("");
            System.out.println(shared + "\t" + notshared);

        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
