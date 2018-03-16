/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.zscoretables;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ZScoreTTestNullDist {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String indir = "/Volumes/iSnackHD/MetaAnalysisFinal/cistrans/2012-05-30-Cistrans/";
            String outdir = indir + "PathwayAnnotation/AllPathwaysVsPCSNP/";
            String distdir = outdir + "dists";
            Gpio.createDir(distdir);
            ArrayList<Double> vals = new ArrayList<Double>();
            int[] realbins = null;
            for (int i = 0; i < 8; i++) {
                String infile = outdir + "Reactome-Plus-TTest.txt";
                if (i > 0) {
                    infile = outdir + "Reactome-Plus-TTest-Permutation" + i + ".txt";
                }
                TextFile f = new TextFile(infile, TextFile.R);
                System.out.println("Loading: " + infile);
                f.readLine();
                String[] elems = f.readLineElems(TextFile.tab);

                while (elems != null) {
                    Double pval = Double.parseDouble(elems[6]);
//                    System.out.println(pval + "\t" + (-Math.log(pval)));
                    double logp = -Math.log(pval);
                    if (logp > 10) {
                        System.out.println(Strings.concat(elems, Strings.tab));
                    }
                    vals.add(logp);
                    elems = f.readLineElems(TextFile.tab);
                }
                f.close();

                // real data distribution
                if (i == 0) {
//                    TextFile out = new TextFile(distdir + "/RealDist.txt", TextFile.W);
                    Collections.sort(vals);
                    int[] bins = convertToBins(vals, 0, 200, 20000);
                    realbins = bins;
//                    for (int q = 0; q < bins.length; q++) {
//                        out.writeln((q * ((double) 200 / 20000)) + "\t" + bins[q]);
//                    }
                    vals = new ArrayList<Double>();
//                    out.close();
                }
            }
            Collections.sort(vals);
            TextFile out = new TextFile(distdir + "/PermDist.txt", TextFile.W);
            Collections.sort(vals);
            int[] bins = convertToBins(vals, 0, 200, 20000);
            for (int q = 0; q < bins.length; q++) {
                double d = (q * ((double) 200 / 20000));
                out.writeln(Math.exp(d) + "\t" + d + "\t" + realbins[q] + "\t" + bins[q]);
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static int[] convertToBins(ArrayList<Double> vals, int min, int max, int numbins) {
        int[] bins = new int[numbins];
        // range = (max - min) / numbins;
        double range = ((double) max - min) / numbins;
//        System.out.println(range);
        for (int i = 0; i < vals.size(); i++) {
            double val = vals.get(i);
            int binno = 0;
            if (val >= max) {
                binno = numbins - 1;
            } else if (val <= min) {
                binno = 0;
            } else {
                binno = (int) Math.floor((val / range));

            }
//            System.out.println(val + "\t" + binno);
            bins[binno]++;
        }
        return bins;
    }
}
