/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import umcg.genetica.containers.Pair;
import umcg.genetica.gwas.Independifier;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class VectorIndependifier {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String vector = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/2013-11-20-MetaAnalysisForPlots-ZScoreComparisons/Vector-CellTypeInteractionZScore-2.txt";
        String dataset = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12v3-1KGImputed/";
        String output = "/Volumes/iSnackHD/AeroFS/CellTypeEQTLDataFiles/GosiaFiles/";
        try {
            Independifier independifier = new Independifier(dataset);

            HashMap<String, Double> neutroSpecific = new HashMap<String, Double>();
            HashMap<String, Double> generic = new HashMap<String, Double>();
            HashMap<String, Double> lymphospecific = new HashMap<String, Double>();

            TextFile tf = new TextFile(vector, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                String snp = elems[0].split("-")[0];
                double z = Double.parseDouble(elems[1]);
                if (z < 0 && Math.abs(z) > 2.605579924) {
                    lymphospecific.put(snp, Math.abs(z));
                } else if (z > 2.605579924) {
                    neutroSpecific.put(snp, z);
                } else {
                    generic.put(snp, Math.abs(z));
                }

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            String[] snps = lymphospecific.keySet().toArray(new String[0]);
            String[] indpendifiedSNPs = independifier.independify(snps, 0.2, 1000000);
            writeMaxEffects(indpendifiedSNPs, lymphospecific, output+"lymphocyteSpecific-r0.2.txt");

            snps = neutroSpecific.keySet().toArray(new String[0]);
            indpendifiedSNPs = independifier.independify(snps, 0.2, 1000000);
            writeMaxEffects(indpendifiedSNPs, neutroSpecific,output+"neutrophilSpecific-r0.2.txt");

            snps = generic.keySet().toArray(new String[0]);
            indpendifiedSNPs = independifier.independify(snps, 0.2, 1000000);
            writeMaxEffects(indpendifiedSNPs, generic,output+"generic-r0.2.txt");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void writeMaxEffects(String[] indpendifiedSNPs, HashMap<String, Double> lymphospecific, String outputfilename) throws IOException {

        TextFile tf = new TextFile(outputfilename, TextFile.W);
        tf.writeln("TopSNP\tSNPFxZScore\tRemainingSNPsInLocus");
        for (String s : indpendifiedSNPs) {
            String[] snps = s.split(";");
            ArrayList<Pair<Double, String>> fxSorter = new ArrayList<Pair<Double, String>>();

            for (String snp : snps) {
                double fx = lymphospecific.get(snp);
                fxSorter.add(new Pair<Double, String>(Math.abs(fx), snp, Pair.SORTBY.LEFT));
            }
            Collections.sort(fxSorter);
            String topFx = fxSorter.get(0).getRight();
            String ln = topFx + "\t" + lymphospecific.get(topFx) + "\t";
            for (int p = 1; p < fxSorter.size(); p++) {
                if (p == 1) {
                    ln += fxSorter.get(p).getRight();
                } else {
                    ln += ";" + fxSorter.get(p).getRight();
                }
            }
            tf.writeln(ln);
        }
        tf.close();

    }

}
