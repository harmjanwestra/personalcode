package nl.harmjanwestra.playground.biogen.freeze2dot1.vcf;


import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

public class LDSc {

    public static void main(String[] args) {
        String file = "D:\\TMP\\ldsc\\ENSG00000099785-mb2dot1ids.sumstats.qualdata.txt";
        LDSc s = new LDSc();
        s.run(file);
    }

    public void run(String file) {
        ArrayList<ArrayList<Double>> l2s = new ArrayList<ArrayList<Double>>();
        ArrayList<ArrayList<Double>> chi2s = new ArrayList<ArrayList<Double>>();


        try {
            TextFile tf = new TextFile(file, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            ArrayList<Pair<Double, Double>> vals = new ArrayList<Pair<Double, Double>>();
            while (elems != null) {
                double z = Double.parseDouble(elems[4]);
                double l2 = Double.parseDouble(elems[elems.length - 1]);
                vals.add(new Pair<Double, Double>(l2, z, Pair.SORTBY.LEFT));
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            Collections.sort(vals);

            int nrvals = vals.size();
            int nrPerBin = vals.size() / 25;
            System.out.println(nrvals);
            System.out.println(nrPerBin);

            for (int i = 0; i < vals.size(); i++) {

            }

        } catch (IOException e) {
            e.printStackTrace();
        }


    }
}
