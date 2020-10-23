package nl.harmjanwestra.playground.biogen.eqtlposthoc;

import it.unimi.dsi.fastutil.Hash;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class FTDCheck {

    public static void main(String[] args) {
        String file = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-01-FTDSNPs\\eur-woENA-woPCA\\eQTLs.txt.gz";
        String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-01-FTDSNPs\\wgsonly\\eQTLs-dscompare.txt";
        String ref = "D:\\tmp\\trans\\2020-05-26-Cortex-EUR-trans-AllGenes-eQTLs-crossMappingEQTLsRemoved-FDR0.05-withDiseaseAnnotation.txt.gz";

        file="D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-01-FTDSNPs\\eur-woENA-woPCA\\eQTLs.txt.gz";
        ref="D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-01-FTDSNPs\\eur-woENA\\eQTLs.txt.gz";
        out="D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-01-FTDSNPs\\compareEURwoENA-vsEURwoENAwoPCA.txt";

        FTDCheck c = new FTDCheck();
        try {
            c.run(file, out, ref);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String file, String outfile, String ref) throws IOException {

        HashMap<String, Double> refqtl = new HashMap<>();

        TextFile tf1 = new TextFile(ref, TextFile.R);
        tf1.readLine();
        String[] elems = tf1.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String gene = elems[4];
            double z = Double.parseDouble(elems[10]);
            refqtl.put(snp + ";" + gene, z);
            elems = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();


        TextFile tf = new TextFile(file, TextFile.R);
        tf.readLine();
        elems = tf.readLineElems(TextFile.tab);

        HashSet<String> genes = new HashSet<String>();
        HashMap<String, HashMap<String, Double>> eqtlspersnp = new HashMap<>();

        while (elems != null) {
            String snp = elems[1];
            String gene = elems[4];
            double z = Double.parseDouble(elems[10]);
            genes.add(gene);
            HashMap<String, Double> eqtls = eqtlspersnp.get(snp);
            if (eqtls == null) {
                eqtls = new HashMap<>();
            }
            eqtls.put(gene, z);
            eqtlspersnp.put(snp, eqtls);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile out = new TextFile(outfile, TextFile.W);
        String header = "Gene";
        ArrayList<String> keys = new ArrayList<String>();
        keys.addAll(eqtlspersnp.keySet());
        for (String key : keys) {
            header += "\t" + key + "-Z\t" + key + "-refz";
        }
        out.writeln(header);

        for (String gene : genes) {
            String ln = gene;
            for (String key : keys) {
                HashMap<String, Double> eqtls = eqtlspersnp.get(key);
                Double z = eqtls.get(gene);
                Double z2 = refqtl.get(key + ";" + gene);
                if (z2 == null) {
                    z2 = 0d;
                }
                ln += "\t" + z + "\t" + z2;
            }
            out.writeln(ln);
        }

        out.close();

    }
}

