package nl.harmjanwestra.playground.biogen.freeze2dot1;


import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashSet;

public class CountEQTLs {


    public static void main(String[] args) {

        String cisefile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\lude\\2020-05-26-Cortex-EUR-cis-eQTLsFDR0.05-ProbeLevel-Iteration1-4-diseaseSNPsAndTopFx-diseaseAnnotation.txt.gz";
        String transefile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\lude\\2020-05-26-Cortex-EUR-trans-eQTLsFDR0.05-withDiseaseAnnotation.txt.gz";
        String output = "";
        String gwaslist = "U:\\IEUGWAS\\2020-06-01-2020-05-03-gwaslist-wgwascatalog-wALS-wMetaBrain.txt.gz";
        String query = "schizo";
        String gwasassocfile = "U:\\IEUGWAS\\2020-06-01-2020-05-03-allTopAssociations-wgwascatalog-wALS-wMetaBrain-MetaBrain2dot1IDs.txt.gz";
        CountEQTLs c = new CountEQTLs();
        try {
            c.countEQTLs(cisefile, transefile, gwaslist, gwasassocfile, query, output);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void countEQTLs(String cisefile, String transefile, String gwaslist, String gwasassocfile, String query, String output) throws IOException {

        HashSet<String> allowedIds = new HashSet<String>();
        TextFile tf = new TextFile(gwaslist, TextFile.R);
        String header = tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (elems[3].toLowerCase().contains(query)) {
                allowedIds.add(elems[0]);
                System.out.println(elems[0] + "\t" + elems[3]);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        HashSet<String> uniqueSNPs = new HashSet<String>();
        tf = new TextFile(gwasassocfile, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            if (allowedIds.contains(elems[0])) {
                uniqueSNPs.add(elems[1]);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        int[] ctrcis = new int[5];

        tf = new TextFile(cisefile, TextFile.R);
        header = tf.readLine();
        elems = tf.readLineElems(TextFile.tab);
        HashSet<String> seenSNPs = new HashSet<String>();
        while (elems != null) {
            String iter = elems[22];
            String iternr = iter.replaceAll("Iteration", "");
            Integer iteri = Integer.parseInt(iternr);
            String[] gwasIds = elems[23].split(";");
            for (String gwasId : gwasIds) {
                if (allowedIds.contains(gwasId)) {
                    ctrcis[iteri]++;
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        int ctrtrans = 0;
        tf = new TextFile(transefile, TextFile.R);
        header = tf.readLine();
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {


            String[] gwasIds = elems[27].split(";");
            for (String gwasId : gwasIds) {
                if (allowedIds.contains(gwasId)) {
                    ctrtrans++;
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(uniqueSNPs.size());
        System.out.println("Trans:\t" + ctrtrans);
        for (int i = 0; i < ctrcis.length; i++) {
            System.out.println("cis:\t" + i + "\t" + ctrcis[i]);
        }
    }

}
