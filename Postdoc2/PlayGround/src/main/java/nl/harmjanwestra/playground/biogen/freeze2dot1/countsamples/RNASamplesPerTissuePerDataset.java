package nl.harmjanwestra.playground.biogen.freeze2dot1.countsamples;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class RNASamplesPerTissuePerDataset {

    public static void main(String[] args) {
        String tissueclass = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-17-tissues\\plots\\2020-02-06-eigenvectors-tissueClassification_residualSampleAssignmentAll.txt";
        String datasetclass = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-17-SamplesBeforeQC\\SampleToDataset.txt";
        String samplelist = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-17-SamplesBeforeQC\\2020-02-06-step4-remove-residual-covariates-pc1_10_residual.txt";
        RNASamplesPerTissuePerDataset s = new RNASamplesPerTissuePerDataset();
        try {
            s.run(datasetclass, tissueclass, samplelist);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String sampleToDsFile, String sampleToTissueFile, String sampleListFile) throws IOException {
        HashMap<String, String> sampleToDs = new HashMap<String, String>();
        HashMap<String, Integer> datasets = new HashMap<String, Integer>();
        ArrayList<String> datasetList = new ArrayList<String>();
        TextFile tf = new TextFile(sampleToDsFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            sampleToDs.put(elems[0], elems[1]);
            if (!datasets.containsKey(elems[1])) {
                datasets.put(elems[1], ctr);
                datasetList.add(elems[1]);
                ctr++;
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        HashMap<String, Integer> tissues = new HashMap<String, Integer>();
        ArrayList<String> tissueList = new ArrayList<String>();
        HashMap<String, String> sampleToTissue = new HashMap<String, String>();
        tf = new TextFile(sampleToTissueFile, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        ctr = 0;
        while (elems != null) {
            sampleToTissue.put(elems[0], elems[1]);
            if (!tissues.containsKey(elems[1])) {
                tissues.put(elems[1], ctr);
                tissueList.add(elems[1]);
                ctr++;
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        int[][] table = new int[datasets.size()][tissues.size()];
        tf = new TextFile(sampleListFile, TextFile.R);
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String sample = elems[0];

            String ds = sampleToDs.get(sample);
            String tissue = sampleToTissue.get(sample);
            if (ds == null) {
                System.out.println(sample + " has no dataset");
            }
            if (tissue == null) {
                System.out.println(sample + " has no tissue");
            }
            if (ds != null && tissue != null) {
                Integer id1 = tissues.get(tissue);
                Integer id2 = datasets.get(ds);
                table[id2][id1]++;
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        String header = "-";
        for (String s : tissueList) {
            header += "\t" + s;
        }
        System.out.println(header);
        for (int d = 0; d < datasetList.size(); d++) {
            String ln = datasetList.get(d);
            for (int t = 0; t < tissueList.size(); t++) {
                ln += "\t" + table[d][t];
            }
            System.out.println(ln);
        }
    }
}
