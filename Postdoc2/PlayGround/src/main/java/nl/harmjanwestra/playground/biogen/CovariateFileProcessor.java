package nl.harmjanwestra.playground.biogen;


import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class CovariateFileProcessor {

    public static void main(String[] args) {

        CovariateFileProcessor cvp = new CovariateFileProcessor();

        String linkfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-GTE-unique.txt-allcombos.txt";
        int from = 3;
        int to = 1;
        String covariatefile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-CategoricalMetaData-ExternalIDs.txt";
        String covariateoutfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-CategoricalMetaData-RNAIDs.txt";

//        try {
//            cvp.replaceSampleIds(linkfile, from, to, covariatefile, covariateoutfile);
//        } catch (IOException e) {
//            e.printStackTrace();
//        }

        linkfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\TargetALS-GTE-unique.txt-allcombos.txt";
        from = 0;
        to = 1;
        covariatefile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\mds\\plinkMds.txt";
        covariateoutfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\mds\\plinkMds-RNAIds.txt";

        try {
            cvp.replaceSampleIdsWithDupReplacements(linkfile, from, to, covariatefile, covariateoutfile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void replaceSampleIds(String linkfile, int from, int to, String covariatefile, String covariateoutfile) throws IOException {
        TextFile tf = new TextFile(linkfile, TextFile.R);
        HashMap<String, String> sampleconv = new HashMap<String, String>();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            sampleconv.put(elems[from], elems[to]);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile covout = new TextFile(covariateoutfile, TextFile.W);
        TextFile covin = new TextFile(covariatefile, TextFile.R);
        covout.writeln(covin.readLine());
        elems = covin.readLineElems(TextFile.tab);
        while (elems != null) {

            String sample = elems[0];
            String sampleout = sampleconv.get(sample);
            if (sampleout == null) {
                System.out.println("could not find sample in link file: " + sample);

            } else {
                elems[0] = sampleout;
            }

            covout.writeln(Strings.concat(elems, Strings.tab));
            elems = covin.readLineElems(TextFile.tab);
        }
        covout.close();

    }

    public void replaceSampleIdsWithDupReplacements(String linkfile, int from, int to, String covariatefile, String covariateoutfile) throws IOException {
        TextFile tf = new TextFile(linkfile, TextFile.R);
        HashMap<String, ArrayList<String>> sampleconv = new HashMap<String, ArrayList<String>>();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            ArrayList<String> list = sampleconv.get(elems[from]);
            if (list == null) {
                list = new ArrayList<>();
            }
            list.add(elems[to]);

            sampleconv.put(elems[from], list);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile covout = new TextFile(covariateoutfile, TextFile.W);
        TextFile covin = new TextFile(covariatefile, TextFile.R);
        covout.writeln(covin.readLine());
        elems = covin.readLineElems(TextFile.tab);
        while (elems != null) {

            String sample = elems[0];
            ArrayList<String> sampleout = sampleconv.get(sample);
            if (sampleout == null) {
                System.out.println("could not find sample in link file: " + sample);

            } else {
                for (String s : sampleout) {
                    elems[0] = s;
                    covout.writeln(Strings.concat(elems, Strings.tab));
                }

            }


            elems = covin.readLineElems(TextFile.tab);
        }
        covout.close();

    }


}
