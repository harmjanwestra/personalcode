package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.*;

public class SplitExpressionDatasets {

    public static void main(String[] args) {


        if (args.length < 3) {
            System.out.println("Usage: input.txt sampleToDs.txt outputdir");
        } else {
            SplitExpressionDatasets exp = new SplitExpressionDatasets();
            String input = args[0];
            String sampleToDs = args[1];
            String output = args[2];
            try {
                exp.run(input, sampleToDs, output);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void run(String input, String sampleToDs, String output) throws IOException {

        System.out.println("Splitting " + input);
        System.out.println("Using: " + sampleToDs);
        System.out.println("Saving output here: " + output);

        Gpio.createDir(output);

        HashMap<String, String> sampleToDsMap = new HashMap<String, String>();
        HashSet<String> datasets = new HashSet<String>();
        HashMap<String, HashSet<String>> dsToSample = new HashMap<>();

        TextFile tf = new TextFile(sampleToDs, TextFile.R);
        String[] elems = tf.readLineElems(Strings.tab);
        while (elems != null) {

            sampleToDsMap.put(elems[0], elems[1]);
            datasets.add(elems[1]);

            HashSet<String> samplelist = dsToSample.get(elems[1]);
            if (samplelist == null) {
                samplelist = new HashSet<>();
            }
            samplelist.add(elems[0]);

            dsToSample.put(elems[1], samplelist);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        System.out.println(datasets.size() + " datasets");

        TextFile[] tfouts = new TextFile[datasets.size()];
        String[][] outputbuffers = new String[datasets.size()][];
        ArrayList<HashMap<String, Integer>> dsSampleMaps = new ArrayList<>();

        ArrayList<String> dslist = new ArrayList<>();


        dslist.addAll(datasets);
        HashMap<String, Integer> dsMap = new HashMap<>();
        TextFile tfi = new TextFile(input, TextFile.R);
        String[] header = tfi.readLineElems(TextFile.tab);

        HashSet<String> allPresentSamples = new HashSet<String>();
        allPresentSamples.addAll(Arrays.asList(header));

        int totalmatched = 0;
        for (int i = 0; i < dslist.size(); i++) {

            dsMap.put(dslist.get(i), i);
            // determine headers
            HashSet<String> samplehash = dsToSample.get(dslist.get(i));

            ArrayList<String> sampleList = new ArrayList<>();
            for (String s : samplehash) {
                if (allPresentSamples.contains(s)) {
                    sampleList.add(s);
                }
            }

            Collections.sort(sampleList);

            HashMap<String, Integer> dsSampleMap = new HashMap<>();
            for (int j = 0; j < sampleList.size(); j++) {
                dsSampleMap.put(sampleList.get(j), j);
            }

            System.out.println(dslist.get(i) + "\t" + sampleList.size() + " samples");
            totalmatched += sampleList.size();
            if (!sampleList.isEmpty()) {
                tfouts[i] = new TextFile(output + dslist.get(i) + ".txt.gz", TextFile.W);
                tfouts[i].writeln("Gene\t" + Strings.concat(sampleList, Strings.tab));
                outputbuffers[i] = new String[dsSampleMap.size()];

            } else {

            }
            dsSampleMaps.add(dsSampleMap);
        }


        System.out.println(header.length + " columns in header");
        System.out.println(totalmatched + " columns matched.");
        if (totalmatched > 0) {
            int[] dsindex = new int[header.length];
            int[] sampleindex = new int[header.length];

            for (int i = 1; i < header.length; i++) {
                String sampleName = header[i];
                String ds = sampleToDsMap.get(sampleName);

                Integer dsId = dsMap.get(ds);
                if (dsId == null) {
                    dsindex[i] = -1;
                    System.out.println("Dataset not found for " + sampleName);
                } else {
                    dsindex[i] = dsId;
                    HashMap<String, Integer> dsSampleMap = dsSampleMaps.get(dsId);
                    Integer sampleId = dsSampleMap.get(sampleName);
                    sampleindex[i] = sampleId;
                }

            }


            String[] elemsi = tfi.readLineElems(TextFile.tab);
            int lctr = 0;
            while (elemsi != null) {

                for (int i = 0; i < outputbuffers.length; i++) {
                    if (outputbuffers[i] != null) {
                        outputbuffers[i][0] = elemsi[0];
                    }
                }

                for (int i = 1; i < elemsi.length; i++) {
                    int dsi = dsindex[i];
                    if (dsi >= 0) {
                        int sai = sampleindex[i];
                        outputbuffers[dsi][sai] = elemsi[i];
                    }
                }

                for (int i = 0; i < outputbuffers.length; i++) {
                    if (tfouts[i] != null) {
                        tfouts[i].writeln(elemsi[0] + "\t" + Strings.concat(outputbuffers[i], Strings.tab));
                    }
                }


                lctr++;
                if (lctr % 100 == 0) {
                    System.out.print("\r" + lctr + " lines parsed.");
                }
                elemsi = tfi.readLineElems(TextFile.tab);
            }
            tfi.close();

            for (int i = 0; i < dslist.size(); i++) {
                if (tfouts[i] != null) {
                    tfouts[i].close();
                }
            }
            System.out.println("");
            System.out.println("Done.");
            System.out.println();
        }

    }

}
